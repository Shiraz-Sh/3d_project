#include <stdio.h>
#include <ctype.h>
#include <math.h>
#include <string.h>

#include "inc_irit/iritprsr.h"
#include "inc_irit/irit_sm.h"
#include "inc_irit/allocate.h"
#include "inc_irit/attribut.h"
#include "inc_irit/bool_lib.h" 
#include "inc_irit/geom_lib.h"
#include "inc_irit/cagd_lib.h"
#include "inc_irit/ip_cnvrt.h"
#include "prsr_loc.h"

/*
 * Create a polygon from 2D points and extrude it into a 3D solid.
 *
 * New behavior:
 * - A new helper `IritPrsrExtrude2DPointsToSolidDir` is provided to accept an
 *   arbitrary 3D extrusion vector (`Dir`). The routine will map the 2D points
 *   into the plane orthogonal to `Dir` (forming a local u/v basis) and then
 *   call the standard extruder with the provided vector.
 *
 * Notes:
 * - This is simple and robust: placing the 2D polygon in the plane
 *   perpendicular to the extrusion direction guarantees the cross-section is
 *   not coplanar with Dir (the extruder requires a non-zero dot product).
 */

IPObjectStruct* IritPrsrExtrude2DPointsToSolidDir(const IrtE2PtStruct* pts,
                                                  int n,
                                                  IrtVecType Dir)
{
    int i;
    IPPolygonStruct* Pl;
    IPVertexStruct* FirstV = NULL, * PrevV = NULL, * V;
    IPObjectStruct* PolyObj = NULL, * Extruded = NULL;
    IrtVecType w, u, vtmp;
    IrtRType dirLen;

    if (pts == NULL || n < 3)
        return NULL;

    /* Compute length of Dir and normalized w = Dir / |Dir| */
    dirLen = sqrt(IRIT_SQR(Dir[0]) + IRIT_SQR(Dir[1]) + IRIT_SQR(Dir[2]));
    if (dirLen <= IRIT_EPS) /* Degenerate direction */
        return NULL;

    w[0] = Dir[0] / dirLen;
    w[1] = Dir[1] / dirLen;
    w[2] = Dir[2] / dirLen;

    /* Choose a stable arbitrary vector not parallel to w */
    if (fabs(w[0]) < 0.9) {
        vtmp[0] = 1.0; vtmp[1] = 0.0; vtmp[2] = 0.0;
    }
    else {
        vtmp[0] = 0.0; vtmp[1] = 1.0; vtmp[2] = 0.0;
    }

    /* u = normalize(cross(vtmp, w)) */
    u[0] = vtmp[1] * w[2] - vtmp[2] * w[1];
    u[1] = vtmp[2] * w[0] - vtmp[0] * w[2];
    u[2] = vtmp[0] * w[1] - vtmp[1] * w[0];
    {
        IrtRType ul = sqrt(IRIT_SQR(u[0]) + IRIT_SQR(u[1]) + IRIT_SQR(u[2]));
        if (ul <= IRIT_EPS) /* extremely unlikely, but guard anyway */
            return NULL;
        u[0] /= ul; u[1] /= ul; u[2] /= ul;
    }

    /* v = cross(w, u) to complete orthonormal basis */
    vtmp[0] = w[1] * u[2] - w[2] * u[1];
    vtmp[1] = w[2] * u[0] - w[0] * u[2];
    vtmp[2] = w[0] * u[1] - w[1] * u[0];

    /* Create polygon and vertex list mapped into plane orthogonal to Dir. */
    Pl = IritPrsrAllocPolygon(0, NULL, NULL);

    for (i = 0; i < n; ++i) {
        /* Allocate a new vertex with no next for now. */
        V = IritPrsrAllocVertex2(NULL);

        /* Map 2D point into 3D in the (u,v) plane:
           P = pts[i].Pt[0] * u + pts[i].Pt[1] * v
           (polygon is centered at origin in this local frame).
        */
        V->Coord[0] = (IrtRType)(pts[i].Pt[0] * u[0] + pts[i].Pt[1] * vtmp[0]);
        V->Coord[1] = (IrtRType)(pts[i].Pt[0] * u[1] + pts[i].Pt[1] * vtmp[1]);
        V->Coord[2] = (IrtRType)(pts[i].Pt[0] * u[2] + pts[i].Pt[1] * vtmp[2]);

        /* Link vertex into circular list. */
        if (FirstV == NULL) {
            FirstV = PrevV = V;
        }
        else {
            PrevV->Pnext = V;
            PrevV = V;
        }
    }

    /* Close loop */
    PrevV->Pnext = FirstV;
    Pl->PVertex = FirstV;

    /* Update polygon plane and cleanup (optional). */
    IritPrsrUpdatePolyPlane(Pl);

    /* Create polygon object. */
    PolyObj = IritPrsrGenPOLYObject(Pl);
    IP_SET_POLYGON_OBJ(PolyObj); /* ensure type flags are set */

    /* Extrude polygon into a solid/mesh along the original (non-normalized) Dir.
       Use the passed vector magnitude so caller controls extrusion length. */
    Extruded = IritGeomPrimGenEXTRUDEObject(PolyObj, Dir, 3);

    /* Free intermediate polygon object (extruder returns new object). */
    IritPrsrFreeObject(PolyObj);

    return Extruded;
}


/*
 * Compute silhouettes for multiple views of a polygonal object.
 *
 * Parameters:
 *   PObj         - polygonal object (input, not modified).
 *   ViewMats     - pointer to an array of IrtHmgnMatType matrices (NumViews entries).
 *   NumViews     - number of view matrices in ViewMats.
 *   GridSize     - grid size for preprocessing (if <= 0 a default of 20 is used).
 *   UsePreprocess- nonzero to preprocess once and query many views, zero to compute directly.
 *
 * Return value:
 *   A LIST IPObjectStruct * whose elements are the silhouette (polyline) objects
 *   computed for each view. NULL on error or if nothing produced.
 *
 * Notes:
 * - If UsePreprocess is true we copy/regularize the input once, build the preprocess
 *   structure and free it when done.
 * - If UsePreprocess is false we copy/regularize the object for each view and call
 *   IritGeomSilExtractSilDirect.
 * - The function assumes polygonal input (IP_IS_POLY_OBJ). Caller is responsible for freeing
 *   the returned LIST object via IritPrsrFreeObject().
 */
IPObjectStruct* IritPrsrGetSilhouettesForViews(const IPObjectStruct* PObj,
    const IrtHmgnMatType* ViewMats,
    int NumViews,
    int GridSize,
    int UsePreprocess)
{
    int i;
    IPObjectStruct* PList = NULL;
    VoidPtr PrepSils = NULL;
    IPObjectStruct* WorkCopy = NULL;
    IPObjectStruct *Head = NULL, *Tail = NULL;
    int resolvedGrid = (GridSize > 0) ? GridSize : 20;

    if (PObj == NULL || ViewMats == NULL || NumViews <= 0)
        return NULL;

    if (!IP_IS_POLY_OBJ((IPObjectStruct*)PObj))
        return NULL;

    /* Preprocess once if requested. */
    if (UsePreprocess) {
        /* Copy the polygonal object and ensure regularization/adjacencies. */
        WorkCopy = IritPrsrCopyObject(NULL, (IPObjectStruct*)PObj, TRUE);
        if (WorkCopy == NULL)
            return NULL;

        IritPrsrOpenPolysToClosed(WorkCopy->U.Pl);
        IritBoolGenAdjacencies(WorkCopy);

        PrepSils = IritGeomSilPreprocessPolys(WorkCopy, resolvedGrid);

        /* We can free the copy after preprocessing because preprocess stores required data. */
        IritPrsrFreeObject(WorkCopy);
        WorkCopy = NULL;

        if (PrepSils == NULL) {
            return NULL;
        }
    }

    /* Collect silhouettes into a simple linked list (Pnext chain). */
    for (i = 0; i < NumViews; ++i) {
        IPObjectStruct* PSil = NULL;

        if (UsePreprocess) {
            PSil = IritGeomSilExtractSil(PrepSils, ViewMats[i]);
        }
        else {
            /* Prepare a regularized per-view copy and compute directly. */
            WorkCopy = IritPrsrCopyObject(NULL, (IPObjectStruct*)PObj, TRUE);
            if (WorkCopy == NULL)
                continue;

            IritPrsrOpenPolysToClosed(WorkCopy->U.Pl);
            IritBoolGenAdjacencies(WorkCopy);

            PSil = IritGeomSilExtractSilDirect(WorkCopy, ViewMats[i]);

            IritPrsrFreeObject(WorkCopy);
            WorkCopy = NULL;
        }

        if (PSil != NULL) {
            /* Append to linked list (use Pnext as next pointer). */
            PSil->Pnext = NULL;
            if (Head == NULL) {
                Head = Tail = PSil;
            }
            else {
                Tail->Pnext = PSil;
                Tail = PSil;
            }
        }
    }

    /* Free preprocess structure if created. */
    if (PrepSils)
        IritGeomSilProprocessFree(PrepSils);

    /* If nothing collected return NULL. */
    if (Head == NULL)
        return NULL;

    /* Convert linked list of objects to a proper LIST object. This function
       will allocate the IPObjectStruct list container and copy the pointers
       from the linked list into it. */
    PList = IritPrsrObjLnkListToListObject(Head);

    return PList;
}


/*
 * Project silhouette polygons into the view's local (u,v) plane.
 *
 * For each polygon in `Sil3D` a polygon object is created whose vertices
 * are the original polygon vertices projected into the plane orthogonal
 * to the view forward vector and placed back in world space:
 *   WorldProjectedPt = Origin + u * x + v * y
 *
 * Returns a LIST IPObjectStruct * whose elements are polygon objects
 * (IP_OBJ_POLY). Caller owns and must free the returned LIST and its members.
 */
IPObjectStruct* IritPrsrProjectSilhouetteToViewPolys(const IPObjectStruct* SilList,
                                                     const IrtHmgnMatType* ViewMats,
                                                     int NumViews)
{
    int vi;

    if (SilList == NULL || ViewMats == NULL || NumViews <= 0)
        return NULL;

    /* Expect a LIST object produced by IritPrsrGetSilhouettesForViews. */
    if (SilList->U.Lst.PObjList == NULL)
        return NULL;

    IPObjectStruct* ResList = IritPrsrGenLISTObject(NULL);
    if (ResList == NULL)
        return NULL;

    /* For each silhouette (list element) use the corresponding view matrix. */
    for (vi = 0; vi < NumViews; ++vi) {
        IPObjectStruct* Sil3D = IritPrsrListObjectGet((IPObjectStruct*)(void*)SilList, vi);
        if (Sil3D == NULL)
            continue;

        /* Extract view forward vector from the 3x3 rotation part (third column). */
        const IrtHmgnMatType* ViewMat = &ViewMats[vi];
        IrtVecType f;
        f[0] = (*ViewMat)[0][2];
        f[1] = (*ViewMat)[1][2];
        f[2] = (*ViewMat)[2][2];

        /* Normalize forward - if degenerate fallback to +Z. */
        {
            IrtRType flen = sqrt(IRIT_SQR(f[0]) + IRIT_SQR(f[1]) + IRIT_SQR(f[2]));
            if (flen <= IRIT_EPS) {
                f[0] = 0.0; f[1] = 0.0; f[2] = 1.0;
                flen = 1.0;
            }
            f[0] /= flen; f[1] /= flen; f[2] /= flen;
        }

        /* We'll collect all projected polygon components for this silhouette
           into a single polygon-chain (HeadPl). After processing all components
           we create one polygon object for the silhouette and append it to ResList. */
        IPPolygonStruct *HeadPl = NULL, *TailPl = NULL;

        /* Project each polygon component of this silhouette using ViewMat.
           Use IritGeomGenProjectionMat to project points onto the plane
           orthogonal to the view forward vector passing through the polygon centroid.
        */
        for (IPPolygonStruct* Pl = Sil3D->U.Pl; Pl != NULL; Pl = Pl->Pnext) {
            IPVertexStruct* V = Pl->PVertex;
            if (V == NULL)
                continue;

            /* Count vertices (handle circular lists). */
            int n = 0;
            IPVertexStruct* cur = V;
            int guard = 0;
            while (cur != NULL) {
                ++n;
                cur = cur->Pnext;
                if (++guard > 200000) { n = 0; break; }
                if (cur == V) break;
            }
            if (n < 1)
                continue;

            /* Collect coordinates and compute centroid. */
            IrtPtType* coords = (IrtPtType*)IritMalloc(sizeof(IrtPtType) * n);
            if (coords == NULL)
                continue;

            cur = V;
            int idx = 0;
            guard = 0;
            IrtRType cx = 0.0, cy = 0.0, cz = 0.0;
            while (cur != NULL && idx < n) {
                coords[idx][0] = cur->Coord[0];
                coords[idx][1] = cur->Coord[1];
                coords[idx][2] = cur->Coord[2];
                cx += coords[idx][0];
                cy += coords[idx][1];
                cz += coords[idx][2];
                ++idx;
                cur = cur->Pnext;
                if (++guard > 200000) break;
                if (cur == V) break;
            }
            if (idx != n) { IritFree(coords); continue; }

            IrtPtType origin;
            origin[0] = cx / n; origin[1] = cy / n; origin[2] = cz / n;

            /* Build projection matrix that projects along direction 'f' onto
               the plane whose normal is 'f' and that passes through 'origin'. */
            IrtPlnType ProjPlane;
            IrtRType EyePos[4];
            IrtHmgnMatType PMat;

            /* Plane equation Ax + By + Cz + D = 0, with (A,B,C) = f and D = -dot(f, origin). */
            ProjPlane[0] = f[0];
            ProjPlane[1] = f[1];
            ProjPlane[2] = f[2];
            ProjPlane[3] = - (f[0] * origin[0] + f[1] * origin[1] + f[2] * origin[2]);

            /* For parallel projection along direction f, set EyePos.w = 0 and EyePos.xyz = f. */
            EyePos[0] = f[0];
            EyePos[1] = f[1];
            EyePos[2] = f[2];
            EyePos[3] = 0.0;

            /* Build projection matrix (ray from EyePos intersects ProjPlane). */
            IritGeomGenProjectionMat(ProjPlane, EyePos, PMat);

            /* Create polygon structure mapped back into world, but lying in view plane. */
            IPPolygonStruct* NewPl = IritPrsrAllocPolygon(0, NULL, NULL);
            IPVertexStruct* FirstV = NULL;
            IPVertexStruct* PrevV = NULL;

            for (int j = 0; j < n; ++j) {
                IPVertexStruct* NV = IritPrsrAllocVertex2(NULL);

                /* Project the world-space vertex into the view plane using the projection matrix. */
                IrtPtType Pproj;
                IritMiscMatMultPtby4by4(Pproj, coords[j], PMat);

                NV->Coord[0] = Pproj[0];
                NV->Coord[1] = Pproj[1];
                NV->Coord[2] = Pproj[2];

                if (FirstV == NULL) {
                    FirstV = PrevV = NV;
                }
                else {
                    PrevV->Pnext = NV;
                    PrevV = NV;
                }
            }

            if (PrevV != NULL) {
                PrevV->Pnext = FirstV;
                NewPl->PVertex = FirstV;
                IritPrsrUpdatePolyPlane(NewPl);

                /* Append NewPl into this silhouette's polygon chain (do not yet create a POLY object). */
                NewPl->Pnext = NULL;
                if (HeadPl == NULL) {
                    HeadPl = TailPl = NewPl;
                }
                else {
                    TailPl->Pnext = NewPl;
                    TailPl = NewPl;
                }
            }
            else {
                IritPrsrFreePolygon(NewPl);
            }

            IritFree(coords);
        } /* per polygon */

        /* If we collected any projected polygon components for this silhouette,
           wrap them into a single polygon object and append to ResList so each
           silhouette yields one 2D object. */
        if (HeadPl != NULL) {
            IPObjectStruct* PolyObj = IritPrsrGenPOLYObject(HeadPl);
            IP_SET_POLYGON_OBJ(PolyObj);
            PolyObj->Pnext = NULL;
            IritPrsrListObjectAppend(ResList, PolyObj);
        }
    } /* per silhouette/list element */

    /* Return NULL if nothing appended. */
    if (ResList->U.Lst.PObjList == NULL || ResList->U.Lst.ListMaxLen == 0) {
        IritPrsrFreeObject(ResList);
        return NULL;
    }

    return ResList;
}

/*
 * Extrude projected view polygons (as returned by IritPrsrProjectSilhouetteToViewPolys)
 * along the view forward direction and return a LIST of resulting solids.
 *
 * Each polygon in ProjPolysList is expected to be a polygon object whose
 * vertices lie in the view plane in world coordinates. For each polygon the
 * function will:
 *  - compute local (u,v) coordinates relative to polygon centroid,
 *  - call IritPrsrExtrude2DPointsToSolidDir with Dir = view_forward * Depth,
 *  - translate the extruded solid to the polygon centroid,
 *  - append the solid to the result LIST.
 */
IPObjectStruct* IritPrsrExtrudeProjectedViewPolys(const IPObjectStruct* ProjPolysList,
    const IrtHmgnMatType ViewMat,
    CagdRType Depth)
{
    if (ProjPolysList == NULL ||
        /* Ensure this is a LIST object with at least one element pointer array. */
        ProjPolysList->U.Lst.PObjList == NULL ||
        Depth == 0.0)
        return NULL;

    /* Extract forward vector from ViewMat (third column). */
    IrtVecType f;
    f[0] = ViewMat[0][2];
    f[1] = ViewMat[1][2];
    f[2] = ViewMat[2][2];
    {
        IrtRType flen = sqrt(IRIT_SQR(f[0]) + IRIT_SQR(f[1]) + IRIT_SQR(f[2]));
        if (flen <= IRIT_EPS) {
            f[0] = 0.0; f[1] = 0.0; f[2] = 1.0;
            flen = 1.0;
        }
        f[0] /= flen; f[1] /= flen; f[2] /= flen;
    }

    IPObjectStruct* ResList = IritPrsrGenLISTObject(NULL);
    if (ResList == NULL)
        return NULL;

    /* Iterate over list elements (each should be a polygon object). */
    int i = 0;
    IPObjectStruct* PObj;
    while ((PObj = IritPrsrListObjectGet((IPObjectStruct*)(void*)ProjPolysList, i++)) != NULL) {
        if (!IP_IS_POLY_OBJ(PObj))
            continue;

        /* For each polygon component in this object */
        for (IPPolygonStruct* Pl = PObj->U.Pl; Pl != NULL; Pl = Pl->Pnext) {
            IPVertexStruct* V = Pl->PVertex;
            if (V == NULL)
                continue;

            /* Count vertices and collect coords */
            int n = 0;
            IPVertexStruct* cur = V;
            int guard = 0;
            while (cur != NULL) {
                ++n;
                cur = cur->Pnext;
                if (++guard > 200000) { n = 0; break; }
                if (cur == V) break;
            }
            if (n < 3)
                continue;

            IrtPtType* coords = (IrtPtType*)IritMalloc(sizeof(IrtPtType) * n);
            if (coords == NULL)
                continue;

            cur = V;
            int idx = 0;
            guard = 0;
            IrtRType cx = 0.0, cy = 0.0, cz = 0.0;
            while (cur != NULL && idx < n) {
                coords[idx][0] = cur->Coord[0];
                coords[idx][1] = cur->Coord[1];
                coords[idx][2] = cur->Coord[2];
                cx += coords[idx][0];
                cy += coords[idx][1];
                cz += coords[idx][2];
                ++idx;
                cur = cur->Pnext;
                if (++guard > 200000) break;
                if (cur == V) break;
            }
            if (idx != n) { IritFree(coords); continue; }

            IrtPtType origin;
            origin[0] = cx / n; origin[1] = cy / n; origin[2] = cz / n;

            /* Build local basis u,v from view forward f (w). */
            IrtVecType w, u, vtmp;
            w[0] = f[0]; w[1] = f[1]; w[2] = f[2];

            IrtVecType arb;
            if (fabs(w[0]) < 0.9) { arb[0] = 1.0; arb[1] = 0.0; arb[2] = 0.0; }
            else { arb[0] = 0.0; arb[1] = 1.0; arb[2] = 0.0; }

            u[0] = arb[1] * w[2] - arb[2] * w[1];
            u[1] = arb[2] * w[0] - arb[0] * w[2];
            u[2] = arb[0] * w[1] - arb[1] * w[0];
            {
                IrtRType ul = sqrt(IRIT_SQR(u[0]) + IRIT_SQR(u[1]) + IRIT_SQR(u[2]));
                if (ul <= IRIT_EPS) { IritFree(coords); continue; }
                u[0] /= ul; u[1] /= ul; u[2] /= ul;
            }

            vtmp[0] = w[1] * u[2] - w[2] * u[1];
            vtmp[1] = w[2] * u[0] - w[0] * u[2];
            vtmp[2] = w[0] * u[1] - w[1] * u[0];

            /* Build 2D pts array relative to origin */
            IrtE2PtStruct* pts2d = (IrtE2PtStruct*)IritMalloc(sizeof(IrtE2PtStruct) * n);
            if (pts2d == NULL) { IritFree(coords); continue; }

            for (int j = 0; j < n; ++j) {
                IrtRType dx = coords[j][0] - origin[0];
                IrtRType dy = coords[j][1] - origin[1];
                IrtRType dz = coords[j][2] - origin[2];

                pts2d[j].Pt[0] = dx * u[0] + dy * u[1] + dz * u[2];
                pts2d[j].Pt[1] = dx * vtmp[0] + dy * vtmp[1] + dz * vtmp[2];
            }

            /* Extrude using view forward scaled by depth */
            IrtVecType Dir;
            Dir[0] = (IrtRType)(w[0] * Depth);
            Dir[1] = (IrtRType)(w[1] * Depth);
            Dir[2] = (IrtRType)(w[2] * Depth);

            IPObjectStruct* Extr = IritPrsrExtrude2DPointsToSolidDir(pts2d, n, Dir);

            IritFree(pts2d);
            IritFree(coords);

            if (Extr == NULL)
                continue;

            /* Translate extruded solid to origin */
            IrtHmgnMatType T;
            IritMiscMatGenUnitMat(T);
            T[0][3] = origin[0];
            T[1][3] = origin[1];
            T[2][3] = origin[2];

            IPObjectStruct* Placed = IritGeomTransformObject(Extr, T);
            IritPrsrFreeObject(Extr);

            if (Placed != NULL) {
                Placed->Pnext = NULL;
                IritPrsrListObjectAppend(ResList, Placed);
            }
        } /* per polygon */
    } /* per list element */

    if (ResList->U.Lst.PObjList == NULL || ResList->U.Lst.ListMaxLen == 0) {
        IritPrsrFreeObject(ResList);
        return NULL;
    }

    return ResList;
}

