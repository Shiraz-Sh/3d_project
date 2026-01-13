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



#include "inc_irit/cagd_lib.h"

IPObjectStruct* IritPrsrExtrude2DPointsToRuledSrf(const IrtE2PtStruct* pts,
    int n,
    IrtVecType Dir)
{
    int i;
    CagdCrvStruct* Crv;
    CagdSrfStruct* Srf;
    CagdVecStruct CagdDir;
    IPObjectStruct* SrfObj;

    if (pts == NULL || n < 2) return NULL;

    /* 1. Create a Linear B-spline Curve (Order 2) from your points. */
    Crv = IritCagdBspCrvNew(n, 2, CAGD_PT_E2_TYPE);

    for (i = 0; i < n; i++) {
        Crv->Points[1][i] = pts[i].Pt[0];
        Crv->Points[2][i] = pts[i].Pt[1];
    }

    /* Use IritCagdBspKnotUniformOpen instead of IritCagdBspKnotUniform. */
    IritCagdBspKnotUniformOpen(n, 2, Crv->KnotVector);

    /* 2. Prepare the extrusion direction. */
    for (i = 0; i < 3; i++) CagdDir.Vec[i] = Dir[i];

    /* 3. Extrude the curve to create a mathematical surface. */
    Srf = IritCagdExtrudeSrf(Crv, &CagdDir);
    IritCagdCrvFree(Crv);

    /* 4. Wrap into an IRIT Surface Object. */
    SrfObj = IritPrsrGenSrfObject("ExtrusionSrf", Srf, NULL);

    return SrfObj;
}
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
    IPObjectStruct* Head = NULL, * Tail = NULL;
    int resolvedGrid = (GridSize > 0) ? GridSize : 20;

    if (PObj == NULL || ViewMats == NULL || NumViews <= 0)
        return NULL;

    if (!IP_IS_POLY_OBJ((IPObjectStruct*)PObj))
        return NULL;

    /* Preprocess once if requested. */
    if (UsePreprocess) {
        /* Copy the polygonal object and ensure regularization/adjacencies. */
        WorkCopy = IritPrsrCopyObject(NULL, (IPObjectStruct*)(void*)PObj, TRUE);
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
            WorkCopy = IritPrsrCopyObject(NULL, (IPObjectStruct*)(void*)PObj, TRUE);
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
 * *****************
 * Change signature: add 'IsPre' flag to let the function name projected polygons accordingly.
 * Replace previous implementation that didn't set ObjName so extrusion step can map each projected
 * polygon to its originating view.
 */
 IPObjectStruct* IritPrsrProjectSilhouetteToViewPolys(const IPObjectStruct* SilList,
         const IrtHmgnMatType* ViewMats,
         int NumViews,
         int IsPre)
 {
     if (SilList == NULL || ViewMats == NULL || NumViews <= 0)
         return NULL;

     IPObjectStruct* ResList = IritPrsrGenLISTObject(NULL);
     if (ResList == NULL)
         return NULL;

     for (int vi = 0; vi < NumViews; ++vi) {
         IPObjectStruct* SilObj = IritPrsrListObjectGet((IPObjectStruct*)(void*)SilList, vi);
         if (SilObj == NULL || !IP_IS_POLY_OBJ(SilObj))
             continue;

         /* Extract forward vector (w) from ViewMat (third column) and normalize. */
         IrtVecType w;
         w[0] = ViewMats[vi][0][2];
         w[1] = ViewMats[vi][1][2];
         w[2] = ViewMats[vi][2][2];
         {
             IrtRType wlen = sqrt(IRIT_SQR(w[0]) + IRIT_SQR(w[1]) + IRIT_SQR(w[2]));
             if (wlen <= IRIT_EPS)
                 continue;
             w[0] /= wlen; w[1] /= wlen; w[2] /= wlen;
         }

         /* Build stable orthonormal basis u,v from w. */
         IrtVecType arb, u, v;
         if (fabs(w[0]) < 0.9) { arb[0] = 1.0; arb[1] = 0.0; arb[2] = 0.0; }
         else { arb[0] = 0.0; arb[1] = 1.0; arb[2] = 0.0; }

         u[0] = arb[1] * w[2] - arb[2] * w[1];
         u[1] = arb[2] * w[0] - arb[0] * w[2];
         u[2] = arb[0] * w[1] - arb[1] * w[0];
         {
             IrtRType ul = sqrt(IRIT_SQR(u[0]) + IRIT_SQR(u[1]) + IRIT_SQR(u[2]));
             if (ul <= IRIT_EPS) continue;
             u[0] /= ul; u[1] /= ul; u[2] /= ul;
         }

         v[0] = w[1] * u[2] - w[2] * u[1];
         v[1] = w[2] * u[0] - w[0] * u[2];
         v[2] = w[0] * u[1] - w[1] * u[0];

         /* For each silhouette polygon, create a projected closed polygon lying on the view plane.
            Implementation uses 4x4 local transform (world -> [u,v,w] local coords with origin)
            and its inverse. We project by setting local z = 0, then transform back. */
         int polyIdx = 0;
         for (IPPolygonStruct* Pl = SilObj->U.Pl; Pl != NULL; Pl = Pl->Pnext, ++polyIdx) {
             IPVertexStruct* V = Pl->PVertex;
             if (V == NULL)
                 continue;

             /* Count vertices. */
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

             /* Collect coords and centroid (origin). */
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

             /* Build world->local matrix MatLocal:
                MatLocal * P_world = [ dot(u,P)-dot(u,origin),
                                       dot(v,P)-dot(v,origin),
                                       dot(w,P)-dot(w,origin) ]
                last row stays unit for homogenous point transform.
             */
             IrtHmgnMatType MatLocal, InvMat;
             IritMiscMatGenUnitMat(MatLocal);

             MatLocal[0][0] = u[0]; MatLocal[0][1] = u[1]; MatLocal[0][2] = u[2];
             MatLocal[0][3] = -(u[0] * origin[0] + u[1] * origin[1] + u[2] * origin[2]);

             MatLocal[1][0] = v[0]; MatLocal[1][1] = v[1]; MatLocal[1][2] = v[2];
             MatLocal[1][3] = -(v[0] * origin[0] + v[1] * origin[1] + v[2] * origin[2]);

             MatLocal[2][0] = w[0]; MatLocal[2][1] = w[1]; MatLocal[2][2] = w[2];
             MatLocal[2][3] = -(w[0] * origin[0] + w[1] * origin[1] + w[2] * origin[2]);

             MatLocal[3][0] = 0.0;  MatLocal[3][1] = 0.0;  MatLocal[3][2] = 0.0; MatLocal[3][3] = 1.0;

             if (!IritMiscMatInverseMatrix(MatLocal, InvMat)) {
                 IritFree(coords);
                 continue;
             }

             /* Create new polygon and vertex list from projected points (z=0 in local coords). */
             IPPolygonStruct* NewPl = IritPrsrAllocPolygon(0, NULL, NULL);
             IPVertexStruct* FirstV = NULL;
             IPVertexStruct* PrevV = NULL;

             for (int j = 0; j < n; ++j) {
                 IrtPtType localP;
                 IritMiscMatMultPtby4by4(localP, coords[j], MatLocal);

                 /* Project into view plane by zeroing local z. */
                 localP[2] = 0.0;

                 /* Map back to world coordinates. */
                 IrtPtType worldProj;
                 IritMiscMatMultPtby4by4(worldProj, localP, InvMat);

                 /* Allocate and link vertex. */
                 IPVertexStruct* NV = IritPrsrAllocVertex2(NULL);
                 NV->Coord[0] = worldProj[0];
                 NV->Coord[1] = worldProj[1];
                 NV->Coord[2] = worldProj[2];

                 if (FirstV == NULL) {
                     FirstV = PrevV = NV;
                 }
                 else {
                     PrevV->Pnext = NV;
                     PrevV = NV;
                 }
             }

             /* Close loop and attach polygon. */
             PrevV->Pnext = FirstV;
             NewPl->PVertex = FirstV;

             IritPrsrUpdatePolyPlane(NewPl);

             IPObjectStruct* PolyObj = IritPrsrGenPOLYObject(NewPl);
             if (PolyObj != NULL) {
                 IP_SET_POLYGON_OBJ(PolyObj);
                 PolyObj->Pnext = NULL;

                 /* Set name encoding view and polygon index so extrusion step knows the view. */
                 char namebuf[256];
                 if (IsPre)
                     snprintf(namebuf, sizeof(namebuf), "proj_pre_view%d_poly%03d", vi, polyIdx);
                 else
                     snprintf(namebuf, sizeof(namebuf), "proj_dir_view%d_poly%03d", vi, polyIdx);

                 PolyObj->ObjName = IritMiscStrdup(namebuf);

                 IritPrsrListObjectAppend(ResList, PolyObj);
             }

             IritFree(coords);
         } /* per polygon */
     } /* per view */

     if (ResList->U.Lst.PObjList == NULL || ResList->U.Lst.ListMaxLen == 0) {
         IritPrsrFreeObject(ResList);
         return NULL;
     }

     return ResList;
 }



/*
 * Extrude projected silhouette polygons (precomputed and direct lists)
 * toward each view's forward direction and save per-polygon OBJ files.
 *
 * Parameters:
 *   ProjListPre   - LIST of 2ds obtained with preprocessing (may be NULL)
 *   ProjListDirect - LIST of 2ds obtained by direct extraction (may be NULL)
 *   Views        - array of view matrices (NumViews entries)
 *   NumViews     - number of view matrices / silhouette list elements
 *   Depth        - extrusion length along the view forward vector
 *
 * Behavior mirrors the logic previously embedded in test.c: for every view,
 * for every polygon in the 2d object, build the polygon's centroid,
 * compute local (u,v) coordinates, call IritPrsrExtrude2DPointsToSolidDir,
 * translate the resulting solid to centroid. return a LIST of the resulting
 * placed solids.
 */
IPObjectStruct* IritPrsrExtrudeSilhouetteListsToViewSolids(IPObjectStruct* ProjListPre,
    IPObjectStruct* ProjListDirect,
    const IrtHmgnMatType* Views,
    int NumViews,
    IrtRType Depth)
{
    if ((ProjListPre == NULL && ProjListDirect == NULL) ||
        Views == NULL || NumViews <= 0 || Depth == 0.0)
        return NULL;

    IPObjectStruct* ResList = IritPrsrGenLISTObject(NULL);
    if (ResList == NULL)
        return NULL;

    /* Helper to process a flattened projected list (elements named "proj_{pre|dir}_view%d_poly%03d"). */
    for (int usePre = 0; usePre < 2; ++usePre) {
        IPObjectStruct* ProjList = usePre ? ProjListPre : ProjListDirect;
        if (ProjList == NULL)
            continue;

        int idx = 0;
        IPObjectStruct* PObj = NULL;
        while ((PObj = IritPrsrListObjectGet(ProjList, idx++)) != NULL) {
            if (!IP_IS_POLY_OBJ(PObj))
                continue;

            /* Determine source view index from ObjName */
            int vi = -1, polyIdx = -1;
            if (PObj->ObjName != NULL) {
                if (usePre)
                    sscanf(PObj->ObjName, "proj_pre_view%d_poly%03d", &vi, &polyIdx);
                else
                    sscanf(PObj->ObjName, "proj_dir_view%d_poly%03d", &vi, &polyIdx);
            }
            /* Fallback: if parsing failed, skip this polygon. */
            if (vi < 0 || vi >= NumViews)
                continue;

            /* Extract view forward vector (third column) and normalize. */
            IrtVecType f;
            f[0] = Views[vi][0][2];
            f[1] = Views[vi][1][2];
            f[2] = Views[vi][2][2];
            {
                IrtRType flen = sqrt(IRIT_SQR(f[0]) + IRIT_SQR(f[1]) + IRIT_SQR(f[2]));
                if (flen <= IRIT_EPS) {
                    f[0] = 0.0; f[1] = 0.0; f[2] = 1.0;
                    flen = 1.0;
                }
                f[0] /= flen; f[1] /= flen; f[2] /= flen;
            }

            /* Build local basis u,v (stable). */
            IrtVecType arb, u, vtmp;
            if (fabs(f[0]) < 0.9) { arb[0] = 1.0; arb[1] = 0.0; arb[2] = 0.0; }
            else { arb[0] = 0.0; arb[1] = 1.0; arb[2] = 0.0; }

            u[0] = arb[1] * f[2] - arb[2] * f[1];
            u[1] = arb[2] * f[0] - arb[0] * f[2];
            u[2] = arb[0] * f[1] - arb[1] * f[0];
            {
                IrtRType ul = sqrt(IRIT_SQR(u[0]) + IRIT_SQR(u[1]) + IRIT_SQR(u[2]));
                if (ul <= IRIT_EPS) continue;
                u[0] /= ul; u[1] /= ul; u[2] /= ul;
            }

            vtmp[0] = f[1] * u[2] - f[2] * u[1];
            vtmp[1] = f[2] * u[0] - f[0] * u[2];
            vtmp[2] = f[0] * u[1] - f[1] * u[0];

            /* For each polygon component inside this projected object (usually one) */
            for (IPPolygonStruct* Pl = PObj->U.Pl; Pl != NULL; Pl = Pl->Pnext) {
                IPVertexStruct* V = Pl->PVertex;
                if (V == NULL)
                    continue;

                /* Count vertices and collect coordinates. */
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
                int idxV = 0;
                guard = 0;
                IrtRType cx = 0.0, cy = 0.0, cz = 0.0;
                while (cur != NULL && idxV < n) {
                    coords[idxV][0] = cur->Coord[0];
                    coords[idxV][1] = cur->Coord[1];
                    coords[idxV][2] = cur->Coord[2];
                    cx += coords[idxV][0];
                    cy += coords[idxV][1];
                    cz += coords[idxV][2];
                    ++idxV;
                    cur = cur->Pnext;
                    if (++guard > 200000) break;
                    if (cur == V) break;
                }
                if (idxV != n) { IritFree(coords); continue; }

                IrtPtType origin;
                origin[0] = cx / n; origin[1] = cy / n; origin[2] = cz / n;

                /* Build 2D points relative to centroid in the view's (u,v) basis. */
                IrtE2PtStruct* pts2d = (IrtE2PtStruct*)IritMalloc(sizeof(IrtE2PtStruct) * n);
                if (pts2d == NULL) { IritFree(coords); continue; }

                for (int j = 0; j < n; ++j) {
                    IrtRType dx = coords[j][0] - origin[0];
                    IrtRType dy = coords[j][1] - origin[1];
                    IrtRType dz = coords[j][2] - origin[2];

                    pts2d[j].Pt[0] = dx * u[0] + dy * u[1] + dz * u[2];
                    pts2d[j].Pt[1] = dx * vtmp[0] + dy * vtmp[1] + dz * vtmp[2];
                }

                /* Extrude in view forward direction scaled by Depth. */
                IrtVecType Dir;
                Dir[0] = (IrtRType)(f[0] * Depth);
                Dir[1] = (IrtRType)(f[1] * Depth);
                Dir[2] = (IrtRType)(f[2] * Depth);

                IPObjectStruct* Extr = IritPrsrExtrude2DPointsToSolidDir(pts2d, n, Dir);

                IritFree(pts2d);
                IritFree(coords);

                if (Extr == NULL)
                    continue;

                /* Translate extruded solid back to polygon centroid (world). */
                IrtHmgnMatType T;
                IritMiscMatGenUnitMat(T);
                T[0][3] = origin[0];
                T[1][3] = origin[1];
                T[2][3] = origin[2];

                IPObjectStruct* Placed = IritGeomTransformObject(Extr, T);
                IritPrsrFreeObject(Extr);

                if (Placed != NULL) {
                    /* keep a meaningful ObjName for saving later */
                    char namebuf[512];
                    if (usePre)
                        snprintf(namebuf, sizeof(namebuf), "sil_pre_view%d_poly%03d_extr", vi, polyIdx);
                    else
                        snprintf(namebuf, sizeof(namebuf), "sil_dir_view%d_poly%03d_extr", vi, polyIdx);

                    Placed->ObjName = IritMiscStrdup(namebuf);
                    Placed->Pnext = NULL;
                    IritPrsrListObjectAppend(ResList, Placed);
                }
            } /* per polygon inside PObj */
        } /* while list elements */
    } /* for pre/direct */

    if (ResList->U.Lst.PObjList == NULL || ResList->U.Lst.ListMaxLen == 0) {
        IritPrsrFreeObject(ResList);
        return NULL;
    }

    return ResList;
}


/* * Helper: Manually generate a LookAt matrix.
 * This replaces the hypothetical 'IritGeomMatGenLookAt'.
 * * Logic:
 * 1. F = Normalize(Center - Eye)
 * 2. S = Normalize(Cross(F, Up))  -> Side/Right vector
 * 3. U = Cross(S, F)              -> Recomputed Up vector orthogonal to F
 * * The resulting matrix maps the world so the camera is at origin,
 * looking down the positive Z axis (IRIT standard).
 */
void GenLookAtMatrix(IrtVecType Eye, IrtVecType Center, IrtVecType Up, IrtHmgnMatType Mat)
{
    IrtVecType F, S, U;

    /* 1. Calculate Forward Vector (F) */
    IRIT_VEC_SUB(F, Center, Eye);
    IRIT_VEC_NORMALIZE(F);

    /* 2. Calculate Side Vector (S) = F x Up */
    IRIT_CROSS_PROD(S, F, Up);
    /* Check for degenerate case where F is parallel to Up */
    if ((IRIT_DOT_PROD(S, S)) < 1e-6) {
        /* Fallback: If looking straight up/down, choose X-axis as side */
        S[0] = 1.0; S[1] = 0.0; S[2] = 0.0;
    }
    IRIT_VEC_NORMALIZE(S);

    /* 3. Calculate True Up Vector (U) = S x F */
    IRIT_CROSS_PROD(U, S, F);
    IRIT_VEC_NORMALIZE(U);

    /* 4. Construct Matrix
       Row 0: S
       Row 1: U
       Row 2: F (or -F depending on coordinate system, IRIT usually +Z view)
       Row 3: Translation
    */
    IritMiscMatGenUnitMat(Mat);

    /* Rotation Part */
    Mat[0][0] = S[0];  Mat[0][1] = S[1];  Mat[0][2] = S[2];
    Mat[1][0] = U[0];  Mat[1][1] = U[1];  Mat[1][2] = U[2];
    /* Note: IRIT often uses +Z as view direction. If standard OpenGL, this row is usually -F */
    Mat[2][0] = F[0];  Mat[2][1] = F[1];  Mat[2][2] = F[2];

    /* Translation Part (-Eye dot Basis) */
    Mat[0][3] = -IRIT_DOT_PROD(S, Eye);
    Mat[1][3] = -IRIT_DOT_PROD(U, Eye);
    Mat[2][3] = -IRIT_DOT_PROD(F, Eye);
}


/* * Helper: Generate a View Matrix from spherical coordinates.
 * Implements the camera position logic from Eq (3) in the paper.
 * r is fixed to 2.0 as per the paper[cite: 171].
 */
void GenViewMatrix(IrtRType phi, IrtRType theta, IrtHmgnMatType Mat)
{
    IrtVecType CamPos, Center, Up;
    IrtRType r = 2.0;

    /* Eq (3): Convert spherical (phi, theta) to Cartesian (x, y, z) */
    CamPos[0] = r * cos(phi) * cos(theta);
    CamPos[1] = r * cos(phi) * sin(theta);
    CamPos[2] = r * sin(phi);

    Center[0] = 0.0; Center[1] = 0.0; Center[2] = 0.0;

    /* Simple Up vector logic (handle singularity at poles) */
    if (fabs(CamPos[0]) < 0.1 && fabs(CamPos[1]) < 0.1) {
        Up[0] = 1.0; Up[1] = 0.0; Up[2] = 0.0;
    }
    else {
        Up[0] = 0.0; Up[1] = 0.0; Up[2] = 1.0;
    }

    /* Use IRIT or standard LookAt logic to create the matrix */
    GenLookAtMatrix(CamPos, Center, Up, Mat);
}

/*
 * Helper: Calculate projected 2D area of a polygon chain.
 * Used as a heuristic: Maximize projected area ~ Maximize feature definition.
 */
IrtRType CalcPolygonArea(IPObjectStruct* PObj)
{
    IrtRType area = 0.0;
    if (PObj == NULL || !IP_IS_POLY_OBJ(PObj)) return 0.0;

    for (IPPolygonStruct* Pl = PObj->U.Pl; Pl != NULL; Pl = Pl->Pnext) {
        /* Standard Green's theorem for polygon area in 2D (XY plane) */
        /* Note: The silhouette extractor usually returns 3D polylines.
           We assume they are projected or we approximate via their major 2D plane.
           For robustness, we can project to the view plane, but for a simple
           heuristic, we summing lengths or bounding box diagonals is also fast.
           Here we implement a simplified 3D polygon area estimation. */

        area += IritGeomPolyOnePolyArea(Pl, TRUE); // Uses IRIT's built-in area function
    }
    return area;
}


/*
 * Main Function: Sampling-Based View Selection
 * * 1. Generates 'NumSamples' viewpoints using a Fibonacci Lattice (uniform sphere).
 * and return this.
 * 
 * 
 * 2. Calls `IritPrsrGetSilhouettesForViews` to batch process them.
 * 3. Evaluates the returned silhouettes and picks the best one.
 * 4. Returns the 3 of the best view and the matrix itself.
 */
IrtHmgnMatType* SelectBestViewSampling(IPObjectStruct* PObj,
    int NumSamples,
    IrtHmgnMatType* ResultMat)
{
    int i, bestIdx = -1;
    IrtRType maxScore = -1.0;
    IrtHmgnMatType* ViewMats = NULL;

    /* 1. Generate Candidate Matrices */
    /* The paper uses Fibonacci grid sampling  */
    ViewMats = (IrtHmgnMatType*)IritMalloc(sizeof(IrtHmgnMatType) * NumSamples);

    IrtRType phi = (sqrt(5.0) - 1.0) / 2.0; /* Golden ratio */

    for (i = 0; i < NumSamples; ++i) {
        IrtRType y = 1 - (i / (IrtRType)(NumSamples - 1)) * 2; /* y goes from 1 to -1 */
        IrtRType radius = sqrt(1 - y * y);

        IrtRType theta = 2 * M_PI * phi * i;

        IrtRType x = cos(theta) * radius;
        IrtRType z = sin(theta) * radius;

        /* Convert this unit vector (x,y,z) on sphere to spherical angles for our helper */
        /* Or directly use LookAt with this vector scaled by R=2.0 */
        IrtVecType CamPos, Center, Up;
        CamPos[0] = x * 2.0; CamPos[1] = y * 2.0; CamPos[2] = z * 2.0;
        Center[0] = 0.0; Center[1] = 0.0; Center[2] = 0.0;
        Up[0] = 0.0; Up[1] = 0.0; Up[2] = 1.0;

        /* Handle pole singularity for Up vector */
        if (fabs(CamPos[0]) < 0.01 && fabs(CamPos[1]) < 0.01) {
            Up[0] = 1.0; Up[1] = 0.0; Up[2] = 0.0;
        }

        GenLookAtMatrix(CamPos, Center, Up, ViewMats[i]);
    }

    /* 2. Batch Process Silhouettes */
    /* We use UsePreprocess = 1 (TRUE) as recommended in code comments for many views */
    printf("Evaluating %d views...\n", NumSamples);
    return ViewMats;
    //IPObjectStruct* SilList = IritPrsrGetSilhouettesForViews(PObj, ViewMats, NumSamples, 20, TRUE);

    //if (SilList == NULL || !IP_IS_OLST_OBJ(SilList)) {
    //    IritFree(ViewMats);
    //    return -1;
    //}

    ///* 3. Evaluate Results */
    ///* SilList is a list of objects, one per view */
    //int viewIdx = 0;
    //IPObjectStruct* SilObj;

    ///* Iterate over the IRIT List Object */
    //for (viewIdx = 0; viewIdx < NumSamples; ++viewIdx) {
    //    SilObj = IritPrsrListObjectGet(SilList, viewIdx);

    //    if (SilObj != NULL) {
    //        /* Heuristic: Maximize projected area */
    //        IrtRType score = CalcPolygonArea(SilObj);

    //        if (score > maxScore) {
    //            maxScore = score;
    //            bestIdx = viewIdx;
    //        }
    //    }
    //}

    ///* 4. Output Result */
    //if (bestIdx >= 0) {
    //    ResultMat = ViewMats[bestIdx]; // Copy struct
    //    printf("Selected View %d with Score: %f\n", bestIdx, maxScore);
    //}

    ///* Cleanup */
    //IritPrsrFreeObject(SilList);
    //IritFree(ViewMats);

    //return bestIdx;
}
