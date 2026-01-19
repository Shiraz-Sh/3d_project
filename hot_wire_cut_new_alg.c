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


/* --------------------------------------------------------------------------
 * Helper: build a stable orthonormal basis (u,v,w) given a view matrix.
 *
 * Inputs:
 *   Mat         - 4x4 view matrix; the "forward" vector is taken from Mat[*][2].
 *   AllowDefault - if nonzero, on degenerate forward vector the helper will
 *                  substitute w=(0,0,1) and continue (fallback). Otherwise
 *                  function returns 0 for failure.
 *
 * Outputs:
 *   u, v, w     - computed orthonormal basis (all three are written).
 *
 * Return:
 *   1 on success, 0 on failure (and u/v/w contents are undefined unless
 *   AllowDefault was used and returned 1).
 *
 * Rationale:
 *   Both `IritPrsrProjectSilhouetteToViewPolys` and
 *   `IritPrsrExtrudeSilhouetteListsToViewSolids` contain identical logic to
 *   extract the view-forward vector and build a stable (u,v) basis. Centralize
 *   that logic here to avoid duplication and to make maintenance simpler.
 * ------------------------------------------------------------------------ */
static int IritPrsrBuildViewBasisFromMat(const IrtHmgnMatType Mat,
    IrtVecType u,
    IrtVecType v,
    IrtVecType w,
    int AllowDefault)
{
    IrtVecType arb;
    IrtRType wlen, ul;

    /* Extract forward vector from third *row* of matrix (Mat[2][0..2]).
       The view/LookAt code in this file stores rotation in rows, so the
       forward vector is on row index 2. */
    w[0] = Mat[2][0];
    w[1] = Mat[2][1];
    w[2] = Mat[2][2];

    wlen = sqrt(IRIT_SQR(w[0]) + IRIT_SQR(w[1]) + IRIT_SQR(w[2]));
    if (wlen <= IRIT_EPS) {
        if (!AllowDefault)
            return 0;
        /* fallback to canonical forward if allowed */
        w[0] = 0.0; w[1] = 0.0; w[2] = 1.0;
    }
    else {
        w[0] /= wlen; w[1] /= wlen; w[2] /= wlen;
    }

    /* Choose stable arbitrary vector not parallel to w */
    if (fabs(w[0]) < 0.9) { arb[0] = 1.0; arb[1] = 0.0; arb[2] = 0.0; }
    else { arb[0] = 0.0; arb[1] = 1.0; arb[2] = 0.0; }

    /* u = normalize(cross(arb, w)) */
    u[0] = arb[1] * w[2] - arb[2] * w[1];
    u[1] = arb[2] * w[0] - arb[0] * w[2];
    u[2] = arb[0] * w[1] - arb[1] * w[0];

    ul = sqrt(IRIT_SQR(u[0]) + IRIT_SQR(u[1]) + IRIT_SQR(u[2]));
    if (ul <= IRIT_EPS) {
        /* extremely unlikely; if allowed try a different arb, else fail. */
        if (!AllowDefault)
            return 0;
        /* try fallback arb */
        arb[0] = 0.0; arb[1] = 1.0; arb[2] = 0.0;
        u[0] = arb[1] * w[2] - arb[2] * w[1];
        u[1] = arb[2] * w[0] - arb[0] * w[2];
        u[2] = arb[0] * w[1] - arb[1] * w[0];
        ul = sqrt(IRIT_SQR(u[0]) + IRIT_SQR(u[1]) + IRIT_SQR(u[2]));
        if (ul <= IRIT_EPS)
            return 0;
    }

    u[0] /= ul; u[1] /= ul; u[2] /= ul;

    /* v = cross(w, u) */
    v[0] = w[1] * u[2] - w[2] * u[1];
    v[1] = w[2] * u[0] - w[0] * u[2];
    v[2] = w[0] * u[1] - w[1] * u[0];

    return 1;
}



/* Helper: collect polygon vertices into an array, remove duplicate closing
 * vertex if it equals the first within IRIT_EPS, and compute centroid.
 *
 * Inputs:
 *   V      - pointer to a polygon vertex circular list (may be NULL).
 *   outN   - pointer to receive number of vertices returned (0 on failure/degenerate).
 *   origin - optional pointer to receive centroid (may be NULL).
 *
 * Returns:
 *   Allocated array of `IrtPtType` with `*outN` entries on success. Caller must IritFree()
 *   the returned array. Returns NULL and sets *outN = 0 on failure or degenerate polygon.
 */
static IrtPtType* CollectPolygonCoordsAndCentroid(const IPVertexStruct* V,
    int* outN,
    IrtPtType origin)
{
    int n = 0;
    *outN = 0;
    if (V == NULL)
        return NULL;

    /* Count vertices (guard against corrupted lists). */
    const IPVertexStruct* cur = V;
    int guard = 0;
    do {
        ++n;
        cur = cur->Pnext;
        if (++guard > 200000) { n = 0; break; }
    } while (cur != NULL && cur != V);

    if (n <= 0)
        return NULL;

    /* Collect coordinates. */
    IrtPtType* coords = (IrtPtType*)IritMalloc(sizeof(IrtPtType) * n);
    if (coords == NULL)
        return NULL;

    cur = V;
    guard = 0;
    int idx = 0;
    while (cur != NULL && idx < n) {
        coords[idx][0] = cur->Coord[0];
        coords[idx][1] = cur->Coord[1];
        coords[idx][2] = cur->Coord[2];
        ++idx;
        cur = cur->Pnext;
        if (++guard > 200000) break;
        if (cur == V) break;
    }
    if (idx != n) {
        IritFree(coords);
        return NULL;
    }

    /* Sanitize: if last equals first (within IRIT_EPS) drop last. */
    if (n > 1) {
        IrtRType dx = coords[n - 1][0] - coords[0][0];
        IrtRType dy = coords[n - 1][1] - coords[0][1];
        IrtRType dz = coords[n - 1][2] - coords[0][2];
        if ((dx * dx + dy * dy + dz * dz) <= IRIT_SQR(IRIT_EPS)) {
            --n;
        }
    }

    /* Need at least 3 distinct vertices to form polygon. */
    if (n < 3) {
        IritFree(coords);
        return NULL;
    }

    /* Compute centroid of the remaining vertices. */
    if (origin != NULL) {
        IrtRType cx = 0.0, cy = 0.0, cz = 0.0;
        for (int i = 0; i < n; ++i) {
            cx += coords[i][0];
            cy += coords[i][1];
            cz += coords[i][2];
        }
        origin[0] = cx / n;
        origin[1] = cy / n;
        origin[2] = cz / n;
    }

    *outN = n;
    return coords;
}


/* ----------------- existing functions (unchanged except for basis usage)  */

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
//IPObjectStruct* IritPrsrProjectSilhouetteToViewPolys(const IPObjectStruct* SilList,
//    const IrtHmgnMatType* ViewMats,
//    int NumViews,
//    int IsPre)
//{
//    if (SilList == NULL || ViewMats == NULL || NumViews <= 0)
//        return NULL;
//
//    IPObjectStruct* ResList = IritPrsrGenLISTObject(NULL);
//    if (ResList == NULL)
//        return NULL;
//
//    for (int vi = 0; vi < NumViews; ++vi) {
//        IPObjectStruct* SilObj = IritPrsrListObjectGet((IPObjectStruct*)(void*)SilList, vi);
//        if (SilObj == NULL || !IP_IS_POLY_OBJ(SilObj))
//            continue;
//
//        /* Build stable orthonormal basis u,v,w from view matrix (no fallback). */
//        IrtVecType u, v, w;
//        if (!IritPrsrBuildViewBasisFromMat(ViewMats[vi], u, v, w, 0))
//            continue;
//
//        /* For each silhouette polygon, create a projected closed polygon lying on the view plane.
//           Implementation uses 4x4 local transform (world -> [u,v,w] local coords with origin)
//           and its inverse. We project by setting local z = 0, then transform back. */
//        int polyIdx = 0;
//        for (IPPolygonStruct* Pl = SilObj->U.Pl; Pl != NULL; Pl = Pl->Pnext, ++polyIdx) {
//            IPVertexStruct* V = Pl->PVertex;
//            if (V == NULL)
//                continue;
//
//            /* Count vertices. */
//            int n = 0;
//            IPVertexStruct* cur = V;
//            int guard = 0;
//            while (cur != NULL) {
//                ++n;
//                cur = cur->Pnext;
//                if (++guard > 200000) { n = 0; break; }
//                if (cur == V) break;
//            }
//            if (n < 3)
//                continue;
//
//            /* Collect coords and centroid (origin). */
//            IrtPtType* coords = (IrtPtType*)IritMalloc(sizeof(IrtPtType) * n);
//            if (coords == NULL)
//                continue;
//
//            cur = V;
//            int idx = 0;
//            guard = 0;
//            IrtRType cx = 0.0, cy = 0.0, cz = 0.0;
//            while (cur != NULL && idx < n) {
//                coords[idx][0] = cur->Coord[0];
//                coords[idx][1] = cur->Coord[1];
//                coords[idx][2] = cur->Coord[2];
//                cx += coords[idx][0];
//                cy += coords[idx][1];
//                cz += coords[idx][2];
//                ++idx;
//                cur = cur->Pnext;
//                if (++guard > 200000) break;
//                if (cur == V) break;
//            }
//            if (idx != n) { IritFree(coords); continue; }
//
//            IrtPtType origin;
//            origin[0] = cx / n; origin[1] = cy / n; origin[2] = cz / n;
//
//            /* Build world->local matrix MatLocal:
//               MatLocal * P_world = [ dot(u,P)-dot(u,origin),
//                                      dot(v,P)-dot(v,origin),
//                                      dot(w,P)-dot(w,origin) ]
//               last row stays unit for homogenous point transform.
//            */
//            IrtHmgnMatType MatLocal, InvMat;
//            IritMiscMatGenUnitMat(MatLocal);
//
//            MatLocal[0][0] = u[0]; MatLocal[0][1] = u[1]; MatLocal[0][2] = u[2];
//            MatLocal[0][3] = -(u[0] * origin[0] + u[1] * origin[1] + u[2] * origin[2]);
//
//            MatLocal[1][0] = v[0]; MatLocal[1][1] = v[1]; MatLocal[1][2] = v[2];
//            MatLocal[1][3] = -(v[0] * origin[0] + v[1] * origin[1] + v[2] * origin[2]);
//
//            MatLocal[2][0] = w[0]; MatLocal[2][1] = w[1]; MatLocal[2][2] = w[2];
//            MatLocal[2][3] = -(w[0] * origin[0] + w[1] * origin[1] + w[2] * origin[2]);
//
//            MatLocal[3][0] = 0.0;  MatLocal[3][1] = 0.0;  MatLocal[3][2] = 0.0; MatLocal[3][3] = 1.0;
//
//            if (!IritMiscMatInverseMatrix(MatLocal, InvMat)) {
//                IritFree(coords);
//                continue;
//            }
//
//            /* Create new polygon and vertex list from projected points (z=0 in local coords). */
//            IPPolygonStruct* NewPl = IritPrsrAllocPolygon(0, NULL, NULL);
//            IPVertexStruct* FirstV = NULL;
//            IPVertexStruct* PrevV = NULL;
//
//            for (int j = 0; j < n; ++j) {
//                IrtPtType localP;
//                IritMiscMatMultPtby4by4(localP, coords[j], MatLocal);
//
//                /* Project into view plane by zeroing local z. */
//                localP[2] = 0.0;
//
//                /* Map back to world coordinates. */
//                IrtPtType worldProj;
//                IritMiscMatMultPtby4by4(worldProj, localP, InvMat);
//
//                /* Allocate and link vertex. */
//                IPVertexStruct* NV = IritPrsrAllocVertex2(NULL);
//                NV->Coord[0] = worldProj[0];
//                NV->Coord[1] = worldProj[1];
//                NV->Coord[2] = worldProj[2];
//
//                if (FirstV == NULL) {
//                    FirstV = PrevV = NV;
//                }
//                else {
//                    PrevV->Pnext = NV;
//                    PrevV = NV;
//                }
//            }
//
//            /* Close loop and attach polygon. */
//            PrevV->Pnext = FirstV;
//            NewPl->PVertex = FirstV;
//
//            IritPrsrUpdatePolyPlane(NewPl);
//
//            IPObjectStruct* PolyObj = IritPrsrGenPOLYObject(NewPl);
//            if (PolyObj != NULL) {
//                IP_SET_POLYGON_OBJ(PolyObj);
//                PolyObj->Pnext = NULL;
//
//                /* Set name encoding view and polygon index so extrusion step knows the view. */
//                char namebuf[256];
//                if (IsPre)
//                    snprintf(namebuf, sizeof(namebuf), "proj_pre_view%d_poly%03d", vi, polyIdx);
//                else
//                    snprintf(namebuf, sizeof(namebuf), "proj_dir_view%d_poly%03d", vi, polyIdx);
//
//                PolyObj->ObjName = IritMiscStrdup(namebuf);
//
//                IritPrsrListObjectAppend(ResList, PolyObj);
//            }
//
//            IritFree(coords);
//        } /* per polygon */
//    } /* per view */
//
//    if (ResList->U.Lst.PObjList == NULL || ResList->U.Lst.ListMaxLen == 0) {
//        IritPrsrFreeObject(ResList);
//        return NULL;
//    }
//
//    return ResList;
//}




/*
 * PerformBooleanUnion - use IRIT native boolean (IritBooleanOR) to union
 * a LIST of polygon objects. Takes ownership of `ListOfPolys`.
 *
 * Behavior:
 *  - If Clist is not a LIST or contains <= 1 element, returns the input unchanged.
 *  - Otherwise, performs sequential union: accum = OR(accum, next).
 *  - On success returns a new LIST object whose elements are the resulting
 *    (possibly multiple) polygon objects produced by the union.
 *  - On failure returns the original `ListOfPolys` (no objects freed).
 *
 * Note: This implementation relies on existing IRIT boolean API found in
 * bool_lib (IritBooleanOR / IritBoolean2D). It uses the IRIT helper
 * routines already present in the code base.
 */
static IPObjectStruct* PerformBooleanUnion(IPObjectStruct* ListOfPolys)
{
    int i, cnt = 0;

    if (ListOfPolys == NULL)
        return NULL;

    /* Only operate on LIST objects. If not a list, nothing to do. */
    if (ListOfPolys->ObjType != IP_OBJ_LIST_OBJ)
        return ListOfPolys;

    /* Count elements in the list. */
    while (IritPrsrListObjectGet(ListOfPolys, cnt) != NULL)
        ++cnt;

    if (cnt <= 1) /* nothing to union */
        return ListOfPolys;

    /* Initialize accumulator with a copy of the first polygon object. */
    IPObjectStruct* first = IritPrsrListObjectGet(ListOfPolys, 0);
    if (first == NULL)
        return ListOfPolys;

    IPObjectStruct* accum = IritPrsrCopyObject(NULL, first, TRUE);
    if (accum == NULL)
        return ListOfPolys;

    /* Sequentially union each subsequent element into accum using IritBooleanOR. */
    for (i = 1; i < cnt; ++i) {
        IPObjectStruct* elem = IritPrsrListObjectGet(ListOfPolys, i);
        if (elem == NULL)
            continue;

        IPObjectStruct* newAccum = IritBooleanOR(accum, elem);
        /* IritBooleanOR returns a newly allocated object (or NULL on failure). */
        IritPrsrFreeObject(accum);
        accum = newAccum;

        if (accum == NULL) {
            /* Failure: leave original list intact. */
            return ListOfPolys;
        }
    }

    /* Convert accum (which is an IP_OBJ_POLY with possibly multiple polygons)
       into a LIST of polygon objects. */
    IPObjectStruct* OutList = IritPrsrGenLISTObject(NULL);
    if (OutList == NULL) {
        IritPrsrFreeObject(accum);
        return ListOfPolys;
    }

    IPPolygonStruct* Pl = accum->U.Pl;
    while (Pl != NULL) {
        IPPolygonStruct* NextPl = Pl->Pnext;
        Pl->Pnext = NULL; /* detach single polygon */

        IPObjectStruct* PolyObj = IritPrsrGenPOLYObject(Pl);
        if (PolyObj != NULL) {
            IP_SET_POLYGON_OBJ(PolyObj);
            PolyObj->Pnext = NULL;
            IritPrsrListObjectAppend(OutList, PolyObj);
        }
        else {
            /* If creation failed, free detached polygon (defensive). */
            IritPrsrFreePolygonList(Pl);
        }

        Pl = NextPl;
    }

    /* Prevent accum destructor from freeing polygons we moved. */
    accum->U.Pl = NULL;
    IritPrsrFreeObject(accum);

    /* Free original input list (we took ownership). */
    IritPrsrFreeObject(ListOfPolys);

    /* If union produced nothing, clean up and return NULL. */
    if (OutList->U.Lst.PObjList == NULL || OutList->U.Lst.ListMaxLen == 0) {
        IritPrsrFreeObject(OutList);
        return NULL;
    }

    return OutList;
}

// Modified projection function: accumulate per-view polygons into TempList,
   //run PerformBooleanUnion on the collected list and append only the union
   //result's elements to ResList. 
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

        /* Build stable orthonormal basis u,v,w from view matrix (no fallback). */
        IrtVecType u, v, w;
        if (!IritPrsrBuildViewBasisFromMat(ViewMats[vi], u, v, w, 0))
            continue;

        /* Temporary list to collect projected polygons for this view.*/
        IPObjectStruct* TempList = IritPrsrGenLISTObject(NULL);
        if (TempList == NULL)
            continue;

        /* For each silhouette polygon, create a projected closed polygon lying on the view plane. */
        int polyIdx = 0;
        for (IPPolygonStruct* Pl = SilObj->U.Pl; Pl != NULL; Pl = Pl->Pnext, ++polyIdx) {
            IPVertexStruct* V = Pl->PVertex;
            if (V == NULL)
                continue;

            /* Use helper to collect, sanitize duplicate closing vertex and compute centroid. */
            IrtPtType origin;
            int n;
            IrtPtType* coords = CollectPolygonCoordsAndCentroid(V, &n, origin);
            if (coords == NULL || n < 3)
                continue;

            /* Build world->local matrix MatLocal and its inverse.*/
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

                /* Append to temporary per-view list (do not append to global ResList yet). */
                IritPrsrListObjectAppend(TempList, PolyObj);
            }

            IritFree(coords);
        } /* per polygon */

        /* Perform 2D union on TempList. This helper takes ownership of TempList
           and returns a list whose elements are the unioned polygons (or the
           original TempList if no union is performed). */
        IPObjectStruct* UnionList = PerformBooleanUnion(TempList);
        /* UnionList now owns the polygon objects. */

        if (UnionList != NULL && IP_IS_OLST_OBJ(UnionList)) {
            /* Move union results into ResList as copies to be safe. */
            int uidx = 0;
            IPObjectStruct* UObj = NULL;
            int outPolyIdx = 0;
            while ((UObj = IritPrsrListObjectGet(UnionList, uidx++)) != NULL) {
                /* Copy polygon object to append into ResList (keeps ownership clear). */
                IPObjectStruct* Copy = IritPrsrCopyObject(NULL, UObj, TRUE);
                if (Copy != NULL) {
                    /* Overwrite name to indicate view and index (makes extruder mapping robust). */
                    char outName[256];
                    if (IsPre)
                        snprintf(outName, sizeof(outName), "proj_pre_view%d_poly%03d", vi, outPolyIdx);
                    else
                        snprintf(outName, sizeof(outName), "proj_dir_view%d_poly%03d", vi, outPolyIdx);
                    Copy->ObjName = IritMiscStrdup(outName);
                    IritPrsrListObjectAppend(ResList, Copy);
                    ++outPolyIdx;
                }
            }

            /* Free the union list container and its internal objects (we copied them). */
            IritPrsrFreeObject(UnionList);
        }
        else {
            /* Union produced nothing - ensure TempList was freed by helper or free it now. */
            if (UnionList != NULL)
                IritPrsrFreeObject(UnionList);
        }
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

    /* Process both pre and direct projected lists. */
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
            if (vi < 0 || vi >= NumViews)
                continue;

            /* Build view basis (allow fallback). */
            IrtVecType u, vtmp, f;
            if (!IritPrsrBuildViewBasisFromMat(Views[vi], u, vtmp, f, 1))
                continue;

            /* Make a working copy of the whole projected polygon object (preserves multiple loops / holes). */
            IPObjectStruct* PolyCopy = IritPrsrCopyObject(NULL, PObj, TRUE);
            if (PolyCopy == NULL)
                continue;

            /* Quick validity check: ensure at least one polygon loop has >= 3 vertices. */
            int hasValidLoop = 0;
            for (IPPolygonStruct* Pl = PolyCopy->U.Pl; Pl != NULL; Pl = Pl->Pnext) {
                if (Pl->PVertex == NULL)
                    continue;
                /* Count vertices in this loop (cheap guard). */
                int vcount = 0;
                IPVertexStruct* V = Pl->PVertex;
                int guard = 0;
                do {
                    ++vcount;
                    V = V->Pnext;
                    if (++guard > 200000) break;
                } while (V != NULL && V != Pl->PVertex);
                if (vcount >= 3) { hasValidLoop = 1; break; }
            }
            if (!hasValidLoop) {
                IritPrsrFreeObject(PolyCopy);
                continue;
            }

            /* Compute extrusion direction from view forward vector scaled by Depth. */
            IrtVecType Dir;
            Dir[0] = (IrtRType)(f[0] * Depth);
            Dir[1] = (IrtRType)(f[1] * Depth);
            Dir[2] = (IrtRType)(f[2] * Depth);

            /* Ensure polygon planes are up-to-date (optional but safe). */
            for (IPPolygonStruct* Pl = PolyCopy->U.Pl; Pl != NULL; Pl = Pl->Pnext) {
                IritPrsrUpdatePolyPlane(Pl);
            }

            /* Extrude the entire polygon object (keeps holes as holes if present). */
            IPObjectStruct* Extr = IritGeomPrimGenEXTRUDEObject(PolyCopy, Dir, 3);

            /* Free the intermediate copy (extruder returns a new object). */
            IritPrsrFreeObject(PolyCopy);

            if (Extr == NULL)
                continue;

            /* Name and append result. */
            char namebuf[512];
            if (usePre)
                snprintf(namebuf, sizeof(namebuf), "sil_pre_view%d_poly%03d_extr", vi, polyIdx);
            else
                snprintf(namebuf, sizeof(namebuf), "sil_dir_view%d_poly%03d_extr", vi, polyIdx);

            Extr->ObjName = IritMiscStrdup(namebuf);
            Extr->Pnext = NULL;
            IritPrsrListObjectAppend(ResList, Extr);
        }
    }

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



/* Replaced SelectBestViewSampling with a genetic sampling-based selector.
   Uses existing IRIT helpers to evaluate silhouette-based fitness.
   Returns an allocated array of view matrices (caller should IritFree it).
*/
//IrtHmgnMatType* SelectBestViewSampling(IPObjectStruct* PObj,
//    int NumSamples,
//    IrtHmgnMatType* ResultMat)
//{
//    /* Parameters adapted from the provided algorithm. */
//    const int N_POP = 30;
//    const int MAX_ITER = 15;
//    const int CONVERGE_ITER = 5;
//    const double RADIUS = 2.0;
//    const double P_UNI = 0.30;
//    const double P_GEO = 0.40;
//    const double P_NCI = 0.30;
//    const double P_GLO = 0.05;
//    const double P_LOC = 0.10;
//    const int K_NEIGHBORS = 10;
//    const double EPS_FITNESS_EQ = 1e-8;
//
//    int i, j;
//
//    if (PObj == NULL || NumSamples <= 0)
//        return NULL;
//
//    /* Allocate candidate view matrices and aux arrays. */
//    IrtHmgnMatType* ViewMats = (IrtHmgnMatType*)IritMalloc(sizeof(IrtHmgnMatType) * NumSamples);
//    if (ViewMats == NULL)
//        return NULL;
//
//    typedef struct {
//        int id;
//        IrtRType x, y, z;
//        IrtRType fitness; /* negative indicates not evaluated yet */
//    } Candidate;
//
//    Candidate* candidates = (Candidate*)IritMalloc(sizeof(Candidate) * NumSamples);
//    if (candidates == NULL) {
//        IritFree(ViewMats);
//        return NULL;
//    }
//
//    /* Generate Fibonacci sphere candidate positions and view matrices. */
//    {
//        IrtRType golden = (sqrt(5.0) - 1.0) / 2.0; /* fractional golden ratio as used elsewhere */
//        for (i = 0; i < NumSamples; ++i) {
//            IrtRType y = 1.0 - (i / (IrtRType)(NumSamples - 1)) * 2.0;
//            IrtRType radius = sqrt(IRIT_MAX(0.0, 1.0 - y * y));
//            IrtRType theta = 2.0 * M_PI * golden * i;
//            IrtRType x = cos(theta) * radius;
//            IrtRType z = sin(theta) * radius;
//
//            /* Scale to radius and create camera position */
//            IrtVecType CamPos;
//            CamPos[0] = (IrtRType)(x * RADIUS);
//            CamPos[1] = (IrtRType)(y * RADIUS);
//            CamPos[2] = (IrtRType)(z * RADIUS);
//
//            /* Choose an Up vector robustly */
//            IrtVecType Center = { 0.0, 0.0, 0.0 }, Up = { 0.0, 0.0, 1.0 };
//            if (fabs(CamPos[0]) < 0.01 && fabs(CamPos[1]) < 0.01) {
//                Up[0] = 1.0; Up[1] = 0.0; Up[2] = 0.0;
//            }
//
//            GenLookAtMatrix(CamPos, Center, Up, ViewMats[i]);
//
//            candidates[i].id = i;
//            candidates[i].x = CamPos[0];
//            candidates[i].y = CamPos[1];
//            candidates[i].z = CamPos[2];
//            candidates[i].fitness = -1.0; /* not evaluated */
//        }
//    }
//
//    /* Utility: evaluate fitness for a candidate index using silhouette area heuristic. */
//    auto evaluate_candidate = [&](int cidx) -> IrtRType {
//        if (cidx < 0 || cidx >= NumSamples)
//            return 0.0;
//
//        if (candidates[cidx].fitness >= 0.0)
//            return candidates[cidx].fitness; /* cached */
//
//        /* Use existing silhouette extractor for a single view and compute area. */
//        IPObjectStruct* SilList = IritPrsrGetSilhouettesForViews(PObj, &ViewMats[cidx], 1, 20, 0);
//        IrtRType score = 0.0;
//        if (SilList != NULL) {
//            IPObjectStruct* SilObj = IritPrsrListObjectGet(SilList, 0);
//            if (SilObj != NULL) {
//                score = CalcPolygonArea(SilObj);
//            }
//            IritPrsrFreeObject(SilList);
//        }
//        /* cache it */
//        candidates[cidx].fitness = score;
//        return score;
//        };
//
//    /* Utility: euclidean distance between two candidates */
//    auto candidate_dist = [&](int a, int b) -> double {
//        double dx = candidates[a].x - candidates[b].x;
//        double dy = candidates[a].y - candidates[b].y;
//        double dz = candidates[a].z - candidates[b].z;
//        return sqrt(dx * dx + dy * dy + dz * dz);
//        };
//
//    /* Utility: nearest candidate index to a given point */
//    auto findNearestCandidate = [&](double cx, double cy, double cz) -> int {
//        int best = -1;
//        double bestD = IRIT_INFNTY;
//        for (i = 0; i < NumSamples; ++i) {
//            double dx = candidates[i].x - cx;
//            double dy = candidates[i].y - cy;
//            double dz = candidates[i].z - cz;
//            double d2 = dx * dx + dy * dy + dz * dz;
//            if (d2 < bestD) {
//                bestD = d2;
//                best = i;
//            }
//        }
//        return best;
//        };
//
//    /* Utility: get K nearest neighbors indices (returns count in out array). */
//    /* Note: K is small, simple insertion sort on K best distances. */
//    auto getKNearest = [&](int targetId, int K, int* outIdx, int outMax) -> int {
//        if (targetId < 0 || targetId >= NumSamples || K <= 0 || outIdx == NULL)
//            return 0;
//        int k = IRIT_MIN(K, outMax);
//        /* pair of (dist, idx) */
//        double* bestD = (double*)IritMalloc(sizeof(double) * k);
//        int* bestI = (int*)IritMalloc(sizeof(int) * k);
//        for (i = 0; i < k; ++i) {
//            bestD[i] = IRIT_INFNTY;
//            bestI[i] = -1;
//        }
//        for (i = 0; i < NumSamples; ++i) {
//            if (i == targetId) continue;
//            double d = candidate_dist(targetId, i);
//            /* insert into best lists if smaller than some entry */
//            for (int j = 0; j < k; ++j) {
//                if (d < bestD[j]) {
//                    /* shift down */
//                    int t;
//                    for (t = k - 1; t > j; --t) {
//                        bestD[t] = bestD[t - 1];
//                        bestI[t] = bestI[t - 1];
//                    }
//                    bestD[j] = d;
//                    bestI[j] = i;
//                    break;
//                }
//            }
//        }
//        int outCount = 0;
//        for (i = 0; i < k; ++i) {
//            if (bestI[i] >= 0) {
//                outIdx[outCount++] = bestI[i];
//            }
//        }
//        IritFree(bestD);
//        IritFree(bestI);
//        return outCount;
//        };
//
//    /* Initialize RNG */
//    srand((unsigned)time(NULL));
//
//    /* -------------- Initialization: Farthest Point Sampling -------------- */
//    int popSize = IRIT_MIN(N_POP, NumSamples);
//    int* population = (int*)IritMalloc(sizeof(int) * popSize);
//    int* selected = (int*)IritMalloc(sizeof(int) * popSize);
//    for (i = 0; i < popSize; ++i) selected[i] = -1;
//
//    /* pick first at random */
//    selected[0] = rand() % NumSamples;
//    for (i = 1; i < popSize; ++i) {
//        double bestMinD = -1.0;
//        int bestIdx = -1;
//        for (j = 0; j < NumSamples; ++j) {
//            int already = 0;
//            for (int s = 0; s < i; ++s) if (selected[s] == j) { already = 1; break; }
//            if (already) continue;
//            double minD = IRIT_INFNTY;
//            for (int s = 0; s < i; ++s) {
//                double d = candidate_dist(j, selected[s]);
//                if (d < minD) minD = d;
//            }
//            if (minD > bestMinD) {
//                bestMinD = minMinD = minD; /* silence potential uninitialized warning pattern */
//                bestIdx = j;
//            }
//        }
//        if (bestIdx < 0)
//            bestIdx = rand() % NumSamples;
//        selected[i] = bestIdx;
//    }
//
//    /* Evaluate initial population fitness and fill population array */
//    for (i = 0; i < popSize; ++i) {
//        population[i] = selected[i];
//        evaluate_candidate(population[i]);
//    }
//
//    IritFree(selected);
//
//    /* ---------- Genetic optimization loop ---------- */
//    int iter = 0;
//    int unchanged = 0;
//    IrtRType bestEver = -1.0;
//
//    while (iter < MAX_ITER && unchanged < CONVERGE_ITER) {
//        /* Generate children until we have popSize children */
//        int* children = (int*)IritMalloc(sizeof(int) * popSize);
//        int childCount = 0;
//
//        /* Precompute population fitness sum for roulette */
//        IrtRType popSum = 0.0;
//        for (i = 0; i < popSize; ++i) {
//            popSum += candidates[population[i]].fitness > 0.0 ? candidates[population[i]].fitness : 0.0;
//        }
//
//        auto roulette_select = [&](void) -> int {
//            if (popSum <= EPS_FITNESS_EQ) {
//                /* fallback to uniform random parent index */
//                return population[rand() % popSize];
//            }
//            double r = ((double)rand() / (double)RAND_MAX) * popSum;
//            double acc = 0.0;
//            for (int pi = 0; pi < popSize; ++pi) {
//                acc += candidates[population[pi]].fitness;
//                if (acc >= r) return population[pi];
//            }
//            return population[popSize - 1];
//            };
//
//        while (childCount < popSize) {
//            /* selection */
//            int p1 = roulette_select();
//            int p2 = roulette_select();
//            /* ensure distinct parents (if possible) */
//            int tries = 0;
//            while (p1 == p2 && tries++ < 8) p2 = roulette_select();
//
//            /* Crossover */
//            double rC = (double)rand() / (double)RAND_MAX;
//            int childIdx = p1;
//            if (rC < P_UNI) {
//                /* uniform: pick either parent */
//                childIdx = (rand() & 1) ? p1 : p2;
//            }
//            else if (rC < P_UNI + P_GEO) {
//                /* geodesic SLERP midpoint on sphere between p1,p2 */
//                double ux1 = candidates[p1].x, uy1 = candidates[p1].y, uz1 = candidates[p1].z;
//                double ux2 = candidates[p2].x, uy2 = candidates[p2].y, uz2 = candidates[p2].z;
//                double dot = (ux1 * ux2 + uy1 * uy2 + uz1 * uz2) / (RADIUS * RADIUS);
//                if (dot > 1.0) dot = 1.0;
//                if (dot < -1.0) dot = -1.0;
//                double omega = acos(dot);
//                if (fabs(omega) < 1e-8) {
//                    childIdx = p1;
//                }
//                else {
//                    double t = 0.5;
//                    double s_omega = sin(omega);
//                    double w1 = sin((1.0 - t) * omega) / s_omega;
//                    double w2 = sin(t * omega) / s_omega;
//                    double cx = w1 * ux1 + w2 * ux2;
//                    double cy = w1 * uy1 + w2 * uy2;
//                    double cz = w1 * uz1 + w2 * uz2;
//                    /* project back to nearest candidate */
//                    childIdx = findNearestCandidate(cx, cy, cz);
//                    if (childIdx < 0) childIdx = p1;
//                }
//            }
//            else {
//                /* neighborhood crossover: pick a neighbor of random parent */
//                int base = (rand() & 1) ? p1 : p2;
//                int neighBuf[K_NEIGHBORS];
//                int ncnt = getKNearest(base, K_NEIGHBORS, neighBuf, K_NEIGHBORS);
//                if (ncnt > 0) {
//                    childIdx = neighBuf[rand() % ncnt];
//                }
//                else {
//                    childIdx = base;
//                }
//            }
//
//            /* Mutation */
//            double rM = (double)rand() / (double)RAND_MAX;
//            if (rM < P_GLO) {
//                childIdx = rand() % NumSamples; /* global reset */
//            }
//            else if (rM < P_GLO + P_LOC) {
//                int neighBuf[K_NEIGHBORS];
//                int ncnt = getKNearest(childIdx, K_NEIGHBORS, neighBuf, K_NEIGHBORS);
//                if (ncnt > 0) childIdx = neighBuf[rand() % ncnt];
//            }
//
//            /* ensure valid */
//            if (childIdx < 0 || childIdx >= NumSamples) childIdx = rand() % NumSamples;
//
//            /* evaluate child fitness if needed */
//            evaluate_candidate(childIdx);
//
//            children[childCount++] = childIdx;
//        }
//
//        /* Combine population + children into a pool, pick top popSize by fitness (elitist) */
//        int poolSize = popSize + childCount;
//        int* pool = (int*)IritMalloc(sizeof(int) * poolSize);
//        for (i = 0; i < popSize; ++i) pool[i] = population[i];
//        for (i = 0; i < childCount; ++i) pool[popSize + i] = children[i];
//
//        /* pick top popSize (selection by repeated max) */
//        int* newPop = (int*)IritMalloc(sizeof(int) * popSize);
//        char* used = (char*)IritMalloc(poolSize);
//        memset(used, 0, poolSize);
//        for (i = 0; i < popSize; ++i) {
//            int bestIdx = -1;
//            IrtRType bestFit = -IRIT_INFNTY;
//            for (j = 0; j < poolSize; ++j) {
//                if (used[j]) continue;
//                IrtRType f = candidates[pool[j]].fitness;
//                if (f > bestFit) {
//                    bestFit = f;
//                    bestIdx = j;
//                }
//            }
//            if (bestIdx < 0) {
//                /* fallback */
//                bestIdx = 0;
//                while (bestIdx < poolSize && used[bestIdx]) ++bestIdx;
//                if (bestIdx >= poolSize) newPop[i] = pool[0];
//                else newPop[i] = pool[bestIdx];
//            }
//            else {
//                newPop[i] = pool[bestIdx];
//                used[bestIdx] = 1;
//            }
//        }
//
//        /* replace population */
//        for (i = 0; i < popSize; ++i) population[i] = newPop[i];
//
//        /* convergence check */
//        /* sort population fitness to find current best quickly */
//        IrtRType currentBest = -1.0;
//        for (i = 0; i < popSize; ++i) {
//            IrtRType f = candidates[population[i]].fitness;
//            if (f > currentBest) currentBest = f;
//        }
//        if (fabs(currentBest - bestEver) < 1e-6) {
//            ++unchanged;
//        }
//        else {
//            unchanged = 0;
//            bestEver = currentBest;
//        }
//
//        /* free temp */
//        IritFree(children);
//        IritFree(pool);
//        IritFree(newPop);
//        IritFree(used);
//
//        ++iter;
//    } /* end GA loop */
//
//    /* Best candidate index is the population member with highest fitness */
//    int bestIdx = population[0];
//    IrtRType bestFit = candidates[bestIdx].fitness;
//    for (i = 1; i < popSize; ++i) {
//        IrtRType f = candidates[population[i]].fitness;
//        if (f > bestFit) { bestFit = f; bestIdx = population[i]; }
//    }
//
//    /* Copy resulting best matrix into ResultMat if provided. */
//    if (ResultMat != NULL && bestIdx >= 0 && bestIdx < NumSamples) {
//        IRIT_GEN_COPY(ResultMat, ViewMats[bestIdx], sizeof(IrtHmgnMatType));
//    }
//
//    /* Clean up temporary arrays but keep ViewMats for caller use. */
//    IritFree(candidates);
//    IritFree(population);
//
//    return ViewMats;
//}
