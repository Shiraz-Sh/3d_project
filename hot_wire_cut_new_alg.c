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

/* Helper: sequential 2D union using IritBoolean2D. Takes ownership of
   ListOfPolys (a LIST of polygon objects). Returns a LIST of polygons
   that represent the unioned result or NULL on failure.
*/
static IPObjectStruct* PerformBooleanUnion2D(IPObjectStruct* ListOfPolys)
{
    int i, cnt = 0;

    if (ListOfPolys == NULL)
        return NULL;

    if (ListOfPolys->ObjType != IP_OBJ_LIST_OBJ)
        return ListOfPolys;

    /* Count elements in the list. */
    while (IritPrsrListObjectGet(ListOfPolys, cnt) != NULL)
        ++cnt;

    if (cnt <= 1) /* nothing to union */
        return ListOfPolys;

    /* Start with a copy of the first element as accumulator (safe ownership). */
    IPObjectStruct* firstObj = IritPrsrListObjectGet(ListOfPolys, 0);
    if (firstObj == NULL)
        return ListOfPolys;

    IPObjectStruct* accumObj = IritPrsrCopyObject(NULL, firstObj, TRUE);
    if (accumObj == NULL)
        return ListOfPolys;

    /* Sequentially union each subsequent element using the 2D boolean routine.
       Note: IritBoolean2D takes IPPolygonStruct* (polygon lists) and a BoolOperType. */
    for (i = 1; i < cnt; ++i) {
        IPObjectStruct* elemObj = IritPrsrListObjectGet(ListOfPolys, i);
        if (elemObj == NULL)
            continue;

        IPPolygonStruct* resPl = IritBoolean2D(accumObj->U.Pl, elemObj->U.Pl, BOOL_OPER_OR);

        /* Free the previous accumulator object (we copied it earlier). */
        IritPrsrFreeObject(accumObj);
        accumObj = NULL;

        if (resPl == NULL) {
            /* Failure: leave original list intact (we didn't touch it). */
            return ListOfPolys;
        }

        /* Wrap the returned polygon list into an IPObjectStruct for the next iteration. */
        accumObj = IritPrsrGenPOLYObject(resPl);
        if (accumObj == NULL) {
            /* Defensive: free polygons if wrapping failed. */
            IritPrsrFreePolygonList(resPl);
            return ListOfPolys;
        }
        IP_SET_POLYGON_OBJ(accumObj);
    }

    /* Convert the accumulated polygon list into a LIST object of polygon objects. */
    IPObjectStruct* OutList = IritPrsrGenLISTObject(NULL);
    if (OutList == NULL) {
        IritPrsrFreeObject(accumObj);
        return ListOfPolys;
    }

    IPPolygonStruct* Pl = accumObj->U.Pl;
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
            IritPrsrFreePolygonList(Pl);
        }

        Pl = NextPl;
    }

    /* Prevent accumObj destructor from freeing polygons we moved. */
    accumObj->U.Pl = NULL;
    IritPrsrFreeObject(accumObj);

    /* Free original input list (we took ownership). */
    IritPrsrFreeObject(ListOfPolys);

    if (OutList->U.Lst.PObjList == NULL || OutList->U.Lst.ListMaxLen == 0) {
        IritPrsrFreeObject(OutList);
        return NULL;
    }

    return OutList;
}

/* Helper: transform a polygon object (IPObjectStruct which is a polygon) by Mat.
   Returns a newly allocated IPObjectStruct (polygon) or NULL. Copies polygons/verts.
*/
static IPObjectStruct* TransformPolyObjectByMat(const IPObjectStruct* SrcObj, const IrtHmgnMatType Mat)
{
    if (SrcObj == NULL || !IP_IS_POLY_OBJ((IPObjectStruct*)SrcObj))
        return NULL;

    IPPolygonStruct* Pl = SrcObj->U.Pl;
    if (Pl == NULL)
        return NULL;

    /* We assume SrcObj may contain multiple polygons, but SrcObj here is one polygon
       object as returned by union list elements. Create a copy of the whole polygon list. */
    IPPolygonStruct* NewPlHead = NULL;
    IPPolygonStruct* NewPlTail = NULL;

    for (; Pl != NULL; Pl = Pl->Pnext) {
        /* Create new polygon and copy transformed vertices. */
        IPPolygonStruct* NewPl = IritPrsrAllocPolygon(0, NULL, NULL);
        IPVertexStruct* FirstV = NULL, * PrevV = NULL;

        IPVertexStruct* V = Pl->PVertex;
        if (V == NULL) {
            IritPrsrFreePolygonList(NewPl);
            continue;
        }

        do {
            IrtPtType localP, worldP;
            IritMiscMatMultPtby4by4(localP, V->Coord, Mat);

            IPVertexStruct* NV = IritPrsrAllocVertex2(NULL);
            NV->Coord[0] = worldP[0] = localP[0];
            NV->Coord[1] = worldP[1] = localP[1];
            NV->Coord[2] = worldP[2] = localP[2];

            if (FirstV == NULL) FirstV = PrevV = NV;
            else { PrevV->Pnext = NV; PrevV = NV; }

            V = V->Pnext;
        } while (V != NULL && V != Pl->PVertex);

        if (PrevV != NULL) {
            PrevV->Pnext = FirstV;
            NewPl->PVertex = FirstV;
            IritPrsrUpdatePolyPlane(NewPl);

            if (NewPlHead == NULL) NewPlHead = NewPlTail = NewPl;
            else {
                NewPlTail->Pnext = NewPl;
                NewPlTail = NewPl;
            }
        }
        else {
            IritPrsrFreePolygonList(NewPl);
        }
    }

    if (NewPlHead == NULL)
        return NULL;

    IPObjectStruct* OutObj = IritPrsrGenPOLYObject(NewPlHead);
    if (OutObj != NULL) {
        IP_SET_POLYGON_OBJ(OutObj);
    }
    else {
        IritPrsrFreePolygonList(NewPlHead);
    }

    return OutObj;
}

/* Corrected IritPrsrProjectSilhouetteToViewPolys:
   - Projects ALL polygons for a view into one common view-local frame (MatLocal with no per-polygon translation),
   - Calls PerformBooleanUnion2D to union flat polygons (Z=0),
   - Transforms unioned polygons back to world coordinates using InvMat.
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

        /* Build stable orthonormal basis u, v, w from view matrix. */
        IrtVecType u, v, w;
        if (!IritPrsrBuildViewBasisFromMat(ViewMats[vi], u, v, w, 0))
            continue;

        /* Build a single MatLocal for the view (world -> view-local).
           We use zero translation so all polygons are placed in the same 2D frame. */
        IrtHmgnMatType MatLocal, InvMat;
        IritMiscMatGenUnitMat(MatLocal);
        MatLocal[0][0] = u[0]; MatLocal[0][1] = u[1]; MatLocal[0][2] = u[2]; MatLocal[0][3] = 0.0;
        MatLocal[1][0] = v[0]; MatLocal[1][1] = v[1]; MatLocal[1][2] = v[2]; MatLocal[1][3] = 0.0;
        MatLocal[2][0] = w[0]; MatLocal[2][1] = w[1]; MatLocal[2][2] = w[2]; MatLocal[2][3] = 0.0;
        MatLocal[3][0] = 0.0;  MatLocal[3][1] = 0.0;  MatLocal[3][2] = 0.0; MatLocal[3][3] = 1.0;

        if (!IritMiscMatInverseMatrix(MatLocal, InvMat)) {
            continue;
        }

        /* Collect all polygons projected into the view-local XY plane (Z=0). */
        IPObjectStruct* LocalPolyList = IritPrsrGenLISTObject(NULL);
        if (LocalPolyList == NULL)
            continue;

        int polyIdx = 0;
        for (IPPolygonStruct* Pl = SilObj->U.Pl; Pl != NULL; Pl = Pl->Pnext, ++polyIdx) {
            IPVertexStruct* V = Pl->PVertex;
            if (V == NULL)
                continue;

            /* Collect coords and skip degenerate polygons. */
            IrtPtType origin;
            int n;
            IrtPtType* coords = CollectPolygonCoordsAndCentroid(V, &n, origin);
            if (coords == NULL || n < 3)
                continue;

            /* Create new polygon in local frame */
            IPPolygonStruct* NewPl = IritPrsrAllocPolygon(0, NULL, NULL);
            IPVertexStruct* FirstV = NULL;
            IPVertexStruct* PrevV = NULL;

            for (int j = 0; j < n; ++j) {
                IrtPtType localP;
                IritMiscMatMultPtby4by4(localP, coords[j], MatLocal);
                localP[2] = 0.0; /* flatten to XY */

                IPVertexStruct* NV = IritPrsrAllocVertex2(NULL);
                NV->Coord[0] = localP[0];
                NV->Coord[1] = localP[1];
                NV->Coord[2] = localP[2];

                if (FirstV == NULL) FirstV = PrevV = NV;
                else { PrevV->Pnext = NV; PrevV = NV; }
            }

            if (PrevV != NULL) {
                PrevV->Pnext = FirstV;
                NewPl->PVertex = FirstV;
                IritPrsrUpdatePolyPlane(NewPl);

                IPObjectStruct* PolyObj = IritPrsrGenPOLYObject(NewPl);
                if (PolyObj != NULL) {
                    IP_SET_POLYGON_OBJ(PolyObj);
                    PolyObj->Pnext = NULL;
                    /* name not necessary for union step, but keep it if you want */
                    IritPrsrListObjectAppend(LocalPolyList, PolyObj);
                }
                else {
                    IritPrsrFreePolygonList(NewPl);
                }
            }
            else {
                IritPrsrFreePolygonList(NewPl);
            }

            IritFree(coords);
        } /* per polygon */

        /* Perform 2D boolean union on the local list */
        IPObjectStruct* UnionLocal = PerformBooleanUnion2D(LocalPolyList);
        if (UnionLocal == NULL) {
            /* nothing produced or failure */
            continue;
        }

        /* Transform union result back to world coordinates by applying InvMat to each union polygon. */
        int uidx = 0;
        IPObjectStruct* UObj = NULL;
        while ((UObj = IritPrsrListObjectGet(UnionLocal, uidx++)) != NULL) {
            /* Transform this union element back to world */
            IPObjectStruct* WorldObj = TransformPolyObjectByMat(UObj, InvMat);
            if (WorldObj == NULL)
                continue;

            /* Name result so extruder can find view index */
            char outName[256];
            if (IsPre)
                snprintf(outName, sizeof(outName), "proj_pre_view%d_poly%03d", vi, uidx - 1);
            else
                snprintf(outName, sizeof(outName), "proj_dir_view%d_poly%03d", vi, uidx - 1);
            WorldObj->ObjName = IritMiscStrdup(outName);
            WorldObj->Pnext = NULL;
            IritPrsrListObjectAppend(ResList, WorldObj);
        }

        IritPrsrFreeObject(UnionLocal);
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



/* -----------------------------------------------------------------------
 * New helper: approximate the 2D outer contour from a view by shrinking
 * a circle of control points toward the unioned projected polygon boundary.
 *
 * Returns a POLY IPObjectStruct* in world coordinates (caller owns).
 *
 * Parameters:
 *   Solid      - polygonal solid to project (IP_OBJ_POLY).
 *   ViewMat    - view matrix to build local XY frame (same convention as other code).
 *   NumCtrl    - number of control points on initial circle (>= 8 recommended).
 *   Iterations - number of shrink iterations (10..200).
 *   ShrinkStep - fraction [0,1] how much to move toward boundary each iter (0.1..0.6).
 *   SmoothW    - Laplacian smoothing weight (0.0..1.0).
 *
 * Notes: This is a robust, simple heuristic — not an exact Newton optimizer,
 * but it respects "stay outside" by snapping inside points to the boundary.
 * ---------------------------------------------------------------------*/
static IrtRType DistPointSegment2D(const IrtPtType p, const IrtPtType a, const IrtPtType b, IrtPtType out)
{
    /* p, a, b are 2D stored in [0],[1]; out is 2D result */
    IrtRType vx = b[0] - a[0];
    IrtRType vy = b[1] - a[1];
    IrtRType wx = p[0] - a[0];
    IrtRType wy = p[1] - a[1];
    IrtRType c1 = vx * wx + vy * wy;
    IrtRType c2 = vx * vx + vy * vy;
    if (c2 <= IRIT_EPS) {
        out[0] = a[0]; out[1] = a[1];
        return sqrt((p[0] - a[0]) * (p[0] - a[0]) + (p[1] - a[1]) * (p[1] - a[1]));
    }
    IrtRType t = c1 / c2;
    if (t < 0.0) t = 0.0;
    else if (t > 1.0) t = 1.0;
    out[0] = a[0] + t * vx;
    out[1] = a[1] + t * vy;
    return sqrt((p[0] - out[0]) * (p[0] - out[0]) + (p[1] - out[1]) * (p[1] - out[1]));
}

/* Find closest point on a polygon list (UnionObj assumed to be LIST of POLY objects
   in the view-local XY plane). Returns TRUE if found and writes closest into out.
   Dist is returned as function value (large on failure). */
static IrtRType FindClosestPointOnPolyList(IPObjectStruct* UnionObj, const IrtPtType p, IrtPtType out)
{
    if (UnionObj == NULL || UnionObj->ObjType != IP_OBJ_LIST_OBJ)
        return IRIT_INFNTY;

    IrtRType best = IRIT_INFNTY;
    IPObjectStruct* Poly = NULL;
    int idx = 0;
    while ((Poly = IritPrsrListObjectGet(UnionObj, idx++)) != NULL) {
        if (!IP_IS_POLY_OBJ(Poly))
            continue;
        for (IPPolygonStruct* Pl = Poly->U.Pl; Pl != NULL; Pl = Pl->Pnext) {
            IPVertexStruct* V = Pl->PVertex;
            if (V == NULL) continue;
            IPVertexStruct* cur = V;
            do {
                IPVertexStruct* nxt = cur->Pnext ? cur->Pnext : Pl->PVertex;
                IrtPtType a, b;
                a[0] = cur->Coord[0]; a[1] = cur->Coord[1];
                b[0] = nxt->Coord[0]; b[1] = nxt->Coord[1];
                IrtPtType cand;
                IrtRType d = DistPointSegment2D(p, a, b, cand);
                if (d < best) {
                    best = d;
                    out[0] = cand[0];
                    out[1] = cand[1];
                }
                cur = cur->Pnext;
            } while (cur != NULL && cur != V);
        }
    }
    return best;
}

/* Ray-casting point-in-polygon for a single polygon loop (2D). */
static int PointInPolyLoop2D(const IrtPtType p, IPPolygonStruct* Pl)
{
    if (Pl == NULL || Pl->PVertex == NULL) return 0;
    int wn = 0;
    IPVertexStruct* V = Pl->PVertex;
    IPVertexStruct* cur = V;
    do {
        IrtRType x1 = cur->Coord[0], y1 = cur->Coord[1];
        IPVertexStruct* nxt = cur->Pnext ? cur->Pnext : V;
        IrtRType x2 = nxt->Coord[0], y2 = nxt->Coord[1];
        if (((y1 <= p[1]) && (y2 > p[1])) || ((y1 > p[1]) && (y2 <= p[1]))) {
            IrtRType vt = (p[1] - y1) / (y2 - y1);
            IrtRType xproj = x1 + vt * (x2 - x1);
            if (xproj < p[0])
                wn ^= 1;
        }
        cur = nxt;
    } while (cur != NULL && cur != V);
    return wn;
}

/* Test point inside any polygon in a LIST (view-local XY). */
static int PointInPolyList2D(IPObjectStruct* UnionObj, const IrtPtType p)
{
    if (UnionObj == NULL || UnionObj->ObjType != IP_OBJ_LIST_OBJ)
        return 0;
    IPObjectStruct* Poly = NULL;
    int idx = 0;
    while ((Poly = IritPrsrListObjectGet(UnionObj, idx++)) != NULL) {
        if (!IP_IS_POLY_OBJ(Poly))
            continue;
        for (IPPolygonStruct* Pl = Poly->U.Pl; Pl != NULL; Pl = Pl->Pnext) {
            if (PointInPolyLoop2D(p, Pl))
                return 1;
        }
    }
    return 0;
}

/* Build union of Solid projected into view-local frame (MatLocal) and return
   LIST of filled polygons in view-local coordinates. Caller owns result.

   Reworked: avoid calling the fragile 2D boolean routines that may crash
   on degenerate inputs. Instead, sanitize per-face loops, gather unique
   projected points and produce a conservative outer contour via convex hull.
*/
static int PtCmpForQSort(const void* a, const void* b)
{
    const IrtPtType* pa = (const IrtPtType*)a;
    const IrtPtType* pb = (const IrtPtType*)b;
    if ((*pa)[0] < (*pb)[0]) return -1;
    if ((*pa)[0] > (*pb)[0]) return 1;
    if ((*pa)[1] < (*pb)[1]) return -1;
    if ((*pa)[1] > (*pb)[1]) return 1;
    return 0;
}

static IrtRType Cross2D(const IrtPtType a, const IrtPtType b, const IrtPtType c)
{
    /* cross of AB x AC in XY */
    return (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]);
}

static IPObjectStruct* BuildProjectedUnionLocal(IPObjectStruct* Solid, const IrtHmgnMatType MatLocal)
{
    if (Solid == NULL || !IP_IS_POLY_OBJ(Solid))
        return NULL;

    /* Collect cleaned projected points across all polygons. */
    const IrtRType eps = 1e-7;
    const IrtRType eps2 = IRIT_SQR(eps);

    /* Temporary dynamic array of points (initial capacity). */
    int cap = 256;
    int npts = 0;
    IrtPtType* pts = (IrtPtType*)IritMalloc(sizeof(IrtPtType) * cap);

    for (IPPolygonStruct* Pl = Solid->U.Pl; Pl != NULL; Pl = Pl->Pnext) {
        IPVertexStruct* V = Pl->PVertex;
        if (V == NULL) continue;

        /* Count vertices */
        int cnt = 0;
        IPVertexStruct* cur = V;
        int guard = 0;
        do {
            ++cnt;
            cur = cur->Pnext;
            if (++guard > 200000) { cnt = 0; break; }
        } while (cur != NULL && cur != V);
        if (cnt <= 0) continue;

        /* Collect transformed coordinates into an array and sanitize consecutive duplicates. */
        IrtPtType* coords = (IrtPtType*)IritMalloc(sizeof(IrtPtType) * cnt);
        cur = V;
        guard = 0;
        int idx = 0;
        while (cur != NULL && idx < cnt) {
            IrtPtType localP;
            IritMiscMatMultPtby4by4(localP, cur->Coord, MatLocal);
            coords[idx][0] = localP[0];
            coords[idx][1] = localP[1];
            coords[idx][2] = localP[2];
            ++idx;
            cur = cur->Pnext;
            if (++guard > 200000) break;
            if (cur == V) break;
        }
        if (idx != cnt) { IritFree(coords); continue; }

        /* Remove near-duplicate consecutive points and trailing duplicate */
        int m = 0;
        for (int i = 0; i < cnt; ++i) {
            if (m == 0) {
                IRIT_PT_COPY(coords[m++], coords[i]);
            }
            else {
                IrtRType dx = coords[i][0] - coords[m - 1][0];
                IrtRType dy = coords[i][1] - coords[m - 1][1];
                if ((dx * dx + dy * dy) > eps2) {
                    IRIT_PT_COPY(coords[m++], coords[i]);
                }
            }
        }
        if (m > 1) {
            IrtRType dx = coords[0][0] - coords[m - 1][0];
            IrtRType dy = coords[0][1] - coords[m - 1][1];
            if ((dx * dx + dy * dy) <= eps2)
                --m;
        }

        /* Skip degenerate loops */
        if (m < 1) { IritFree(coords); continue; }

        /* Append distinct points into global pts array (we'll unique them later). */
        for (int k = 0; k < m; ++k) {
            /* ensure capacity */
            if (npts + 1 > cap) {
                int newcap = cap * 2;
                IrtPtType* newpts = (IrtPtType*)IritMalloc(sizeof(IrtPtType) * newcap);
                memcpy(newpts, pts, sizeof(IrtPtType) * npts);
                IritFree(pts);
                pts = newpts;
                cap = newcap;
            }
            IRIT_PT_COPY(pts[npts++], coords[k]);
        }

        IritFree(coords);
    }

    if (npts < 3) {
        IritFree(pts);
        return NULL;
    }

    /* Sort points and unique them (by XY) */
    qsort(pts, npts, sizeof(IrtPtType), PtCmpForQSort);
    int uniq = 0;
    for (int i = 0; i < npts; ++i) {
        if (uniq == 0) IRIT_PT_COPY(pts[uniq++], pts[i]);
        else {
            IrtRType dx = pts[i][0] - pts[uniq - 1][0];
            IrtRType dy = pts[i][1] - pts[uniq - 1][1];
            if ((dx * dx + dy * dy) > eps2) {
                IRIT_PT_COPY(pts[uniq++], pts[i]);
            }
        }
    }
    npts = uniq;

    if (npts < 3) {
        IritFree(pts);
        return NULL;
    }

    /* Monotone chain convex hull (2D) */
    IrtPtType* hull = (IrtPtType*)IritMalloc(sizeof(IrtPtType) * npts * 2);
    int k = 0;
    /* lower hull */
    for (int i = 0; i < npts; ++i) {
        while (k >= 2 && Cross2D(hull[k - 2], hull[k - 1], pts[i]) <= 0.0) --k;
        IRIT_PT_COPY(hull[k++], pts[i]);
    }
    /* upper hull */
    int t = k + 1;
    for (int i = npts - 2; i >= 0; --i) {
        while (k >= t && Cross2D(hull[k - 2], hull[k - 1], pts[i]) <= 0.0) --k;
        IRIT_PT_COPY(hull[k++], pts[i]);
    }
    if (k <= 2) { IritFree(pts); IritFree(hull); return NULL; }
    int hull_n = k - 1; /* last point repeats first */

    /* Build POLY object from hull (view-local XY, Z = 0) */
    IPPolygonStruct* NewPl = IritPrsrAllocPolygon(0, NULL, NULL);
    IPVertexStruct* FirstV = NULL, * PrevV = NULL;
    for (int i = 0; i < hull_n; ++i) {
        IPVertexStruct* NV = IritPrsrAllocVertex2(NULL);
        NV->Coord[0] = hull[i][0];
        NV->Coord[1] = hull[i][1];
        NV->Coord[2] = 0.0;
        if (FirstV == NULL) FirstV = PrevV = NV;
        else { PrevV->Pnext = NV; PrevV = NV; }
    }
    if (PrevV != NULL) {
        PrevV->Pnext = FirstV;
        NewPl->PVertex = FirstV;
        IritPrsrUpdatePolyPlane(NewPl);
    }
    else {
        IritPrsrFreePolygonList(NewPl);
        IritFree(pts);
        IritFree(hull);
        return NULL;
    }

    IritFree(pts);
    IritFree(hull);

    IPObjectStruct* PolyObj = IritPrsrGenPOLYObject(NewPl);
    if (PolyObj == NULL) {
        IritPrsrFreePolygonList(NewPl);
        return NULL;
    }
    IP_SET_POLYGON_OBJ(PolyObj);
    PolyObj->Pnext = NULL;

    /* Put single polygon into a LIST object to match callers' expectations. */
    IPObjectStruct* OutList = IritPrsrGenLISTObject(NULL);
    if (OutList == NULL) {
        IritPrsrFreeObject(PolyObj);
        return NULL;
    }
    IritPrsrListObjectAppend(OutList, PolyObj);

    return OutList;
}


/* Main approximation driver:
   Projects Solid, computes union, shrinks a circle into an outer contour,
   returns a world-space POLY object approximating the outer contour. */
IPObjectStruct* IritPrsrApproxBSplineContourFromSolidView(IPObjectStruct* Solid,
    const IrtHmgnMatType ViewMat,
    int NumCtrl,
    int Iterations,
    IrtRType ShrinkStep,
    IrtRType SmoothW)
{
    if (Solid == NULL || !IP_IS_POLY_OBJ(Solid) || NumCtrl < 6 || Iterations <= 0)
        return NULL;

    /* Build view-local frame like ProjectSilhouetteToViewPolys */
    IrtVecType u, v, w;
    if (!IritPrsrBuildViewBasisFromMat(ViewMat, u, v, w, 1))
        return NULL;

    IrtHmgnMatType MatLocal, InvMat;
    IritMiscMatGenUnitMat(MatLocal);
    MatLocal[0][0] = u[0]; MatLocal[0][1] = u[1]; MatLocal[0][2] = u[2]; MatLocal[0][3] = 0.0;
    MatLocal[1][0] = v[0]; MatLocal[1][1] = v[1]; MatLocal[1][2] = v[2]; MatLocal[1][3] = 0.0;
    MatLocal[2][0] = w[0]; MatLocal[2][1] = w[1]; MatLocal[2][2] = w[2]; MatLocal[2][3] = 0.0;
    MatLocal[3][0] = 0.0;  MatLocal[3][1] = 0.0;  MatLocal[3][2] = 0.0; MatLocal[3][3] = 1.0;

    if (!IritMiscMatInverseMatrix(MatLocal, InvMat))
        return NULL;

    /* Build union of projected faces in local XY. InvMat maps World -> Local. */
    IPObjectStruct* UnionLocal = BuildProjectedUnionLocal(Solid, InvMat);
    if (UnionLocal == NULL) {
        return NULL;
    }

    /* Compute bounding box of union polygons */
    IrtRType minx = IRIT_INFNTY, miny = IRIT_INFNTY, maxx = -IRIT_INFNTY, maxy = -IRIT_INFNTY;
    IPObjectStruct* Poly = NULL;
    int pidx = 0;
    while ((Poly = IritPrsrListObjectGet(UnionLocal, pidx++)) != NULL) {
        for (IPPolygonStruct* Pl = Poly->U.Pl; Pl != NULL; Pl = Pl->Pnext) {
            IPVertexStruct* V = Pl->PVertex;
            if (V == NULL) continue;
            IPVertexStruct* cur = V;
            do {
                if (cur->Coord[0] < minx) minx = cur->Coord[0];
                if (cur->Coord[1] < miny) miny = cur->Coord[1];
                if (cur->Coord[0] > maxx) maxx = cur->Coord[0];
                if (cur->Coord[1] > maxy) maxy = cur->Coord[1];
                cur = cur->Pnext;
            } while (cur != NULL && cur != V);
        }
    }
    if (minx > maxx || miny > maxy) {
        IritPrsrFreeObject(UnionLocal);
        return NULL;
    }

    IrtRType cx = 0.5 * (minx + maxx);
    IrtRType cy = 0.5 * (miny + maxy);
    IrtRType radius = 0.5 * IRIT_MAX(maxx - minx, maxy - miny) * 1.2;

    /* Initialize control points on circle in view-local XY */
    IrtPtType* ctrl = (IrtPtType*)IritMalloc(sizeof(IrtPtType) * NumCtrl);
    for (int i = 0; i < NumCtrl; ++i) {
        IrtRType ang = 2.0 * M_PI * i / NumCtrl;
        ctrl[i][0] = cx + radius * cos(ang);
        ctrl[i][1] = cy + radius * sin(ang);
        ctrl[i][2] = 0.0;
    }

    /* Iteratively shrink + smooth */
    IrtPtType* nextctrl = (IrtPtType*)IritMalloc(sizeof(IrtPtType) * NumCtrl);
    for (int it = 0; it < Iterations; ++it) {
        /* For each control point find closest boundary point and move toward it */
        for (int i = 0; i < NumCtrl; ++i) {
            IrtPtType cand;
            IrtRType d = FindClosestPointOnPolyList(UnionLocal, ctrl[i], cand);
            if (d >= IRIT_INFNTY / 2.0) {
                /* nothing found, keep place */
                IRIT_PT_COPY(nextctrl[i], ctrl[i]);
            }
            else {
                /* vector toward boundary */
                IrtRType vx = cand[0] - ctrl[i][0];
                IrtRType vy = cand[1] - ctrl[i][1];
                IrtRType len = sqrt(vx * vx + vy * vy);
                if (len <= IRIT_EPS) {
                    IRIT_PT_COPY(nextctrl[i], ctrl[i]);
                }
                else {
                    IrtRType step = ShrinkStep * d;
                    nextctrl[i][0] = ctrl[i][0] + (vx / len) * step;
                    nextctrl[i][1] = ctrl[i][1] + (vy / len) * step;
                    nextctrl[i][2] = 0.0;
                }
            }
        }

        /* Laplacian smoothing (preserve overall shrink) */
        for (int i = 0; i < NumCtrl; ++i) {
            int im1 = (i - 1 + NumCtrl) % NumCtrl;
            int ip1 = (i + 1) % NumCtrl;
            ctrl[i][0] = (nextctrl[im1][0] + nextctrl[ip1][0] + SmoothW * nextctrl[i][0]) / (2.0 + SmoothW);
            ctrl[i][1] = (nextctrl[im1][1] + nextctrl[ip1][1] + SmoothW * nextctrl[i][1]) / (2.0 + SmoothW);
            ctrl[i][2] = 0.0;
        }

        /* Enforce "stay outside" by snapping any inside control to closest boundary */
        for (int i = 0; i < NumCtrl; ++i) {
            if (PointInPolyList2D(UnionLocal, ctrl[i])) {
                IrtPtType snap;
                IrtRType d = FindClosestPointOnPolyList(UnionLocal, ctrl[i], snap);
                if (d < IRIT_INFNTY / 2.0) {
                    ctrl[i][0] = snap[0];
                    ctrl[i][1] = snap[1];
                }
            }
        }
    }

    IritFree(nextctrl);

    /* Transform final control polygon back to world coords using MatLocal (Local -> World) */
    IPPolygonStruct* FinalPl = IritPrsrAllocPolygon(0, NULL, NULL);
    IPVertexStruct* FirstV = NULL;
    IPVertexStruct* PrevV = NULL;
    for (int i = 0; i < NumCtrl; ++i) {
        IrtPtType worldP;
        IrtPtType localP;
        localP[0] = ctrl[i][0]; localP[1] = ctrl[i][1]; localP[2] = 0.0;
        IritMiscMatMultPtby4by4(worldP, localP, MatLocal);
        IPVertexStruct* NV = IritPrsrAllocVertex2(NULL);
        NV->Coord[0] = worldP[0];
        NV->Coord[1] = worldP[1];
        NV->Coord[2] = worldP[2];
        if (FirstV == NULL) FirstV = PrevV = NV;
        else { PrevV->Pnext = NV; PrevV = NV; }
    }
    if (PrevV != NULL) {
        PrevV->Pnext = FirstV;
        FinalPl->PVertex = FirstV;
        IritPrsrUpdatePolyPlane(FinalPl);
    }
    else {
        IritPrsrFreePolygonList(FinalPl);
        FinalPl = NULL;
    }

    IritFree(ctrl);
    IritPrsrFreeObject(UnionLocal);

    if (FinalPl == NULL) return NULL;

    IPObjectStruct* FinalObj = IritPrsrGenPOLYObject(FinalPl);
    if (FinalObj != NULL) IP_SET_POLYGON_OBJ(FinalObj);
    else IritPrsrFreePolygonList(FinalPl);

    return FinalObj;
}




IrtHmgnMatType* SelectBestViewSampling(IPObjectStruct* PObj,
    int NumSamples,
    IrtHmgnMatType* ResultMat)
{
    if (PObj == NULL || NumSamples <= 0)
        return NULL;

    NumSamples = NumSamples * 2; /* Generate more candidates than needed for better selection */

    IrtHmgnMatType* ViewMats = (IrtHmgnMatType*)IritMalloc(sizeof(IrtHmgnMatType) * NumSamples);
    if (ViewMats == NULL)
        return NULL;

    /* We'll include canonical axes first to make small sample counts predictable. */
    int next = 0;
    if (NumSamples >= 1) {
        IrtVecType CamPos, Center, Up;
        Center[0] = Center[1] = Center[2] = 0.0;

        /* +X */
        if (next < NumSamples) {
            CamPos[0] = 2.0; CamPos[1] = 0.0; CamPos[2] = 0.0;
            Up[0] = 0.0; Up[1] = 0.0; Up[2] = 1.0;
            GenLookAtMatrix(CamPos, Center, Up, ViewMats[next++]);
        }

        /* +Y */
        if (next < NumSamples) {
            CamPos[0] = 0.0; CamPos[1] = 2.0; CamPos[2] = 0.0;
            Up[0] = 0.0; Up[1] = 0.0; Up[2] = 1.0;
            GenLookAtMatrix(CamPos, Center, Up, ViewMats[next++]);
        }

        /* +Z (use Y as up to avoid pole singularity) */
        if (next < NumSamples) {
            CamPos[0] = 0.0; CamPos[1] = 0.0; CamPos[2] = 2.0;
            Up[0] = 0.0; Up[1] = 1.0; Up[2] = 0.0;
            GenLookAtMatrix(CamPos, Center, Up, ViewMats[next++]);
        }
    }

    /* 1. Generate remaining candidate view matrices (Fibonacci lattice). */
    IrtRType phi = (sqrt(5.0) - 1.0) / 2.0; /* Golden ratio */
    for (int s = 0; next < NumSamples && s < NumSamples * 4; ++s) {
        IrtRType y = 1.0 - (s / (IrtRType)(NumSamples - 1)) * 2.0; /* y in [1,-1] */
        IrtRType radius = sqrt(IRIT_MAX(0.0, 1.0 - y * y));
        IrtRType theta = 2.0 * M_PI * phi * s;
        IrtRType x = cos(theta) * radius;
        IrtRType z = sin(theta) * radius;

        /* Skip points near poles if they duplicate canonical axes already added.
           Use small epsilon to avoid exact duplicates. */
        const IrtRType eps = 1e-6;
        if ((fabs(x - 1.0) < eps && fabs(y) < eps && fabs(z) < eps) ||
            (fabs(y - 1.0) < eps && fabs(x) < eps && fabs(z) < eps) ||
            (fabs(z - 1.0) < eps && fabs(x) < eps && fabs(y) < eps)) {
            continue;
        }

        IrtVecType CamPos, Center, Up;
        CamPos[0] = x * 2.0; CamPos[1] = y * 2.0; CamPos[2] = z * 2.0;
        Center[0] = Center[1] = Center[2] = 0.0;
        Up[0] = 0.0; Up[1] = 0.0; Up[2] = 1.0;
        if (fabs(CamPos[0]) < 0.01 && fabs(CamPos[1]) < 0.01) {
            Up[0] = 1.0; Up[1] = 0.0; Up[2] = 0.0;
        }
        GenLookAtMatrix(CamPos, Center, Up, ViewMats[next++]);
    }

    /* 2. Score each candidate by constructing an approximated outer contour. */
    IrtRType* scores = (IrtRType*)IritMalloc(sizeof(IrtRType) * NumSamples);
    if (scores == NULL) {
        IritFree(ViewMats);
        return NULL;
    }

    /* Parameters for approximation - tune for performance/quality. */
    const int NumCtrl = 32;
    const int Iterations = 80;
    const IrtRType ShrinkStep = 0.25;
    const IrtRType SmoothW = 0.3;

    for (int i = 0; i < NumSamples; ++i) {
        scores[i] = 0.0;

        /* Reject near-vertical views (where wxy_len < 0.1) since the horizontal wire
           cannot reach the required vertical slope. */
        IrtRType wxy_len = sqrt(ViewMats[i][2][0] * ViewMats[i][2][0] + ViewMats[i][2][1] * ViewMats[i][2][1]);
        if (wxy_len < 0.1) {
            continue;
        }

        IPObjectStruct* Contour = IritPrsrApproxBSplineContourFromSolidView(PObj,
            ViewMats[i],
            NumCtrl,
            Iterations,
            ShrinkStep,
            SmoothW);
        if (Contour != NULL) {
            scores[i] = CalcPolygonArea(Contour);
            IritPrsrFreeObject(Contour);
        }
    }

    /* 3. Reorder ViewMats in-place by descending score (simple selection sort). */
    for (int i = 0; i < NumSamples - 1; ++i) {
        int best = i;
        for (int j = i + 1; j < NumSamples; ++j) {
            if (scores[j] > scores[best])
                best = j;
        }
        if (best != i) {
            IrtHmgnMatType tmpMat;
            IrtRType tmpScore = scores[i];

            memcpy(tmpMat, ViewMats[i], sizeof(IrtHmgnMatType));
            memcpy(ViewMats[i], ViewMats[best], sizeof(IrtHmgnMatType));
            memcpy(ViewMats[best], tmpMat, sizeof(IrtHmgnMatType));

            scores[i] = scores[best];
            scores[best] = tmpScore;
        }
    }

    if (ResultMat != NULL) {
        memcpy(ResultMat, ViewMats, sizeof(IrtHmgnMatType) * NumSamples);
    }

    IritFree(scores);
    return ViewMats;
}




/* -----------------------------------------------------------------------
 * IritPrsrHWCBuildSilhouetteRuledSrf
 *
 * Build a ruled surface suitable for 2D hot-wire silhouette cutting.
 *
 * The hot-wire machine traces the silhouette contour in machine XZ space
 * while keeping both wire endpoints (front/back clamps) at the same XZ
 * position, separated only in Y by FoamDepth. This requires a surface
 * whose V-isoparametric curves are the silhouette path:
 *
 *   S(u, vMin) = silhouette at Y = -FoamDepth/2   (front clamp path)
 *   S(u, vMax) = silhouette at Y = +FoamDepth/2   (back  clamp path)
 *
 * Steps:
 *   1. Project each contour vertex into view-local 2D (ignoring world Z).
 *   2. Compute AABB and scale uniformly to fit FoamWidth x FoamHeight
 *      with a 10% safety margin, centered at the machine origin.
 *   3. Build an E3 B-spline curve: X=local_u_scaled, Z=local_v_scaled+FoamH/2,
 *      Y=-FoamDepth/2 (front face relative to foam center).
 *   4. Extrude the curve along (0, FoamDepth, 0) to get S(u,v).
 *   5. Tag with IRIT_ATTR_ID_Dir = CAGD_CONST_U_DIR so IritPrsrHWCCreatePath
 *      samples in V direction (silhouette), giving both clamps the full
 *      silhouette path in sync.
 *
 * Parameters:
 *   Contour - polygon object from IritPrsrApproxBSplineContourFromSolidView.
 *   ViewMat - view matrix (row 2 = forward/view direction).
 *   Params  - HWC machine parameters; uses FoamWidth, FoamHeight, FoamDepth.
 *
 * Return: IPObjectStruct* (surface) on success, NULL on failure.
 *         Caller owns the returned object and must free it.
 * --------------------------------------------------------------------- */
IPObjectStruct* IritPrsrHWCBuildSilhouetteRuledSrf(
    const IPObjectStruct* Contour,
    const IrtHmgnMatType         ViewMat,
    const IritPrsrHWCDataStruct* Params)
{
    int k, n, guard;
    IrtRType minx, miny, maxx, maxy, scale_x, scale_z, scale, cx, cy;
    IrtVecType u_vec, v_vec, w_vec;
    IrtHmgnMatType MatLocal;
    IrtPtType* pts2d;
    IPVertexStruct* V, * cur;
    CagdCrvStruct* Crv;
    CagdVecStruct YDir;
    CagdSrfStruct* Srf;
    IPObjectStruct* SrfObj;

    if (Contour == NULL || !IP_IS_POLY_OBJ((IPObjectStruct*)Contour) ||
        Contour->U.Pl == NULL || Params == NULL)
        return NULL;

    /* 1. Build view-local orthonormal basis (u_vec=right, v_vec=up, w_vec=fwd). */
    if (!IritPrsrBuildViewBasisFromMat(ViewMat, u_vec, v_vec, w_vec, 1))
        return NULL;

    /* MatLocal maps world -> view-local (rows = basis vectors). */
    IritMiscMatGenUnitMat(MatLocal);
    MatLocal[0][0] = u_vec[0]; MatLocal[0][1] = u_vec[1];
    MatLocal[0][2] = u_vec[2]; MatLocal[0][3] = 0.0;
    MatLocal[1][0] = v_vec[0]; MatLocal[1][1] = v_vec[1];
    MatLocal[1][2] = v_vec[2]; MatLocal[1][3] = 0.0;
    MatLocal[2][0] = w_vec[0]; MatLocal[2][1] = w_vec[1];
    MatLocal[2][2] = w_vec[2]; MatLocal[2][3] = 0.0;
    MatLocal[3][0] = 0.0;  MatLocal[3][1] = 0.0;
    MatLocal[3][2] = 0.0;  MatLocal[3][3] = 1.0;

    /* 2. Count and validate polygon vertices. */
    V = Contour->U.Pl->PVertex;
    if (V == NULL) return NULL;

    n = 0;
    cur = V;
    guard = 0;
    do {
        ++n;
        cur = cur->Pnext;
        if (++guard > 200000) { n = 0; break; }
    } while (cur != NULL && cur != V);

    if (n < 3) return NULL;

    /* 3. Project each vertex to view-local 2D, ignoring world Z.
     *   localP[0] = view-right component -> machine X
     *   localP[1] = view-up   component  -> machine Z (vertical) */
    pts2d = (IrtPtType*)IritMalloc(sizeof(IrtPtType) * n);
    cur = V;
    guard = 0;
    for (k = 0; k < n; ++k) {
        /* Compute true dot products to project world coordinates onto the view plane. */
        pts2d[k][0] = cur->Coord[0] * u_vec[0] + cur->Coord[1] * u_vec[1] + cur->Coord[2] * u_vec[2];
        pts2d[k][1] = cur->Coord[0] * v_vec[0] + cur->Coord[1] * v_vec[1] + cur->Coord[2] * v_vec[2];
        pts2d[k][2] = 0.0;
        cur = cur->Pnext;
        if (cur == NULL || cur == V) break;
        if (++guard > 200000) break;
    }

    /* 4. Compute AABB. */
    minx = maxx = pts2d[0][0];
    miny = maxy = pts2d[0][1];
    for (k = 1; k < n; ++k) {
        if (pts2d[k][0] < minx) minx = pts2d[k][0];
        if (pts2d[k][0] > maxx) maxx = pts2d[k][0];
        if (pts2d[k][1] < miny) miny = pts2d[k][1];
        if (pts2d[k][1] > maxy) maxy = pts2d[k][1];
    }

    if ((maxx - minx) < IRIT_EPS || (maxy - miny) < IRIT_EPS) {
        IritFree(pts2d);
        return NULL;
    }

    /* 5. Uniform scale: fit silhouette into 90% of FoamWidth x FoamHeight. */
    scale_x = (Params->FoamWidth * 0.9) / (maxx - minx);
    scale_z = (Params->FoamHeight * 0.9) / (maxy - miny);
    scale = IRIT_MIN(scale_x, scale_z);
    cx = 0.5 * (minx + maxx);
    cy = 0.5 * (miny + maxy);

    /* 6. Validate the view is horizontally cuttable.
     *
     *  The hot-wire machine spans the wire along the horizontal Y axis between two clamps.
     *  Therefore it can only cleanly sweep along views where the horizontal length
     *  of the view vector is non-negligible. We reject completely vertical views. */
    IrtRType wxy_len = sqrt(w_vec[0] * w_vec[0] + w_vec[1] * w_vec[1]);
    if (wxy_len < 0.1) {
        IritFree(pts2d);
        return NULL;
    }

    /* 7. Build silhouette curve mapped exactly into 3D World space.
     *    The 2D contour is originally in the (u_vec, v_vec) plane of the view.
     *    We rebuild the exact 3D orientation.
     *
     *    The extruded surface must be aligned with w_vec in 3D.
     *    The start of the extrusion will be -w_vec * FoamDepth/2.
     *    We also shift the entire cut upwards to be centered in the foam block. */
    Crv = IritCagdBspCrvNew(n, 2, CAGD_PT_E3_TYPE);
    if (Crv == NULL) {
        IritFree(pts2d);
        return NULL;
    }

    for (k = 0; k < n; ++k) {
        IrtRType sx = (pts2d[k][0] - cx) * scale;
        IrtRType sy = (pts2d[k][1] - cy) * scale;

        /* Map back to 3D world coordinates. */
        IrtRType wx = sx * u_vec[0] + sy * v_vec[0];
        IrtRType wy = sx * u_vec[1] + sy * v_vec[1];
        IrtRType wz = sx * u_vec[2] + sy * v_vec[2];

        /* To bypass a critical math bug in the legacy hot_wire_cut.c logic,
         * we must pre-extrapolate the Z-coordinate. The legacy code computes:
         *   a_buggy = Crv2_Z + slope * Crv2_Y_rot
         *   z_buggy = a_buggy - slope * MACHINE_MAX_Y (where MACHINE_MAX_Y is 390.0)
         * But the correct math for extrapolating to Y = +/- 195.0 is:
         *   a_target = Crv2_Z + slope * (195.0 - Crv2_Y_rot)
         * By setting a_buggy = a_target, we can algebraically derive the exact Z shift
         * needed per-point to cancel the bug out completely.
         *   z_shift = slope * (195.0 - 2 * Crv2_Y_rot)
         * After substituting Crv2_Y_rot = -wz * slope + wxy_len * FoamDepth / 2, we get:
         */
        IrtRType slope = w_vec[2] / wxy_len;
        IrtRType z_shift = slope * (195.0 - wxy_len * Params->FoamDepth) + 2.0 * wz * slope * slope;

        Crv->Points[1][k] = wx - w_vec[0] * Params->FoamDepth * 0.5;
        Crv->Points[2][k] = wy - w_vec[1] * Params->FoamDepth * 0.5;
        Crv->Points[3][k] = wz - w_vec[2] * Params->FoamDepth * 0.5 + Params->FoamHeight * 0.5 + z_shift;
    }
    IritCagdBspKnotUniformOpen(n, 2, Crv->KnotVector);

    IritFree(pts2d);

    /* 8. Extrude exactly along the true 3D view direction by FoamDepth.
     *    This creates a cylinder (surface) whose rulings are exactly w_vec.
     *    When IritPrsrHWCProcessRulingPair processes this surface, it computes
     *    Angle = atan2(Dx, Dy) = atan2(w_vec[0], w_vec[1]), and rotates the
     *    clamps accordingly to perfectly align the Y-axis of the machine with w_vec. */
    YDir.Vec[0] = w_vec[0] * Params->FoamDepth;
    YDir.Vec[1] = w_vec[1] * Params->FoamDepth;
    YDir.Vec[2] = w_vec[2] * Params->FoamDepth;

    Srf = IritCagdExtrudeSrf(Crv, &YDir);
    IritCagdCrvFree(Crv);
    if (Srf == NULL) return NULL;

    /* 8. Wrap into IPObjectStruct and tag ruling direction as U.
     *   With CAGD_CONST_U_DIR the HWC generator sets:
     *     RulingDirection   = U  (front-to-back)
     *     SamplingDirection = V  (along silhouette)
     *   => Crv1 = S(u,vMin) and Crv2 = S(u,vMax), so both clamps trace
     *      the full silhouette contour in sync. */
    SrfObj = IritPrsrGenSrfObject("sil_ruled_srf", Srf, NULL);
    if (SrfObj != NULL) {
        IritMiscAttrIDSetObjectIntAttrib(SrfObj, IRIT_ATTR_ID_Dir,
            CAGD_CONST_U_DIR);
    }
    else {
        IritCagdSrfFree(Srf);
    }

    return SrfObj;
}
