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


/**
 * @brief Build a stable orthonormal basis (u,v,w) from a view matrix.
 *
 * Short intro:
 *   Extracts the view 'forward' vector from the given 4x4 matrix (row 2),
 *   normalizes it as `w`, picks an arbitrary vector not parallel to `w` and
 *   produces `u = normalize(cross(arb, w))` and `v = cross(w, u)`.
 *
 * Parameters:
 *   Mat - input 4x4 view matrix (rotation stored in rows).
 *   u, v, w - output 3-vectors (must point to 3-element arrays).
 *   AllowDefault - if non-zero, fall back to canonical directions on degenerate inputs.
 *
 * Returns:
 *   1 on success, 0 on failure (unless AllowDefault allowed a fallback).
 *
 * Notes:
 *   - Function is static (internal). Handles degenerate forwards when allowed.
 */
static int IritPrsrBuildViewBasisFromMat(const IrtHmgnMatType Mat,
    IrtVecType u,
    IrtVecType v,
    IrtVecType w,
    int AllowDefault)
{
    IrtVecType arb;
    IrtRType wlen, ul;

    /* Extract forward vector from third *row* of matrix (Mat[2][0..2]). */
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


/* ----------------- existing functions (unchanged except for basis usage)  */

/**
 * @brief Create a polygon from 2D points (in a local plane) and extrude into a 3D solid.
 *
 * Short intro:
 *   Places the 2D polygon into the plane orthogonal to `Dir`, constructs a closed
 *   polygon, and extrudes it along the original `Dir` vector (magnitude controls extrusion length).
 *
 * Parameters:
 *   pts - array of 2D points (IrtE2PtStruct) in XY local coordinates.
 *   n   - number of points (must be >= 3).
 *   Dir - 3D extrusion vector (magnitude is extrusion length).
 *
 * Returns:
 *   New IPObjectStruct* (solid) on success or NULL on failure.
 *
 * Notes:
 *   - Caller owns returned object and must free it.
 *   - Returns NULL on degenerate Dir or invalid input.
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


/**
 * @brief Extrude projected silhouette polygons into solids for each view.
 *
 * Short intro:
 *   For each polygon in the supplied projected lists (preprocessed and direct),
 *   compute extrusion direction using the view forward vector scaled by Depth,
 *   extrude and return a LIST of placed solids.
 *
 * Parameters:
 *   ProjListPre    - LIST of polys from preprocess (may be NULL)
 *   ProjListDirect - LIST of polys extracted directly (may be NULL)
 *   Views          - array of view matrices (NumViews entries)
 *   NumViews       - number of view matrices
 *   Depth          - extrusion scalar (applied along view forward)
 *
 * Returns:
 *   LIST IPObjectStruct* of extruded solids (caller owns) or NULL on error.
 *
 * Notes:
 *   - Each returned solid has ObjName like "sil_dir_view%d_poly%03d_extr".
 *   - Polygons with invalid ObjName/view index are skipped.
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



/**
 * @brief Generate a LookAt view matrix (camera at Eye looking at Center).
 *
 * Short intro:
 *   Builds an IRIT-style view matrix that places the camera at `Eye`,
 *   looking at `Center` with `Up` as a hint. Matrix rotation stored in rows.
 *
 * Parameters:
 *   Eye, Center, Up - input 3-vectors.
 *   Mat - output 4x4 matrix (written).
 *
 * Notes:
 *   - Handles degenerate Up/F direction by choosing a fallback side vector.
 *   - IRIT convention: positive Z is the view direction row (row 2).
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


/**
 * @brief Compute a heuristic projected area for a polygon object.
 *
 * Short intro:
 *   Sums IRIT's per-loop area estimate across polygon loops. Used as a
 *   heuristic score for view ranking (larger projected contour area preferred).
 *
 * Parameters:
 *   PObj - polygon object (IP_OBJ_POLY) to measure.
 *
 * Returns:
 *   Floating area measure (0.0 for NULL or non-polygon).
 *
 * Notes:
 *   - Not a strict projection-based area computation; a fast heuristic.
 */
IrtRType CalcPolygonArea(IPObjectStruct* PObj)
{
    IrtRType area = 0.0;
    if (PObj == NULL || !IP_IS_POLY_OBJ(PObj)) return 0.0;

    for (IPPolygonStruct* Pl = PObj->U.Pl; Pl != NULL; Pl = Pl->Pnext) {
        area += IritGeomPolyOnePolyArea(Pl, TRUE); // Uses IRIT's built-in area function
    }
    return area;
}


/**
 * @brief Distance from point p to segment [a,b] in 2D and closest point.
 *
 * Short intro:
 *   Projects p onto segment ab (clamped) and returns Euclidean distance;
 *   writes the closest point in 'out'.
 *
 * Parameters:
 *   p - query 2D point.
 *   a,b - segment endpoints (2D stored in index 0,1).
 *   out - output closest 2D point.
 *
 * Returns:
 *   Distance from p to the closest point on segment ab.
 */
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

/**
 * @brief Find closest point on any edge of a LIST of polygon objects.
 *
 * Short intro:
 *   Iterates the provided LIST of polygon objects (view-local XY) and computes
 *   the closest point on any polygon edge to the query point p.
 *
 * Parameters:
 *   UnionObj - LIST object containing POLY objects.
 *   p - query 2D point.
 *   out - closest point (written on success).
 *
 * Returns:
 *   Smallest distance found, or IRIT_INFNTY on failure/no polygons.
 */
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

/**
 * @brief Ray-casting point-in-polygon test for a single polygon loop.
 *
 * Short intro:
 *   Determines whether the 2D point p is inside the given polygon loop using
 *   a simple crossing test.
 *
 * Parameters:
 *   p - query 2D point.
 *   Pl - polygon loop (uses Coord[0], Coord[1]).
 *
 * Returns:
 *   1 if inside, 0 otherwise.
 */
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

/**
 * @brief Test whether a point is inside any polygon in a LIST (2D).
 *
 * Short intro:
 *   Iterates all POLY objects in the supplied LIST and returns true when the
 *   point lies inside any loop (using PointInPolyLoop2D).
 *
 * Parameters:
 *   UnionObj - LIST of POLY objects.
 *   p - query 2D point.
 *
 * Returns:
 *   1 if inside any polygon loop, 0 otherwise.
 */
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

/**
 * @brief Comparison function for qsort to sort 2D points by (x,y).
 *
 * Parameters:
 *   a, b - pointers to IrtPtType.
 *
 * Returns:
 *   -1, 0, 1 for ordering by x then y.
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

/**
 * @brief 2D cross product (signed area) used by convex hull routines.
 *
 * Parameters:
 *   a, b, c - 2D points (IrtPtType) where cross = (b-a)x(c-a) in XY.
 *
 * Returns:
 *   Signed cross value (>0 left turn, <0 right turn).
 */
static IrtRType Cross2D(const IrtPtType a, const IrtPtType b, const IrtPtType c)
{
    /* cross of AB x AC in XY */
    return (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]);
}

/**
 * @brief Project a polygonal Solid into a view-local frame and build a conservative outer contour.
 *
 * Short intro:
 *   Projects all polygon vertices of `Solid` using `MatLocal` (rotation-only)
 *   collects unique projected points, computes a monotone-chain convex hull and
 *   returns a LIST containing a single POLY (the hull) in view-local XY (Z=0).
 *
 * Parameters:
 *   Solid - polygonal object (IP_OBJ_POLY).
 *   MatLocal - world -> view-local rotation matrix (zero translation).
 *
 * Returns:
 *   LIST IPObjectStruct* containing one POLY (outer hull) on success, NULL on failure.
 *
 * Notes:
 *   - Conservative approach: produces outer convex hull rather than exact union.
 *   - Caller owns returned LIST and must free it.
 */
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


/**
 * @brief Approximate the outer contour of a Solid as seen from a view using a shrinking control circle.
 *
 * Short intro:
 *   Projects the solid into a view-local frame, builds a conservative union hull,
 *   then iteratively shrinks and smooths a circular control polygon toward the hull
 *   to produce a smooth outer contour. The final contour is transformed back to world coords.
 *
 * Parameters:
 *   Solid - input polygonal solid (IP_OBJ_POLY).
 *   ViewMat - view matrix defining projection.
 *   NumCtrl - number of control points on initial circle (>= 6 recommended).
 *   Iterations - number of shrink iterations.
 *   ShrinkStep - fraction of distance toward boundary to move per iteration.
 *   SmoothW - Laplacian smoothing weight.
 *
 * Returns:
 *   POLY IPObjectStruct* in world coordinates approximating the outer contour (caller owns), or NULL.
 *
 * Notes:
 *   - Ensures "stay outside" by snapping inside control points to the boundary.
 *   - Uses BuildProjectedUnionLocal to avoid fragile 2D boolean operations.
 */
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

    /* Build union of projected faces in local XY */
    IPObjectStruct* UnionLocal = BuildProjectedUnionLocal(Solid, MatLocal);
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

    /* Transform final control polygon back to world coords using InvMat, build closed polygon object */
    IPPolygonStruct* FinalPl = IritPrsrAllocPolygon(0, NULL, NULL);
    IPVertexStruct* FirstV = NULL;
    IPVertexStruct* PrevV = NULL;
    for (int i = 0; i < NumCtrl; ++i) {
        IrtPtType worldP;
        IrtPtType localP;
        localP[0] = ctrl[i][0]; localP[1] = ctrl[i][1]; localP[2] = 0.0;
        IritMiscMatMultPtby4by4(worldP, localP, InvMat);
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



/**
 * @brief Select best view sampling for a polygonal object by scoring candidate views.
 *
 * Short intro:
 *   Generates candidate camera positions (canonical axes + Fibonacci lattice),
 *   approximates outer contours for each candidate, scores them by projected area,
 *   sorts candidates by descending score and returns an allocated array of view matrices.
 *
 * Parameters:
 *   PObj - polygonal object to sample.
 *   NumSamples - desired number of samples (function internally doubles this to generate more candidates).
 *   ResultMat - optional preallocated buffer to copy resulting matrices into (may be NULL).
 *
 * Returns:
 *   Newly allocated array of IrtHmgnMatType (size = NumSamples*2 internally) ordered by descending score.
 *   Caller must free with IritFree().
 *
 * Notes:
 *   - Expensive: each candidate runs a contour approximation.
 *   - The function multiplies NumSamples by two internally to improve selection.
 */
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



