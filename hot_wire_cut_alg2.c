/*****************************************************************************
* Generates 5-Axis Hot-Wire Cutter tool paths using Visual Hull Extrusion.   *
******************************************************************************
* (C) Gershon Elber, Technion, Israel Institute of Technology                *
******************************************************************************
* Written by:  Shiraz Shmulman and Dani Rifkind       Ver 1.0, May 2026      *
******************************************************************************
* ALGORITHM OVERVIEW:                                                        *
*                                                                            *
* This module implements a novel 5-axis hot-wire cutting algorithm based on  *
* silhouette volume intersection. It evaluates optimal viewing directions    *
* for a given 3D object (in xy plane), extracts the 2D projected silhouette  *
* from each view, and extrudes that silhouette into a 3D ruled surface       *
* aligned with the true 3D viewing direction. By sequencing these views      *
* and rotating the foam on the machine's rotary axis (B), the machine cuts   *
* out the visual hull of the target 3D object.                               *
*                                                                            *
* CONVENTIONS:                                                               *
*                                                                            *
* 1. All functions have the prefix of HWC                                    *
*                                                                            *
* MACHINE COORDINATE SYSTEM & INFORMATION:                                   *
*                                                                            *
* - Origin (0,0,0): Back-Right-Bottom corner of the machine volume.          *
*                                                                            *
* - Axes interpretation:                                                     *
*   X : Left clamp horizontal position.                                      *
*   Y : Right clamp horizontal position.                                     *
*   Z : Left clamp vertical position.                                        *
*   A : Right clamp vertical position.                                       *
*   B : Foam rotation angle (rotary table) in degrees.                       *
*                                                                            *
* - Physical dimensions (in mm):                                             *
*   MACHINE_MAX_X = 390.0                                                    *
*   MACHINE_MAX_Y = 390.0                                                    *
*   MACHINE_MAX_Z = 340.0                                                    *
*   Default Foam Size = 135.0 (W) x 135.0 (D) x 135.0 (H) mm                 *
*                                                                            *
* NOTE: All coordinates and dimensions generated are in millimeters (mm).    *
*****************************************************************************/

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
#define HWCSIL_GRID_RES 1024

typedef struct IritPrsrHWCEdgeStruct {
    IrtPtType Pt1;
    IrtPtType Pt2;
} IritPrsrHWCEdgeStruct;

static int HWCIritPrsrBuildViewBasisFromMat(const IrtHmgnMatType Mat,
                                            IrtVecType u,
                                            IrtVecType v,
                                            IrtVecType w,
                                            int AllowDefault);

static void HWCGenLookAtMatrix(IrtVecType Eye,
                               IrtVecType Center,
                               IrtVecType Up,
                               IrtHmgnMatType Mat);

static IrtRType HWCCalcPolygonArea(IPObjectStruct *PObj);

static IPObjectStruct *HWCBuildProjectUnionLocal(IPObjectStruct *Solid,
                                                 const IrtHmgnMatType MatLocal);

static IPObjectStruct *HWCIritPrsrApproxBSplineContourFromSolidView(
    IPObjectStruct *Solid,
    const IrtHmgnMatType ViewMat,
    int NumCtrl);

static IrtHmgnMatType *HWCSelectBestViewSampling(IPObjectStruct *PObj,
                                                 int NumSamples,
                                                 IrtHmgnMatType *ResultMat);

static IrtHmgnMatType *HWCSelectViewsSampling(int NumSamples);

static IritPrsrHWCEdgeStruct *HWCFindBoundaryEdges(const IPObjectStruct *Solid,
                                                   int *OutNumEdges);

static IrtRType HWCDistPointSegment2D(IrtPtType P,
                                      IrtPtType A,
                                      IrtPtType B);

static IPObjectStruct *HWCBuildSilhouetteRuledSrf(const IPObjectStruct *Contour,
                                                  const IPObjectStruct *Solid,
                                                  const IrtHmgnMatType ViewMat,
                                                  const IritPrsrHWCDataStruct *Params);

static void HWCCombineGCodeFiles(const char const **GcodeFiles,
                                 int NumFiles,
                                 const char *OutPath,
                                 const IritPrsrHWCDataStruct *Params);

static void HWCConvertPolylinesToPolygons(IPObjectStruct *PObj);

IrtRType HWCGenerateGCodeFromObj(IPObjectStruct *RawModel,
                                 const char *OutputGCodePath,
                                 IrtRType NumViews,
                                 IrtRType OutputITDType,
                                 IrtRType ITDLength,
                                 const char *fileName,
                                 int uniformCut);

/*****************************************************************************
* AUXILIARY:                                                                 *
*                                                                            *
* Auxiliary function to extract orthonormal basis from view matrix.          *
*****************************************************************************/
static int HWCIritPrsrBuildViewBasisFromMat(const IrtHmgnMatType Mat,
                                            IrtVecType u,
                                            IrtVecType v,
                                            IrtVecType w,
                                            int AllowDefault)
{
    IrtVecType arb;
    IrtRType wlen, ul;

    /* Extract forward vector from third *row* of matrix (Mat[2][0..2]). */
    /* The view/LookAt code in this file stores rotation in rows, so the */
    /* forward vector is on row index 2.                                 */
    w[0] = Mat[2][0];
    w[1] = Mat[2][1];
    w[2] = Mat[2][2];

    wlen = sqrt(IRIT_SQR(w[0]) + IRIT_SQR(w[1]) + IRIT_SQR(w[2]));
    if (wlen <= IRIT_EPS) {
        if (!AllowDefault)
            return 0;
        /* Fallback to canonical forward if allowed. */
        w[0] = 0.0;
        w[1] = 0.0;
        w[2] = 1.0;
    }
    else {
        w[0] /= wlen;
        w[1] /= wlen;
        w[2] /= wlen;
    }

    /* Choose stable arbitrary vector not parallel to w. */
    if (fabs(w[0]) < 0.9) {
        arb[0] = 1.0;
        arb[1] = 0.0;
        arb[2] = 0.0;
    }
    else {
        arb[0] = 0.0;
        arb[1] = 1.0;
        arb[2] = 0.0;
    }

    /* u = normalize(cross(arb, w)). */
    u[0] = arb[1] * w[2] - arb[2] * w[1];
    u[1] = arb[2] * w[0] - arb[0] * w[2];
    u[2] = arb[0] * w[1] - arb[1] * w[0];

    ul = sqrt(IRIT_SQR(u[0]) + IRIT_SQR(u[1]) + IRIT_SQR(u[2]));
    if (ul <= IRIT_EPS) {
        /* Extremely unlikely; if allowed try a different arb, else fail. */
        if (!AllowDefault)
            return 0;
        /* Try fallback arb. */
        arb[0] = 0.0;
        arb[1] = 1.0;
        arb[2] = 0.0;
        u[0] = arb[1] * w[2] - arb[2] * w[1];
        u[1] = arb[2] * w[0] - arb[0] * w[2];
        u[2] = arb[0] * w[1] - arb[1] * w[0];
        ul = sqrt(IRIT_SQR(u[0]) + IRIT_SQR(u[1]) + IRIT_SQR(u[2]));
        if (ul <= IRIT_EPS)
            return 0;
    }

    u[0] /= ul;
    u[1] /= ul;
    u[2] /= ul;

    /* v = cross(w, u). */
    v[0] = w[1] * u[2] - w[2] * u[1];
    v[1] = w[2] * u[0] - w[0] * u[2];
    v[2] = w[0] * u[1] - w[1] * u[0];

    return 1;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
* Generate a LookAt matrix from eye position, center, and up vector.         *
* Constructs a view matrix mapping the world so the camera is at origin,     *
* looking down the positive Z axis (IRIT standard).                          *
*                                                                            *
* PARAMETERS:                                                                *
* Eye, Center, Up: Vectors defining the camera pose.                         *
* Mat:             4x4 homogeneous transformation matrix.                    *
*                                                                            *
* RETURN VALUE:                                                              *
* Void.                                                                      *
*                                                                            *
* KEYWORDS:                                                                  *
* LookAt, view matrix, camera.                                               *
*****************************************************************************/
static void HWCGenLookAtMatrix(IrtVecType Eye,
                               IrtVecType Center,
                               IrtVecType Up,
                               IrtHmgnMatType Mat)
{
    IrtVecType F, S, U;

    /* Calculate Forward Vector (F). */
    IRIT_VEC_SUB(F, Center, Eye);
    IRIT_VEC_NORMALIZE(F);

    /* Calculate Side Vector (S) = F x Up. */
    IRIT_CROSS_PROD(S, F, Up);
    /* Check for degenerate case where F is parallel to Up. */
    if ((IRIT_DOT_PROD(S, S)) < 1e-6) {
        /* Fallback: If looking straight up/down, choose X-axis as side. */
        S[0] = 1.0;
        S[1] = 0.0;
        S[2] = 0.0;
    }
    IRIT_VEC_NORMALIZE(S);

    /* Calculate True Up Vector (U) = S x F. */
    IRIT_CROSS_PROD(U, S, F);
    IRIT_VEC_NORMALIZE(U);

    /* Construct Matrix. Row 0: S, Row 1: U, Row 2: F. */
    IritMiscMatGenUnitMat(Mat);

    /* Rotation Part. */
    Mat[0][0] = S[0];
    Mat[0][1] = S[1];
    Mat[0][2] = S[2];
    Mat[1][0] = U[0];
    Mat[1][1] = U[1];
    Mat[1][2] = U[2];
    /* IRIT uses +Z as view direction. */
    Mat[2][0] = F[0];
    Mat[2][1] = F[1];
    Mat[2][2] = F[2];

    /* Translation Part (-Eye dot Basis). */
    Mat[0][3] = -IRIT_DOT_PROD(S, Eye);
    Mat[1][3] = -IRIT_DOT_PROD(U, Eye);
    Mat[2][3] = -IRIT_DOT_PROD(F, Eye);
}

/*****************************************************************************
* DESCRIPTION:                                                               *
* Calculate projected 2D area of a polygon chain.                            *
* Used as a heuristic: maximize projected area ~ maximize feature definition.*
*                                                                            *
* PARAMETERS:                                                                *
* PObj: Polygon object.                                                      *
*                                                                            *
* RETURN VALUE:                                                              *
* IrtRType: Estimated polygon area.                                          *
*****************************************************************************/
static IrtRType HWCCalcPolygonArea(IPObjectStruct *PObj)
{
    IrtRType
        area = 0.0;
    IPPolygonStruct *Pl;

    if (PObj == NULL || !IP_IS_POLY_OBJ(PObj))
        return 0.0;

    for (Pl = PObj -> U.Pl; Pl != NULL; Pl = Pl -> Pnext) {
        /* Uses IRIT's built-in area function. */
        area += IritGeomPolyOnePolyArea(Pl, TRUE);
    }
    return area;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
* Build a rasterized projection of polygonal solid into view-local 2D frame  *
* and extract boundary contour using Moore neighborhood tracing. Robust      *
* method for computing silhouette by scanline rasterization and contour      *
* extraction.                                                                *
*                                                                            *
* Projects all polygons of input solid into view-local coordinate system,    *
* rasterizes them into a 2D binary grid using scanline fill algorithm,       *
* then traces the filled region boundary using Moore neighborhood algorithm. *
* Returns contour as a high-resolution 2D polygon in view-local XY plane.    *
*                                                                            *
* PARAMETERS:                                                                *
* Solid: IPObjectStruct pointer to valid polygonal object (must pass         *
*        IP_IS_POLY_OBJ check). Multiple loops allowed. Input not            *
*        modified. NULL returns NULL.                                        *
* MatLocal: 4x4 transformation matrix mapping world coordinates to           *
*           view-local frame. Typically built from orthonormal basis         *
*           (u=right, v=up, w=forward). Singular matrices return NULL.       *
*                                                                            *
* RETURN VALUE:                                                              *
* IPObjectStruct *: LIST object containing single POLY with boundary contour *
*                   in view-local XY coordinates (Z=0), or NULL on failure   *
*                   (NULL input, degenerate projection, allocation failure,  *
*                   or empty boundary). Caller owns returned object and must *
*                   free via IritPrsrFreeObject().                           *
*                                                                            *
* ALGORITHM:                                                                 *
* 1. Project all polygon vertices to local XY, compute bounding box with 5%  *
*    margin to prevent boundary clipping.                                    *
* 2. Allocate SIL_GRID_RES x SIL_GRID_RES (1024x1024) binary raster grid.    *
* 3. For each input polygon, execute scanline fill:                          *
*    - Map polygon vertices to grid coordinates.                             *
*    - For each scanline intersecting polygon Y-range, compute               *
       X-intersections.                                                      *
*    - Fill horizontal spans between paired intersections (parity rule).     *
* 4. Find first filled pixel (leftmost, then topmost).                       *
* 5. Moore neighborhood contour tracing: starting from filled pixel, follow  *
*    filled-to-empty boundary using 8-neighbor search, always turning 135°   *
*    left to maintain consistent winding. Traces until returning to start.   *
* 6. Convert traced pixel path back to view-local floating-point coordinates *
*    using inverse grid-to-world mapping.                                    *
* 7. Create POLY with traced vertices, wrap in LIST, return.                 *
*                                                                            *
* PERFORMANCE:                                                               *
* - O(V*max_scanlines) for rasterization (V=vertex count).                   *
* - O(perimeter) for contour tracing (proportional to boundary length).      *
* - Memory: 1MB grid + temporary arrays. Suitable for interactive use.       *
*                                                                            *
* SEE ALSO:                                                                  *
* HWCIritPrsrApproxBSplineContourFromSolidView (higher-level wrapper),       *
* HWCIritPrsrBuildViewBasisFromMat (constructs MatLocal).                    *
*****************************************************************************/
static IPObjectStruct *HWCBuildProjectUnionLocal(IPObjectStruct *Solid,
    const IrtHmgnMatType MatLocal)
{
    IrtRType
        minx = IRIT_INFNTY,
        miny = IRIT_INFNTY,
        maxx = -IRIT_INFNTY,
        maxy = -IRIT_INFNTY, 
        dx, dy, margin;
    unsigned char *grid;
    int *px, *py, *ints, x, y, i, j,
        startX = -1,
        startY = -1;
    IPPolygonStruct *Pl;

    if (Solid == NULL || !IP_IS_POLY_OBJ(Solid))
        return NULL;

    /* 1. Project points to local XY, collect them to compute bounding box. */
    for (Pl = Solid->U.Pl; Pl != NULL; Pl = Pl->Pnext) {
        IPVertexStruct *cur,
            * V = Pl->PVertex;
        int guard = 0;

        if (V == NULL)
            continue;
        cur = V;
        do {
            IrtPtType localP;

            IritMiscMatMultPtby4by4(localP, cur->Coord, MatLocal);
            if (localP[0] < minx)
                minx = localP[0];
            if (localP[1] < miny)
                miny = localP[1];
            if (localP[0] > maxx)
                maxx = localP[0];
            if (localP[1] > maxy)
                maxy = localP[1];

            cur = cur->Pnext;
            if (++guard > 200000)
                break;
        } while (cur != NULL && cur != V);
    }

    if (minx > maxx)
        return NULL;

    /* Add a small 5% margin to prevent tracing out of bounds. */
    dx = maxx - minx;
    dy = maxy - miny;
    margin = IRIT_MAX(dx, dy) * 0.05;
    if (margin < 1e-5)
        margin = 1e-5;
    minx -= margin;
    maxx += margin;
    miny -= margin;
    maxy += margin;
    dx = maxx - minx;
    dy = maxy - miny;

    /* 2. Allocate and clear 2D grid. */
    grid = (unsigned char*) IritMalloc(HWCSIL_GRID_RES * HWCSIL_GRID_RES);
    memset(grid, 0, HWCSIL_GRID_RES * HWCSIL_GRID_RES);

    /* Assume max polygon vertex count <= 4096. */
    px = (int*) IritMalloc(sizeof(int) * 4096);
    py = (int*) IritMalloc(sizeof(int) * 4096);
    ints = (int*) IritMalloc(sizeof(int) * 4096);

    /* 3. Scanline fill each projected polygon. */
    for (Pl = Solid -> U.Pl; Pl != NULL; Pl = Pl -> Pnext) {
        IPVertexStruct *cur,
            * V = Pl -> PVertex;
        int pminy, pmaxy,
            npts = 0,
            guard = 0;

        if (V == NULL)
            continue;

        cur = V;
        do {
            IrtPtType localP;
            int px_val, py_val;

            IritMiscMatMultPtby4by4(localP, cur -> Coord, MatLocal);

            px_val = (int)((localP[0] - minx) / dx * (HWCSIL_GRID_RES - 1));
            py_val = (int)((localP[1] - miny) / dy * (HWCSIL_GRID_RES - 1));
            if (px_val < 0)
                px_val = 0;
            if (px_val >= HWCSIL_GRID_RES)
                px_val = HWCSIL_GRID_RES - 1;
            if (py_val < 0)
                py_val = 0;
            if (py_val >= HWCSIL_GRID_RES)
                py_val = HWCSIL_GRID_RES - 1;

            if (npts < 4096) {
                px[npts] = px_val;
                py[npts] = py_val;
                npts++;
            }
            cur = cur -> Pnext;
            if (++guard > 200000)
                break;
        } while (cur != NULL && cur != V);

        if (npts < 3)
            continue;

        /* Scanline logic */
        pminy = HWCSIL_GRID_RES;
        pmaxy = -1;
        for (i = 0; i < npts; ++i) {
            if (py[i] < pminy)
                pminy = py[i];
            if (py[i] > pmaxy)
                pmaxy = py[i];
        }

        for (y = pminy; y <= pmaxy; ++y) {
            int a, b,
                count = 0;

            for (i = 0; i < npts; ++i) {
                int j_idx = (i + 1) % npts,
                    y1 = py[i],
                    y2 = py[j_idx],
                    x1 = px[i],
                    x2 = px[j_idx];

                if (y1 == y2)
                    continue;

                if ((y1 <= y && y < y2) || (y2 <= y && y < y1)) {
                    int x_int = x1 + (y - y1) * (x2 - x1) / (y2 - y1);
                    if (count < 4096) {
                        ints[count++] = x_int;
                    }
                }
            }

            /* Sort intersections. */
            for (a = 0; a < count - 1; ++a) {
                for (b = a + 1; b < count; ++b) {
                    if (ints[a] > ints[b]) {
                        int tmp = ints[a];

                        ints[a] = ints[b];
                        ints[b] = tmp;
                    }
                }
            }

            /* Fill pairs. */
            for (a = 0; a < count; a += 2) {
                if (a + 1 < count) {
                    int x_start = ints[a];
                    int x_end = ints[a + 1];
                    int x_val;

                    if (x_start < 0)
                        x_start = 0;
                    if (x_end >= HWCSIL_GRID_RES)
                        x_end = HWCSIL_GRID_RES - 1;
                    for (x_val = x_start; x_val <= x_end; ++x_val) {
                        grid[y * HWCSIL_GRID_RES + x_val] = 1;
                    }
                }
            }
        }
    }

    IritFree(px);
    IritFree(py);
    IritFree(ints);

    /* 4. Find starting pixel for Moore Neighborhood Tracing. */
    for (y = 0; y < HWCSIL_GRID_RES; ++y) {
        for (x = 0; x < HWCSIL_GRID_RES; ++x) {
            if (grid[y * HWCSIL_GRID_RES + x]) {
                startX = x;
                startY = y;
                break;
            }
        }
        if (startX != -1) break;
    }

    if (startX == -1) {
        IritFree(grid);
        return NULL;
    }

    /* 5. Moore Neighborhood Tracing. */
    {
        int m_dx[8] = { 0, 1, 1, 1, 0, -1, -1, -1 },
            m_dy[8] = { -1, -1, 0, 1, 1, 1, 0, -1 },
            max_path = HWCSIL_GRID_RES * HWCSIL_GRID_RES,
            path_len = 0,
            cx = startX,
            cy = startY,
            first_step = 1,
            dir = 2, /* Pretend we moved East (2), so left is North (0). */
            *pathX = (int*) IritMalloc(sizeof(int) * max_path),
            *pathY = (int*) IritMalloc(sizeof(int) * max_path);

        do {
            int next_dir, i,
                found = 0;

            if (path_len < max_path) {
                pathX[path_len] = cx;
                pathY[path_len] = cy;
                path_len++;
            }
            else break;

            for (i = 0; i < 8; ++i) {
                int nx, ny;
                /* Turn 135 deg left to ensure we trace outside boundary. */
                next_dir = (dir + 5 + i) % 8;
                nx = cx + m_dx[next_dir];
                ny = cy + m_dy[next_dir];

                if (nx >= 0 && nx < HWCSIL_GRID_RES && ny >= 0 && ny < HWCSIL_GRID_RES) {
                    if (grid[ny * HWCSIL_GRID_RES + nx]) {
                        cx = nx;
                        cy = ny;
                        dir = next_dir;
                        found = 1;
                        break;
                    }
                }
            }

            if (!found) break; /* Isolated pixel. */

            if (!first_step && cx == startX && cy == startY) {
                break;
            }
            first_step = 0;
        } while (1);

        IritFree(grid);

        if (path_len < 3) {
            IritFree(pathX);
            IritFree(pathY);
            return NULL;
        }

        /* 6. Convert traced pixel path back to view-local XY. */
        {
            IPPolygonStruct *NewPl = IritPrsrAllocPolygon(0, NULL, NULL);
            IPVertexStruct *FirstV = NULL,
                *PrevV = NULL;
            int i_path;
            IPObjectStruct *PolyObj,
                *OutList;

            for (i_path = 0; i_path < path_len; ++i_path) {
                IPVertexStruct *NV = IritPrsrAllocVertex2(NULL);
                NV -> Coord[0] = minx + pathX[i_path] * dx / (HWCSIL_GRID_RES - 1);
                NV -> Coord[1] = miny + pathY[i_path] * dy / (HWCSIL_GRID_RES - 1);
                NV -> Coord[2] = 0.0;
                if (FirstV == NULL) FirstV = PrevV = NV;
                else { PrevV -> Pnext = NV; PrevV = NV; }
            }
            IritFree(pathX);
            IritFree(pathY);

            if (PrevV != NULL) {
                PrevV -> Pnext = FirstV;
                NewPl -> PVertex = FirstV;
                IritPrsrUpdatePolyPlane(NewPl);
            }
            else {
                IritPrsrFreePolygonList(NewPl);
                return NULL;
            }

            PolyObj = IritPrsrGenPOLYObject(NewPl);
            if (PolyObj == NULL) {
                IritPrsrFreePolygonList(NewPl);
                return NULL;
            }
            IP_SET_POLYGON_OBJ(PolyObj);
            PolyObj -> Pnext = NULL;

            OutList = IritPrsrGenLISTObject(NULL);
            if (OutList == NULL) {
                IritPrsrFreeObject(PolyObj);
                return NULL;
            }
            IritPrsrListObjectAppend(OutList, PolyObj);

            return OutList;
        }
    }
}

/*****************************************************************************
* DESCRIPTION:                                                               *
* Approximate visible silhouette contour of polygonal solid from given view. *
* Computes a smooth feature-rich 2D contour polygon in world coordinates     *
* using high-resolution rasterization and Laplacian smoothing.               *
* Contour closely approximates the visual silhouette boundary as seen        *
* from the given viewpoint.                                                  *
*                                                                            *
* Algorithm:                                                                 *
* 1: Extracts orthonormal basis (u, v, w) from view matrix.                  *
* 2: Builds view-local transformation mapping world to local XY plane.       *
* 3: Rasterizes solid into 1024×1024 grid, extracts boundary via Moore       *
*    neighborhood tracing, producing high-resolution pixel-level contour.    *
* 4: Applies Laplacian smoothing to remove pixel stair-stepping artifacts.   *
* 5: Subsamples points to requested count (NumCtrl).                         *
* 6: Transforms back to world coordinates using inverse view matrix.         *
* 7: Returns closed polygon in world frame, ready for further processing.    *
*                                                                            *
* PARAMETERS:                                                                *
* Solid: IPObjectStruct pointer to valid polygonal object (checked           *
*        via IP_IS_POLY_OBJ). Multiple polygons/loops allowed.               *
*        Input not modified. NULL or non-polygon returns NULL.               *
* ViewMat: 4x4 homogeneous view matrix. Row 2 (Mat[2][0..2])                 *
*          interpreted as forward direction vector (normalized               *
*          internally). Defines camera frame via orthonormal basis.          *
* NumCtrl: Target number of control points in final contour (positive        *
*          integer). Actual count may vary due to subsampling step           *
*          calculation. Larger values preserve detail; typical 32-128.       *
*                                                                            *
* RETURN VALUE:                                                              *
* IPObjectStruct *: Newly allocated POLY object in world coordinates with    *
*                   closed contour polygon, or NULL on failure (NULL input,  *
*                   invalid matrix, degenerate projection, allocation        *
*                   failure or insufficient boundary points). Caller owns    *
*                   returned object and must free via IritPrsrFreeObject().  *
*                                                                            *
* PERFORMANCE:                                                               *
* - Time: O(V + G) where V = total vertices, G = grid area (1M operations).  *
*   Typically <100ms for complex solids on modern hardware.                  *
* - Memory: O(G) for grid (1MB) + temporary contour arrays. Total ~2-5MB.    *
* - Bottleneck: scanline rasterization; proportional to polygon count and    *
*   projected area.                                                          *
*                                                                            *
* SEE ALSO:                                                                  *
* HWCBuildProjectUnionLocal (low-level grid rasterization),                  *
* HWCIritPrsrBuildViewBasisFromMat (view frame extraction),                  *
* HWCBuildSilhouetteRuledSrf (contour to ruled surface),                     *
* HWCSelectBestViewSampling (multi-view contour scoring).                    *
*****************************************************************************/
static IPObjectStruct *HWCIritPrsrApproxBSplineContourFromSolidView(
                                                 IPObjectStruct *Solid,
                                                 const IrtHmgnMatType ViewMat,
                                                 int NumCtrl)
{
    IrtVecType u, v, w;
    IrtHmgnMatType MatLocal, InvMat;
    IPObjectStruct *UnionLocal, *Poly, *FinalObj;
    IPPolygonStruct *Pl, *FinalPl;
    IPVertexStruct *V, *cur,
        *FirstV = NULL,
        *PrevV = NULL;
    IrtPtType *pts,
        *tmp_pts;
    int i, it, step,
        npts = 0,
        smooth_iters = 5;

    if (Solid == NULL || !IP_IS_POLY_OBJ(Solid))
        return NULL;

    if (!HWCIritPrsrBuildViewBasisFromMat(ViewMat, u, v, w, 1))
        return NULL;

    IritMiscMatGenUnitMat(MatLocal);
    MatLocal[0][0] = u[0]; MatLocal[0][1] = u[1];
    MatLocal[0][2] = u[2]; MatLocal[0][3] = 0.0;
    MatLocal[1][0] = v[0]; MatLocal[1][1] = v[1];
    MatLocal[1][2] = v[2]; MatLocal[1][3] = 0.0;
    MatLocal[2][0] = w[0]; MatLocal[2][1] = w[1];
    MatLocal[2][2] = w[2]; MatLocal[2][3] = 0.0;
    MatLocal[3][0] = 0.0;  MatLocal[3][1] = 0.0;
    MatLocal[3][2] = 0.0; MatLocal[3][3] = 1.0;

    if (!IritMiscMatInverseMatrix(MatLocal, InvMat))
        return NULL;

    /* Build union of projected faces in local XY
    using robust grid rasterization. */
    UnionLocal = HWCBuildProjectUnionLocal(Solid, InvMat);
    if (UnionLocal == NULL) 
        return NULL;
    

    /* UnionLocal contains a single high-resolution pixel-traced boundary.
       Extract it, smooth it lightly to remove pixel steps,
       and return in world coords. */
    Poly = IritPrsrListObjectGet(UnionLocal, 0);
    if (Poly == NULL || !IP_IS_POLY_OBJ(Poly) || Poly -> U.Pl == NULL) {
        IritPrsrFreeObject(UnionLocal);
        return NULL;
    }

    Pl = Poly -> U.Pl;
    V = Pl -> PVertex;
    cur = V;
    do {
        npts++;
        cur = cur -> Pnext;
    } while (cur != NULL && cur != V);

    pts = (IrtPtType*) IritMalloc(sizeof(IrtPtType) * npts);
    tmp_pts = (IrtPtType*) IritMalloc(sizeof(IrtPtType) * npts);

    cur = V;
    for (i = 0; i < npts; ++i) {
        pts[i][0] = cur -> Coord[0];
        pts[i][1] = cur -> Coord[1];
        pts[i][2] = 0.0;
        cur = cur -> Pnext;
    }

    /* Laplacian smoothing to remove pixel stair-stepping artifacts. */
    for (it = 0; it < smooth_iters; ++it) {
        for (i = 0; i < npts; ++i) {
            int im1 = (i - 1 + npts) % npts,
                ip1 = (i + 1) % npts;
            tmp_pts[i][0] = (pts[im1][0] + pts[ip1][0] + pts[i][0]) / 3.0;
            tmp_pts[i][1] = (pts[im1][1] + pts[ip1][1] + pts[i][1]) / 3.0;
            tmp_pts[i][2] = 0.0;
        }
        memcpy(pts, tmp_pts, sizeof(IrtPtType) * npts);
    }

    /* Subsample points. */
    step = npts / NumCtrl;
    if (step < 1) step = 1;

    FinalPl = IritPrsrAllocPolygon(0, NULL, NULL);

    for (i = 0; i < npts; i += step) {
        IrtPtType worldP;
        IPVertexStruct *NV;
        IritMiscMatMultPtby4by4(worldP, pts[i], MatLocal);
        NV = IritPrsrAllocVertex2(NULL);
        NV -> Coord[0] = worldP[0];
        NV -> Coord[1] = worldP[1];
        NV -> Coord[2] = worldP[2];
        if (FirstV == NULL)
            FirstV = PrevV = NV;
        else
            PrevV -> Pnext = NV; PrevV = NV;
    }

    IritFree(pts);
    IritFree(tmp_pts);
    IritPrsrFreeObject(UnionLocal);

    if (PrevV != NULL) {
        PrevV -> Pnext = FirstV;
        FinalPl -> PVertex = FirstV;
        IritPrsrUpdatePolyPlane(FinalPl);
    }
    else {
        IritPrsrFreePolygonList(FinalPl);
        return NULL;
    }

    FinalObj = IritPrsrGenPOLYObject(FinalPl);
    if (FinalObj != NULL) 
        IP_SET_POLYGON_OBJ(FinalObj);
    else 
        IritPrsrFreePolygonList(FinalPl);

    return FinalObj;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
* Select optimal viewing directions for hot-wire silhouette cutting using    *
* area-based scoring over a constrained circular sampling in the XY plane.   *
*                                                                            *
* Candidate views are generated uniformly along a circle centered at the     *
* origin, with camera positions defined by angular parameter theta in        *
* [0, 2π).                                                                   *
* For each candidate view, a silhouette contour is computed and scored based *
* on its projected 2D area. Larger silhouettes are preferred as they         *
* typically preserve more geometric features during cutting.                 *
*                                                                            *
* The function returns all generated views sorted implicitly by their        *
* computed scores (sorting stage assumed to follow). Only views that satisfy *
* basic geometric feasibility (non-vertical direction) are considered.       *
*                                                                            *
* It assumes object is centered at origin; no automatic centering is         *
* performed.                                                                 *
*                                                                            *
* PARAMETERS:                                                                *
* PObj: IPObjectStruct pointer to a polygonal object. Must satisfy           *
*       IP_IS_POLY_OBJ. The object is not modified. NULL returns NULL.       *
* NumSamples: Number of desired output samples (positive integer).           *
*             Internally doubled (2*NumSamples) to improve sampling density  *
*             and robustness of scoring.                                     *
* ResultMat: Optional output buffer. If non-NULL, receives copy of generated *
*            view matrices (size must be 2*NumSamples). If NULL, no copy     *
*            is performed.                                                   *
*                                                                            *
* RETURN VALUE:                                                              *
* IrtHmgnMatType *: Newly allocated array of 2*NumSamples view matrices      *
*                   corresponding to uniformly sampled viewpoints along      *
*                   a circle in the XY plane. Each view is associated with   *
*                   a silhouette-based score (not returned explicitly).      *
*                   NULL is returned on failure (invalid input or            *
*                   allocation). Caller must free using IritFree().          *
*                                                                            *
* SEE ALSO:                                                                  *
* HWCIritPrsrApproxBSplineContourFromSolidView, HWCGenLookAtMatrix,          *
* HWCCalcPolygonArea.                                                        *
*****************************************************************************/
static IrtHmgnMatType *HWCSelectBestViewSampling(IPObjectStruct *PObj,
                                                 int NumSamples,
                                                 IrtHmgnMatType *ResultMat)
{
    int next, s, i, j, best;
    IrtHmgnMatType *ViewMats;
    IrtRType *scores;
    const int NumCtrl = 128;

    if (PObj == NULL || NumSamples <= 0)
        return NULL;

    NumSamples = NumSamples * 2;

    ViewMats = (IrtHmgnMatType*) IritMalloc(sizeof(IrtHmgnMatType) * NumSamples);
    if (ViewMats == NULL)
        return NULL;

    next = 0;

    /* Generate candidate view matrices (Uniform circle in XY plane). */
    for (s = 0; s < NumSamples; ++s) {
        IrtRType theta = 2.0 * 3.14159265358979323846 * s / NumSamples;
        IrtRType r = 2.0;
        IrtVecType CamPos, Center, Up;

        CamPos[0] = cos(theta) * r;
        CamPos[1] = sin(theta) * r;
        CamPos[2] = 0.0; /* Fixed Z value, no variation in height. */

        Center[0] = 0.0;
        Center[1] = 0.0;
        Center[2] = 0.0;

        Up[0] = 0.0;
        Up[1] = 0.0;
        Up[2] = 1.0;

        HWCGenLookAtMatrix(CamPos, Center, Up, ViewMats[next++]);
    }

    scores = (IrtRType*) IritMalloc(sizeof(IrtRType) * NumSamples);
    if (scores == NULL) {
        IritFree(ViewMats);
        return NULL;
    }

    for (i = 0; i < NumSamples; ++i) {
        IPObjectStruct *Contour;
        IrtRType wxy_len;

        scores[i] = 0.0;

        /* Reject near-vertical views (where wxy_len < 0.1) since the horizontal
           wire cannot reach the required vertical slope. */
        wxy_len = sqrt(ViewMats[i][2][0] * ViewMats[i][2][0] + ViewMats[i][2][1] *
            ViewMats[i][2][1]);
        if (wxy_len < 0.1) {
            continue;
        }

        Contour = HWCIritPrsrApproxBSplineContourFromSolidView(PObj,
            ViewMats[i],
            NumCtrl);
        if (Contour != NULL) {
            scores[i] = HWCCalcPolygonArea(Contour);
            IritPrsrFreeObject(Contour);
        }
    }

    /* Reorder ViewMats in-place by descending score (simple selection sort). */
    for (i = 0; i < NumSamples - 1; ++i) {
        best = i;
        for (j = i + 1; j < NumSamples; ++j) {
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

    /* Build final view list: first 5 are fixed cardinal + diagonal views,
       then fill remaining slots with best diverse scored views. */
    {
        int RequestedViews = NumSamples / 2;
        int num_selected = 0;
        IrtHmgnMatType *SelectedMats = (IrtHmgnMatType*)
            IritMalloc(sizeof(IrtHmgnMatType) * NumSamples);
        /* 0°=+X, 45°=diagonal, 90°=+Y, 135°=other diagonal */
        IrtRType fixed_angles[4] = { 0.0, 45.0, 90.0, 135.0 };

        for (i = 0; i < 4 && num_selected < RequestedViews; ++i) {
            IrtRType theta = fixed_angles[i] * 3.14159265358979323846 / 180.0;
            IrtVecType CamPos, Center, Up;
            CamPos[0] = cos(theta) * 2.0;
            CamPos[1] = sin(theta) * 2.0;
            CamPos[2] = 0.0;
            Center[0] = 0.0; Center[1] = 0.0; Center[2] = 0.0;
            Up[0] = 0.0; Up[1] = 0.0; Up[2] = 1.0;
            HWCGenLookAtMatrix(CamPos, Center, Up,
                SelectedMats[num_selected++]);
        }

        /* Fill remaining slots with best-scored views that are not too
           similar (within 15 degrees) to any already-selected view. */
        for (i = 0; i < NumSamples && num_selected < RequestedViews; ++i) {
            int is_similar = 0;
            for (j = 0; j < num_selected; ++j) {
                IrtRType dot =
                    ViewMats[i][2][0] * SelectedMats[j][2][0] +
                    ViewMats[i][2][1] * SelectedMats[j][2][1] +
                    ViewMats[i][2][2] * SelectedMats[j][2][2];
                if (IRIT_FABS(dot) >
                    cos(15.0 * 3.14159265358979323846 / 180.0)) {
                    is_similar = 1;
                    break;
                }
            }
            if (!is_similar) {
                memcpy(SelectedMats[num_selected++], ViewMats[i],
                    sizeof(IrtHmgnMatType));
            }
        }

        /* If still not enough, relax the threshold and accept any
           non-duplicate view. */
        for (i = 0; i < NumSamples && num_selected < RequestedViews; ++i) {
            int is_exact = 0;
            for (j = 0; j < num_selected; ++j) {
                IrtRType dot =
                    ViewMats[i][2][0] * SelectedMats[j][2][0] +
                    ViewMats[i][2][1] * SelectedMats[j][2][1] +
                    ViewMats[i][2][2] * SelectedMats[j][2][2];
                if (IRIT_FABS(dot) > 0.999) {
                    is_exact = 1;
                    break;
                }
            }
            if (!is_exact) {
                memcpy(SelectedMats[num_selected++], ViewMats[i],
                    sizeof(IrtHmgnMatType));
            }
        }

        memcpy(ViewMats, SelectedMats, sizeof(IrtHmgnMatType) * num_selected);
        IritFree(SelectedMats);
    }

    if (ResultMat != NULL) {
        memcpy(ResultMat, ViewMats, sizeof(IrtHmgnMatType) * NumSamples);
    }

    IritFree(scores);
    return ViewMats;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Select NumSamples uniformly spaced view directions in the XY plane.      *
*   Views start at 0 degrees and are spaced 180/NumSamples degrees apart,    *
*   covering X, Y, and diagonal directions without symmetric duplicates      *
*   (since a view at angle A produces the same silhouette as A + 180).       *
*                                                                            *
* PARAMETERS:                                                                *
*   NumSamples:  Number of desired view directions (positive integer).       *
*                                                                            *
* RETURN VALUE:                                                              *
*   IrtHmgnMatType *: Newly allocated array of NumSamples view matrices.     *
*                     NULL on failure. Caller must free using IritFree().    *
*                                                                            *
* SEE ALSO:                                                                  *
*   HWCSelectBestViewSampling, HWCGenerateGCodeFromObj                       *
*****************************************************************************/
static IrtHmgnMatType *HWCSelectViewsSampling(int NumSamples)
{
    int i;
    IrtRType step;
    IrtHmgnMatType *ViewMats;

    if (NumSamples <= 0)
        return NULL;

    ViewMats = (IrtHmgnMatType*) IritMalloc(sizeof(IrtHmgnMatType) * NumSamples);
    if (ViewMats == NULL)
        return NULL;

    /* Space views by 180/N degrees so that no two views are mirrors
       of each other (a view at angle A covers the same silhouette
       as A + 180). */
    step = 180.0 / NumSamples;

    for (i = 0; i < NumSamples; ++i) {
        IrtRType theta = (step * i) * 3.14159265358979323846 / 180.0;
        IrtVecType CamPos, Center, Up;

        CamPos[0] = cos(theta) * 2.0;
        CamPos[1] = sin(theta) * 2.0;
        CamPos[2] = 0.0;

        Center[0] = 0.0;
        Center[1] = 0.0;
        Center[2] = 0.0;

        Up[0] = 0.0;
        Up[1] = 0.0;
        Up[2] = 1.0;

        HWCGenLookAtMatrix(CamPos, Center, Up, ViewMats[i]);
    }

    return ViewMats;
}

/*****************************************************************************
* AUXILIARY:								                                 *
* Auxiliary function to find missing boundary edges in an open solid.	     *
*****************************************************************************/
static IritPrsrHWCEdgeStruct *HWCFindBoundaryEdges(const IPObjectStruct *Solid,
                                                   int *OutNumEdges)
{
    int i, j, k, StackTop, CurrentLoopCount,
        MaxEdges = 0,
        NumEdges = 0,
        NumRawBoundaryEdges = 0,
        BestLoopCount = 0;
    IrtRType CurrentLoopLength, CurrentLoopAvgZ,
        BestLoopLength = -1.0,
        BestLoopAvgZ = IRIT_INFNTY,
        EPS = 1e-4;
    int *EdgeCounts, *Visited, *Stack, *CurrentLoopIndices;
    IPPolygonStruct *Pl;
    IPVertexStruct *V, *VNext;
    IritPrsrHWCEdgeStruct *AllEdges,
        *RawBoundaryEdges = NULL,
        *BestLoopEdges = NULL;

    if (OutNumEdges)
        *OutNumEdges = 0;
    if (!Solid || !IP_IS_POLY_OBJ(Solid))
        return NULL;

    /* Count total possible edges. */
    for (Pl = Solid -> U.Pl; Pl != NULL; Pl = Pl -> Pnext) {
        V = Pl -> PVertex;
        if (!V)
            continue;
        do {
            MaxEdges++;
            V = V -> Pnext;
        } while (V != NULL && V != Pl -> PVertex);
    }

    if (MaxEdges == 0)
        return NULL;

    AllEdges = (IritPrsrHWCEdgeStruct*)
        IritMalloc(sizeof(IritPrsrHWCEdgeStruct) * MaxEdges);
    EdgeCounts = (int*) IritMalloc(sizeof(int) * MaxEdges);

    /* Collect all edges. */
    for (Pl = Solid -> U.Pl; Pl != NULL; Pl = Pl -> Pnext) {
        V = Pl -> PVertex;
        if (!V)
            continue;
        do {
            VNext = V -> Pnext ? V -> Pnext : Pl -> PVertex;
            IRIT_PT_COPY(AllEdges[NumEdges].Pt1, V -> Coord);
            IRIT_PT_COPY(AllEdges[NumEdges].Pt2, VNext -> Coord);
            EdgeCounts[NumEdges] = 1;
            NumEdges++;
            V = V -> Pnext;
        } while (V != NULL && V != Pl -> PVertex);
    }

    /* O(N^2) robust edge matching. */
    for (i = 0; i < NumEdges; i++) {
        if (EdgeCounts[i] == 0)
            continue;

        for (j = i + 1; j < NumEdges; j++) {
            int MatchFwd = 1,
                MatchRev = 1;

            if (EdgeCounts[j] == 0)
                continue;

            for (k = 0; k < 3; k++) {
                if (IRIT_FABS(AllEdges[i].Pt1[k] - AllEdges[j].Pt1[k]) > EPS ||
                    IRIT_FABS(AllEdges[i].Pt2[k] - AllEdges[j].Pt2[k]) > EPS) {
                    MatchFwd = 0;
                }
                if (IRIT_FABS(AllEdges[i].Pt1[k] - AllEdges[j].Pt2[k]) > EPS ||
                    IRIT_FABS(AllEdges[i].Pt2[k] - AllEdges[j].Pt1[k]) > EPS) {
                    MatchRev = 0;
                }
            }

            if (MatchFwd || MatchRev) {
                EdgeCounts[i]++;
                EdgeCounts[j] = 0;
            }
        }
    }

    RawBoundaryEdges = (IritPrsrHWCEdgeStruct*)
        IritMalloc(sizeof(IritPrsrHWCEdgeStruct) * NumEdges);
    for (i = 0; i < NumEdges; i++) {
        if (EdgeCounts[i] == 1) {
            RawBoundaryEdges[NumRawBoundaryEdges++] = AllEdges[i];
        }
    }

    IritFree(AllEdges);
    IritFree(EdgeCounts);

    if (NumRawBoundaryEdges == 0) {
        IritFree(RawBoundaryEdges);
        return NULL;
    }

    /* Group boundary edges into loops using an adjacency graph and DFS. */
    Visited = (int*) IritMalloc(sizeof(int) * NumRawBoundaryEdges);
    for (i = 0; i < NumRawBoundaryEdges; i++)
        Visited[i] = 0;

    Stack = (int*) IritMalloc(sizeof(int) * NumRawBoundaryEdges);
    CurrentLoopIndices = (int*) IritMalloc(sizeof(int) * NumRawBoundaryEdges);

    for (i = 0; i < NumRawBoundaryEdges; i++) {
        if (Visited[i])
            continue;

        CurrentLoopLength = 0.0;
        CurrentLoopCount = 0;
        StackTop = 0;

        Stack[StackTop++] = i;
        Visited[i] = 1;

        while (StackTop > 0) {
            int Curr = Stack[--StackTop];
            IrtRType Len;

            CurrentLoopIndices[CurrentLoopCount++] = Curr;
            Len = sqrt(IRIT_SQR(RawBoundaryEdges[Curr].Pt1[0] -
                RawBoundaryEdges[Curr].Pt2[0]) +
                IRIT_SQR(RawBoundaryEdges[Curr].Pt1[1] -
                    RawBoundaryEdges[Curr].Pt2[1]) +
                IRIT_SQR(RawBoundaryEdges[Curr].Pt1[2] -
                    RawBoundaryEdges[Curr].Pt2[2]));
            CurrentLoopLength += Len;

            for (j = 0; j < NumRawBoundaryEdges; j++) {
                int Match11 = 1, Match12 = 1, Match21 = 1, Match22 = 1;

                if (Visited[j])
                    continue;

                for (k = 0; k < 3; k++) {
                    if (IRIT_FABS(RawBoundaryEdges[Curr].Pt1[k] -
                        RawBoundaryEdges[j].Pt1[k]) > EPS) Match11 = 0;
                    if (IRIT_FABS(RawBoundaryEdges[Curr].Pt1[k] -
                        RawBoundaryEdges[j].Pt2[k]) > EPS) Match12 = 0;
                    if (IRIT_FABS(RawBoundaryEdges[Curr].Pt2[k] -
                        RawBoundaryEdges[j].Pt1[k]) > EPS) Match21 = 0;
                    if (IRIT_FABS(RawBoundaryEdges[Curr].Pt2[k] -
                        RawBoundaryEdges[j].Pt2[k]) > EPS) Match22 = 0;
                }

                if (Match11 || Match12 || Match21 || Match22) {
                    Visited[j] = 1;
                    Stack[StackTop++] = j;
                }
            }
        }

        /* Compute average Z of this loop's edge midpoints. */
        CurrentLoopAvgZ = 0.0;
        for (j = 0; j < CurrentLoopCount; j++) {
            int ei = CurrentLoopIndices[j];
            CurrentLoopAvgZ += (RawBoundaryEdges[ei].Pt1[2] +
                                RawBoundaryEdges[ei].Pt2[2]) * 0.5;
        }
        if (CurrentLoopCount > 0)
            CurrentLoopAvgZ /= CurrentLoopCount;

        /* Always pick the loop with the lowest average Z, so the bottom
           opening is chosen over the top one regardless of size. */
        if (CurrentLoopAvgZ < BestLoopAvgZ) {
            BestLoopLength = CurrentLoopLength;
            BestLoopAvgZ = CurrentLoopAvgZ;
            if (BestLoopEdges)
                IritFree(BestLoopEdges);
            BestLoopEdges = (IritPrsrHWCEdgeStruct*)
                IritMalloc(sizeof(IritPrsrHWCEdgeStruct) * CurrentLoopCount);
            for (j = 0; j < CurrentLoopCount; j++) {
                BestLoopEdges[j] = RawBoundaryEdges[CurrentLoopIndices[j]];
            }
            BestLoopCount = CurrentLoopCount;
        }
    }

    IritFree(Stack);
    IritFree(CurrentLoopIndices);
    IritFree(Visited);
    IritFree(RawBoundaryEdges);

    if (OutNumEdges)
        *OutNumEdges = BestLoopCount;

    return BestLoopEdges;
}

/*****************************************************************************
* AUXILIARY:								                                 *
* Auxiliary function to calculate point to segment distance in 2D.	         *
*****************************************************************************/
static IrtRType HWCDistPointSegment2D(IrtPtType P,
                                      IrtPtType A,
                                      IrtPtType B)
{
    IrtRType
        L2 = IRIT_SQR(A[0] - B[0]) + IRIT_SQR(A[1] - B[1]);

    if (L2 < 1e-10)
        return sqrt(IRIT_SQR(P[0] - A[0]) + IRIT_SQR(P[1] - A[1]));
    else {
        IrtRType
            T = ((P[0] - A[0]) * (B[0] - A[0]) +
                (P[1] - A[1]) * (B[1] - A[1])) / L2;
        IrtPtType Proj;

        if (T < 0.0)
            T = 0.0;
        else if (T > 1.0)
            T = 1.0;

        Proj[0] = A[0] + T * (B[0] - A[0]);
        Proj[1] = A[1] + T * (B[1] - A[1]);

        return sqrt(IRIT_SQR(P[0] - Proj[0]) + IRIT_SQR(P[1] - Proj[1]));
    }
}

/*****************************************************************************
* DESCRIPTION:                                                               *
* Build a ruled surface suitable for 2D hot-wire silhouette cutting.         *
*                                                                            *
* The hot-wire machine traces the silhouette contour in machine XZ space     *
* while keeping both wire endpoints (front/back clamps) at the same XZ       *
* position, separated only in Y by FoamDepth. This requires a surface        *
* whose V-isoparametric curves are the silhouette path:                      *
*                                                                            *
*   S(u, vMin) = silhouette at Y = -FoamDepth/2   (front clamp path)         *
*   S(u, vMax) = silhouette at Y = +FoamDepth/2   (back  clamp path)         *
*                                                                            *
* Steps:                                                                     *
*   1. Project each contour vertex into view-local 2D (ignoring world Z).    *
*   2. Compute AABB and scale uniformly to fit FoamWidth x FoamHeight        *
*      with a 10% safety margin, centered at the machine origin.             *
*   3. Build an E3 B-spline curve: X=local_u_scaled,                         *
*                                  Z=local_v_scaled+FoamH/2                  *
*      Y=-FoamDepth/2 (front face relative to foam center).                  *
*   4. Extrude the curve along (0, FoamDepth, 0) to get S(u,v).              *
*   5. Tag with IRIT_ATTR_ID_Dir = CAGD_CONST_U_DIR so IritPrsrHWCCreatePath *
*      samples in V direction (silhouette), giving both clamps the full      *
*      silhouette path in sync.                                              *
*                                                                            *
* PARAMETERS:                                                                *
*   Contour: Polygon object from                                             *
*            HWCIritPrsrApproxBSplineContourFromSolidView.                   *
*   Solid: The original 3D polygonal solid used to find missing boundaries.  *
*   ViewMat: View matrix (row 2 = forward/view direction).                   *
*   Params: HWC machine parameters; uses FoamWidth, FoamHeight, FoamDepth.   *
*                                                                            *
* RETRUN VALUE:                                                              *
* IPObjectStruct*: Surface on success, NULL on failure.                      *
*         Caller owns the returned object and must free it.                  *
*****************************************************************************/
static IPObjectStruct *HWCBuildSilhouetteRuledSrf(const IPObjectStruct *Contour,
                                                  const IPObjectStruct *Solid,
                                                  const IrtHmgnMatType ViewMat,
                                                  const IritPrsrHWCDataStruct *Params)
{
    int k, n, guard,
        best_start = 0,
        best_len = 0,
        num_boundary_edges = 0,
        top_start = 0,
        top_len = n;
    IrtRType minx, miny, maxx, maxy, scale_x, scale_z, scale, cx, cy;
    IrtVecType u_vec, v_vec, w_vec;
    IrtHmgnMatType MatLocal;
    IrtPtType *pts2d;
    IPVertexStruct *V, *cur;
    CagdCrvStruct *Crv;
    CagdVecStruct YDir;
    CagdSrfStruct *Srf;
    IPObjectStruct *SrfObj;

    if (Contour == NULL || !IP_IS_POLY_OBJ((IPObjectStruct*) Contour) ||
        Contour -> U.Pl == NULL || Params == NULL)
        return NULL;

    /* 1. Build view-local orthonormal basis
       (u_vec=right, v_vec=up, w_vec=fwd). */
    if (!HWCIritPrsrBuildViewBasisFromMat(ViewMat, u_vec, v_vec, w_vec, 1))
        return NULL;

    /* MatLocal maps world  ->  view-local (rows = basis vectors). */
    IritMiscMatGenUnitMat(MatLocal);
    MatLocal[0][0] = u_vec[0];
    MatLocal[0][1] = u_vec[1];
    MatLocal[0][2] = u_vec[2];
    MatLocal[0][3] = 0.0;
    MatLocal[1][0] = v_vec[0];
    MatLocal[1][1] = v_vec[1];
    MatLocal[1][2] = v_vec[2];
    MatLocal[1][3] = 0.0;
    MatLocal[2][0] = w_vec[0];
    MatLocal[2][1] = w_vec[1];
    MatLocal[2][2] = w_vec[2];
    MatLocal[2][3] = 0.0;
    MatLocal[3][0] = 0.0;
    MatLocal[3][1] = 0.0;
    MatLocal[3][2] = 0.0;
    MatLocal[3][3] = 1.0;

    /* 2. Count and validate polygon vertices. */
    V = Contour -> U.Pl -> PVertex;
    if (V == NULL) return NULL;

    n = 0;
    cur = V;
    guard = 0;
    do {
        ++n;
        cur = cur -> Pnext;
        if (++guard > 200000) { n = 0; break; }
    } while (cur != NULL && cur != V);

    if (n < 3) return NULL;

    /* 3. Project each vertex to view-local 2D, ignoring world Z.
     *   localP[0] = view-right component -> machine X
     *   localP[1] = view-up   component -> machine Z (vertical) */
    pts2d = (IrtPtType*) IritMalloc(sizeof(IrtPtType) * n);
    cur = V;
    guard = 0;
    for (k = 0; k < n; ++k) {
        /* Compute true dot products to project world coordinates onto the view plane. */
        pts2d[k][0] = cur -> Coord[0] * u_vec[0] + cur -> Coord[1] * u_vec[1]
            + cur -> Coord[2] * u_vec[2];
        pts2d[k][1] = cur -> Coord[0] * v_vec[0] + cur -> Coord[1] * v_vec[1]
            + cur -> Coord[2] * v_vec[2];
        pts2d[k][2] = 0.0;
        cur = cur -> Pnext;
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
    scale_x = (Params -> FoamWidth * 0.9) / (maxx - minx);
    scale_z = (Params -> FoamHeight * 0.9) / (maxy - miny);
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

    /* 7. Robust Bottom Line Removal based on exact distance to missing edges.
     *    We project the 3D boundary edges to 2D line segments.
     *    Then we measure the distance from every point on the contour to these segments.
     *    Points very close to ANY boundary segment are identified as part
     *    of the missing bottom. We then keep the longest contiguous
     *    sequence of points that are NOT missing. */
    IritPrsrHWCEdgeStruct *boundary_edges = HWCFindBoundaryEdges(Solid, &num_boundary_edges);

    if (boundary_edges != NULL && num_boundary_edges > 0) {
        int b,
            *is_boundary = (int*) IritMalloc(sizeof(int) * n);
        IrtRType 
            diag = sqrt((maxx - minx) * (maxx - minx) + (maxy - miny) * (maxy - miny)),
            dist_thresh = diag * 0.005; /* 0.5% of diagonal as tolerance */

        for (k = 0; k < n; ++k) {
            is_boundary[k] = 0;
            for (b = 0; b < num_boundary_edges; ++b) {
                IrtPtType proj_p1, proj_p2;
                proj_p1[0] = boundary_edges[b].Pt1[0] * u_vec[0] + boundary_edges[b].Pt1[1]
                    * u_vec[1] + boundary_edges[b].Pt1[2] * u_vec[2];
                proj_p1[1] = boundary_edges[b].Pt1[0] * v_vec[0] + boundary_edges[b].Pt1[1]
                    * v_vec[1] + boundary_edges[b].Pt1[2] * v_vec[2];
                proj_p2[0] = boundary_edges[b].Pt2[0] * u_vec[0] + boundary_edges[b].Pt2[1]
                    * u_vec[1] + boundary_edges[b].Pt2[2] * u_vec[2];
                proj_p2[1] = boundary_edges[b].Pt2[0] * v_vec[0] + boundary_edges[b].Pt2[1]
                    * v_vec[1] + boundary_edges[b].Pt2[2] * v_vec[2];

                if (HWCDistPointSegment2D(pts2d[k], proj_p1, proj_p2) < dist_thresh) {
                    is_boundary[k] = 1;
                    break;
                }
            }
        }
        IritFree(boundary_edges);

        /* Find the longest contiguous sequence of is_boundary == 0 (wrapping around). */
        for (k = 0; k < n; ++k) {
            if (is_boundary[k] == 0) {
                int len = 0;
                int idx = k;
                while (len < n && is_boundary[idx] == 0) {
                    len++;
                    idx = (idx + 1) % n;
                }
                if (len > best_len) {
                    best_len = len;
                    best_start = k;
                }
                /* Fast forward k to the end of this non-boundary sequence. */
                /* But only up to n-1 to avoid infinite loops if the whole thing wraps. */
                k += len - 1;
            }
        }

        IritFree(is_boundary);

        if (best_len > 0) {
            top_start = best_start;
            top_len = best_len;
        }
    }
    else {
        /* Fallback if no boundary edges found: remove bottom 10%. */
        IrtRType
            bottom_thresh = miny + (maxy - miny) * 0.10,
            min_x_bottom = maxx + 1.0,
            max_x_bottom = minx - 1.0;
        int left_idx = -1,
            right_idx = -1;

        for (k = 0; k < n; ++k) {
            if (pts2d[k][1] <= bottom_thresh) {
                if (pts2d[k][0] < min_x_bottom) {
                    min_x_bottom = pts2d[k][0];
                    left_idx = k;
                }
                if (pts2d[k][0] > max_x_bottom) {
                    max_x_bottom = pts2d[k][0];
                    right_idx = k;
                }
            }
        }

        if (left_idx != -1 && right_idx != -1 && left_idx != right_idx &&
            (max_x_bottom - min_x_bottom) > (maxx - minx) * 0.10) {

            /* Path 1: left_idx to right_idx (incrementing). */
            IrtRType 
                max_y_path1 = miny;
            int len1 = (right_idx - left_idx + n) % n;
            for (k = 0; k <= len1; ++k) {
                int idx = (left_idx + k) % n;
                if (pts2d[idx][1] > max_y_path1) 
                    max_y_path1 = pts2d[idx][1];
            }

            /* Path 2: right_idx to left_idx (incrementing). */
            IrtRType 
                max_y_path2 = miny;
            int len2 = (left_idx - right_idx + n) % n;
            for (k = 0; k <= len2; ++k) {
                int idx = (right_idx + k) % n;
                if (pts2d[idx][1] > max_y_path2) 
                    max_y_path2 = pts2d[idx][1];
            }

            if (max_y_path1 > max_y_path2) {
                top_start = left_idx;
                top_len = len1 + 1;
            }
            else {
                top_start = right_idx;
                top_len = len2 + 1;
            }
        }
    }

    /* Extract the top path. */
    IrtPtType *top_pts = (IrtPtType*) IritMalloc(sizeof(IrtPtType) * top_len);
    for (k = 0; k < top_len; ++k) {
        IRIT_PT_COPY(top_pts[k], pts2d[(top_start + k) % n]);
    }

    /* 8. Build silhouette curve mapped exactly into 3D World space.
     *    The 2D contour is originally in the (u_vec, v_vec) plane of the view.
     *    We rebuild the exact 3D orientation.
     *
     *    The extruded surface must be aligned with w_vec in 3D.
     *    The start of the extrusion will be -w_vec * FoamDepth/2.
     *    We also shift the entire cut upwards to be centered in the foam block. */
    Crv = IritCagdBspCrvNew(top_len, 2, CAGD_PT_E3_TYPE);
    if (Crv == NULL) {
        IritFree(pts2d);
        IritFree(top_pts);
        return NULL;
    }

    for (k = 0; k < top_len; ++k) {
        IrtRType sx = (top_pts[k][0] - cx) * scale;
        IrtRType sy = (top_pts[k][1] - cy) * scale;

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
        IrtRType 
            slope = w_vec[2] / wxy_len,
            z_shift = slope * (195.0 - wxy_len * Params -> FoamDepth) +
            2.0 * wz * slope * slope;

        Crv -> Points[1][k] = wx - w_vec[0] * Params -> FoamDepth * 0.5;
        Crv -> Points[2][k] = wy - w_vec[1] * Params -> FoamDepth * 0.5;
        Crv -> Points[3][k] = wz - w_vec[2] * Params -> FoamDepth * 0.5 +
            Params -> FoamHeight * 0.5 + z_shift;
    }
    IritCagdBspKnotUniformOpen(top_len, 2, Crv -> KnotVector);

    IritFree(pts2d);
    IritFree(top_pts);

    /* 8. Extrude exactly along the true 3D view direction by FoamDepth.
     *    This creates a cylinder (surface) whose rulings are exactly w_vec.
     *    When IritPrsrHWCProcessRulingPair processes this surface, it computes
     *    Angle = atan2(Dx, Dy) = atan2(w_vec[0], w_vec[1]), and rotates the
     *    clamps accordingly to perfectly align the Y-axis of the machine with w_vec. */
    YDir.Vec[0] = w_vec[0] * Params -> FoamDepth;
    YDir.Vec[1] = w_vec[1] * Params -> FoamDepth;
    YDir.Vec[2] = w_vec[2] * Params -> FoamDepth;

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

/*****************************************************************************
* AUXILIARY:                                                                 *
*                                                                            *
* Auxiliary function to combine per-view GCode files into one.               *
*****************************************************************************/
static void HWCCombineGCodeFiles(const char const **GcodeFiles,
                                 int NumFiles,
                                 const char *OutPath,
                                 const IritPrsrHWCDataStruct *Params)
{
#define IRIT_MAX_LINE_LEN 512
    FILE *fout, *fin;
    int fi,
        f = 300,
        safe_entry_lines_to_modify = 0;
    double
        bOffset = 0.0,
        lastB = 0.0,
        x = 0.0,
        y = 0.0,
        z = 0.0,
        a = 0.0,
        b = 0.0,
        min_z_seen = 99999.0,
        ExtraSafeZClearance = 150.0;
    char line[IRIT_MAX_LINE_LEN];

    fout = fopen(OutPath, "w");
    if (fout == NULL) {
        fprintf(stderr, "HWCCombineGCodeFiles: cannot open '%s'.\n", OutPath);
        return;
    }
#ifdef DEBUG
    fprintf(fout, "; Combined GCode - %d views\n\n", NumFiles);
#endif /* DEBUG */

    for (fi = 0; fi < NumFiles; ++fi) {
        lastB = bOffset;

        fin = fopen(GcodeFiles[fi], "r");
        if (fin == NULL) {
            fprintf(stderr, "HWCCombineGCodeFiles: cannot open '%s'.\n",
                GcodeFiles[fi]);
            continue;
        }
#ifdef DEBUG
        fprintf(fout, "; === View %d start ===\n", fi);
#endif /* DEBUG */

        while (fgets(line, sizeof(line), fin) != NULL) {
            if (strncmp(line, ";Performing safe entry movement", 31) == 0) {
                safe_entry_lines_to_modify = 2;
            }

            if (strncmp(line, "; Performing vertical safe lift to safe Z height.", 49) == 0) {
                safe_entry_lines_to_modify = 1;
            }

            if (strncmp(line, "G1 ", 3) == 0 &&
                strstr(line, "B") != NULL &&
                strstr(line, "X") != NULL) {
                if (sscanf(line, "G1 X%lf Y%lf Z%lf A%lf B%lf F%d",
                    &x, &y, &z, &a, &b, &f) == 6) {
                    /* B from the per-view GCode is the ABSOLUTE rotation
                       angle for this view (each view resets _TotalFoamRotation
                       to 0).  Adding bOffset blindly treated it as incremental,
                       making the machine accumulate extra full rotations.
                       Instead, compute the shortest angular path. */
                    {
                        IrtRType diff = fmod(b - lastB, 360.0);
                        if (diff > 180.0) 
                            diff -= 360.0;
                        if (diff < -180.0) 
                            diff += 360.0;
                        lastB = lastB + diff;
                    }
                    if (z < min_z_seen) 
                        min_z_seen = z;

                    /* Override any safe Z exits to provide massive clearance for taller foam blocks. */
                    if (Params != NULL && z >= Params -> FoamHeight) 
                        z = Params -> FoamHeight + ExtraSafeZClearance;
                    
                    if (Params != NULL && a >= Params -> FoamHeight) 
                        a = Params -> FoamHeight + ExtraSafeZClearance;
                    
                    if (safe_entry_lines_to_modify > 0) {
                        f = 350;
                        safe_entry_lines_to_modify--;
                    }

                    fprintf(fout,
                        "G1 X%.3f Y%.3f Z%.3f A%.3f B%.3f F%d\n",
                        x, y, z, a, lastB, f);
                    continue;
                }
            }
            /* Skip per-view Go Home commands; the combined file
               emits its own G28 at the very end. */
            if (strncmp(line, "G28", 3) == 0)
                continue;
            fputs(line, fout);
        }

        fclose(fin);
        bOffset = lastB;
#ifdef DEBUG
        fprintf(fout, "; === View %d end ===\n\n", fi);
#endif /* DEBUG */
    }

    if (Params != NULL) {
        fprintf(fout, "\n; === Slice off base ===\n");
        /* Move safely above the foam. */
        z = Params -> FoamHeight + ExtraSafeZClearance;
        a = Params -> FoamHeight + ExtraSafeZClearance;
        fprintf(fout, "G1 X%.3f Y%.3f Z%.3f A%.3f B%.3f F300\n", x, y, z, a, lastB);
        /* Move to the starting side of the machine (safely clear of foam). */
        x = 20.0;
        y = 20.0;
        fprintf(fout, "G1 X%.3f Y%.3f Z%.3f A%.3f B%.3f F300\n", x, y, z, a, lastB);
        /* Drop down to the slice height (1mm below the lowest point of the object). */
        z = min_z_seen - 1.0;
        if (z < Params -> MinimalHeight) 
            z = Params -> MinimalHeight;
        a = z;
        fprintf(fout, "G1 X%.3f Y%.3f Z%.3f A%.3f B%.3f F150\n", x, y, z, a, lastB);
        /* Slice entirely through the foam to the other side. */
        x = 450.0;
        y = 450.0;
        fprintf(fout, "G1 X%.3f Y%.3f Z%.3f A%.3f B%.3f F150\n", x, y, z, a, lastB);
        /* Pull back up to safe height. */
        z = Params -> FoamHeight + ExtraSafeZClearance;
        a = Params -> FoamHeight + ExtraSafeZClearance;
        fprintf(fout, "G1 X%.3f Y%.3f Z%.3f A%.3f B%.3f F300\n", x, y, z, a, lastB);
    }

    fprintf(fout, "\n; === Return Home Synchronously ===\n");
    if (Params != NULL) {
        /* Move XY and B to 0 while maintaining the high Z clearance to avoid collisions. */
        fprintf(fout, "G1 X0.000 Y0.000 Z%.3f A%.3f B0.000 F300\n", z, a);
    }
    /* Finally drop Z back to 0 synchronously to avoid snapping the wire. */
    fprintf(fout, "G1 X0.000 Y0.000 Z0.000 A0.000 B0.000 F300\n");

    /* Emit G28 only to ensure firmware resets coordinate systems if necessary,
       but since we are already at 0,0,0,0,0 it won't physically move sequentially. */
    fprintf(fout, "G28; Go Home\n");
    fclose(fout);
#ifdef DEBUG
    printf("Combined GCode written to: %s\n", OutPath);
#endif // DEBUG
#undef IRIT_MAX_LINE_LEN
}

/*****************************************************************************
* AUXILIARY:                                                                 *
*   Recursive function to convert polylines to polygons in a given object.   *
*   It traverses lists and converts matching polyline objects in-place.      *
*                                                                            *
* PARAMETERS:                                                                *
*   PObj:       The IPObjectStruct to process (can be a list or polyline).   *
*****************************************************************************/
static void HWCConvertPolylinesToPolygons(IPObjectStruct *PObj) {
    if (PObj == NULL) 
        return;
    if (IP_IS_OLST_OBJ(PObj)) {
        IPObjectStruct *PList;
        int i = 0;
        while ((PList = IritPrsrListObjectGet(PObj, i++)) != NULL) {
            HWCConvertPolylinesToPolygons(PList);
        }
    }
    else if (IP_IS_POLY_OBJ(PObj) && IP_IS_POLYLINE_OBJ(PObj)) 
        IP_SET_POLYGON_OBJ(PObj);
    
    if (PObj -> Pnext != NULL) 
        HWCConvertPolylinesToPolygons(PObj -> Pnext);
    
}

/*****************************************************************************
* DESCRIPTION:                                                               M
* Main public entry point. Encapsulates entire hot-wire cutting pipeline.    M
*                                                                            *
* PARAMETERS:                                                                M
* RawModel:         Input 3D model object.                                   M
* OutputGCodePath:  Path to output GCode file.                               M
* NumViews:         Number of view directions to generate.                   M
* OutputITDType:    0 to skip ITD, 1 for polylines, 2 for polygons.          M
* ITDLength:        Wire path visualization length (0.0 for default).        M
* fileName:         Name of the output ITD file.                             M
* uniformCut:       View selection mode. 0 to use score-based best-view      M
*                   selection (HWCSelectBestViewSampling), 1 to use          M
*                   uniformly spaced views at 180/N degree intervals         M
*                   (HWCSelectViewsSampling).                                M
*                                                                            *
* RETURN VALUE:                                                              M
*   IrtRType: 1 on success, 0 on failure.                                    M
*                                                                            *
* SEE ALSO:                                                                  M
*   IritPrsrHWCCreatePath, HWCBuildSilhouetteRuledSrf,                       M
*   HWCSelectBestViewSampling, HWCSelectViewsSampling                        M
*                                                                            *
* KEYWORDS:                                                                  M
*   HWCGenerateGCodeFromObj, hot-wire, GCode, pipeline.                      M
*****************************************************************************/
IrtRType HWCGenerateGCodeFromObj(IPObjectStruct *RawModel,
                                 const char *OutputGCodePath,
                                 IrtRType NumViews,
                                 IrtRType OutputITDType,
                                 IrtRType ITDLength,
                                 const char *fileName,
                                 int uniformCut)
{
    IPObjectStruct *Contour, *RuledSrf, *SimObj, *AllSimObjs,
        *Solid = NULL;
    IrtHmgnMatType
        *Views = NULL;
    IritPrsrHWCDataStruct HWCParams;
    char **gcodeFiles, gcodeFile[512], contourFile[512];
    int vi, i,
        gcodeCount = 0,
        result = 0;

    if (RawModel == NULL || OutputGCodePath == NULL || NumViews <= 0) {
        fprintf(stderr,
            "HWCGenerateGCodeFromObj: invalid parameters.\n");
        return 0;
    }

    gcodeFiles = NULL;
    AllSimObjs = IritPrsrGenLISTObject(NULL);

    /* Convert any FreeForm surfaces to polygonal mesh. */
    Solid = IritPrsrConvertFreeFormHierachy(RawModel, &IritPrsrFFCState,
        FALSE, FALSE);
    RawModel = NULL;

    Solid = IritPrsrFlattenTree(Solid);
#ifdef DEBUG
    printf("Model ready. Type = %d.\n", Solid -> ObjType);
    fflush(stdout);
#endif /* DEBUG */

    /* Save a copy as OBJ for visual inspection. */
#ifdef DEBUG
    printf("Saving solid.obj...\n"); fflush(stdout);
    if (IritPrsrOBJSaveFile(Solid, "solid.obj", FALSE, 2, 0))
        printf("Wrote solid.obj.\n");
    fflush(stdout);
#endif /* DEBUG */

    /* 2. Select best view directions. */
#ifdef DEBUG
    printf("Selecting %d best view directions...\n", NumViews);
    fflush(stdout);
#endif /* DEBUG */

    if (uniformCut) 
        Views = HWCSelectViewsSampling(NumViews);
    else 
        Views = HWCSelectBestViewSampling(Solid, NumViews, NULL);
    
    if (Views == NULL) {
        fprintf(stderr, "View sampling failed.\n");

        /* Clean up. */
        if (gcodeFiles != NULL) {
            for (i = 0; i < gcodeCount; ++i) {
                if (gcodeFiles[i] != NULL) {
                    IritFree(gcodeFiles[i]);
                }
            }
            IritFree(gcodeFiles);
        }

        if (Solid != NULL)
            IritPrsrFreeObject(Solid);

        if (Views != NULL) 
            IritFree(Views);
    }
#ifdef DEBUG 
    printf("View sampling done.\n");
    fflush(stdout);
#endif /* DEBUG */

    /* 3. Set up HWC machine parameters. */
    IritPrsrHWCSetDfltParams(&HWCParams);
    HWCParams.MinimalHeight = 0.1;
    HWCParams.RuledApproxDir = CAGD_CONST_U_DIR;
    HWCParams.PieceWiseRuledApproximation = 0.01;
    if (ITDLength > 0.0) 
        HWCParams.SimulationArcLength = ITDLength;

    /* 4. Allocate temporary GCode file list. */
    gcodeFiles = (char** )IritMalloc(sizeof(char*) * NumViews);
    if (gcodeFiles == NULL) {
        fprintf(stderr, "Failed to allocate GCode file list.\n");

        /* Clean up. */
        if (gcodeFiles != NULL) {
            for (i = 0; i < gcodeCount; ++i) {
                if (gcodeFiles[i] != NULL) {
                    IritFree(gcodeFiles[i]);
                }
            }
            IritFree(gcodeFiles);
        }

        if (Solid != NULL) {
            IritPrsrFreeObject(Solid);
        }

        if (Views != NULL) {
            IritFree(Views);
        }
    }

    /* 5. Per-view: contour -> ruled surface -> GCode. */
    for (vi = 0; vi < NumViews; ++vi) {
        Contour = NULL;
        RuledSrf = NULL;
        SimObj = NULL;

        /* 5a. Approximate silhouette contour. */
#ifdef DEBUG
        printf("  View %d: Computing silhouette contour...\n", vi);
        fflush(stdout);
#endif /* DEBUG */

        Contour = HWCIritPrsrApproxBSplineContourFromSolidView(Solid,
            Views[vi],
            128);

        if (Contour == NULL) {
            printf("View %d: contour approximation failed, skipping.\n", vi);
            continue;
        }

        /* Save the raw 2D contour polygon for inspection. */
#ifdef DEBUG
        snprintf(contourFile, sizeof(contourFile),
            "contour_view%d.obj", vi);
        fflush(stdout);
        if (!IritPrsrOBJSaveFile(Contour, contourFile, FALSE, 0, 0))
            fprintf(stderr, "Failed saving %s.\n", contourFile);
        else
            printf("Wrote %s\n", contourFile);
#endif /* DEBUG */

        /* 5b. Build ruled surface. */
#ifdef DEBUG
        printf("  View %d: Building ruled surface...\n", vi);
        fflush(stdout);
#endif /* DEBUG */
        RuledSrf = HWCBuildSilhouetteRuledSrf(Contour, Solid, Views[vi],
            &HWCParams);
        IritPrsrFreeObject(Contour);
        Contour = NULL;

        if (RuledSrf == NULL) {
            printf("View %d: failed to build ruled surface, skipping.\n",
                vi);
            continue;
        }

        /* 5c. Generate GCode for this view. */
#ifdef DEBUG
        printf("  View %d: Generating GCode...\n", vi);
        fflush(stdout);
#endif
        snprintf(gcodeFile, sizeof(gcodeFile),
            "view%d_silhouette.gcode", vi);

        SimObj = IritPrsrHWCCreatePath(RuledSrf, NULL, NULL, gcodeFile,
            &HWCParams);
        IritPrsrFreeObject(RuledSrf);
        RuledSrf = NULL;

        if (SimObj != NULL) {
#ifdef DEBUG
            printf("View %d: GCode written to %s\n", vi, gcodeFile);
#endif /* DEBUG */

            /* Append to AllSimObjs to combine all paths. */
            IritPrsrListObjectAppend(AllSimObjs, SimObj);
            SimObj = NULL;

            /* Store filename for combiner. */
            gcodeFiles[gcodeCount] = IritMiscStrdup(gcodeFile);
            ++gcodeCount;
        }
        else {
            printf("View %d: IritPrsrHWCCreatePath failed, skipping.\n",
                vi);
        }
    }

    /* 6. Combine all per-view GCode files. */
    if (gcodeCount > 0) {
#ifdef DEBUG
        printf("Combining %d GCode files...\n", gcodeCount);
        fflush(stdout);
#endif /* DEBUG */
        HWCCombineGCodeFiles((const char const**)gcodeFiles, gcodeCount,
            OutputGCodePath, &HWCParams);
        result = 1;
    }
    else {
        fprintf(stderr, "No per-view GCode was produced.\n");
        result = 0;
    }

    /* Save all combined paths to a single ITD file. */
    if (AllSimObjs != NULL) {
        if (OutputITDType > 0) {
            if (OutputITDType == 2) {
                HWCConvertPolylinesToPolygons(AllSimObjs);
            }
            IritPrsrPutObjectToFile3(fileName, AllSimObjs, 0);
#ifdef DEBUG
            printf("All path geometries written to %s\n", fileName);
#endif /* DEBUG */
        }

        IritPrsrFreeObject(AllSimObjs);
        AllSimObjs = NULL;
    }

    /* 7. Clean up intermediate filenames and delete files in release mode. */
    if (gcodeFiles != NULL) {
        for (i = 0; i < gcodeCount; ++i) {
            if (gcodeFiles[i] != NULL) {
#ifndef DEBUG
                remove(gcodeFiles[i]);
#endif /* DEBUG */
                IritFree(gcodeFiles[i]);
            }
        }
        IritFree(gcodeFiles);
    }


    return result;
}


