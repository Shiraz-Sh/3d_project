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
//#define DEBUG_HOT_WIRE_CUT_ALG2

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

static IrtRType HWCCalcPolygonArea(IritPrsrObjectStruct *PObj);

static IritPrsrObjectStruct *HWCBuildProjectUnionLocal(
                                                 IritPrsrObjectStruct *Solid,
                                                 IrtHmgnMatType MatLocal);

static IritPrsrObjectStruct *HWCIritPrsrApproxBSplineContourFromSolidView(
                                                  IritPrsrObjectStruct *Solid,
						  const IrtHmgnMatType ViewMat,
						  int NumCtrl);

static IrtHmgnMatType *HWCSelectBestViewSampling(IritPrsrObjectStruct *PObj,
                                                 int NumSamples,
                                                 IrtHmgnMatType *ResultMat);

static IrtHmgnMatType *HWCSelectViewsSampling(int NumSamples);

static IritPrsrHWCEdgeStruct *HWCFindBoundaryEdges(
					     const IritPrsrObjectStruct *Solid,
					     int *OutNumEdges);

static IrtRType HWCDistPointSegment2D(IrtPtType P,
                                      IrtPtType A,
                                      IrtPtType B);

static IritPrsrObjectStruct *HWCBuildSilhouetteRuledSrf(
					  const IritPrsrObjectStruct *Contour,
					  const IritPrsrObjectStruct *Solid,
					  const IrtHmgnMatType ViewMat,
					  const IritPrsrHWCDataStruct *Params);

static void HWCCombineGCodeFiles(const char *const *GcodeFiles,
                                 int NumFiles,
                                 const char *OutPath,
                                 const IritPrsrHWCDataStruct *Params);

static void HWCConvertPolylinesToPolygons(IritPrsrObjectStruct *PObj);

IrtRType HWCGenerateGCodeFromObj(IritPrsrObjectStruct *RawModel,
                                 const char *OutputGCodePath,
                                 int NumViews,
                                 int OutputITDType,
                                 int RulingExtent,
                                 const char *outputITDFileName,
                                 int uniformCut,
                                 IrtRType FineNess);

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
    IrtVecType Arb;
    IrtRType wLen, Ul;

    /* Extract forward vector from third *row* of matrix (Mat[2][0..2]). */
    /* The view/LookAt code in this file stores rotation in rows, so the */
    /* forward vector is on row index 2.                                 */
    w[0] = Mat[2][0];
    w[1] = Mat[2][1];
    w[2] = Mat[2][2];

    wLen = sqrt(IRIT_SQR(w[0]) + IRIT_SQR(w[1]) + IRIT_SQR(w[2]));
    if (wLen <= IRIT_EPS) {
        if (!AllowDefault)
            return 0;
        /* Fallback to canonical forward if allowed. */
        w[0] = 0.0;
        w[1] = 0.0;
        w[2] = 1.0;
    }
    else {
        w[0] /= wLen;
        w[1] /= wLen;
        w[2] /= wLen;
    }

    /* Choose stable arbitrary vector not parallel to w. */
    if (fabs(w[0]) < 0.9) {
        Arb[0] = 1.0;
        Arb[1] = 0.0;
        Arb[2] = 0.0;
    }
    else {
        Arb[0] = 0.0;
        Arb[1] = 1.0;
        Arb[2] = 0.0;
    }

    /* u = normalize(cross(arb, w)). */
    u[0] = Arb[1] * w[2] - Arb[2] * w[1];
    u[1] = Arb[2] * w[0] - Arb[0] * w[2];
    u[2] = Arb[0] * w[1] - Arb[1] * w[0];

    Ul = sqrt(IRIT_SQR(u[0]) + IRIT_SQR(u[1]) + IRIT_SQR(u[2]));
    if (Ul <= IRIT_EPS) {
        /* Extremely unlikely; if allowed try a different arb, else fail. */
        if (!AllowDefault)
            return 0;
        /* Try fallback arb. */
        Arb[0] = 0.0;
        Arb[1] = 1.0;
        Arb[2] = 0.0;
        u[0] = Arb[1] * w[2] - Arb[2] * w[1];
        u[1] = Arb[2] * w[0] - Arb[0] * w[2];
        u[2] = Arb[0] * w[1] - Arb[1] * w[0];
        Ul = sqrt(IRIT_SQR(u[0]) + IRIT_SQR(u[1]) + IRIT_SQR(u[2]));
        if (Ul <= IRIT_EPS)
            return 0;
    }

    u[0] /= Ul;
    u[1] /= Ul;
    u[2] /= Ul;

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
static IrtRType HWCCalcPolygonArea(IritPrsrObjectStruct *PObj)
{
    IrtRType
        Area = 0.0;
    IritPrsrPolygonStruct *Pl;

    if (PObj == NULL || !IRIT_PRSR_IS_POLY_OBJ(PObj))
        return 0.0;

    for (Pl = PObj -> U.Pl; Pl != NULL; Pl = Pl -> Pnext) {
        /* Uses IRIT's built-in area function. */
        Area += IritGeomPolyOnePolyArea(Pl, TRUE);
    }
    return Area;
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
* Solid: IritPrsrObjectStruct pointer to valid polygonal object (must pass   *
*        IRIT_PRSR_IS_POLY_OBJ check). Multiple loops allowed. Input not     *
*        modified. NULL returns NULL.                                        *
* MatLocal: 4x4 transformation matrix mapping world coordinates to           *
*           view-local frame. Typically built from orthonormal basis         *
*           (u=right, v=up, w=forward). Singular matrices return NULL.       *
*                                                                            *
* RETURN VALUE:                                                              *
* IritPrsrObjectStruct *: LIST object containing single POLY with boundary   *
*               contour in view-local XY coordinates (Z=0), or NULL on       *
*               failure (NULL input, degenerate projection, allocation       *
*               failure,  or empty boundary). Caller owns returned object    *
*               and must free via IritPrsrFreeObject().                      *
*                                                                            *
* ALGORITHM:                                                                 *
* 1. Project all polygon vertices to local XY, compute bounding box with 5%  *
*    margin to prevent boundary clipping.                                    *
* 2. Allocate SIL_GRID_RES x SIL_GRID_RES (1024x1024) binary raster grid.    *
* 3. For each input polygon, execute scanline fill:                          *
*    - Map polygon vertices to grid coordinates.                             *
*    - For each scanline intersecting polygon Y-range, compute               *
*      X-intersections.                                                      *
*    - Fill horizontal spans between paired intersections (parity rule).     *
* 4. Find first filled pixel (leftmost, then topmost).                       *
* 5. Moore neighborhood contour tracing: starting from filled pixel, follow  *
*    filled-to-empty boundary using 8-neighbor search, always turning 135    *
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
static IritPrsrObjectStruct *HWCBuildProjectUnionLocal(
                                                  IritPrsrObjectStruct *Solid,
						  IrtHmgnMatType MatLocal)
{
    IrtRType Dx, Dy, Margin,
        minX = IRIT_INFNTY,
        minY = IRIT_INFNTY,
        maxX = -IRIT_INFNTY,
        maxY = -IRIT_INFNTY;
    unsigned char *Grid;
    int *Px, *Py, *Ints, x, y, i, j,
        startX = -1,
        startY = -1;
    IritPrsrPolygonStruct *Pl;

    if (Solid == NULL || !IRIT_PRSR_IS_POLY_OBJ(Solid))
        return NULL;

    /* 1. Project points to local XY, collect them to compute bounding box. */
    for (Pl = Solid -> U.Pl; Pl != NULL; Pl = Pl -> Pnext) {
        IritPrsrVertexStruct *Cur,
            * V = Pl->PVertex;
        int Guard = 0;

        if (V == NULL)
            continue;
        Cur = V;
        do {
            IrtPtType localP;

            IritMiscMatMultPtby4by4(localP, Cur->Coord, MatLocal);
            if (localP[0] < minX)
                minX = localP[0];
            if (localP[1] < minY)
                minY = localP[1];
            if (localP[0] > maxX)
                maxX = localP[0];
            if (localP[1] > maxY)
                maxY = localP[1];

            Cur = Cur->Pnext;
            if (++Guard > 200000)
                break;
        } while (Cur != NULL && Cur != V);
    }

    if (minX > maxX)
        return NULL;

    /* Add a small 5% margin to prevent tracing out of bounds. */
    Dx = maxX - minX;
    Dy = maxY - minY;
    Margin = IRIT_MAX(Dx, Dy) * 0.05;
    if (Margin < 1e-5)
        Margin = 1e-5;
    minX -= Margin;
    maxX += Margin;
    minY -= Margin;
    maxY += Margin;
    Dx = maxX - minX;
    Dy = maxY - minY;

    /* 2. Allocate and clear 2D grid. */
    Grid = (unsigned char *) IritMalloc(HWCSIL_GRID_RES * HWCSIL_GRID_RES);
    memset(Grid, 0, HWCSIL_GRID_RES * HWCSIL_GRID_RES);

    /* Assume max polygon vertex count <= 4096. */
    Px = (int*) IritMalloc(sizeof(int) * 4096);
    Py = (int*) IritMalloc(sizeof(int) * 4096);
    Ints = (int*) IritMalloc(sizeof(int) * 4096);

    /* 3. Scanline fill each projected polygon. */
    for (Pl = Solid -> U.Pl; Pl != NULL; Pl = Pl -> Pnext) {
        IritPrsrVertexStruct *Cur,
            *V = Pl -> PVertex;
        int pminY, pmaxY,
            nPts = 0,
            Guard = 0;

        if (V == NULL)
            continue;

        Cur = V;
        do {
            IrtPtType localP;
            int PxVal, PyVal;

            IritMiscMatMultPtby4by4(localP, Cur -> Coord, MatLocal);

            PxVal = (int) ((localP[0] - minX) / Dx * (HWCSIL_GRID_RES - 1));
            PyVal = (int) ((localP[1] - minY) / Dy * (HWCSIL_GRID_RES - 1));
            if (PxVal < 0)
                PxVal = 0;
            if (PxVal >= HWCSIL_GRID_RES)
                PxVal = HWCSIL_GRID_RES - 1;
            if (PyVal < 0)
                PyVal = 0;
            if (PyVal >= HWCSIL_GRID_RES)
                PyVal = HWCSIL_GRID_RES - 1;

            if (nPts < 4096) {
                Px[nPts] = PxVal;
                Py[nPts] = PyVal;
                nPts++;
            }
            Cur = Cur -> Pnext;
            if (++Guard > 200000)
                break;
        } while (Cur != NULL && Cur != V);

        if (nPts < 3)
            continue;

        /* Scanline logic */
        pminY = HWCSIL_GRID_RES;
        pmaxY = -1;
        for (i = 0; i < nPts; ++i) {
            if (Py[i] < pminY)
                pminY = Py[i];
            if (Py[i] > pmaxY)
                pmaxY = Py[i];
        }

        for (y = pminY; y <= pmaxY; ++y) {
            int a, b,
                count = 0;

            for (i = 0; i < nPts; ++i) {
                int jIDx = (i + 1) % nPts,
                    y1 = Py[i],
                    y2 = Py[jIDx],
                    x1 = Px[i],
                    x2 = Px[jIDx];

                if (y1 == y2)
                    continue;

                if ((y1 <= y && y < y2) || (y2 <= y && y < y1)) {
                    int xInt = x1 + (y - y1) * (x2 - x1) / (y2 - y1);
                    if (count < 4096) {
                        Ints[count++] = xInt;
                    }
                }
            }

            /* Sort intersections. */
            for (a = 0; a < count - 1; ++a) {
                for (b = a + 1; b < count; ++b) {
                    if (Ints[a] > Ints[b]) {
                        int tmp = Ints[a];

                        Ints[a] = Ints[b];
                        Ints[b] = tmp;
                    }
                }
            }

            /* Fill pairs. */
            for (a = 0; a < count; a += 2) {
                if (a + 1 < count) {
                    int xStart = Ints[a];
                    int xEnd = Ints[a + 1];
                    int xVal;

                    if (xStart < 0)
                        xStart = 0;
                    if (xEnd >= HWCSIL_GRID_RES)
                        xEnd = HWCSIL_GRID_RES - 1;
                    for (xVal = xStart; xVal <= xEnd; ++xVal) {
                        Grid[y * HWCSIL_GRID_RES + xVal] = 1;
                    }
                }
            }
        }
    }

    IritFree(Px);
    IritFree(Py);
    IritFree(Ints);

    /* 4. Find starting pixel for Moore Neighborhood Tracing. */
    for (y = 0; y < HWCSIL_GRID_RES; ++y) {
        for (x = 0; x < HWCSIL_GRID_RES; ++x) {
            if (Grid[y * HWCSIL_GRID_RES + x]) {
                startX = x;
                startY = y;
                break;
            }
        }
        if (startX != -1) break;
    }

    if (startX == -1) {
        IritFree(Grid);
        return NULL;
    }

    /* 5. Moore Neighborhood Tracing. */
    {
        int xDx[8] = { 0, 1, 1, 1, 0, -1, -1, -1 },
            xDy[8] = { -1, -1, 0, 1, 1, 1, 0, -1 },
            maxPath = HWCSIL_GRID_RES * HWCSIL_GRID_RES,
            pathLen = 0,
            Cx = startX,
            Cy = startY,
            firstStep = 1,
            Dir = 2,   /* Pretend we moved East (2), so left is North (0). */
            *pathX = (int*) IritMalloc(sizeof(int) * maxPath),
            *pathY = (int*) IritMalloc(sizeof(int) * maxPath);

        do {
            int nextDir, i,
                Found = 0;

            if (pathLen < maxPath) {
                pathX[pathLen] = Cx;
                pathY[pathLen] = Cy;
                pathLen++;
            }
            else break;

            for (i = 0; i < 8; ++i) {
                int Nx, Ny; 

                /* Turn 135 deg left to ensure we trace outside boundary. */
                nextDir = (Dir + 5 + i) % 8;
                Nx = Cx + xDx[nextDir];
                Ny = Cy + xDy[nextDir];

                if (Nx >= 0 &&
		    Nx < HWCSIL_GRID_RES && Ny >= 0 &&
		    Ny < HWCSIL_GRID_RES) {
                    if (Grid[Ny * HWCSIL_GRID_RES + Nx]) {
                        Cx = Nx;
                        Cy = Ny;
                        Dir = nextDir;
                        Found = 1;
                        break;
                    }
                }
            }

            if (!Found) break; /* Isolated pixel. */

            if (!firstStep && Cx == startX && Cy == startY) {
                break;
            }
            firstStep = 0;
        } while (1);

        IritFree(Grid);

        if (pathLen < 3) {
            IritFree(pathX);
            IritFree(pathY);
            return NULL;
        }

        /* 6. Convert traced pixel path back to view-local XY. */
        {
            IritPrsrPolygonStruct
                *NewPl = IritPrsrAllocPolygon(0, NULL, NULL);
            IritPrsrVertexStruct
                *FirstV = NULL,
                *PrevV = NULL;
            int iPath;
            IritPrsrObjectStruct *PolyObj,
                *OutList;

            for (iPath = 0; iPath < pathLen; ++iPath) {
                IritPrsrVertexStruct
                    *NV = IritPrsrAllocVertex2(NULL);

                NV -> Coord[0] = minX +
		                 pathX[iPath] * Dx / (HWCSIL_GRID_RES - 1);
                NV -> Coord[1] = minY +
		                 pathY[iPath] * Dy / (HWCSIL_GRID_RES - 1);
                NV -> Coord[2] = 0.0;
                if (FirstV == NULL)
                    FirstV = PrevV = NV;
                else {
                    PrevV -> Pnext = NV; PrevV = NV;
                }
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
            IRIT_PRSR_SET_POLYGON_OBJ(PolyObj);
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
* 3: Rasterizes solid into 1024×1024 grid, extracts boundary via Moore   *
*    neighborhood tracing, producing high-resolution pixel-level contour.    *
* 4: Applies Laplacian smoothing to remove pixel stair-stepping artifacts.   *
* 5: Subsamples points to requested count (NumCtrl).                         *
* 6: Transforms back to world coordinates using inverse view matrix.         *
* 7: Returns closed polygon in world frame, ready for further processing.    *
*                                                                            *
* PARAMETERS:                                                                *
* Solid: IritPrsrObjectStruct pointer to valid polygonal object (checked     *
*        via IRIT_PRSR_IS_POLY_OBJ). Multiple polygons/loops allowed.        *
*        Input not modified. NULL or non-polygon returns NULL.               *
* ViewMat: 4x4 homogeneous view matrix. Row 2 (Mat[2][0..2])                 *
*          interpreted as forward direction vector (normalized               *
*          internally). Defines camera frame via orthonormal basis.          *
* NumCtrl: Target number of control points in final contour (positive        *
*          integer). Actual count may vary due to subsampling step           *
*          calculation. Larger values preserve detail; typical 32-128.       *
*                                                                            *
* RETURN VALUE:                                                              *
* IritPrsrObjectStruct *: Newly allocated POLY object in world coordinates   *
*               with closed contour polygon, or NULL on failure (NULL input, *
*               invalid matrix, degenerate projection, allocation            *
*               failure or insufficient boundary points). Caller owns        *
*               returned object and must free via IritPrsrFreeObject().      *
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
static IritPrsrObjectStruct *HWCIritPrsrApproxBSplineContourFromSolidView(
                                                 IritPrsrObjectStruct *Solid,
                                                 const IrtHmgnMatType ViewMat,
                                                 int NumCtrl)
{
    IrtVecType u, v, w;
    IrtHmgnMatType MatLocal, InvMat;
    IritPrsrObjectStruct *UnionLocal, *Poly, *FinalObj;
    IritPrsrPolygonStruct *Pl, *FinalPl;
    IritPrsrVertexStruct *V, *Cur,
        *FirstV = NULL,
        *PrevV = NULL;
    IrtPtType *Pts, *tmpPts;
    int i, It, Step,
        nPts = 0,
        smoothIters = 5;

    if (Solid == NULL || !IRIT_PRSR_IS_POLY_OBJ(Solid))
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

    /* Build union of projected faces in local XY using robust grid         */
    /* rasterization.							    */
    UnionLocal = HWCBuildProjectUnionLocal(Solid, InvMat);
    if (UnionLocal == NULL) 
        return NULL;

    /* UnionLocal contains a single high-resolution pixel-traced boundary.  */
    /* Extract it, smooth it lightly to remove pixel steps, and return in   */
    /* world coords.							    */
    Poly = IritPrsrListObjectGet(UnionLocal, 0);
    if (Poly == NULL || !IRIT_PRSR_IS_POLY_OBJ(Poly) || Poly -> U.Pl == NULL) {
        IritPrsrFreeObject(UnionLocal);
        return NULL;
    }

    Pl = Poly -> U.Pl;
    V = Pl -> PVertex;
    Cur = V;
    do {
        nPts++;
        Cur = Cur -> Pnext;
    } while (Cur != NULL && Cur != V);

    Pts = (IrtPtType *) IritMalloc(sizeof(IrtPtType) * nPts);
    tmpPts = (IrtPtType *) IritMalloc(sizeof(IrtPtType) * nPts);

    Cur = V;
    for (i = 0; i < nPts; ++i) {
        Pts[i][0] = Cur -> Coord[0];
        Pts[i][1] = Cur -> Coord[1];
        Pts[i][2] = 0.0;
        Cur = Cur -> Pnext;
    }

    /* Laplacian smoothing to remove pixel stair-stepping artifacts. */
    for (It = 0; It < smoothIters; ++It) {
        for (i = 0; i < nPts; ++i) {
            int Im1 = (i - 1 + nPts) % nPts,
                Ip1 = (i + 1) % nPts;

            tmpPts[i][0] = (Pts[Im1][0] + Pts[Ip1][0] + Pts[i][0]) / 3.0;
            tmpPts[i][1] = (Pts[Im1][1] + Pts[Ip1][1] + Pts[i][1]) / 3.0;
            tmpPts[i][2] = 0.0;
        }
        memcpy(Pts, tmpPts, sizeof(IrtPtType) * nPts);
    }

    /* Subsample points. */
    Step = nPts / NumCtrl;
    if (Step < 1) Step = 1;

    FinalPl = IritPrsrAllocPolygon(0, NULL, NULL);

    for (i = 0; i < nPts; i += Step) {
        IrtPtType worldP;
        IritPrsrVertexStruct *NV;

        IritMiscMatMultPtby4by4(worldP, Pts[i], MatLocal);
        NV = IritPrsrAllocVertex2(NULL);
        NV -> Coord[0] = worldP[0];
        NV -> Coord[1] = worldP[1];
        NV -> Coord[2] = worldP[2];
        if (FirstV == NULL)
            FirstV = PrevV = NV;
        else
            PrevV -> Pnext = NV; PrevV = NV;
    }

    IritFree(Pts);
    IritFree(tmpPts);
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
        IRIT_PRSR_SET_POLYGON_OBJ(FinalObj);
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
* [0, 2π).                                                               *
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
* PObj: IritPrsrObjectStruct pointer to a polygonal object. Must satisfy     *
*       IRIT_PRSR_IS_POLY_OBJ. The object is not modified.                   *
*                           NULL returns NULL.                               *
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
static IrtHmgnMatType *HWCSelectBestViewSampling(IritPrsrObjectStruct *PObj,
                                                 int NumSamples,
                                                 IrtHmgnMatType *ResultMat)
{
    int Next, s, i, j, Best;
    IrtHmgnMatType *ViewMats;
    IrtRType *Scores;
    const int 
        NumCtrl = 128;

    if (PObj == NULL || NumSamples <= 0)
        return NULL;

    NumSamples = NumSamples * 2;

    ViewMats = (IrtHmgnMatType *)
			      IritMalloc(sizeof(IrtHmgnMatType) * NumSamples);
    if (ViewMats == NULL)
        return NULL;

    Next = 0;

    /* Generate candidate view matrices (Uniform circle in XY plane). */
    for (s = 0; s < NumSamples; ++s) {
        IrtRType
            Theta = 2.0 * 3.14159265358979323846 * s / NumSamples;
        IrtRType r = 2.0;
        IrtVecType CamPos, Center, Up;

        CamPos[0] = cos(Theta) * r;
        CamPos[1] = sin(Theta) * r;
        CamPos[2] = 0.0;          /* Fixed Z value, no variation in height. */

        Center[0] = 0.0;
        Center[1] = 0.0;
        Center[2] = 0.0;

        Up[0] = 0.0;
        Up[1] = 0.0;
        Up[2] = 1.0;

        HWCGenLookAtMatrix(CamPos, Center, Up, ViewMats[Next++]);
    }

    Scores = (IrtRType*) IritMalloc(sizeof(IrtRType) * NumSamples);
    if (Scores == NULL) {
        IritFree(ViewMats);
        return NULL;
    }

    for (i = 0; i < NumSamples; ++i) {
        IritPrsrObjectStruct *Contour;
        IrtRType wxyLen;

        Scores[i] = 0.0;

        /* Reject near-vertical views (where wxyLen < 0.1) since the       */
	/* horizontal wire cannot reach the required vertical slope.        */
        wxyLen = sqrt(ViewMats[i][2][0] * ViewMats[i][2][0] +
		       ViewMats[i][2][1] * ViewMats[i][2][1]);
        if (wxyLen < 0.1) {
            continue;
        }

        Contour = HWCIritPrsrApproxBSplineContourFromSolidView(PObj,
							       ViewMats[i],
							       NumCtrl);
        if (Contour != NULL) {
            Scores[i] = HWCCalcPolygonArea(Contour);
            IritPrsrFreeObject(Contour);
        }
    }

    /* Reorder ViewMats in-place by descending score (simple selection      */
    /* sort).							            */
    for (i = 0; i < NumSamples - 1; ++i) {
        Best = i;
        for (j = i + 1; j < NumSamples; ++j) {
            if (Scores[j] > Scores[Best])
                Best = j;
        }
        if (Best != i) {
            IrtHmgnMatType tmpMat;
            IrtRType tmpScore = Scores[i];

            memcpy(tmpMat, ViewMats[i], sizeof(IrtHmgnMatType));
            memcpy(ViewMats[i], ViewMats[Best], sizeof(IrtHmgnMatType));
            memcpy(ViewMats[Best], tmpMat, sizeof(IrtHmgnMatType));

            Scores[i] = Scores[Best];
            Scores[Best] = tmpScore;
        }
    }

    /* Build final view list: first 5 are fixed cardinal + diagonal views,  */
    /* then fill remaining slots with best diverse scored views.            */
    {
        int RequestedViews = NumSamples / 2,
            numSelected = 0;
        IrtHmgnMatType
            *SelectedMats = (IrtHmgnMatType *)
                               IritMalloc(sizeof(IrtHmgnMatType) * NumSamples);

        /* 0 = X, 45 deg = diagonal, 90 deg =+Y, 135 deg = other diagonal. */
        IrtRType fixedAngles[4] = { 0.0, 45.0, 90.0, 135.0 };

        for (i = 0; i < 4 && numSelected < RequestedViews; ++i) {
            IrtRType Theta = fixedAngles[i] * 3.14159265358979323846 / 180.0;
            IrtVecType CamPos, Center, Up;

            CamPos[0] = cos(Theta) * 2.0;
            CamPos[1] = sin(Theta) * 2.0;
            CamPos[2] = 0.0;
            Center[0] = 0.0; Center[1] = 0.0; Center[2] = 0.0;
            Up[0] = 0.0; Up[1] = 0.0; Up[2] = 1.0;
            HWCGenLookAtMatrix(CamPos, Center, Up,
                SelectedMats[numSelected++]);
        }

        /* Fill remaining slots with best-scored views that are not too     */
	/* similar (within 15 degrees) to any already-selected view.	    */
        for (i = 0; i < NumSamples && numSelected < RequestedViews; ++i) {
            int isSimilar = 0;

            for (j = 0; j < numSelected; ++j) {
                IrtRType Dot =
                    ViewMats[i][2][0] * SelectedMats[j][2][0] +
                    ViewMats[i][2][1] * SelectedMats[j][2][1] +
                    ViewMats[i][2][2] * SelectedMats[j][2][2];

                if (IRIT_FABS(Dot) >
                    cos(15.0 * 3.14159265358979323846 / 180.0)) {
                    isSimilar = 1;
                    break;
                }
            }
            if (!isSimilar) {
                memcpy(SelectedMats[numSelected++], ViewMats[i],
                    sizeof(IrtHmgnMatType));
            }
        }

        /* If still not enough, relax the threshold and accept any          */
	/* non-duplicate view.						    */
        for (i = 0; i < NumSamples && numSelected < RequestedViews; ++i) {
            int isExact = 0;

            for (j = 0; j < numSelected; ++j) {
                IrtRType 
                    Dot = ViewMats[i][2][0] * SelectedMats[j][2][0] +
                          ViewMats[i][2][1] * SelectedMats[j][2][1] +
                          ViewMats[i][2][2] * SelectedMats[j][2][2];

                if (IRIT_FABS(Dot) > 0.999) {
                    isExact = 1;
                    break;
                }
            }
            if (!isExact) {
                memcpy(SelectedMats[numSelected++], ViewMats[i],
                    sizeof(IrtHmgnMatType));
            }
        }

        memcpy(ViewMats, SelectedMats, sizeof(IrtHmgnMatType) * numSelected);
        IritFree(SelectedMats);
    }

    if (ResultMat != NULL) {
        memcpy(ResultMat, ViewMats, sizeof(IrtHmgnMatType) * NumSamples);
    }

    IritFree(Scores);
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
    IrtRType Step;
    IrtHmgnMatType *ViewMats;

    if (NumSamples <= 0)
        return NULL;

    ViewMats = (IrtHmgnMatType *)
			      IritMalloc(sizeof(IrtHmgnMatType) * NumSamples);
    if (ViewMats == NULL)
        return NULL;

    /* Space views by 180/N degrees so that no two views are mirrors of     */
    /* each other (a view at angle A covers the same silhouette as A + 180).*/
    Step = 180.0 / NumSamples;

    for (i = 0; i < NumSamples; ++i) {
        IrtRType
	    Theta = (Step * i) * 3.14159265358979323846 / 180.0;
        IrtVecType CamPos, Center, Up;

        CamPos[0] = cos(Theta) * 2.0;
        CamPos[1] = sin(Theta) * 2.0;
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
* AUXILIARY:								     *
* Auxiliary function to find missing boundary edges in an open solid.	     *
*****************************************************************************/
static IritPrsrHWCEdgeStruct *HWCFindBoundaryEdges(
					    const IritPrsrObjectStruct *Solid,
					    int *OutNumEdges)
{
    int i, j, k, StackTop, CurrentLoopCount, *EdgeCounts, *Visited, *Stack,
        *CurrentLoopIndices,
        MaxEdges = 0,
        NumEdges = 0,
        NumRawBoundaryEdges = 0,
        BestLoopCount = 0;
    IrtRType 
        CurrentLoopLength, CurrentLoopAvgZ,
        BestLoopLength = -1.0,
        BestLoopAvgZ = IRIT_INFNTY,
        EPS = 1e-4;
    IritPrsrPolygonStruct *Pl;
    IritPrsrVertexStruct *V, *VNext;
    IritPrsrHWCEdgeStruct *AllEdges,
        *RawBoundaryEdges = NULL,
        *BestLoopEdges = NULL;

    if (OutNumEdges)
        *OutNumEdges = 0;
    if (!Solid || !IRIT_PRSR_IS_POLY_OBJ(Solid))
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

    AllEdges = (IritPrsrHWCEdgeStruct *)
                        IritMalloc(sizeof(IritPrsrHWCEdgeStruct) * MaxEdges);
    EdgeCounts = (int *) IritMalloc(sizeof(int) * MaxEdges);

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
                int Match11 = 1,
		    Match12 = 1,
		    Match21 = 1,
		    Match22 = 1;

                if (Visited[j])
                    continue;

                for (k = 0; k < 3; k++) {
                    if (IRIT_FABS(RawBoundaryEdges[Curr].Pt1[k] -
                        RawBoundaryEdges[j].Pt1[k]) > EPS)
		        Match11 = 0;
                    if (IRIT_FABS(RawBoundaryEdges[Curr].Pt1[k] -
                        RawBoundaryEdges[j].Pt2[k]) > EPS)
		        Match12 = 0;
                    if (IRIT_FABS(RawBoundaryEdges[Curr].Pt2[k] -
                        RawBoundaryEdges[j].Pt1[k]) > EPS)
		        Match21 = 0;
                    if (IRIT_FABS(RawBoundaryEdges[Curr].Pt2[k] -
                        RawBoundaryEdges[j].Pt2[k]) > EPS)
		        Match22 = 0;
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
            int Ei = CurrentLoopIndices[j];

            CurrentLoopAvgZ += (RawBoundaryEdges[Ei].Pt1[2] +
                                RawBoundaryEdges[Ei].Pt2[2]) * 0.5;
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
* AUXILIARY:							             *
* Auxiliary function to calculate point to segment distance in 2D.           *
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
*   5. Tag with IRIT_ATTR_ID_Dir = IRIT_CAGD_CONST_U_DIR                     *
*                            so IritPrsrHWCCreatePath                        *
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
* IritPrsrObjectStruct*: Surface on success, NULL on failure.                *
*         Caller owns the returned object and must free it.                  *
*****************************************************************************/
static IritPrsrObjectStruct *HWCBuildSilhouetteRuledSrf(
                                          const IritPrsrObjectStruct *Contour,
					  const IritPrsrObjectStruct *Solid,
					  const IrtHmgnMatType ViewMat,
					  const IritPrsrHWCDataStruct *Params)
{
    int k, n, Guard,
        bestStart = 0,
        bestLen = 0,
        numBoundaryEdges = 0,
        topStart = 0,
        topLen = 0;
    IrtRType minX, minY, maxX, maxY, scaleX, scaleZ, Scale, Cx, Cy;
    IrtVecType uVec, vVec, wVec;
    IrtHmgnMatType MatLocal;
    IrtPtType *Pts2d;
    IritPrsrVertexStruct *V, *Cur;
    IritCagdCrvStruct *Crv;
    IritCagdVecStruct YDir;
    IritCagdSrfStruct *Srf;
    IritPrsrObjectStruct *SrfObj;

    if (Contour == NULL ||
	!IRIT_PRSR_IS_POLY_OBJ((IritPrsrObjectStruct *) Contour) ||
        Contour -> U.Pl == NULL || Params == NULL)
        return NULL;

    /* 1. Build view-local orthonormal basis				    */
    /* (uVec=right, vVec=up, wVec=fwd).				    */
    if (!HWCIritPrsrBuildViewBasisFromMat(ViewMat, uVec, vVec, wVec, 1))
        return NULL;

    /* MatLocal maps world  ->  view-local (rows = basis vectors). */
    IritMiscMatGenUnitMat(MatLocal);
    MatLocal[0][0] = uVec[0];
    MatLocal[0][1] = uVec[1];
    MatLocal[0][2] = uVec[2];
    MatLocal[0][3] = 0.0;
    MatLocal[1][0] = vVec[0];
    MatLocal[1][1] = vVec[1];
    MatLocal[1][2] = vVec[2];
    MatLocal[1][3] = 0.0;
    MatLocal[2][0] = wVec[0];
    MatLocal[2][1] = wVec[1];
    MatLocal[2][2] = wVec[2];
    MatLocal[2][3] = 0.0;
    MatLocal[3][0] = 0.0;
    MatLocal[3][1] = 0.0;
    MatLocal[3][2] = 0.0;
    MatLocal[3][3] = 1.0;

    /* 2. Count and validate polygon vertices. */
    V = Contour -> U.Pl -> PVertex;
    if (V == NULL) return NULL;

    n = 0;
    Cur = V;
    Guard = 0;
    do {
        ++n;
        Cur = Cur -> Pnext;
        if (++Guard > 200000) { n = 0; break; }
    } while (Cur != NULL && Cur != V);

    if (n < 3) return NULL;

    topLen = n;

    /* 3. Project each vertex to view-local 2D, ignoring world Z.	    */
    /*   localP[0] = view-right component -> machine X			    */
    /*   localP[1] = view-up   component -> machine Z (vertical).	    */
    Pts2d = (IrtPtType*) IritMalloc(sizeof(IrtPtType) * n);
    Cur = V;
    Guard = 0;
    for (k = 0; k < n; ++k) {
        /* Compute true dot products to project world coordinates onto the  */
        /* view plane.							    */
        Pts2d[k][0] = Cur -> Coord[0] * uVec[0] + Cur -> Coord[1] * uVec[1]
            + Cur -> Coord[2] * uVec[2];
        Pts2d[k][1] = Cur -> Coord[0] * vVec[0] + Cur -> Coord[1] * vVec[1]
            + Cur -> Coord[2] * vVec[2];
        Pts2d[k][2] = 0.0;
        Cur = Cur -> Pnext;
        if (Cur == NULL || Cur == V) 
            break;
        if (++Guard > 200000) 
            break;
    }

    /* 4. Compute AABB. */
    minX = maxX = Pts2d[0][0];
    minY = maxY = Pts2d[0][1];
    for (k = 1; k < n; ++k) {
        if (Pts2d[k][0] < minX) 
            minX = Pts2d[k][0];
        if (Pts2d[k][0] > maxX) 
            maxX = Pts2d[k][0];
        if (Pts2d[k][1] < minY) 
            minY = Pts2d[k][1];
        if (Pts2d[k][1] > maxY) 
            maxY = Pts2d[k][1];
    }

    if ((maxX - minX) < IRIT_EPS || (maxY - minY) < IRIT_EPS) {
        IritFree(Pts2d);
        return NULL;
    }

    /* 5. Uniform scale: fit silhouette into 90% of FoamWidth x FoamHeight. */
    scaleX = (Params -> FoamWidth * 0.9) / (maxX - minX);
    scaleZ = (Params -> FoamHeight * 0.9) / (maxY - minY);
    Scale = IRIT_MIN(scaleX, scaleZ);
    Cx = 0.5 * (minX + maxX);
    Cy = 0.5 * (minY + maxY);

    /* 6. Validate the view is horizontally cuttable.			    */
    /*									    */
    /*  The hot-wire machine spans the wire along the horizontal Y axis     */
    /* between two clamps.						    */
    /*  Therefore it can only cleanly sweep along views where the           */
    /* horizontal length of the view vector is non-negligible. We reject    */
    /* completely vertical views.					    */
    IrtRType wxyLen = sqrt(wVec[0] * wVec[0] + wVec[1] * wVec[1]);

    if (wxyLen < 0.1) {
        IritFree(Pts2d);
        return NULL;
    }

    /* 7. Robust Bottom Line Removal based on exact distance to missing     */
    /*    edges.							    */
    /*    We project the 3D boundary edges to 2D line segments.		    */
    /*    Then we measure the distance from every point on the contour to   */
    /*    these segments.						    */
    /*    Points very close to ANY boundary segment are identified as part  */
    /*    of the missing bottom. We then keep the longest contiguous	    */
    /*    sequence of points that are NOT missing.			    */
    IritPrsrHWCEdgeStruct
        *boundaryEdges = HWCFindBoundaryEdges(Solid, &numBoundaryEdges);

    if (boundaryEdges != NULL && numBoundaryEdges > 0) {
        int b,
            *isBoundary = (int*) IritMalloc(sizeof(int) * n);
        IrtRType 
            Diag = sqrt((maxX - minX) * (maxX - minX) +
			(maxY - minY) * (maxY - minY)),
            distThresh = Diag * 0.005;   /* 0.5% of diagonal as tolerance. */

        for (k = 0; k < n; ++k) {
            isBoundary[k] = 0;
            for (b = 0; b < numBoundaryEdges; ++b) {
                IrtPtType projP1, projP2;

                projP1[0] = boundaryEdges[b].Pt1[0] * uVec[0] +
		             boundaryEdges[b].Pt1[1] * uVec[1] +
		             boundaryEdges[b].Pt1[2] * uVec[2];
                projP1[1] = boundaryEdges[b].Pt1[0] * vVec[0] + boundaryEdges[b].Pt1[1]
                    * vVec[1] + boundaryEdges[b].Pt1[2] * vVec[2];
                projP2[0] = boundaryEdges[b].Pt2[0] * uVec[0] + boundaryEdges[b].Pt2[1]
                    * uVec[1] + boundaryEdges[b].Pt2[2] * uVec[2];
                projP2[1] = boundaryEdges[b].Pt2[0] * vVec[0] + boundaryEdges[b].Pt2[1]
                    * vVec[1] + boundaryEdges[b].Pt2[2] * vVec[2];

                if (HWCDistPointSegment2D(Pts2d[k], projP1, projP2) < distThresh) {
                    isBoundary[k] = 1;
                    break;
                }
            }
        }
        IritFree(boundaryEdges);

        /* Find the longest contiguous sequence of isBoundary == 0 (wrapping around). */
        for (k = 0; k < n; ++k) {
            if (isBoundary[k] == 0) {
                int Len = 0;
                int iDx = k;
                while (Len < n && isBoundary[iDx] == 0) {
                    Len++;
                    iDx = (iDx + 1) % n;
                }
                if (Len > bestLen) {
                    bestLen = Len;
                    bestStart = k;
                }
                /* Fast forward k to the end of this non-boundary sequence. */
                /* But only up to n-1 to avoid infinite loops if the whole thing wraps. */
                k += Len - 1;
            }
        }

        IritFree(isBoundary);

        if (bestLen > 0) {
            topStart = bestStart;
            topLen = bestLen;
        }
    }
    else {
        /* Fallback if no boundary edges found: remove bottom 10%. */
        IrtRType
            bottomThresh = minY + (maxY - minY) * 0.10,
            minXBottom = maxX + 1.0,
            maxXBottom = minX - 1.0;
        int leftIDx = -1,
            rightIDx = -1;

        for (k = 0; k < n; ++k) {
            if (Pts2d[k][1] <= bottomThresh) {
                if (Pts2d[k][0] < minXBottom) {
                    minXBottom = Pts2d[k][0];
                    leftIDx = k;
                }
                if (Pts2d[k][0] > maxXBottom) {
                    maxXBottom = Pts2d[k][0];
                    rightIDx = k;
                }
            }
        }

        if (leftIDx != -1 && rightIDx != -1 && leftIDx != rightIDx &&
            (maxXBottom - minXBottom) > (maxX - minX) * 0.10) {

            /* Path 1: leftIdx to rightIdx (incrementing). */
            IrtRType 
                maxYPath1 = minY;
            int Len1 = (rightIDx - leftIDx + n) % n;
            for (k = 0; k <= Len1; ++k) {
                int iDx = (leftIDx + k) % n;
                if (Pts2d[iDx][1] > maxYPath1) 
                    maxYPath1 = Pts2d[iDx][1];
            }

            /* Path 2: rightIdx to leftIdx (incrementing). */
            IrtRType 
                maxYPath2 = minY;
            int Len2 = (leftIDx - rightIDx + n) % n;
            for (k = 0; k <= Len2; ++k) {
                int iDx = (rightIDx + k) % n;
                if (Pts2d[iDx][1] > maxYPath2) 
                    maxYPath2 = Pts2d[iDx][1];
            }

            if (maxYPath1 > maxYPath2) {
                topStart = leftIDx;
                topLen = Len1 + 1;
            }
            else {
                topStart = rightIDx;
                topLen = Len2 + 1;
            }
        }
    }

    /* Extract the top path. */
    IrtPtType 
        *top_Pts = (IrtPtType*) IritMalloc(sizeof(IrtPtType) * topLen);

    for (k = 0; k < topLen; ++k) {
        IRIT_PT_COPY(top_Pts[k], Pts2d[(topStart + k) % n]);
    }

    /* 8. Build silhouette curve mapped exactly into 3D World space.
     *    The 2D contour is originally in the (uVec, vVec) plane of the view.
     *    We rebuild the exact 3D orientation.
     *
     *    The extruded surface must be aligned with wVec in 3D.
     *    The start of the extrusion will be -wVec * FoamDepth/2.
     *    We also shift the entire cut upwards to be centered in the foam block. */
    Crv = IritCagdBspCrvNew(topLen, 2, IRIT_CAGD_PT_E3_TYPE);
    if (Crv == NULL) {
        IritFree(Pts2d);
        IritFree(top_Pts);
        return NULL;
    }

    for (k = 0; k < topLen; ++k) {
        IrtRType 
            Sx = (top_Pts[k][0] - Cx) * Scale,
            Sy = (top_Pts[k][1] - Cy) * Scale;

        /* Map back to 3D world coordinates. */
        IrtRType 
            Wx = Sx * uVec[0] + Sy * vVec[0],
            Wy = Sx * uVec[1] + Sy * vVec[1],
            Wz = Sx * uVec[2] + Sy * vVec[2];

        /* To bypass a critical math bug in the legacy hot_wire_cut.c logic,
         * we must pre-extrapolate the Z-coordinate. The legacy code computes:
         *   aBuggy = Crv2_Z + slope * Crv2_Y_rot
         *   zBuggy = aBuggy - slope * MACHINE_MAX_Y (where MACHINE_MAX_Y is 390.0)
         * But the correct math for extrapolating to Y = +/- 195.0 is:
         *   aTarget = Crv2_Z + slope * (195.0 - Crv2_Y_rot)
         * By setting aBuggy = aTarget, we can algebraically derive the exact Z shift
         * needed per-point to cancel the bug out completely.
         *   zShift = slope * (195.0 - 2 * Crv2_Y_rot)
         * After substituting Crv2_Y_rot = -wz * slope + wxyLen * FoamDepth / 2, we get:
         */
        IrtRType 
            Slope = wVec[2] / wxyLen,
            ZShift = Slope * (195.0 - wxyLen * Params -> FoamDepth) +
            2.0 * Wz * Slope * Slope;

        Crv -> Points[1][k] = Wx - wVec[0] * Params -> FoamDepth * 0.5;
        Crv -> Points[2][k] = Wy - wVec[1] * Params -> FoamDepth * 0.5;
        Crv -> Points[3][k] = Wz - wVec[2] * Params -> FoamDepth * 0.5 +
            Params -> FoamHeight * 0.5 + ZShift;
    }
    IritCagdBspKnotUniformOpen(topLen, 2, Crv -> KnotVector);

    IritFree(Pts2d);
    IritFree(top_Pts);

    /* 8. Extrude exactly along the true 3D view direction by FoamDepth.
     *    This creates a cylinder (surface) whose rulings are exactly wVec.
     *    When IritPrsrHWCProcessRulingPair processes this surface, it computes
     *    Angle = atan2(Dx, Dy) = atan2(wVec[0], wVec[1]), and rotates the
     *    clamps accordingly to perfectly align the Y-axis of the machine with wVec. */
    YDir.Vec[0] = wVec[0] * Params -> FoamDepth;
    YDir.Vec[1] = wVec[1] * Params -> FoamDepth;
    YDir.Vec[2] = wVec[2] * Params -> FoamDepth;

    Srf = IritCagdExtrudeSrf(Crv, &YDir);
    IritCagdCrvFree(Crv);
    if (Srf == NULL) return NULL;

    /* 8. Wrap into IritPrsrObjectStruct and tag ruling direction as U.
     *   With IRIT_CAGD_CONST_U_DIR the HWC generator sets:
     *     RulingDirection   = U  (front-to-back)
     *     SamplingDirection = V  (along silhouette)
     *   => Crv1 = S(u,vMin) and Crv2 = S(u,vMax), so both clamps trace
     *      the full silhouette contour in sync. */
    SrfObj = IritPrsrGenSrfObject("sil_ruled_srf", Srf, NULL);
    if (SrfObj != NULL) {
        IritMiscAttrIDSetObjectIntAttrib(SrfObj, IRIT_ATTR_ID_Dir,
            IRIT_CAGD_CONST_U_DIR);
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
static void HWCCombineGCodeFiles(const char *const *GcodeFiles,
                                 int NumFiles,
                                 const char *OutPath,
                                 const IritPrsrHWCDataStruct *Params)
{
#define IRIT_MAX_LINE_LEN 512
    FILE *Fout, *Fin;
    int Fi,
        f = 300,
        safeEntryLinesToModify = 0;
    double
        bOffset = 0.0,
        lastB = 0.0,
        x = 0.0,
        y = 0.0,
        z = 0.0,
        a = 0.0,
        b = 0.0,
        minZSeen = 99999.0,
        ExtraSafeZClearance = 150.0;
    char Line[IRIT_MAX_LINE_LEN];

    Fout = fopen(OutPath, "w");
    if (Fout == NULL) {
        fprintf(stderr, "HWCCombineGCodeFiles: cannot open '%s'.\n", OutPath);
        return;
    }
#ifdef DEBUG_HOT_WIRE_CUT_ALG2
    fprintf(Fout, "; Combined GCode - %d views\n\n", NumFiles);
#endif /* DEBUG_HOT_WIRE_CUT_ALG2 */

    for (Fi = 0; Fi < NumFiles; ++Fi) {
        lastB = bOffset;

        Fin = fopen(GcodeFiles[Fi], "r");
        if (Fin == NULL) {
            fprintf(stderr, "HWCCombineGCodeFiles: cannot open '%s'.\n",
                GcodeFiles[Fi]);
            continue;
        }
#ifdef DEBUG_HOT_WIRE_CUT_ALG2
        fprintf(Fout, "; === View %d start ===\n", Fi);
#endif /* DEBUG_HOT_WIRE_CUT_ALG2 */

        while (fgets(Line, sizeof(Line), Fin) != NULL) {
            if (strncmp(Line, ";Performing safe entry movement", 31) == 0) {
                safeEntryLinesToModify = 2;
            }

            if (strncmp(Line, "; Performing vertical safe lift to safe Z height.", 49) == 0) {
                safeEntryLinesToModify = 1;
            }

            if (strncmp(Line, "G1 ", 3) == 0 &&
                strstr(Line, "B") != NULL &&
                strstr(Line, "X") != NULL) {
                if (sscanf(Line, "G1 X%lf Y%lf Z%lf A%lf B%lf F%d",
                    &x, &y, &z, &a, &b, &f) == 6) {
                    /* B from the per-view GCode is the ABSOLUTE rotation
                       angle for this view (each view resets _TotalFoamRotation
                       to 0).  Adding bOffset blindly treated it as incremental,
                       making the machine accumulate extra full rotations.
                       Instead, compute the shortest angular path. */
                    {
                        IrtRType Diff = fmod(b - lastB, 360.0);
                        if (Diff > 180.0) 
                            Diff -= 360.0;
                        if (Diff < -180.0) 
                            Diff += 360.0;
                        lastB = lastB + Diff;
                    }
                    if (z < minZSeen) 
                        minZSeen = z;

                    /* Override any safe Z exits to provide massive clearance for taller foam blocks. */
                    if (Params != NULL && z >= Params -> FoamHeight) 
                        z = Params -> FoamHeight + ExtraSafeZClearance;
                    
                    if (Params != NULL && a >= Params -> FoamHeight) 
                        a = Params -> FoamHeight + ExtraSafeZClearance;
                    
                    if (safeEntryLinesToModify > 0) {
                        f = 350;
                        safeEntryLinesToModify--;
                    }

                    fprintf(Fout,
                        "G1 X%.3f Y%.3f Z%.3f A%.3f B%.3f F%d\n",
                        x, y, z, a, lastB, f);
                    continue;
                }
            }
            /* Skip per-view Go Home commands; the combined file
               emits its own G28 at the very end. */
            if (strncmp(Line, "G28", 3) == 0)
                continue;
            fputs(Line, Fout);
        }

        fclose(Fin);
        bOffset = lastB;
#ifdef DEBUG_HOT_WIRE_CUT_ALG2
        fprintf(Fout, "; === View %d end ===\n\n", Fi);
#endif /* DEBUG_HOT_WIRE_CUT_ALG2 */
    }

    if (Params != NULL) {
        fprintf(Fout, "\n; === Slice off base ===\n");
        /* Move safely above the foam. */
        z = Params -> FoamHeight + ExtraSafeZClearance;
        a = Params -> FoamHeight + ExtraSafeZClearance;
        fprintf(Fout, "G1 X%.3f Y%.3f Z%.3f A%.3f B%.3f F300\n", x, y, z, a, lastB);
        /* Move to the starting side of the machine (safely clear of foam). */
        x = 20.0;
        y = 20.0;
        fprintf(Fout, "G1 X%.3f Y%.3f Z%.3f A%.3f B%.3f F300\n", x, y, z, a, lastB);
        /* Drop down to the slice height (1mm below the lowest point of the object). */
        z = minZSeen - 1.0;
        if (z < Params -> MinimalHeight) 
            z = Params -> MinimalHeight;
        a = z;
        fprintf(Fout, "G1 X%.3f Y%.3f Z%.3f A%.3f B%.3f F150\n", x, y, z, a, lastB);
        /* Slice entirely through the foam to the other side. */
        x = 450.0;
        y = 450.0;
        fprintf(Fout, "G1 X%.3f Y%.3f Z%.3f A%.3f B%.3f F150\n", x, y, z, a, lastB);
        /* Pull back up to safe height. */
        z = Params -> FoamHeight + ExtraSafeZClearance;
        a = Params -> FoamHeight + ExtraSafeZClearance;
        fprintf(Fout, "G1 X%.3f Y%.3f Z%.3f A%.3f B%.3f F300\n", x, y, z, a, lastB);
    }

    fprintf(Fout, "\n; === Return Home Synchronously ===\n");
    if (Params != NULL) {
        /* Move XY and B to 0 while maintaining the high Z clearance to avoid collisions. */
        fprintf(Fout, "G1 X0.000 Y0.000 Z%.3f A%.3f B0.000 F300\n", z, a);
    }
    /* Finally drop Z back to 0 synchronously to avoid snapping the wire. */
    fprintf(Fout, "G1 X0.000 Y0.000 Z0.000 A0.000 B0.000 F300\n");

    /* Emit G28 only to ensure firmware resets coordinate systems if necessary,
       but since we are already at 0,0,0,0,0 it won't physically move sequentially. */
    fprintf(Fout, "G28; Go Home\n");
    fclose(Fout);
#ifdef DEBUG_HOT_WIRE_CUT_ALG2
    printf("Combined GCode written to: %s\n", OutPath);
#endif // DEBUG_HOT_WIRE_CUT_ALG2
#undef IRIT_MAX_LINE_LEN
}

/*****************************************************************************
* AUXILIARY:                                                                 *
*   Recursive function to convert polylines to polygons in a given object.   *
*   It traverses lists and converts matching polyline objects in-place.      *
*                                                                            *
* PARAMETERS:                                                                *
*   PObj:       The IritPrsrObjectStruct to                                  * 
*                                   process (can be a list or polyline).     *
*****************************************************************************/
static void HWCConvertPolylinesToPolygons(IritPrsrObjectStruct *PObj) {
    if (PObj == NULL) 
        return;
    if (IRIT_PRSR_IS_OLST_OBJ(PObj)) {
        IritPrsrObjectStruct *PList; 
        int i = 0;
        while ((PList = IritPrsrListObjectGet(PObj, i++)) != NULL) {
            HWCConvertPolylinesToPolygons(PList);
        }
    }
    else if (IRIT_PRSR_IS_POLY_OBJ(PObj) && IRIT_PRSR_IS_POLYLINE_OBJ(PObj)) 
        IRIT_PRSR_SET_POLYGON_OBJ(PObj);
    
    if (PObj -> Pnext != NULL) 
        HWCConvertPolylinesToPolygons(PObj -> Pnext);
    
}

/*****************************************************************************
* DESCRIPTION:                                                               M
* Main public entry point. Encapsulates entire hot-wire cutting pipeline.    M
*                                                                            *
* PARAMETERS:                                                                M
* RawModel:         Input 3D model object. A something something @Shiraz     M
* OutputGCodePath:  Path to output GCode file.                               M
* NumViews:         Number of view directions to generate.                   M
* OutputITDType:    0 to skip ITD, 1 for polylines, 2 for polygons.          M
* RulingExtent:        Wire path visualization length (0.0 for default).     M
* outputITDFileName:         Name of the output ITD file.                    M
* uniformCut:       View selection mode. 0 to use score-based best-view      M
*                   selection (HWCSelectBestViewSampling), 1 to use          M
*                   uniformly spaced views at 180/N degree intervals         M
*                   (HWCSelectViewsSampling).                                M
* FineNess:         Polygon tessellation fineness for freeform-to-polygon    M
*                   conversion. If 0 (or negative), uses the default         M
*                   IRIT_GRAP_DEFAULT_POLYGON_FINENESS. Higher values        M
*                   produce denser polygon meshes.                           M
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
IrtRType HWCGenerateGCodeFromObj(IritPrsrObjectStruct *RawModel,
                                 const char *OutputGCodePath,
                                 int NumViews,
                                 int OutputITDType,
                                 int RulingExtent,
                                 const char *outputITDFileName,
                                 int uniformCut,
                                 IrtRType FineNess)
{
    IritPrsrObjectStruct 
        *Contour, *RuledSrf, *SimObj, *AllSimObjs,
        *Solid = NULL;
    IrtHmgnMatType
        *Views = NULL;
    IritPrsrHWCDataStruct 
        HWCParams;
    char **gcodeFiles, gcodeFile[512], contourFile[512];
    int Vi, i,
        gcodeCount = 0,
        Result = 0;

    if (RawModel == NULL || OutputGCodePath == NULL || NumViews <= 0) {
        fprintf(stderr,
            "HWCGenerateGCodeFromObj: invalid parameters.\n");
        return 0;
    }

    gcodeFiles = NULL;
    AllSimObjs = IritPrsrGenLISTObject(NULL);

    /* Convert any FreeForm surfaces to polygonal mesh. */
    {
        IrtRType OldFineNess = IritPrsrFFCState.FineNess;
        if (FineNess > 0.0)
            IritPrsrFFCState.FineNess = FineNess;
        Solid = IritPrsrConvertFreeFormHierachy(RawModel, &IritPrsrFFCState,
                                                    FALSE, FALSE);
        IritPrsrFFCState.FineNess = OldFineNess;
    }
    RawModel = NULL;

    Solid = IritPrsrFlattenTree(Solid);
#ifdef DEBUG_HOT_WIRE_CUT_ALG2
    printf("Model ready. Type = %d.\n", Solid -> ObjType);
    fflush(stdout);
#endif /* DEBUG_HOT_WIRE_CUT_ALG2 */

    /* Save a copy as OBJ for visual inspection. */
#ifdef DEBUG_HOT_WIRE_CUT_ALG2
    printf("Saving solid.obj...\n"); fflush(stdout);
    if (IritPrsrOBJSaveFile(Solid, "solid.obj", FALSE, 2, 0))
        printf("Wrote solid.obj.\n");
    fflush(stdout);
#endif /* DEBUG_HOT_WIRE_CUT_ALG2 */

    /* 2. Select best view directions. */
#ifdef DEBUG_HOT_WIRE_CUT_ALG2
    printf("Selecting %d best view directions...\n", NumViews);
    fflush(stdout);
#endif /* DEBUG_HOT_WIRE_CUT_ALG2 */

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
#ifdef DEBUG_HOT_WIRE_CUT_ALG2 
    printf("View sampling done.\n");
    fflush(stdout);
#endif /* DEBUG_HOT_WIRE_CUT_ALG2 */

    /* 3. Set up HWC machine parameters. */
    IritPrsrHWCSetDfltParams(&HWCParams);
    HWCParams.MinimalHeight = 0.1;
    HWCParams.RuledApproxDir = IRIT_CAGD_CONST_U_DIR;
    HWCParams.PieceWiseRuledApproximation = 0.01;
    if (RulingExtent > 0.0) 
        HWCParams.SimulationArcLength = RulingExtent;

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
    for (Vi = 0; Vi < NumViews; ++Vi) {
        Contour = NULL;
        RuledSrf = NULL;
        SimObj = NULL;

        /* 5a. Approximate silhouette contour. */
#ifdef DEBUG_HOT_WIRE_CUT_ALG2
        printf("  View %d: Computing silhouette contour...\n", Vi);
        fflush(stdout);
#endif /* DEBUG_HOT_WIRE_CUT_ALG2 */

        Contour = HWCIritPrsrApproxBSplineContourFromSolidView(Solid,
                                                                Views[Vi], 128);

        if (Contour == NULL) {
            printf("View %d: contour approximation failed, skipping.\n", Vi);
            continue;
        }

        /* Save the raw 2D contour polygon for inspection. */
#ifdef DEBUG_HOT_WIRE_CUT_ALG2
        sprint(contourFile, sizeof(contourFile),
            "contour_view%d.obj", Vi);
        fflush(stdout);
        if (!IritPrsrOBJSaveFile(Contour, contourFile, FALSE, 0, 0))
            fprintf(stderr, "Failed saving %s.\n", contourFile);
        else
            printf("Wrote %s\n", contourFile);
#endif /* DEBUG_HOT_WIRE_CUT_ALG2 */

        /* 5b. Build ruled surface. */
#ifdef DEBUG_HOT_WIRE_CUT_ALG2
        printf("  View %d: Building ruled surface...\n", Vi);
        fflush(stdout);
#endif /* DEBUG_HOT_WIRE_CUT_ALG2 */
        RuledSrf = HWCBuildSilhouetteRuledSrf(Contour, Solid, Views[Vi],
                                                 &HWCParams);
        IritPrsrFreeObject(Contour);
        Contour = NULL;

        if (RuledSrf == NULL) {
            printf("View %d: failed to build ruled surface, skipping.\n",
                Vi);
            continue;
        }

        /* 5c. Generate GCode for this view. */
#ifdef DEBUG_HOT_WIRE_CUT_ALG2
        printf("  View %d: Generating GCode...\n", Vi);
        fflush(stdout);
#endif
        sprint(gcodeFile, sizeof(gcodeFile),
            "view%d_silhouette.gcode", Vi);

        SimObj = IritPrsrHWCCreatePath(RuledSrf, NULL, NULL, gcodeFile,
                                                &HWCParams);
        IritPrsrFreeObject(RuledSrf);
        RuledSrf = NULL;

        if (SimObj != NULL) {
#ifdef DEBUG_HOT_WIRE_CUT_ALG2
            printf("View %d: GCode written to %s\n", Vi, gcodeFile);
#endif /* DEBUG_HOT_WIRE_CUT_ALG2 */

            /* Append to AllSimObjs to combine all paths. */
            IritPrsrListObjectAppend(AllSimObjs, SimObj);
            SimObj = NULL;

            /* Store filename for combiner. */
            gcodeFiles[gcodeCount] = IritMiscStrdup(gcodeFile);
            ++gcodeCount;
        }
        else {
            printf("View %d: IritPrsrHWCCreatePath failed, skipping.\n",
                Vi);
        }
    }

    /* 6. Combine all per-view GCode files. */
    if (gcodeCount > 0) {
#ifdef DEBUG_HOT_WIRE_CUT_ALG2
        printf("Combining %d GCode files...\n", gcodeCount);
        fflush(stdout);
#endif /* DEBUG_HOT_WIRE_CUT_ALG2 */
        HWCCombineGCodeFiles((const char *const *)gcodeFiles, gcodeCount,
            OutputGCodePath, &HWCParams);
        Result = 1;
    }
    else {
        fprintf(stderr, "No per-view GCode was produced.\n");
        Result = 0;
    }

    /* Save all combined paths to a single ITD file. */
    if (AllSimObjs != NULL) {
        if (OutputITDType > 0) {
            if (OutputITDType == 2) {
                HWCConvertPolylinesToPolygons(AllSimObjs);
            }
            IritPrsrPutObjectToFile3(outputITDFileName, AllSimObjs, 0);
#ifdef DEBUG_HOT_WIRE_CUT_ALG2
            printf("All path geometries written to %s\n", outputITDFileName);
#endif /* DEBUG_HOT_WIRE_CUT_ALG2 */
        }

        IritPrsrFreeObject(AllSimObjs);
        AllSimObjs = NULL;
    }

    /* 7. Clean up intermediate filenames and delete files in release mode. */
    if (gcodeFiles != NULL) {
        for (i = 0; i < gcodeCount; ++i) {
            if (gcodeFiles[i] != NULL) {
#ifndef DEBUG_HOT_WIRE_CUT_ALG2
                remove(gcodeFiles[i]);
#endif /* DEBUG_HOT_WIRE_CUT_ALG2 */
                IritFree(gcodeFiles[i]);
            }
        }
        IritFree(gcodeFiles);
    }


    return Result;
}


