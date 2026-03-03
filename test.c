#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "inc_irit/irit_sm.h"
#include "inc_irit/iritprsr.h"
#include "inc_irit/allocate.h"
#include "inc_irit/attribut.h"
#include "inc_irit/bool_lib.h"
#include "inc_irit/geom_lib.h"
#include "inc_irit/misc_lib.h"
#include "inc_irit/ip_cnvrt.h" /* Added: prototype for IritPrsrOpenPolysToClosed */
#include "../prsr_lib/prsr_loc.h"
/* For IritPrsrHWCDataStruct and friends */

/* Ensure M_PI is available on MSVC builds that require _USE_MATH_DEFINES. */
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif


/* Prototypes for the helpers implemented in prsr_lib/hot_wire_cut_new_alg.c */
IPObjectStruct* IritPrsrExtrude2DPointsToSolidDir(const IrtE2PtStruct* pts,
    int n,
    IrtVecType Dir);

/* Approximate outer B-spline-like contour from a view (implemented in hot_wire_cut_new_alg.c) */
IPObjectStruct* IritPrsrApproxBSplineContourFromSolidView(IPObjectStruct* Solid,
    const IrtHmgnMatType ViewMat,
    int NumCtrl,
    int Iterations,
    IrtRType ShrinkStep,
    IrtRType SmoothW);


/* Now operate on projected polygon lists (ProjPre / ProjDir) */
IPObjectStruct* IritPrsrExtrudeSilhouetteListsToViewSolids(IPObjectStruct* ProjListPre,
    IPObjectStruct* ProjListDirect,
    const IrtHmgnMatType* Views,
    int NumViews,
    IrtRType Depth);

/* GenLookAtMatrix (defined in hot_wire_cut_new_alg.c) - declare for test usage */
void GenLookAtMatrix(IrtVecType Eye, IrtVecType Center, IrtVecType Up, IrtHmgnMatType Mat);

/* Already in hot_wire_cut.c - ensure prototype visible to the compiler. */
void IritPrsrHWCSetDfltParams(IritPrsrHWCDataStruct* Data);

IPObjectStruct* IritPrsrHWCCreatePath(const IPObjectStruct* ModelMainPart,
    const IPObjectStruct* ModelTopCover,
    const IPObjectStruct* ModelBottomCover,
    const char* GCodeOutputFilePath,
    const IritPrsrHWCDataStruct* Params);

/* New/changed helper exported from hot_wire_cut_new_alg.c */
IrtHmgnMatType* SelectBestViewSampling(IPObjectStruct* PObj,
    int NumSamples,
    IrtHmgnMatType* ResultMat);


/* Count elements inside a LIST object using the public accessor. */
static int CountListObjects(IPObjectStruct* ListObj)
{
    int cnt = 0;
    if (ListObj == NULL)
        return 0;

    while (IritPrsrListObjectGet(ListObj, cnt) != NULL)
        ++cnt;

    return cnt;
}

#define NUM_PTS 12

// Main with only direct approach.

int main(void)
{
    /* Build a 2D half-dome outline and create a simple test solid.
       The half-dome is a semicircular arc (y >= 0) plus a straight diameter
       segment between the arc endpoints to close the polygon. */
    int nArc = NUM_PTS;                          /* samples along the semicircle */
    int nBase = IRIT_MAX(1, NUM_PTS / 8);       /* number of points along the flat base (excluding arc endpoints) */
    int NumPts = nArc + nBase;
    IrtE2PtStruct pts[NUM_PTS * 2]; /* plenty of room */

    IrtRType R = 1.0; /* dome radius */
    /* Sample semicircle from angle = 0..pi (right -> left) */
    for (int i = 0; i < nArc; ++i) {
        IrtRType ang = M_PI * (IrtRType)i / (IrtRType)(nArc - 1);
        pts[i].Pt[0] = R * cos(ang);   /* x */
        pts[i].Pt[1] = R * sin(ang);   /* y (>=0) */
    }
    /* Add base (diameter) points from left endpoint back to right endpoint,
       excluding endpoints to avoid duplicate consecutive points. */
    for (int j = 0; j < nBase; ++j) {
        IrtRType t = (j + 1) / (IrtRType)(nBase + 1); /* in (0,1) */
        /* map t along x from -R to +R, y = 0 */
        pts[nArc + j].Pt[0] = -R + t * (2.0 * R);
        pts[nArc + j].Pt[1] = 0.0;
    }

    IrtVecType extrDir;
    IRIT_PT_RESET(extrDir);
    extrDir[2] = 1.0;

    IPObjectStruct* Solid = IritPrsrExtrude2DPointsToSolidDir(pts, NumPts, extrDir);

    if (Solid == NULL) {
        fprintf(stderr, "Failed to create test solid.\n");
        return 1;
    }

    if (IritPrsrOBJSaveFile(Solid, "solid.obj", FALSE, 2, 0))
        printf("Wrote full model to solid.obj (triangulated).\n");
    else
        fprintf(stderr, "Failed to write solid.obj\n");


        /* Choose views using the new helper in hot_wire_cut_new_alg.c.
       This replaces the manual construction of three canonical views. */
    const int NumViews = 3;
    IrtHmgnMatType* Views = SelectBestViewSampling(Solid, NumViews, NULL);
    if (Views == NULL) {
        fprintf(stderr, "SelectBestViewSampling failed to provide view matrices.\n");
        IritPrsrFreeObject(Solid);
        return 1;
    }

    /* Choose three canonical viewpoints: +X, +Y, +Z. */
    //const int NumViews = 3;
    //IrtHmgnMatType Views[3];
    //{
    //    IrtVecType Eye, Center, Up;

    //    Center[0] = 0.0; Center[1] = 0.0; Center[2] = 0.0;

    //    /* +X */
    //    Eye[0] = 2.0; Eye[1] = 0.0; Eye[2] = 0.0;
    //    Up[0] = 0.0;  Up[1] = 0.0;  Up[2] = 1.0;
    //    GenLookAtMatrix(Eye, Center, Up, Views[0]);

    //    /* +Y */
    //    Eye[0] = 0.0; Eye[1] = 2.0; Eye[2] = 0.0;
    //    Up[0] = 0.0;  Up[1] = 0.0;  Up[2] = 1.0;
    //    GenLookAtMatrix(Eye, Center, Up, Views[1]);

    //    /* +Z (top) - use Y as Up to avoid singularity for strict top) */
    //    Eye[0] = 0.0; Eye[1] = 0.0; Eye[2] = 2.0;
    //    Up[0] = 0.0;  Up[1] = 1.0;  Up[2] = 0.0;
    //    GenLookAtMatrix(Eye, Center, Up, Views[2]);
    //}

    /* ---------- REPLACEMENT: Use B-spline contour approximation per view ----------
       For each view, compute a robust outer contour polygon using the new
       IritPrsrApproxBSplineContourFromSolidView helper, collect them in ProjDir
       (world-space polygons), and continue with extrusion as before.
    */

    IPObjectStruct* ProjDir = IritPrsrGenLISTObject(NULL);
    if (ProjDir == NULL) {
        fprintf(stderr, "Failed allocating ProjDir list\n");
        IritPrsrFreeObject(Solid);
        IritFree(Views);
        return 1;
    }

    int created = 0;
    for (int vi = 0; vi < NumViews; ++vi) {
        /* Approximate an outer contour for this view.
           Parameters: NumCtrl=32, Iterations=80, ShrinkStep=0.25, SmoothW=0.3 */
        IPObjectStruct* Contour = IritPrsrApproxBSplineContourFromSolidView(Solid, Views[vi], 32, 80, 0.25, 0.3);
        if (Contour == NULL) {
            printf("View %d: contour approximation failed\n", vi);
            continue;
        }

        /* Name and append to list so later code can find view index if needed. */
        char cname[256];
        snprintf(cname, sizeof(cname), "proj_dir_view%d_poly%03d", vi, 0);
        Contour->ObjName = IritMiscStrdup(cname);
        Contour->Pnext = NULL;
        IritPrsrListObjectAppend(ProjDir, Contour);
        ++created;

        /* Save raw contour poly for inspection. */
        char fname[512];
        snprintf(fname, sizeof(fname), "contour_view%d.obj", vi);
        if (!IritPrsrOBJSaveFile(Contour, fname, FALSE, 0, 0))
            fprintf(stderr, "Failed saving %s\n", fname);
        else
            printf("Wrote %s\n", fname);
    }

    if (created == 0) {
        printf("No contours produced from approximation.\n");
        IritPrsrFreeObject(ProjDir);
        IritPrsrFreeObject(Solid);
        IritFree(Views);
        return 0;
    }

    /* Extrude projected polygons in their original view directions and collect solids. */
    IPObjectStruct* SolList = IritPrsrExtrudeSilhouetteListsToViewSolids(NULL, ProjDir, Views, NumViews, 0.5);
    if (SolList == NULL) {
        printf("No solids generated from projected contours.\n");
    }
    else {
        /* Save each generated solid using its ObjName (set by the extruder). */
        char fname[1024];
        int sidx = 0;
        IPObjectStruct* S = NULL;
        while ((S = IritPrsrListObjectGet(SolList, sidx)) != NULL) {
            const char* base = S->ObjName ? S->ObjName : "dir_extr_poly";
            snprintf(fname, sizeof(fname), "%s.obj", base);
            if (!IritPrsrOBJSaveFile(S, fname, FALSE, 2, 0))
                fprintf(stderr, "Failed saving %s\n", fname);
            else
                printf("Wrote %s\n", fname);
            ++sidx;
        }
    }

    /* 6. GCode Generation Block (Only added here at the end). */
    if (SolList != NULL) {
        IritPrsrHWCDataStruct HWCParams;
        IritPrsrHWCSetDfltParams(&HWCParams);
        /* Adjust parameters to clear the machine bed and avoid Z bound errors. */
        HWCParams.MinimalHeight = 20.0;
        HWCParams.RuledApproxDir = CAGD_CONST_V_DIR;

        int sidx = 0;
        IPObjectStruct* S = NULL;
        while ((S = IritPrsrListObjectGet(SolList, sidx++)) != NULL) {
            int vi, pi;
            /* Map the extruded solid back to its source projection polygon. */
            if (S->ObjName && sscanf(S->ObjName, "sil_dir_view%d_poly%03d_extr", &vi, &pi) == 2) {
                char search[256], gcname[1024];
                snprintf(search, sizeof(search), "proj_dir_view%d_poly%03d", vi, pi);
                IPObjectStruct* P = NULL;
                for (int j = 0; (P = IritPrsrListObjectGet(ProjDir, j)) != NULL; j++) {
                    if (P->ObjName && strcmp(P->ObjName, search) == 0) break;
                }

                if (P && IP_IS_POLY_OBJ(P)) {
                    /* Create a surface representation required for HWC ruling calculations. */
                    int n = IritPrsrVrtxListLen(P->U.Pl->PVertex);
                    CagdCrvStruct* Crv = IritCagdBspCrvNew(n, 2, CAGD_PT_E3_TYPE);
                    IPVertexStruct* V = P->U.Pl->PVertex;
                    for (int k = 0; k < n; k++, V = V->Pnext) {
                        Crv->Points[1][k] = V->Coord[0];
                        Crv->Points[2][k] = V->Coord[1];
                        /* Lift the object by 100mm to move it into positive machine Z space. */
                        Crv->Points[3][k] = V->Coord[2] + 100.0;
                    }
                    IritCagdBspKnotUniformOpen(n, 2, Crv->KnotVector);

                    CagdVecStruct Dir;
                    IrtRType Depth = 10.0; /* Extrusion depth for HWC simulation. */
                    Dir.Vec[0] = Views[vi][2][0] * Depth;
                    Dir.Vec[1] = Views[vi][2][1] * Depth;
                    Dir.Vec[2] = Views[vi][2][2] * Depth;

                    CagdSrfStruct* Srf = IritCagdExtrudeSrf(Crv, &Dir);
                    IPObjectStruct* TempSrfObj = IritPrsrGenSrfObject("temp", Srf, NULL);
                    IritMiscAttrIDSetObjectIntAttrib(TempSrfObj, IRIT_ATTR_ID_Dir, CAGD_CONST_V_DIR);

                    snprintf(gcname, sizeof(gcname), "%s.gcode", S->ObjName);
                    IPObjectStruct* Sim = IritPrsrHWCCreatePath(TempSrfObj, NULL, NULL, gcname, &HWCParams);
                    if (Sim) {
                        printf("GCode successfully written to: %s\n", gcname);
                        IritPrsrFreeObject(Sim);
                    }
                    else {
                        printf("This is some wierd object I can't cut\n");
                    }

                    IritPrsrFreeObject(TempSrfObj);
                    IritCagdCrvFree(Crv);
                    /* DO NOT free Views here - it is used later/only freed once after all processing. */
                }
                else {
                    printf("This is some wierd object I can't cut\n");
                }

            }
        }
    }

    /* Cleanup allocated IRIT objects. */

    /* Cleanup */
    IritPrsrFreeObject(ProjDir);
    IritPrsrFreeObject(Solid);
    IritPrsrFreeObject(SolList);
    IritFree(Views);

    printf("Contour-approximation silhouette test completed.\n");
    return 0;
}
