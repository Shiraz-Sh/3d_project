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

/* Prototypes for the helpers implemented in prsr_lib/hot_wire_cut_new_alg.c */
IPObjectStruct* IritPrsrExtrude2DPointsToSolidDir(const IrtE2PtStruct* pts,
    int n,
    IrtVecType Dir);

IPObjectStruct* IritPrsrGetSilhouettesForViews(const IPObjectStruct* PObj,
    const IrtHmgnMatType* ViewMats,
    int NumViews,
    int GridSize,
    int UsePreprocess);

/* Project-only helper (new) - updated signature includes IsPre flag */
IPObjectStruct* IritPrsrProjectSilhouetteToViewPolys(const IPObjectStruct* SilList,
    const IrtHmgnMatType* ViewMats,
    int NumViews,
    int IsPre);

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
    /* Build a 2D outline and create a simple test solid. */
    int NumPts = NUM_PTS;
    IrtE2PtStruct pts[NUM_PTS];
    for (int i = 0; i < NumPts; ++i) {
        IrtRType r = (i % 2 == 0) ? 1.0 : 0.5;
        IrtRType ang = 2.0 * M_PI * i / NUM_PTS;
        pts[i].Pt[0] = r * cos(ang);
        pts[i].Pt[1] = r * sin(ang);
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

    /* Build three canonical view matrices: +X, +Y, +Z */
    IrtHmgnMatType Views[3];
    {
        IrtVecType Eye, Center, Up;

        Eye[0] = 2.0; Eye[1] = 0.0; Eye[2] = 0.0;
        Center[0] = 0.0; Center[1] = 0.0; Center[2] = 0.0;
        Up[0] = 0.0; Up[1] = 0.0; Up[2] = 1.0;
        GenLookAtMatrix(Eye, Center, Up, Views[0]);

        Eye[0] = 0.0; Eye[1] = 2.0; Eye[2] = 0.0;
        GenLookAtMatrix(Eye, Center, Up, Views[1]);

        Eye[0] = 0.0; Eye[1] = 0.0; Eye[2] = 2.0;
        Up[0] = 0.0; Up[1] = 1.0; Up[2] = 0.0;
        GenLookAtMatrix(Eye, Center, Up, Views[2]);
    }

    const int NumViews = 3;

    /* ---------- DIRECT-ONLY SILHOUETTE WORKFLOW (no preprocessing) ---------- */

    IPObjectStruct* SilListDirect = IritPrsrGetSilhouettesForViews(Solid, Views, NumViews, 20, 0);
    if (SilListDirect == NULL) {
        fprintf(stderr, "IritPrsrGetSilhouettesForViews (direct) returned NULL\n");
        IritPrsrFreeObject(Solid);
        return 1;
    }
    printf("Direct mode: silhouettes returned = %d (expected %d)\n", CountListObjects(SilListDirect), NumViews);

    /* Save each per-view silhouette to its own OBJ file (so you get the raw 3D silhouette polylines). */
    {
        char fname[512];
        for (int vi = 0; vi < NumViews; ++vi) {
            IPObjectStruct* SV = IritPrsrListObjectGet(SilListDirect, vi);
            if (SV == NULL)
                continue;
            snprintf(fname, sizeof(fname), "sil_direct_view%d.obj", vi);
            if (!IritPrsrOBJSaveFile(SV, fname, FALSE, 0, 0))
                fprintf(stderr, "Failed saving %s\n", fname);
            else
                printf("Wrote %s\n", fname);
        }
    }

    /* Project direct silhouettes (flatten to 2D polygons per view). */
    IPObjectStruct* ProjDir = IritPrsrProjectSilhouetteToViewPolys(SilListDirect, Views, NumViews, 0);
    if (ProjDir == NULL) {
        printf("No projected polygons produced from direct silhouettes.\n");
        IritPrsrFreeObject(SilListDirect);
        IritPrsrFreeObject(Solid);
        return 0;
    }

    printf("Projected (direct) polygons = %d\n", CountListObjects(ProjDir));
    /* Save each projected 2D polygon object to its own OBJ file. */
    {
        char fname[512];
        int pidx = 0;
        IPObjectStruct* P = NULL;
        while ((P = IritPrsrListObjectGet(ProjDir, pidx)) != NULL) {
            /* If ObjName was set by projector, use it; otherwise use index-based name. */
            const char *base = P->ObjName ? P->ObjName : "dir_proj_poly";
            snprintf(fname, sizeof(fname), "%s.obj", base);
            if (!IritPrsrOBJSaveFile(P, fname, FALSE, 0, 0))
                fprintf(stderr, "Failed saving %s\n", fname);
            else
                printf("Wrote %s\n", fname);
            ++pidx;
        }
    }

    /* Extrude projected polygons in their original view directions and collect solids. */
    IPObjectStruct* SolList = IritPrsrExtrudeSilhouetteListsToViewSolids(NULL, ProjDir, Views, NumViews, 0.5);
    if (SolList == NULL) {
        printf("No solids generated from projected polygons.\n");
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

    /* Cleanup */
    IritPrsrFreeObject(ProjDir);
    IritPrsrFreeObject(SilListDirect);
    IritPrsrFreeObject(Solid);
    IritPrsrFreeObject(SolList);

    printf("Direct-only silhouette test completed.\n");
    return 0;
}


// Main with testing direct and prprocess approaches.

//int main(void)
//{
//    /* Build a more interesting 2D outline: a 12-point star/octagon and extrude it. */
//    int NumPts = NUM_PTS;
//    IrtE2PtStruct pts[NUM_PTS];
//    for (int i = 0; i < NumPts; ++i) {
//        /* alternating radii to produce a star-like outline */
//        IrtRType r = (i % 2 == 0) ? 1.0 : 0.5;
//        IrtRType ang = 2.0 * M_PI * i / NUM_PTS;
//        pts[i].Pt[0] = r * cos(ang);
//        pts[i].Pt[1] = r * sin(ang);
//    }
//
//    /* Extrude the 2D outline to a polygonal solid along Z (axis = 2).
//       Use the vector-based API to allow arbitrary directions: */
//    IrtVecType extrDir;
//    IRIT_PT_RESET(extrDir);
//    extrDir[2] = 1.0; /* length 1.0 along Z */
//    IPObjectStruct* Solid = IritPrsrExtrude2DPointsToSolidDir(pts, NumPts, extrDir);
//    if (Solid == NULL) {
//        fprintf(stderr, "Failed to create test solid.\n");
//        return 1;
//    }
//
//    /* Write the entire 3D model (Solid) to OBJ as triangles. */
//    if (IritPrsrOBJSaveFile(Solid, "solid.obj", FALSE, 2, 0)) {
//        printf("Wrote full model to solid.obj (triangulated).\n");
//    }
//    else {
//        fprintf(stderr, "Failed to write solid.obj\n");
//    }
//
//    /* Prepare three view matrices that look along +X, +Y and +Z axes.
//       We construct LookAt matrices with camera positioned on the +axis,
//       looking at the origin. */
//    IrtHmgnMatType Views[3];
//    {
//        IrtVecType Eye, Center, Up;
//
//        /* View 0: look from +X toward origin. */
//        Eye[0] = 2.0; Eye[1] = 0.0; Eye[2] = 0.0;
//        Center[0] = 0.0; Center[1] = 0.0; Center[2] = 0.0;
//        Up[0] = 0.0; Up[1] = 0.0; Up[2] = 1.0; /* Z up */
//        GenLookAtMatrix(Eye, Center, Up, Views[0]);
//
//        /* View 1: look from +Y toward origin. */
//        Eye[0] = 0.0; Eye[1] = 2.0; Eye[2] = 0.0;
//        Center[0] = 0.0; Center[1] = 0.0; Center[2] = 0.0;
//        Up[0] = 0.0; Up[1] = 0.0; Up[2] = 1.0; /* Z up */
//        GenLookAtMatrix(Eye, Center, Up, Views[1]);
//
//        /* View 2: look from +Z toward origin (top). Use Y as up to avoid singularity. */
//        Eye[0] = 0.0; Eye[1] = 0.0; Eye[2] = 2.0;
//        Center[0] = 0.0; Center[1] = 0.0; Center[2] = 0.0;
//        Up[0] = 0.0; Up[1] = 1.0; Up[2] = 0.0; /* Y up for top view */
//        GenLookAtMatrix(Eye, Center, Up, Views[2]);
//    }
//
//    const int NumViews = 3;
//
//    /* Test with preprocessing (should be faster for manyviews). */
//    IPObjectStruct* SilListPre = IritPrsrGetSilhouettesForViews(Solid, Views, NumViews, 20, 1);
//    if (SilListPre == NULL) {
//        fprintf(stderr, "IritPrsrGetSilhouettesForViews (preprocess) returned NULL\n");
//        IritPrsrFreeObject(Solid);
//        return 1;
//    }
//    int cntPre = CountListObjects(SilListPre);
//    printf("Preprocess mode: silhouettes returned = %d (expected %d)\n", cntPre, NumViews);
//
//    /* Test without preprocessing (direct per-view computation). */
//    IPObjectStruct* SilListDirect = IritPrsrGetSilhouettesForViews(Solid, Views, NumViews, 20, 0);
//    if (SilListDirect == NULL) {
//        fprintf(stderr, "IritPrsrGetSilhouettesForViews (direct) returned NULL\n");
//        IritPrsrFreeObject(SilListPre);
//        IritPrsrFreeObject(Solid);
//        return 1;
//    }
//    int cntDirect = CountListObjects(SilListDirect);
//    printf("Direct mode: silhouettes returned = %d (expected %d)\n", cntDirect, NumViews);
//
//    /* Basic validation. */
//    if (cntPre != NumViews || cntDirect != NumViews) {
//        fprintf(stderr, "Unexpected number of silhouette objects returned.\n");
//        IritPrsrFreeObject(SilListPre);
//        IritPrsrFreeObject(SilListDirect);
//        IritPrsrFreeObject(Solid);
//        return 1;
//    }
//
//    /* Optionally inspect/print polygon counts per silhouette (simple diagnostic). */
//    for (int i = 0; i < NumViews; ++i) {
//        IPObjectStruct* S = IritPrsrListObjectGet(SilListPre, i);
//        int numPl = 0;
//        if (S && IP_IS_POLY_OBJ(S)) {
//            for (IPPolygonStruct* Pl = S->U.Pl; Pl != NULL; Pl = Pl->Pnext)
//                ++numPl;
//        }
//        printf("Pre silhouette %d polygons: %d\n", i, numPl);
//    }
//
//    /* Also write the LIST objects as an OBJ (grouped) using the same exporter. */
//    IritPrsrOBJSaveFile(SilListPre, "sil_pre.obj", FALSE, 0, 0);
//    IritPrsrOBJSaveFile(SilListDirect, "sil_direct.obj", FALSE, 0, 0);
//
//    /* --- Save each silhouette element separately as its own OBJ file --- */
//    {
//        char fname[512];
//        IPObjectStruct *S = NULL;
//        int idx = 0;
//
//        /* Save preprocessed silhouettes individually */
//        idx = 0;
//        while ((S = IritPrsrListObjectGet(SilListPre, idx)) != NULL) {
//            snprintf(fname, sizeof(fname), "sil_pre_view%d_sil%d.obj", idx / 1, idx); /* idx corresponds to view index order */
//            if (!IritPrsrOBJSaveFile(S, fname, FALSE, 0, 0))
//                fprintf(stderr, "Failed saving %s\n", fname);
//            else
//                printf("Wrote %s\n", fname);
//            ++idx;
//        }
//
//        /* Save direct silhouettes individually */
//        idx = 0;
//        while ((S = IritPrsrListObjectGet(SilListDirect, idx)) != NULL) {
//            snprintf(fname, sizeof(fname), "sil_direct_view%d_sil%d.obj", idx / 1, idx);
//            if (!IritPrsrOBJSaveFile(S, fname, FALSE, 0, 0))
//                fprintf(stderr, "Failed saving %s\n", fname);
//            else
//                printf("Wrote %s\n", fname);
//            ++idx;
//        }
//    }
//
//   
//    /* --- Save projected 2D polygons (use new API that accepts whole silhouette LIST + Views array) --- */
//    IPObjectStruct* ProjPre = NULL;
//    IPObjectStruct* ProjDir = NULL;
//    {
//        char fname[512];
//
//        /* Project all preprocessed silhouettes using the Views array.
//           The new helper takes the silhouette LIST, the Views array and NumViews.
//           It returns a LIST of polygon objects (flattened across silhouettes). */
//        ProjPre = IritPrsrProjectSilhouetteToViewPolys(SilListPre, Views, NumViews, 1);
//        if (ProjPre != NULL) {
//            int pidx = 0;
//            IPObjectStruct* P = NULL;
//            int total = CountListObjects(ProjPre);
//            printf("Projected (pre) polygons = %d\n", total);
//            while ((P = IritPrsrListObjectGet(ProjPre, pidx)) != NULL) {
//                snprintf(fname, sizeof(fname), "sil_pre_proj_poly%03d.obj", pidx);
//                if (!IritPrsrOBJSaveFile(P, fname, FALSE, 0, 0))
//                    fprintf(stderr, "Failed saving %s\n", fname);
//                else
//                    printf("Wrote %s\n", fname);
//                ++pidx;
//            }
//        }
//
//        /* Project all direct silhouettes similarly. */
//        ProjDir = IritPrsrProjectSilhouetteToViewPolys(SilListDirect, Views, NumViews, 0);
//        if (ProjDir != NULL) {
//            int pidx = 0;
//            IPObjectStruct* P = NULL;
//            int total = CountListObjects(ProjDir);
//            printf("Projected (direct) polygons = %d\n", total);
//            while ((P = IritPrsrListObjectGet(ProjDir, pidx)) != NULL) {
//                snprintf(fname, sizeof(fname), "sil_direct_proj_poly%03d.obj", pidx);
//                if (!IritPrsrOBJSaveFile(P, fname, FALSE, 0, 0))
//                    fprintf(stderr, "Failed saving %s\n", fname);
//                else
//                    printf("Wrote %s\n", fname);
//                ++pidx;
//            }
//        }
//    }
//
//    /* ---  extrude 2D polygons in the direction of the view point --- */
//    IPObjectStruct* SolList = NULL;
//    {
//        /* NOTE: pass projected polygon lists (ProjPre / ProjDir) to the extruder */
//        SolList = IritPrsrExtrudeSilhouetteListsToViewSolids(ProjPre, ProjDir, Views, NumViews, 0.5);
//
//        if (SolList != NULL) {
//            /* Save each generated solid using its ObjName (set by the function). */
//            char fname[1024];
//            int sidx = 0;
//            IPObjectStruct* S = NULL;
//            while ((S = IritPrsrListObjectGet(SolList, sidx)) != NULL) {
//                const char* base = S->ObjName ? S->ObjName : "extr_poly";
//                snprintf(fname, sizeof(fname), "%s.obj", base);
//                if (!IritPrsrOBJSaveFile(S, fname, FALSE, 2, 0))
//                    fprintf(stderr, "Failed saving %s\n", fname);
//                else
//                    printf("Wrote %s\n", fname);
//                ++sidx;
//            }
//            /* Keep SolList around for later use (do not free here). */
//        }
//        else {
//            printf("No solids generated from projected polygons.\n");
//        }
//    }
//
//    /* Free allocated objects. */
//    IritPrsrFreeObject(ProjPre);
//    IritPrsrFreeObject(ProjDir);
//    IritPrsrFreeObject(SilListPre);
//    IritPrsrFreeObject(SilListDirect);
//    IritPrsrFreeObject(Solid);
//
//    printf("Silhouette tests completed successfully. Files written: sil_pre.obj, sil_direct.obj, solid.obj and per-view solids\n");
//    return 0;
//}
