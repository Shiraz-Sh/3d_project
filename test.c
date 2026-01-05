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

/* Project-only helper (new) */
IPObjectStruct* IritPrsrProjectSilhouetteToViewPolys(const IPObjectStruct* SilList,
    const IrtHmgnMatType* ViewMats,
    int NumViews);


/* Already in hot_wire_cut.c - ensure prototype visible to the compiler. */
void IritPrsrHWCSetDfltParams(IritPrsrHWCDataStruct* Data);

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

int main(void)
{
    /* Build a more interesting 2D outline: a 12-point star/octagon and extrude it. */
    int NumPts = NUM_PTS;
    IrtE2PtStruct pts[NUM_PTS];
    for (int i = 0; i < NumPts; ++i) {
        /* alternating radii to produce a star-like outline */
        IrtRType r = (i % 2 == 0) ? 1.0 : 0.5;
        IrtRType ang = 2.0 * M_PI * i / NUM_PTS;
        pts[i].Pt[0] = r * cos(ang);
        pts[i].Pt[1] = r * sin(ang);
    }

    /* Extrude the 2D outline to a polygonal solid along Z (axis = 2).
       Use the vector-based API to allow arbitrary directions: */
    IrtVecType extrDir;
    IRIT_PT_RESET(extrDir);
    extrDir[2] = 1.0; /* length 1.0 along Z */
    IPObjectStruct* Solid = IritPrsrExtrude2DPointsToSolidDir(pts, NumPts, extrDir);
    if (Solid == NULL) {
        fprintf(stderr, "Failed to create test solid.\n");
        return 1;
    }

    /* Write the entire 3D model (Solid) to OBJ as triangles. */
    if (IritPrsrOBJSaveFile(Solid, "solid.obj", FALSE, 2, 0)) {
        printf("Wrote full model to solid.obj (triangulated).\n");
    }
    else {
        fprintf(stderr, "Failed to write solid.obj\n");
    }

    /* Prepare three view matrices:
       - Views[0]: identity (front)
       - Views[1]: rotated 45 degrees about Y (angled)
       - Views[2]: top-down view that looks onto the XY plane (rotate -180 deg about X)
    */
    IrtHmgnMatType Views[3];
    IritMiscMatGenUnitMat(Views[0]);
    IritMiscMatGenMatRotY1(IRIT_DEG2RAD(45.0), Views[1]);
    IritMiscMatGenMatRotX1(IRIT_DEG2RAD(-180.0), Views[2]);

    const int NumViews = 3;

    /* Test with preprocessing (should be faster for manyviews). */
    IPObjectStruct* SilListPre = IritPrsrGetSilhouettesForViews(Solid, Views, NumViews, 20, 1);
    if (SilListPre == NULL) {
        fprintf(stderr, "IritPrsrGetSilhouettesForViews (preprocess) returned NULL\n");
        IritPrsrFreeObject(Solid);
        return 1;
    }
    int cntPre = CountListObjects(SilListPre);
    printf("Preprocess mode: silhouettes returned = %d (expected %d)\n", cntPre, NumViews);

    /* Test without preprocessing (direct per-view computation). */
    IPObjectStruct* SilListDirect = IritPrsrGetSilhouettesForViews(Solid, Views, NumViews, 20, 0);
    if (SilListDirect == NULL) {
        fprintf(stderr, "IritPrsrGetSilhouettesForViews (direct) returned NULL\n");
        IritPrsrFreeObject(SilListPre);
        IritPrsrFreeObject(Solid);
        return 1;
    }
    int cntDirect = CountListObjects(SilListDirect);
    printf("Direct mode: silhouettes returned = %d (expected %d)\n", cntDirect, NumViews);

    /* Basic validation. */
    if (cntPre != NumViews || cntDirect != NumViews) {
        fprintf(stderr, "Unexpected number of silhouette objects returned.\n");
        IritPrsrFreeObject(SilListPre);
        IritPrsrFreeObject(SilListDirect);
        IritPrsrFreeObject(Solid);
        return 1;
    }

    /* Optionally inspect/print polygon counts per silhouette (simple diagnostic). */
    for (int i = 0; i < NumViews; ++i) {
        IPObjectStruct* S = IritPrsrListObjectGet(SilListPre, i);
        int numPl = 0;
        if (S && IP_IS_POLY_OBJ(S)) {
            for (IPPolygonStruct* Pl = S->U.Pl; Pl != NULL; Pl = Pl->Pnext)
                ++numPl;
        }
        printf("Pre silhouette %d polygons: %d\n", i, numPl);
    }

    /* Also write the LIST objects as an OBJ (grouped) using the same exporter. */
    IritPrsrOBJSaveFile(SilListPre, "sil_pre.obj", FALSE, 0, 0);
    IritPrsrOBJSaveFile(SilListDirect, "sil_direct.obj", FALSE, 0, 0);

    /* --- Save each silhouette element separately as its own OBJ file --- */
    {
        char fname[512];
        IPObjectStruct *S = NULL;
        int idx = 0;

        /* Save preprocessed silhouettes individually */
        idx = 0;
        while ((S = IritPrsrListObjectGet(SilListPre, idx)) != NULL) {
            snprintf(fname, sizeof(fname), "sil_pre_view%d_sil%d.obj", idx / 1, idx); /* idx corresponds to view index order */
            if (!IritPrsrOBJSaveFile(S, fname, FALSE, 0, 0))
                fprintf(stderr, "Failed saving %s\n", fname);
            else
                printf("Wrote %s\n", fname);
            ++idx;
        }

        /* Save direct silhouettes individually */
        idx = 0;
        while ((S = IritPrsrListObjectGet(SilListDirect, idx)) != NULL) {
            snprintf(fname, sizeof(fname), "sil_direct_view%d_sil%d.obj", idx / 1, idx);
            if (!IritPrsrOBJSaveFile(S, fname, FALSE, 0, 0))
                fprintf(stderr, "Failed saving %s\n", fname);
            else
                printf("Wrote %s\n", fname);
            ++idx;
        }
    }

    /* --- Save projected 2D polygons (use new API that accepts whole silhouette LIST + Views array) --- */
    {
        char fname[512];

        /* Project all preprocessed silhouettes using the Views array.
           The new helper takes the silhouette LIST, the Views array and NumViews.
           It returns a LIST of polygon objects (flattened across silhouettes). */
        IPObjectStruct* ProjPre = IritPrsrProjectSilhouetteToViewPolys(SilListPre, Views, NumViews);
        if (ProjPre != NULL) {
            int pidx = 0;
            IPObjectStruct* P = NULL;
            int total = CountListObjects(ProjPre);
            printf("Projected (pre) polygons = %d\n", total);
            while ((P = IritPrsrListObjectGet(ProjPre, pidx)) != NULL) {
                snprintf(fname, sizeof(fname), "sil_pre_proj_poly%03d.obj", pidx);
                if (!IritPrsrOBJSaveFile(P, fname, FALSE, 0, 0))
                    fprintf(stderr, "Failed saving %s\n", fname);
                else
                    printf("Wrote %s\n", fname);
                ++pidx;
            }
            IritPrsrFreeObject(ProjPre);
        }

        /* Project all direct silhouettes similarly. */
        IPObjectStruct* ProjDir = IritPrsrProjectSilhouetteToViewPolys(SilListDirect, Views, NumViews);
        if (ProjDir != NULL) {
            int pidx = 0;
            IPObjectStruct* P = NULL;
            int total = CountListObjects(ProjDir);
            printf("Projected (direct) polygons = %d\n", total);
            while ((P = IritPrsrListObjectGet(ProjDir, pidx)) != NULL) {
                snprintf(fname, sizeof(fname), "sil_direct_proj_poly%03d.obj", pidx);
                if (!IritPrsrOBJSaveFile(P, fname, FALSE, 0, 0))
                    fprintf(stderr, "Failed saving %s\n", fname);
                else
                    printf("Wrote %s\n", fname);
                ++pidx;
            }
            IritPrsrFreeObject(ProjDir);
        }
    }

    /* Free allocated objects. */
    IritPrsrFreeObject(SilListPre);
    IritPrsrFreeObject(SilListDirect);
    IritPrsrFreeObject(Solid);

    printf("Silhouette tests completed successfully. Files written: sil_pre.obj, sil_direct.obj, solid.obj and per-view solids\n");
    return 0;
}