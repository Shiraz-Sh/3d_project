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
#include "inc_irit/ip_cnvrt.h"
#include "../prsr_lib/prsr_loc.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ---------------------------------------------------------------------------
 * Prototypes for functions implemented in prsr_lib/hot_wire_cut_new_alg.c
 * ------------------------------------------------------------------------- */
IPObjectStruct* IritPrsrGetSilhouettesForViews(const IPObjectStruct* PObj,
    const IrtHmgnMatType* ViewMats,
    int NumViews,
    int GridSize,
    int UsePreprocess);

IPObjectStruct* IritPrsrApproxBSplineContourFromSolidView(IPObjectStruct* Solid,
    const IrtHmgnMatType ViewMat,
    int NumCtrl,
    int Iterations,
    IrtRType ShrinkStep,
    IrtRType SmoothW);

/* Build the silhouette-to-ruled-surface converter (new in hot_wire_cut_new_alg.c). */
IPObjectStruct* IritPrsrHWCBuildSilhouetteRuledSrf(
    const IPObjectStruct* Contour,
    const IrtHmgnMatType         ViewMat,
    const IritPrsrHWCDataStruct* Params);

IrtHmgnMatType* SelectBestViewSampling(IPObjectStruct* PObj,
    int NumSamples,
    IrtHmgnMatType* ResultMat);

void GenLookAtMatrix(IrtVecType Eye,
    IrtVecType Center,
    IrtVecType Up,
    IrtHmgnMatType Mat);

/* Prototypes from hot_wire_cut.c */
void IritPrsrHWCSetDfltParams(IritPrsrHWCDataStruct* Data);

IPObjectStruct* IritPrsrHWCCreatePath(const IPObjectStruct* ModelMainPart,
    const IPObjectStruct* ModelTopCover,
    const IPObjectStruct* ModelBottomCover,
    const char* GCodeOutputFilePath,
    const IritPrsrHWCDataStruct* Params);

/* ---------------------------------------------------------------------------
 * CombineGCodeFiles
 *
 * Concatenate per-view GCode files into one combined file, adjusting B
 * (rotation) values so rotation is continuous across views.
 *
 * For each line of the form "G1 X... Y... Z... A... B<val> F..." the B
 * value is offset by `bOffset` before writing. After reading each file the
 * function scans for the last B value written and uses it as the offset
 * base for the next file.
 *
 * Parameters:
 *   GcodeFiles  - array of file path strings, length NumFiles.
 *   NumFiles    - number of files to combine.
 *   OutPath     - path for the combined output file.
 * ------------------------------------------------------------------------- */
static void CombineGCodeFiles(const char* const* GcodeFiles,
    int NumFiles,
    const char* OutPath)
{
#define _MAX_LINE 512
    FILE* fout;
    int fi;
    double bOffset = 0.0; /* cumulative B offset */

    fout = fopen(OutPath, "w");
    if (fout == NULL) {
        fprintf(stderr, "CombineGCodeFiles: cannot open '%s'\n", OutPath);
        return;
    }

    fprintf(fout, "; Combined GCode - %d views\n\n", NumFiles);

    for (fi = 0; fi < NumFiles; ++fi) {
        FILE* fin;
        char line[_MAX_LINE];
        double lastB = bOffset; /* track last B seen in this view */

        fin = fopen(GcodeFiles[fi], "r");
        if (fin == NULL) {
            fprintf(stderr, "CombineGCodeFiles: cannot open '%s'\n",
                GcodeFiles[fi]);
            continue;
        }

        fprintf(fout, "; === View %d start ===\n", fi);

        while (fgets(line, sizeof(line), fin) != NULL) {
            /* Check for G1 motion lines that carry a B value. */
            if (strncmp(line, "G1 ", 3) == 0 &&
                strstr(line, "B") != NULL &&
                strstr(line, "X") != NULL) {
                /* Parse the B value from the line. */
                double x, y, z, a, b;
                int    f;
                if (sscanf(line, "G1 X%lf Y%lf Z%lf A%lf B%lf F%d",
                    &x, &y, &z, &a, &b, &f) == 6) {
                    lastB = bOffset + b;
                    fprintf(fout, "G1 X%.3f Y%.3f Z%.3f A%.3f B%.3f F%d\n",
                        x, y, z, a, lastB, f);
                    continue;
                }
            }
            /* Pass all other lines through unchanged. */
            fputs(line, fout);
        }

        fclose(fin);

        /* The next view's B values are offset by the last B we wrote. */
        bOffset = lastB;

        fprintf(fout, "; === View %d end ===\n\n", fi);
    }

    /* Final home command. */
    fprintf(fout, "\nG28; Go Home\n");
    fclose(fout);
    printf("Combined GCode written to: %s\n", OutPath);
#undef _MAX_LINE
}

/* ---------------------------------------------------------------------------
 * CountListObjects - count elements in a LIST object.
 * ------------------------------------------------------------------------- */
static int CountListObjects(IPObjectStruct* ListObj)
{
    int cnt = 0;
    if (ListObj == NULL) return 0;
    while (IritPrsrListObjectGet(ListObj, cnt) != NULL)
        ++cnt;
    return cnt;
}

/* ===========================================================================
 * main
 * =========================================================================== */
int main(int argc, char** argv)
{
    IPObjectStruct* RawModel;
    IPObjectStruct* Solid;
    IrtHmgnMatType* Views;
    const int NumViews = 3;
    IritPrsrHWCDataStruct HWCParams;
    char** gcodeFiles;
    int gcodeCount = 0;
    int vi, i;
    const char* inputFile = NULL;

    if (argc < 2) {
        fprintf(stderr, "Usage: %s <input_model_file.itd>\n", argv[0]);
        return 1;
    }
    inputFile = argv[1];

    /* ------------------------------------------------------------------
     * 1. Load 3D model from file.
     *    IritPrsrGetObjects2 reads any IRIT-supported format (.itd, .obj,
     *    .igs, .stl) by file extension.
     * ------------------------------------------------------------------ */
    printf("Loading '%s'...\n", inputFile); fflush(stdout);
    RawModel = IritPrsrGetObjects2(inputFile);
    if (RawModel == NULL) {
        fprintf(stderr, "Failed to load model from '%s'.\n"
            "Make sure the file is in the working directory.\n",
            inputFile);
        return 1;
    }
    printf("File read OK. Type = %d\n", RawModel->ObjType); fflush(stdout);

    /* Convert any FreeForm surfaces (like the Dome.itd surface of revolution)
       into a Polygonal mesh, because the contour algorithms purely operate on polygons. */
    Solid = IritPrsrConvertFreeFormHierachy(RawModel, &IritPrsrFFCState, FALSE, FALSE);
    RawModel = NULL;

    Solid = IritPrsrFlattenTree(Solid);
    /* FlattenTree consumes its argument and frees list headers internally if applicable.
       Solid now represents the flattened polygonal geometry. */

    printf("Model ready. Type = %d  (5 = polygon, 6 = list).\n",
        Solid->ObjType);
    fflush(stdout);

    /* Save a copy as OBJ for visual inspection. */
    printf("Saving solid.obj...\n"); fflush(stdout);
    if (IritPrsrOBJSaveFile(Solid, "solid.obj", FALSE, 2, 0))
        printf("Wrote solid.obj\n");
    fflush(stdout);

    /* ------------------------------------------------------------------
     * 2. Select best view directions.
     * ------------------------------------------------------------------ */
    printf("Selecting %d best view directions...\n", NumViews); fflush(stdout);
    Views = SelectBestViewSampling(Solid, NumViews, NULL);
    if (Views == NULL) {
        fprintf(stderr, "SelectBestViewSampling failed.\n");
        IritPrsrFreeObject(Solid);
        return 1;
    }
    printf("View sampling done.\n"); fflush(stdout);

    /* ------------------------------------------------------------------
     * 3. Set up HWC machine parameters once (shared across all views).
     * ------------------------------------------------------------------ */
    IritPrsrHWCSetDfltParams(&HWCParams);
    HWCParams.MinimalHeight = 20.0;  /* mm above machine bed */
    HWCParams.RuledApproxDir = CAGD_CONST_U_DIR; /* ruling = U */
    HWCParams.PieceWiseRuledApproximation = 0.01; /* near-exact for ruled */

    /* ------------------------------------------------------------------
     * 4. Per-view: compute silhouette contour, build correct ruled surface,
     *    generate GCode.
     * ------------------------------------------------------------------ */
     /* Track per-view GCode filenames for the combiner. */
    gcodeFiles = (char**)IritMalloc(sizeof(char*) * NumViews);

    for (vi = 0; vi < NumViews; ++vi) {
        char contourFile[512], gcodeFile[512];

        /* 4a. Approximate outer silhouette contour from this view. */
        IPObjectStruct* Contour = IritPrsrApproxBSplineContourFromSolidView(
            Solid, Views[vi], 32, 80, 0.25, 0.3);

        if (Contour == NULL) {
            printf("View %d: contour approximation failed, skipping.\n", vi);
            continue;
        }

        /* Save the raw 2D contour polygon for inspection. */
        snprintf(contourFile, sizeof(contourFile), "contour_view%d.obj", vi);
        if (!IritPrsrOBJSaveFile(Contour, contourFile, FALSE, 0, 0))
            fprintf(stderr, "Failed saving %s\n", contourFile);
        else
            printf("Wrote %s\n", contourFile);

        /* 4b. Build the ruled surface for the HWC machine.
         *
         *  IritPrsrHWCBuildSilhouetteRuledSrf:
         *    - Projects contour to view-local 2D (ignores world Z).
         *    - Maps 2D silhouette uniformly into FoamWidth x FoamHeight space.
         *    - Builds extrusion surface along Y (FoamDepth).
         *    - Tags surface with CAGD_CONST_U_DIR so HWC samples in V,
         *      making both clamps trace the same silhouette shape in sync.
         */
        IPObjectStruct* RuledSrf = IritPrsrHWCBuildSilhouetteRuledSrf(
            Contour, Views[vi], &HWCParams);

        IritPrsrFreeObject(Contour);

        if (RuledSrf == NULL) {
            printf("View %d: failed to build ruled surface, skipping.\n", vi);
            continue;
        }

        /* 4c. Generate GCode for this view's silhouette cut. */
        snprintf(gcodeFile, sizeof(gcodeFile), "view%d_silhouette.gcode", vi);

        IPObjectStruct* SimObj = IritPrsrHWCCreatePath(RuledSrf, NULL, NULL,
            gcodeFile, &HWCParams);
        IritPrsrFreeObject(RuledSrf);

        if (SimObj != NULL) {
            printf("View %d: GCode written to %s\n", vi, gcodeFile);
            IritPrsrFreeObject(SimObj);

            /* Store filename for the combiner. */
            gcodeFiles[gcodeCount] = IritMiscStrdup(gcodeFile);
            ++gcodeCount;
        }
        else {
            printf("View %d: IritPrsrHWCCreatePath returned NULL "
                "(cannot cut this silhouette).\n", vi);
        }
    }

    /* ------------------------------------------------------------------
     * 5. Combine all per-view GCode files into one continuous flow file.
     *    B rotation values are accumulated across views for continuity.
     * ------------------------------------------------------------------ */
    if (gcodeCount > 0) {
        CombineGCodeFiles((const char* const*)gcodeFiles, gcodeCount,
            "combined_silhouette.gcode");
    }
    else {
        printf("No per-view GCode was produced; combined file not written.\n");
    }

    /* ------------------------------------------------------------------
     * 6. Cleanup.
     * ------------------------------------------------------------------ */
    for (i = 0; i < gcodeCount; ++i)
        IritFree(gcodeFiles[i]);
    IritFree(gcodeFiles);
    IritPrsrFreeObject(Solid);
    IritFree(Views);

    printf("Silhouette cut test completed.\n");
    return 0;
}
