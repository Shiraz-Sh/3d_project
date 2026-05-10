#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "inc_irit/irit_sm.h"
#include "inc_irit/iritprsr.h"

/* Prototype for main GCode generation function from hot_wire_cut_new_alg.c */
int IritPrsrHWCGenerateGCodeFromObj(IPObjectStruct *RawModel,
    const char* OutputGCodePath,
    int NumViews,
    int OutputITDType,
    double ITDLength);

/* ===========================================================================
 * main
 * 
 * Simple entry point: Takes an .itd file and generates combined GCode.
 * =========================================================================== */
int main(int argc, char** argv)
{
    const char* inputFile = NULL;
    const char* outputFile = "combined_silhouette.gcode";
    const int NumViews = 10;
    int result;

    if (argc < 2) {
        fprintf(stderr, "Usage: %s <input_model_file.itd>\n", argv[0]);
        fprintf(stderr, "Output will be saved as: %s\n", outputFile);
        return 1;
    }

    inputFile = argv[1];

    /* Call the main GCode generation function */
    printf("========================================\n");
    printf("Hot-Wire Cut GCode Generator\n");
    printf("========================================\n\n");

    IPObjectStruct *RawModel = IritPrsrGetObjects2(inputFile);
    if (RawModel == NULL) {
        fprintf(stderr, "Failed to load model from '%s'.\n", inputFile);
        return 1;
    }

    result = IritPrsrHWCGenerateGCodeFromObj(RawModel, outputFile, NumViews, 1, 0.0, "all_views_paths.itd");

    if (result) {
        printf("\n========================================\n");
        printf("SUCCESS: GCode generated and saved to:\n  %s\n", outputFile);
        printf("========================================\n");
        return 0;
    }
    else {
        printf("\n========================================\n");
        printf("FAILED: Could not generate GCode\n");
        printf("========================================\n");
        return 1;
    }
}
