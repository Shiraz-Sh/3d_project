#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "inc_irit/irit_sm.h"
#include "inc_irit/iritprsr.h"

/* Prototype for main GCode generation function from hot_wire_cut_new_alg.c */
int IritPrsrHWCGenerateGCodeFromFile(const char* InputModelPath,
    const char* OutputGCodePath,
    int NumViews);

/* ===========================================================================
 * main
 * 
 * Simple entry point: Takes an .itd file and generates combined GCode.
 * =========================================================================== */
int main(int argc, char** argv)
{
    const char* inputFile = NULL;
    const char* outputFile = "combined_silhouette.gcode";
    const int NumViews = 5;
    int result;

    if (argc < 2) {
        fprintf(stderr, "Usage: %s <input_model_file.itd>\n", argv[0]);
        fprintf(stderr, "Output will be saved as: %s\n", outputFile);
        return 1;
    }

    inputFile = argv[1];

    /* Call the main GCode generation function */
    #ifdef DEBUG
    printf("========================================\n");
    printf("Hot-Wire Cut GCode Generator\n");
    printf("========================================\n\n");
    #endif /* DEBUG */

    result = IritPrsrHWCGenerateGCodeFromFile(inputFile, outputFile, NumViews);

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
