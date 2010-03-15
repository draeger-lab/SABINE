/****************************************************************************
 *
 * Globals.h
 *
 * Header file for Globals definition.
 ****************************************************************************/

#include <stdio.h>



#define TRUE 1
#define FALSE 0

extern int COMPUTE_KERNEL_MATRIX;
extern int USE_WEIGHTED_MATRIX;
extern int ONLY_EXACT_MISMATCHES;
extern int PRINT_KMER_HISTOGRAM;
extern int PRINT_SOURCE_WEIGHTS;

extern int QUIET_MODE;

extern double **kernelMatrix;

extern int numTotalRecords;
extern struct inputRecord **globalRecords;


extern int *globalSequence;
extern int globalSequenceLength;
extern struct intPointerVector *exampleInputSequences;


extern double MAX_MISMATCH_WEIGHT;
extern int KERNEL_OUTPUT_SCALING;

extern int NUM_INPUT_SYMBOLS;
extern int NUM_OUTPUT_SYMBOLS;
extern int MAX_MISMATCHES;
extern int MIN_SIGNAL;
extern int MAX_DEPTH;
extern int VERBOSE;

extern char INPUT_SYMBOL_FILENAME[61];
extern char OUTPUT_SYMBOL_FILENAME[61];
extern char MUTATION_MATRIX_FILENAME[61];

extern struct node *staticNewNode;

extern int inSequence;

extern double **mutationMatrix;

extern struct symbolTable *globInputSymbols;
extern struct symbolTable *globOutputSymbols;


extern int sequenceInput;

extern int *currentSequence;

extern int *currentSequence;

extern int *sourceOffsets;
extern int NUM_SOURCES;


void readGlobalsFile(char *);
void printGlobals();
void setDefaultGlobals();
void freeGlobals();
