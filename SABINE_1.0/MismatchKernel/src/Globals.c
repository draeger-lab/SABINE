/****************************************************************************
 *
 * Globals.c
 *
 * Implementation for Globals variables.
 ****************************************************************************/

#include <stdio.h>
#include <string.h>

#include "includes.h"

#include <stdlib.h>

// This file declares all of the global variables and has functions
// for reading global variables files.  The file also has functions
// for reading and writing the global variables files to disk.


/******************* external global variables **********************/
// These variables are defined in the parameter file.

int COMPUTE_KERNEL_MATRIX;
int USE_WEIGHTED_MATRIX;
int ONLY_EXACT_MISMATCHES;
int PRINT_KMER_HISTOGRAM;
int PRINT_SOURCE_WEIGHTS;

int QUIET_MODE;

double **kernelMatrix;


int numTotalRecords=0;
struct inputRecord **globalRecords;



int *currentSequence;

int maxSignalCount=0;
int *maxSignalSequence;


int *sourceOffsets;


// This is the minimum signal for printing.
int MIN_SIGNAL;

// This is the number of source files.
int NUM_SOURCES;

// This is the number of symbols in the input (output) alphabets.
int NUM_INPUT_SYMBOLS;
int NUM_OUTPUT_SYMBOLS;


// The number of Phis at a node.  This is always MAX_SKIP+2
int MAX_MISMATCHES;

// The maximum depth of a tree.
int MAX_DEPTH;

// Variable to control verbose output.
int VERBOSE;



// The input and output filenames.
char INPUT_SYMBOL_FILENAME[61];
char OUTPUT_SYMBOL_FILENAME[61];
char MUTATION_MATRIX_FILENAME[61];

/****************** internel global variables *********************/
// These variables are internal global variables.

// These variables store the input and output symbolTables.
struct symbolTable *globInputSymbols;
struct symbolTable *globOutputSymbols;

double **mutationMatrix;
double MAX_MISMATCH_WEIGHT;
int KERNEL_OUTPUT_SCALING;

int *globalSequence;
int globalSequenceLength;


// freeGlobals() frees memory used by the global variables.
void freeGlobals() {

  //  free(backGroundProbs);
  //free(SPStarLogs);
  free(sourceOffsets);
}

// setDefaultGlobals() sets default values for the globals.
void setDefaultGlobals() {
  
  NUM_INPUT_SYMBOLS=-1;
  NUM_OUTPUT_SYMBOLS=-1;
  MAX_MISMATCHES=1;
  MAX_MISMATCH_WEIGHT=1000;
  KERNEL_OUTPUT_SCALING=1;
  NUM_SOURCES=0;
  MIN_SIGNAL=1;
  MAX_DEPTH=5;
  VERBOSE=0;
  QUIET_MODE=1;
  COMPUTE_KERNEL_MATRIX=0;
  USE_WEIGHTED_MATRIX=0;
  ONLY_EXACT_MISMATCHES=0;
  PRINT_KMER_HISTOGRAM=0;
  PRINT_SOURCE_WEIGHTS=0;


  strcpy(INPUT_SYMBOL_FILENAME,"../data/Amino.txt");
  strcpy(OUTPUT_SYMBOL_FILENAME,"../data/IDs.txt");
  strcpy(MUTATION_MATRIX_FILENAME,"../data/substitution-matrix.txt");


}


// matchTag() is a function that compares strings used in processing
// globals files.
int matchTag(char *tag, char *input) {
  if (strcmp(tag,input)==0)
    return 1;
  else
    return 0;
}


// readGlobalsFile() is a function for reading a parameter files that
// sets all of the global variables.
void readGlobalsFile(char *fileName) {
  
  char line[101];
  int input;
  double dInput;
  char tag[101];
  char sInput[101];
  
  FILE *fp;

  if ((fp = fopen(fileName,"r")) == NULL) {
    fprintf(stderr, "Error opening Globals File: %s\n",fileName);
    exit(1);
 }

  while (fgets(line,100,fp)!=NULL) {

    if (strlen(line)>1) {
      if (line[0]!='#') {
	sscanf(line,"%s%d",tag,&input);
	if (matchTag("NUM_INPUT_SYMBOLS",tag)) {
	  NUM_INPUT_SYMBOLS=input;
	}
	if (matchTag("NUM_OUTPUT_SYMBOLS",tag)) {
	  NUM_OUTPUT_SYMBOLS=input;
	}

	if (matchTag("MAX_MISMATCHES",tag)) {
	  MAX_MISMATCHES=input;
	}
	if (matchTag("NUM_SOURCES",tag)) {
	  NUM_SOURCES=input;
	}
	if (matchTag("MIN_SIGNAL",tag)) {
	  MIN_SIGNAL=input;
	}

	if (matchTag("MAX_DEPTH",tag)) {
	  MAX_DEPTH=input;
	}
	if (matchTag("KERNEL_OUTPUT_SCALING",tag)) {
	  KERNEL_OUTPUT_SCALING=input;
	}
	
	if (matchTag("MAX_MISMATCH_WEIGHT",tag)) {
	  MAX_MISMATCH_WEIGHT=dInput;
	}

	sscanf(line,"%s%s",tag,sInput);
	if (matchTag("INPUT_SYMBOL_FILENAME",tag)) {
	  strcpy(INPUT_SYMBOL_FILENAME,sInput);
	}
	if (matchTag("OUTPUT_SYMBOL_FILENAME",tag)) {
	  strcpy(OUTPUT_SYMBOL_FILENAME,sInput);
	}
	if (matchTag("MUTATION_MATRIX_FILENAME",tag)) {
	  strcpy(MUTATION_MATRIX_FILENAME,sInput);
	}
	if (matchTag("VERBOSE",tag)) {
	  if (matchTag("TRUE",sInput)) {
	    VERBOSE=1;
	  } else {
	    VERBOSE=0;
	  }
	}
	if (matchTag("QUIET_MODE",tag)) {
	  if (matchTag("TRUE",sInput)) {
	    QUIET_MODE=1;
	  } else {
	    QUIET_MODE=0;
	  }
	}
	if (matchTag("COMPUTE_KERNEL_MATRIX",tag)) {
	  if (matchTag("TRUE",sInput)) {
	    COMPUTE_KERNEL_MATRIX=1;
	  } else {
	    COMPUTE_KERNEL_MATRIX=0;
	  }
	}
	if (matchTag("USE_WEIGHTED_MATRIX",tag)) {
	  if (matchTag("TRUE",sInput)) {
	    USE_WEIGHTED_MATRIX=1;
	  } else {
	    USE_WEIGHTED_MATRIX=0;
	  }
	}
	if (matchTag("ONLY_EXACT_MISMATCHES",tag)) {
	  if (matchTag("TRUE",sInput)) {
	    ONLY_EXACT_MISMATCHES=1;
	  } else {
	    ONLY_EXACT_MISMATCHES=0;
	  }
	}
	if (matchTag("PRINT_KMER_HISTOGRAM",tag)) {
	  if (matchTag("TRUE",sInput)) {
	    PRINT_KMER_HISTOGRAM=1;
	  } else {
	    PRINT_KMER_HISTOGRAM=0;
	  }
	}
	if (matchTag("PRINT_SOURCE_WEIGHTS",tag)) {
	  if (matchTag("TRUE",sInput)) {
	    PRINT_SOURCE_WEIGHTS=1;
	  } else {
	    PRINT_SOURCE_WEIGHTS=0;
	  }
	}



      }
    }
  }
  fclose(fp);
}

// printGlobals() prints out the values of the global variables to
// stdout.
void printGlobals() {

  printf("Global Variables:\n");
  printf("NUM_INPUT_SYMBOLS=%d\n",NUM_INPUT_SYMBOLS);
  printf("NUM_OUTPUT_SYMBOLS=%d\n",NUM_OUTPUT_SYMBOLS);
  printf("MAX_MISMATCHES=%d\n",MAX_MISMATCHES);
  printf("NUM_SOURCES=%d\n",NUM_SOURCES);
  printf("MIN_SIGNAL=%d\n",MIN_SIGNAL);
  printf("MAX_DEPTH=%d\n",MAX_DEPTH);
  if (VERBOSE)
    printf("VERBOSE=TRUE\n");
  else
    printf("VERBOSE=FALSE\n");
  if (QUIET_MODE)
    printf("QUIET_MODE=TRUE\n");
  else
    printf("QUIET_MODE=FALSE\n");
  printf("INPUT_SYMBOL_FILENAME=%s\n",INPUT_SYMBOL_FILENAME);
  printf("OUTPUT_SYMBOL_FILENAME=%s\n",OUTPUT_SYMBOL_FILENAME);
  printf("MUTATION_MATRIX_FILENAME=%s\n",MUTATION_MATRIX_FILENAME);


}

