/*****************************************************************************


string-kernel

A program that computes the string kernel

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

//Include for getopt
#include <unistd.h>


#include "includes.h"



void printHelpInfo() {

  fprintf(stderr,"Usage: string-kernel -OPTIONS FILENAME\n");
  fprintf(stderr,"Valid Options:\n");
  fprintf(stderr,"-h\tHelp\t\tPrint this help.\n");
  fprintf(stderr,"-v\tVerbose\t\tEnter Verbose Mode.\n");
  fprintf(stderr,"-s\tSequence\tRead in Sequence File.\n");
  fprintf(stderr,"-i FILE\tInput Symbols\tInput Symbol FILE.\n");
  fprintf(stderr,"-o FILE\tOutput Symbols\tOutput Symbol FILE.\n");
  fprintf(stderr,"-g FILE\tGlobals\t\tRead Globals FILE.\n");
  fprintf(stderr,"-L NUM\tLength\t\tSet the Kmer length to NUM.\n");
  fprintf(stderr,"-D NUM\tMismatches\t\tSet the Maximum Mismatches to NUM.\n");
  fprintf(stderr,"-N NUM\tMin Kmers\tSet the Minimum Number of Kmers to NUM.\n");
  fprintf(stderr,"-K \tKernel Matrix\tToggle Computation of Kernel Matrix.\n");  
  fprintf(stderr,"-H \tKmer Histogram\tToggle Computation of Kmer Histogram.\n");

}


// createStaticVars() allocates space for the static variables.
void createStaticVars() {
  int i;

  currentSequence=calloc(MAX_DEPTH,sizeof(int));


  for (i=0; i<MAX_DEPTH; i++) {
    currentSequence[i]=0;
  }


  NUM_SOURCES=NUM_OUTPUT_SYMBOLS;


}


void setSourceOffsets() {

  int i;

  sourceOffsets=calloc(NUM_OUTPUT_SYMBOLS, sizeof(int));
  for (i=0; i<NUM_OUTPUT_SYMBOLS; i++) {
    sourceOffsets[i]=(int)(globOutputSymbols->updateAmounts[i]);
    if (VERBOSE)
      printf("Source : %s : %d\n",lookUpSymbol(globOutputSymbols,i),sourceOffsets[i]);
  }
}

void incrementKernelMatrix(struct inputRecord **currentInstances, 
			   int numInstances, int *currentMismatchCounts,
			   double *instanceWeights) {

  int i,j;
  for (i=0; i<numInstances; i++) {
    //printf("instanceWeight[%d]=%f\n",i,instanceWeights[i]);
    for (j=0; j<numInstances; j++) {
      if ((!ONLY_EXACT_MISMATCHES) 
	  || ((currentMismatchCounts[i]==MAX_MISMATCHES)
	       && (currentMismatchCounts[j]==MAX_MISMATCHES))) {
	if (USE_WEIGHTED_MATRIX)
	  kernelMatrix[currentInstances[i]->source][currentInstances[j]->source]+=exp(-(instanceWeights[i]+instanceWeights[j]));
	else
	  kernelMatrix[currentInstances[i]->source][currentInstances[j]->source]++;
      }
    }
  }
}


void processLeaf(struct inputRecord **currentInstances, 
		 int numInstances, int *currentMismatchCounts,
		 double *instanceWeights) {
  
  int i;
  double sourceWeight;

  if (COMPUTE_KERNEL_MATRIX) {
    incrementKernelMatrix(currentInstances, numInstances, 
			  currentMismatchCounts, instanceWeights);
  }
  if (PRINT_KMER_HISTOGRAM) {
    printf("KMER: ");
    for (i=0; i<MAX_DEPTH; i++) {
      printf("%s ", lookUpSymbol(globInputSymbols,currentSequence[i]));
    }
    printf(": %d ",numInstances);
    if (PRINT_SOURCE_WEIGHTS) {
      sourceWeight=0.0;
      for (i=0; i<numInstances; i++) {
	sourceWeight+=globOutputSymbols->updateAmounts[currentInstances[i]->source];
      }
      printf(" %f ",sourceWeight);
    }
    printf("\n");
  }
}
void traverse(struct inputRecord **currentInstances, 
		    int *currentMismatchCounts, 
		    double *currentMismatchWeights,
		    int numInstances,
		    int currentDepth) {

  int i;
  int currentSymbol;
  struct inputRecord **newInstances;
  int *newMismatchCounts;
  double *newMismatchWeights;
  int numNewInstances;
  //int *sourceInstances,numSources;


  if (currentDepth<=4 || VERBOSE) {
    if (!QUIET_MODE) {
      printf("At Node: ");
      for (i=0; i<currentDepth; i++) {
	printf("%s ", lookUpSymbol(globInputSymbols,currentSequence[i]));
      }
      printf(": %d\n",numInstances);
    }
  }
  //fflush(stdout);

  if (currentDepth==MAX_DEPTH) {
    processLeaf(currentInstances, numInstances, 
		currentMismatchCounts, currentMismatchWeights);
  } else {

    newInstances=calloc(numInstances,sizeof(struct inputRecord *));
    newMismatchCounts=calloc(numInstances,sizeof(int));
    newMismatchWeights=calloc(numInstances,sizeof(double));


    for (currentSymbol=0; currentSymbol<NUM_INPUT_SYMBOLS; currentSymbol++) {
      //printf("Current Symbol = %d\n",currentSymbol);
      numNewInstances=0;
      currentSequence[currentDepth]=currentSymbol;


      for (i=0; i<numInstances; i++) {
	if (currentInstances[i]->inputs[currentDepth]!=currentSymbol) {
	  //printf("No Match\n");
	  if ((currentMismatchCounts[i]<MAX_MISMATCHES) ) {
	    // && (currentMismatchWeights[i]+
            // mutationMatrix[currentInstances[i]
            //->inputs[currentDepth]][currentSymbol] < MAX_MISMATCH_WEIGHT)) {
	    newInstances[numNewInstances]=currentInstances[i];
	    newMismatchCounts[numNewInstances]=currentMismatchCounts[i]+1;
	    newMismatchWeights[numNewInstances]=currentMismatchWeights[i]+mutationMatrix[currentInstances[i]->inputs[currentDepth]][currentSymbol];
	    numNewInstances++;
	  }
	} else {
	  //printf("Match %d\n",i);
	  newInstances[numNewInstances]=currentInstances[i];
	  newMismatchCounts[numNewInstances]=currentMismatchCounts[i];
	  newMismatchWeights[numNewInstances]=currentMismatchWeights[i]+mutationMatrix[currentInstances[i]->inputs[currentDepth]][currentSymbol];
	  numNewInstances++;
	}
      }
      if (numNewInstances>=MIN_SIGNAL) {
	traverse(newInstances,newMismatchCounts,
		       newMismatchWeights,numNewInstances,
		       currentDepth+1);
      }
    
    }
    free(newInstances);
    free(newMismatchCounts);
    free(newMismatchWeights);

  }
}


void readMutationMatrix(char *fileName) {

  FILE *fp;
  char line[1000001];
 
  char *token;
  double dvalue;

  int i,j;

  mutationMatrix=calloc(NUM_INPUT_SYMBOLS,sizeof(double *));
  for (i=0; i<NUM_INPUT_SYMBOLS; i++) {
    mutationMatrix[i]=calloc(NUM_INPUT_SYMBOLS,sizeof(double));
  }

  if ((fp = fopen(fileName,"r"))==NULL) {
    fprintf(stderr, "Invalid Mutation Matrix File:%s\n",fileName);
    exit(1);
  }


  for (i=0; i<NUM_INPUT_SYMBOLS; i++) {

    fgets(line,100000,fp);
      
      //printf("Line = %s\n",line);

    token=strtok(line," \t\n");
    for (j=0; j<NUM_INPUT_SYMBOLS; j++) {


      //printf("Token = %s\n",token);
      sscanf(token,"%lf",&dvalue);
      mutationMatrix[i][j]=-log(dvalue);

      //printf("mutationMatrix[%d][%d]=%f\n",i,j,mutationMatrix[i][j]);

      //mutationMatrix[i][j]=dvalue;
      token=strtok(NULL," \t\n");
      
    }
  }
}

void initializeKernelMatrix() {
  int i,j;

  kernelMatrix=calloc(NUM_OUTPUT_SYMBOLS,sizeof(double *));
  for (j=0; j<NUM_OUTPUT_SYMBOLS; j++) {
    kernelMatrix[j]=calloc(NUM_OUTPUT_SYMBOLS,sizeof(double));
  }
  for (i=0; i<NUM_OUTPUT_SYMBOLS; i++) {
    for (j=0; j<NUM_OUTPUT_SYMBOLS; j++) {
      kernelMatrix[i][j]=0;
    }
  }
}
  
int
main (int argc, char **argv)
{



  //char *inputFile;

  char **sequenceFiles;

  // Command Line Flags.
  int iflag = 0;//Input symbol file flag
  char *ivalue = NULL;//Input symbol file filename.
  int oflag = 0;//Output symbol file flag
  char *ovalue = NULL;//Output symbol file filename.
  int gflag = 0;//Globals flag.
  char *gvalue = NULL;//Globals filename.
  int vflag = 0;//Verbose Mode
  int qflag = 0; // Quiet mode
  int hflag = 0;// Print help.
  int Lflag = 0; // Length Command Line
  int Hflag = 0; // print histogram
  char *Lvalue = NULL;
  int Dflag = 0; // Mismatch Command Line
  char *Dvalue = NULL;
  int Nflag = 0;
  char *Nvalue = NULL;
  int Kflag = 0; //Compute String Kernel

  int c,i,j;

  int *mismatchCounts;
  int numInputFiles;

  double *mismatchWeights;

  globInputSymbols=NULL;
  globOutputSymbols=NULL;


  
  opterr = 0;
  
  while ((c = getopt (argc, argv, "L:D:N:vhqKHg:i:o:")) != -1)
    switch (c)
      {
      case 'i':
	iflag = 1;
	ivalue = optarg;
	break;
      case 'o':
	oflag = 1;
	ovalue = optarg;
	break;
      case 'K':
	Kflag=1;
	break;
      case 'g':
	gflag=1;
	gvalue = optarg;
	break;
      case 'v':
	vflag=1;
	break;
      case 'q':
	qflag=1;
	break;
      case 'h':
	hflag=1;
	break;
      case 'L':
	Lflag=1;
	Lvalue = optarg;
	break;
      case 'D':
	Dflag=1;
	Dvalue = optarg;
	break;
      case 'N':
	Nflag=1;
	Nvalue = optarg;
	break;
      case 'H':
	Hflag=1;
	break;
      case '?':
	fprintf (stderr, "Unknown option `-%c'.\n", optopt);
	return 1;
      default:
	abort ();
      }

  // Print Help Info and exit.
  if (hflag) {
    printHelpInfo();
    exit(0);
  }
  
  //The input files are the first arguments.
  if (optind < argc) {
    sequenceFiles=argv+optind;
    numInputFiles=argc-optind;
  } else {
    fprintf(stderr, "Error: no input file.\n");
    printHelpInfo();
    exit(1);
  }



  // First set the globals to the defaults;
  setDefaultGlobals();


  // Read in the Globals File if given.
  if (gflag) {
    readGlobalsFile(gvalue);
  }


  if (Lflag) {
    sscanf(Lvalue,"%d",&MAX_DEPTH);
  }
  if (Dflag) {
    sscanf(Dvalue,"%d",&MAX_MISMATCHES);
  }
  if (Nflag) {
    sscanf(Nvalue,"%d",&MIN_SIGNAL);
  }


  // Read in the symbols files if given.
  if (iflag) {
    strcpy(INPUT_SYMBOL_FILENAME,ivalue);
  }
  if (oflag) {
    strcpy(OUTPUT_SYMBOL_FILENAME,ovalue);
  }



  // Toggle Verbose mode
  if (vflag) {
    if (VERBOSE)
      VERBOSE=0;
    else {
      VERBOSE=TRUE;
      QUIET_MODE=0;
    }
  }

  // Toggle Quiet Mode
  if (qflag) {
    if (QUIET_MODE) {
      QUIET_MODE=0;
    } else {
      QUIET_MODE=1;
      VERBOSE=0;
    }
  }


  //If the symbol tables are not loaded in from the saved tree file
  // we need to load in the syumbol tables from a file.
  if ((globInputSymbols==NULL) || (globOutputSymbols==NULL)) {

    // Load in the symbols files and exit if inconsistent.
    globInputSymbols=loadSymbols(INPUT_SYMBOL_FILENAME);
    // Default to the number of symbols in the symbol file.
    if (NUM_INPUT_SYMBOLS==-1)
      NUM_INPUT_SYMBOLS=globInputSymbols->size;
    
    
    if (globInputSymbols->size!=NUM_INPUT_SYMBOLS) {
      fprintf(stderr,"Invalid Number of Symbols in Input Symbol File.\n");
      fprintf(stderr,"Global NUM_INPUT_SYMBOLS=%d\n",NUM_INPUT_SYMBOLS);
      fprintf(stderr,"Symbols in %s=%d\n",ivalue,globInputSymbols->size);
      exit(1);
    }
    globOutputSymbols=loadSymbols(OUTPUT_SYMBOL_FILENAME);
    // Default to the number of symbols in the symbol file.
    if (NUM_OUTPUT_SYMBOLS==-1)
    NUM_OUTPUT_SYMBOLS=globOutputSymbols->size;
    
    if (globOutputSymbols->size!=NUM_OUTPUT_SYMBOLS) {
      fprintf(stderr,"Invalid Number of Symbols in Output Symbol File.\n");
      fprintf(stderr,"Global NUM_OUTPUT_SYMBOLS=%d\n",NUM_OUTPUT_SYMBOLS);
      fprintf(stderr,"Symbols in %s=%d\n",ovalue,globOutputSymbols->size);
      exit(1);
    }
  }



  // set global variables.
  createStaticVars();
  //setBackgroundProbs();
  setSourceOffsets();
  
  // Print out the globals if verbose mode is true.
  if (VERBOSE==TRUE)
    printGlobals();


  // compute string kernel
  if (Kflag) {
    if (COMPUTE_KERNEL_MATRIX)
      COMPUTE_KERNEL_MATRIX=0;
    else
      COMPUTE_KERNEL_MATRIX=1;
  }
  if (Hflag) {
    if (PRINT_KMER_HISTOGRAM)
      PRINT_KMER_HISTOGRAM=0;
    else
      PRINT_KMER_HISTOGRAM=1;
  }

  globalRecords=readMultipleSequences(sequenceFiles, 
				      numInputFiles, 
				      globInputSymbols, 
				      globOutputSymbols);
  if (!QUIET_MODE) {
    printf("Num Total Instances = %d\n",numTotalRecords);
    printf("(LENGTH, MAX_MISMATCHES, SEQUENCES) = (%d,%d,%d)\n",MAX_DEPTH,MAX_MISMATCHES,NUM_OUTPUT_SYMBOLS);
  }
  
  readMutationMatrix(MUTATION_MATRIX_FILENAME);

  mismatchCounts=(int *)calloc(numTotalRecords,sizeof(int));
  mismatchWeights=(double *)calloc(numTotalRecords,sizeof(double));
  for (j=0; j<numTotalRecords; j++) {
    mismatchCounts[j]=0;
    mismatchWeights[j]=0.0;
  }

  if (COMPUTE_KERNEL_MATRIX) {
    initializeKernelMatrix();
  }

  traverse(globalRecords,mismatchCounts, mismatchWeights, numTotalRecords, 0);

  if (COMPUTE_KERNEL_MATRIX) {
    //printf("Kernel Matrix:\n");
    
    printf("example\t");
    for (i=0; i<NUM_OUTPUT_SYMBOLS; i++) {
      printf("%s\t",lookUpSymbol(globOutputSymbols,i));
    }
    printf("\n");

    //for (i=0; i<NUM_OUTPUT_SYMBOLS; i++) {
    //  printf("JUNK\t");
    //}
    //printf("\n");

    for (i=0; i<NUM_OUTPUT_SYMBOLS; i++) {
      printf("%s\t",lookUpSymbol(globOutputSymbols,i));
      for (j=0; j<NUM_OUTPUT_SYMBOLS; j++) {
	if (USE_WEIGHTED_MATRIX)
	  printf("%.2f\t",kernelMatrix[i][j]*KERNEL_OUTPUT_SCALING);
	else
	  printf("%4d\t",(int)kernelMatrix[i][j]);
      }
      printf("\n");
    }
  }
  
  

  // Free the data before ending.  This is uncommented when debugging
  // and looking for memory leaks.

  for (i=0; i<numTotalRecords; i++) {
    free(globalRecords[i]);
  }
  free(globalRecords);
  free(mismatchCounts);


  if (COMPUTE_KERNEL_MATRIX) {
    for (i=0; i<NUM_OUTPUT_SYMBOLS; i++) {
      free(kernelMatrix[i]);
    }
    free(kernelMatrix);
  }

  freeGlobals();
  freeSymbolTable(globInputSymbols);
  freeSymbolTable(globOutputSymbols);
  free(globalSequence);
  free(currentSequence);


  



  return 0;
}



