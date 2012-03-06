/****************************************************************************
 *
 * Input.c
 *
 * Source file for handling Input Files.
 ****************************************************************************/


#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

#include "includes.h"

// This file handles I/O for the SST package.  
// It handles printing out input records and predictions.  
// It also handles reading in data from files or standard input.

/**************** Printing Input Records and Predictions *********************/


// printInputRecordInputs() prints out the input sequence
// corresponding to an input record.
void printInputRecordInputs(struct inputRecord *inputRec, 
		      struct symbolTable *inputSymbols, 
		      struct symbolTable *outputSymbols,
		      char *additionalComment) {

  int i;


  for (i=0; i<inputRec->length; i++) {
    printf("%s ", lookUpSymbol(inputSymbols,inputRec->inputs[i]));
  }
  printf("\t%s%s",additionalComment,inputRec->comment);
}

// pintInputRecordOutput() prints out the output symbol of an input record.
void printInputRecordSource(struct inputRecord *inputRec, 
		      struct symbolTable *inputSymbols, 
		      struct symbolTable *outputSymbols,
		      char *additionalComment) {

  printf("%s: ",lookUpSymbol(outputSymbols,inputRec->source));
}




/************** Reading Sequences *****************************************/

int getSequenceLength(FILE *fp) {

  int length=0;
  char tempBuffer[61];

  while (fscanf(fp,"%s",tempBuffer)!=EOF) {
    if (tempBuffer[0]!='>')
      length++;
  }
  return length;
}

int getNumSources(FILE *fp) {

  int sources=0;
  char tempBuffer[61];

  while (fscanf(fp,"%s",tempBuffer)!=EOF) {
    if (tempBuffer[0]=='>')
      sources++;
  }
  return sources;
}


/************** Reading Multiple Sequences **********************************/

// TODO RETURN NUM RECORDS THROUGH POINTER AND CHEANGE IN GLOABL SO IT DEOS NOT CRASH.

struct inputRecord **readInputMultSequence(FILE *fp,
					   struct symbolTable *inputSymbols,
					   struct symbolTable *outputSymbols, 
					   int sequenceLength,
					   int sequenceNumber) {

  int i,j;
  int numRecords=sequenceLength;
  int *sequence=(int *)calloc(sequenceLength,sizeof(int));
  int *sources=(int *)calloc(sequenceLength,sizeof(int));
  int currentSource=-1;
  int position;
  int *newSequence;
  char tempSymbols[61];
  struct inputRecord **records
    =(struct inputRecord **)calloc(numRecords,sizeof(struct inputRecord *));
  

  for (i =0; i<sequenceLength; i++) {

    fscanf(fp,"%s",tempSymbols);
    if (tempSymbols[0]!='>') {
      if (currentSource==-1) {
	fprintf(stderr,"Invalid Input File\n");
	exit(1);
      }
      sequence[i]=mapSymbol(inputSymbols,tempSymbols);
      sources[i]=currentSource;
    } else {
      currentSource=mapSymbol(outputSymbols,tempSymbols+1);
      numRecords-=MAX_DEPTH-1;
      if (i>MAX_DEPTH) {
	for (j=0; j<MAX_DEPTH-1; j++) {
	  sources[i-j-1]=-1;
	}
      }
      i--;
    }
  }
  if (i>MAX_DEPTH) {
    for (j=0; j<MAX_DEPTH-1; j++) {
      sources[i-j-1]=-1;
    }
  }
  

  

  newSequence=(int *)calloc(sequenceLength+globalSequenceLength,sizeof(int));
  for (i=0; i<globalSequenceLength; i++) {
    newSequence[i]=globalSequence[i];
  }
  for (i=globalSequenceLength; i<globalSequenceLength+sequenceLength; i++) {
    newSequence[i]=sequence[i-globalSequenceLength];
  }
  if (globalSequenceLength>0) {
    free(globalSequence);
  }
  globalSequence=newSequence;

  j=0;
  position=0;
  for (i=0; i<numRecords; i++) {
    while (sources[j]==-1) {
      j++;
      position=0;
    }
    position++;

    records[i]=(struct inputRecord *)malloc(sizeof(struct inputRecord));
    records[i]->source=sources[j];
    records[i]->inputs=&sequence[j];
    records[i]->position=position;
    records[i]->length=MAX_DEPTH;
    records[i]->comment=NULL;
    j++;
  }
  

  globalSequenceLength+=sequenceLength;


  return records;
}

struct inputRecord **readMultipleSequences(char **fileNames,
					   int numFiles,
					   struct symbolTable *inputSymbols, 
					   struct symbolTable *outputSymbols) {
  
  int sequenceLength=0;
  struct inputRecord **records,**newRecords,**oldRecords;
  int i,j;
  int numSources,numRecords;
  FILE *fp;

  numTotalRecords=0;

  globalSequenceLength=0;



  for (j=0; j<numFiles; j++) {

    if (VERBOSE)
      printf("Current File Name = %s\n",fileNames[j]);

    if ((fp = fopen(fileNames[j],"r"))==NULL) {
      fprintf(stderr, "Invalid Input Sequence File:%s\n",fileNames[j]);
      exit(1);
    }
    
    sequenceLength=getSequenceLength(fp);
    fclose(fp);

    if ((fp = fopen(fileNames[j],"r"))==NULL) {
      fprintf(stderr, "Invalid Input Sequence File:%s\n",fileNames[j]);
      exit(1);
    }
    
    numSources=getNumSources(fp);
    fclose(fp);

    if ((fp = fopen(fileNames[j],"r"))==NULL) {
      fprintf(stderr, "Invalid Input Sequence File:%s\n",fileNames[j]);
      exit(1);
    }

    numRecords=sequenceLength-(numSources*(MAX_DEPTH-1));
    records=readInputMultSequence(fp,inputSymbols,outputSymbols,sequenceLength,j);

    oldRecords=newRecords;
    newRecords=calloc(numTotalRecords+numRecords,sizeof(struct inputRecord *));
    for (i=0; i<numTotalRecords; i++) {
      newRecords[i]=oldRecords[i];
    }
    if (numTotalRecords>0)
      free(oldRecords);

    for (i=numTotalRecords; i<numTotalRecords+numRecords; i++) {
      
      newRecords[i]=records[i-numTotalRecords];
      if (VERBOSE && i % 1000==0) {
	printf("Reading Example %d\n",i);
      }
    }
    free(records);
    numTotalRecords+=numRecords;
  }
  return newRecords;
} 



