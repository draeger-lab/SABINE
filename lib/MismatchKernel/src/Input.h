/****************************************************************************
 *
 * Input.h
 *
 * Header file for handling Input Files.
 ****************************************************************************/



struct inputRecord {

  int source;
  int *inputs;
  int position;
  char length;
  char *comment;
};



struct inputRecord **readMultipleSequences(char **fileNames,
					   int numFiles,
					   struct symbolTable *inputSymbols, 
					   struct symbolTable *outputSymbols);


//void readSequence(struct node *currentTree, 
//		  char *fileName, 
//		  struct symbolTable *inputSymbols, 
//		  struct symbolTable *outputSymbols,
//		  int update);


