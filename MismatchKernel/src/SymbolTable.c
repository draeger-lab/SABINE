/***************************************************************************
 *
 * SymbolTable.
 * 
 * Implementation of SymbolTables
 *
 ***************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>


#include "includes.h"

// This file handles the reading in of symbol files as well as reading
// and writing them to disk.  The file also handles the mapping of
// symbols to integers and vice versa.


//This method gets the number of symbols in a symbol file.
//This is done by counting the symbols in the file.
int getSymbolTableSize(char *fileName) {

  int length=0;
  char line[MAXLINE+1];

  FILE *fp;

  if ((fp = fopen(fileName,"r")) == NULL) {
    fprintf(stderr, "Error opening Symbol File: %s\n",fileName);
    exit(1);
  }

  while (fgets(line,MAXLINE-1,fp)!=NULL) {

    if (strlen(line)>1) {
      length++;
    } else {
      fprintf(stderr, "Invalid Symbol File Entry: %s\n",line);
      fclose(fp);
      exit(1);
    }
  }
  fclose(fp);
  return length;
}


// charHashValue() computes a hashvalue for a string.
int charHashValue(char *string) {

  int i=0;
  int total=0;
  int length=strlen(string);


  for (i=0; i<length; i++) {
    total=((total*50)+string[i]-'a'+1)%0xfffffff;
  }
  return abs(total);
}

// loadSymbols() loads a symbol table from a symol file.  The second
// column in the symbol file gives the update amount of the symbol
// used if WEIGHTED_COUNTS is true.
struct symbolTable *loadSymbols(char *fileName) {
  
 
  int i;
  struct symbolTable *symbols;
  char line[MAXLINE+1];
  char *token;
  char *newChar;
  double updateAmount;
  
  FILE *fp;

  int length = getSymbolTableSize(fileName);


  // Error checking.
  if (length==0) {
    fprintf(stderr, "Invalid Symbol File:%s\n",fileName);
    exit(1);
    return NULL;
  }
  if ((fp = fopen(fileName,"r"))==NULL) {
    fprintf(stderr, "Invalid Symbol File:%s\n",fileName);
    exit(1);
    return NULL;
  }


  symbols = (struct symbolTable *) malloc(sizeof(struct symbolTable));

  symbols->symbols=newHashTable(1,CHAR_HASH);
  symbols->numbers=newHashTable(1,INT_HASH);

  symbols->updateAmounts=calloc(length+1,sizeof(double));

  //Put unknown symbol in symbol table.
  strcpy(line,UNKNOWN_SYMBOL_SYM);
  newChar=(char *) calloc(strlen(line)+1,sizeof(char));
  strcpy(newChar,line);
  putHashValueC(symbols->symbols,length,newChar);
  putHashValueI(symbols->numbers,charHashValue(newChar),length);
  symbols->updateAmounts[length]=1.0;

  for (i=0; i<length; i++) {
    fgets(line,MAXLINE,fp);


    token=strtok(line," \t\n");


    newChar=(char *) calloc(strlen(token)+1,sizeof(char));
    strcpy(newChar,token);

    
    putHashValueC(symbols->symbols,i,newChar);
    putHashValueI(symbols->numbers,charHashValue(newChar),i);

    token=strtok(NULL," \t\n");
    if (token!=NULL) {
      sscanf(token,"%lf",&updateAmount);
      symbols->updateAmounts[i]=updateAmount;
    } else {
      symbols->updateAmounts[i]=1.0;
    }
  }
  symbols->size=length;
  fclose(fp);

  return symbols;
}


// mapSymbol() converts a string to an integer using the symbolTable.
// If the sybol is not found, it exits the program.
int mapSymbol(struct symbolTable *symbols, char *symbol) {

  int value=getHashValueI(symbols->numbers,charHashValue(symbol));
  
  if (value!=-1)
    return value;
  else {
    fprintf(stderr,"Unknown Symbol: %s\n",symbol);
    exit(1);
  }
  return -1;
}

// lookUpSymbol() converts an integer to its appropriate string.
char *lookUpSymbol(struct symbolTable *symbols, int value) {

  return getHashValueC(symbols->symbols,value);
}

// printSymbolTable() prints the symbol table (for debugging purposes)
void printSymbolTable(struct symbolTable *symbols) {

  printHashTable(symbols->symbols);

}

// freeSymbolTable() frees the memory for a symbol table.
void freeSymbolTable(struct symbolTable *symbols) {

  freeData(symbols->symbols->data,
	   symbols->symbols->numBuckets,
	   symbols->symbols->type,
	   1);
  free(symbols->symbols);
  freeData(symbols->numbers->data,
	   symbols->numbers->numBuckets,
	   symbols->numbers->type,
	   1);
  free(symbols->numbers);
  free(symbols->updateAmounts);
  free(symbols);

}

// readSymbolTable() reads a symbol table from a file
struct symbolTable  *readSymbolTable(FILE *fp) {

  struct symbolTable *newSymbolTable;

  newSymbolTable = (struct symbolTable *) malloc(sizeof(struct symbolTable));
  fread(newSymbolTable,sizeof(struct symbolTable),1,fp);
  newSymbolTable->symbols=newHashTable(1,CHAR_HASH);
  newSymbolTable->numbers=newHashTable(1,INT_HASH);

  readHashTableI(newSymbolTable->numbers, fp);
  readHashTableC(newSymbolTable->symbols, fp);

  newSymbolTable->updateAmounts=calloc(newSymbolTable->size+1,sizeof(double));

  fread(newSymbolTable->updateAmounts,sizeof(double),newSymbolTable->size+1,fp);
  return newSymbolTable;
}

// writeSymbolTable() writes a symbol table to a file.
void writeSymbolTable(struct symbolTable *currentSymbolTable, FILE *fp) {

  fwrite(currentSymbolTable, sizeof(struct symbolTable),1,fp);
  writeHashTableI(currentSymbolTable->numbers,fp);
  writeHashTableC(currentSymbolTable->symbols,fp);
  fwrite(currentSymbolTable->updateAmounts,sizeof(double),currentSymbolTable->size+1,fp);

}
