/***************************************************************************
 *
 * SymbolTable.h
 * 
 * Global Variables for SymbolTables
 *
 ***************************************************************************/

#define MAXLINE 100

#define UNKNOWN_SYMBOL_NUM NUM_INPUT_SYMBOLS
#define UNKNOWN_SYMBOL_SYM "UKS"

struct symbolTable {
  struct hashTable *symbols;
  struct hashTable *numbers;
  int size;
  double *updateAmounts;
};



struct symbolTable *loadSymbols(char *);
int mapSymbol(struct symbolTable *, char *);
char *lookUpSymbol(struct symbolTable *, int);
void printSymbolTable(struct symbolTable *);
void freeSymbolTable(struct symbolTable *);

struct symbolTable *readSymbolTable(FILE *);
void writeSymbolTable(struct symbolTable *, FILE *);
