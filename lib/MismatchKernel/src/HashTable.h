/***************************************************************************
 *
 * HashTable.h
 * 
 * Global Variables for HashTables
 *
 ***************************************************************************/


#include <stdio.h>

#define INITIAL_BUCKETS 2

#define NUMBUCKETS 10

#define INT_HASH 0
#define DOUBLE_HASH 2
#define CHAR_HASH 3



struct hashRecord {

  int key;
  union {
    void *ptr;
    int num;
    double dnum;
  } value;
  struct hashRecord *next;
};

struct hashTable {

  struct hashRecord **data;
  unsigned int numBuckets;
  unsigned int numRecords : 28;
  unsigned int type : 4;
};

void removeHashValue(struct hashTable *, int);

void initHashTable(struct hashTable *, int, int);
void initHashTableN(struct hashTable *);
void initHashTableI(struct hashTable *);
void initHashTableD(struct hashTable *);

struct hashTable *newHashTable();
void *getHashValue(struct hashTable *, int); 
void putHashValue(struct hashTable *, int, void *);

struct hashRecord *getHashRecord(struct hashTable *, int); 

void copyHashTableSameData(struct hashTable *,struct hashTable *);

struct hashTable *newHashTableI();
int getHashValueI(struct hashTable *, int); 
void putHashValueI(struct hashTable *, int, int);
void incHashValueI(struct hashTable *, int);
int getHashTotalI(struct hashTable *);

struct hashTable *newHashTableD();
double getHashValueD(struct hashTable *, int); 
void putHashValueD(struct hashTable *, int, double);
void incHashValueD(struct hashTable *, int, double);



struct hashTable *newHashTableC();
char *getHashValueC(struct hashTable *, int); 
void putHashValueC(struct hashTable *, int, char *);



void printHashTable(struct hashTable *);
struct hashRecord *newHashRecord(int, void *, struct hashRecord *);

void writeHashTableI(struct hashTable *, FILE *);
void readHashTableI(struct hashTable *, FILE *);

void writeHashTableD(struct hashTable *, FILE *);
void readHashTableD(struct hashTable *, FILE *);

void writeHashTableC(struct hashTable *, FILE *);
void readHashTableC(struct hashTable *, FILE *);

void growHashTable(struct hashTable *);
void freeData(struct hashRecord **, int, int, int);
