/***************************************************************************
 *
 * HashTable.c
 * 
 * Implementation of HashTable
 *
 ***************************************************************************/


#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "includes.h"

// This file implements hash tables which are used everywhere in the
// system.  The basic hashtable is the same.  Different hash table
// types store values of different types and have different interfaces
// to them.

#define INIT_BUCKETS 1

// newHashRecord() creates and allocates memory for a new hashRecord
// where the value is a void pointer.
struct hashRecord *newHashRecord(int key, void *value, struct hashRecord *next) {

  struct hashRecord *newRecord = (struct hashRecord *) malloc(sizeof(struct hashRecord));
  newRecord->key=key;
  newRecord->value.ptr=value;
  newRecord->next=next;
  return(newRecord);
}
    
// newHashRecordI() creates and allocates memory for a new hashRecord
// where the value is an integer.
struct hashRecord *newHashRecordI(int key, int value, struct hashRecord *next) {

  struct hashRecord *newRecord = (struct hashRecord *) malloc(sizeof(struct hashRecord));
  newRecord->key=key;
  newRecord->value.num=value;
  newRecord->next=next;
  return(newRecord);
}

// newHashRecordD() creates and allocates memory for a new hashRecord
// where the value is a double.
struct hashRecord *newHashRecordD(int key, double value, struct hashRecord *next) {

  struct hashRecord *newRecord = (struct hashRecord *) malloc(sizeof(struct hashRecord));
  newRecord->key=key;
  newRecord->value.dnum=value;
  newRecord->next=next;
  return(newRecord);
}

// getHashRecord() is a function that returns the hashrecord
// corresponding to a key or it returns null if there is no record
// corresponding to the key.
struct hashRecord *getHashRecord(struct hashTable *hash, int key) 
{
  struct hashRecord *currentBucket=hash->data[key%(hash->numBuckets)];


  while (currentBucket != NULL) {

    if (currentBucket->key==key) {
      return(currentBucket);
    } else {
      currentBucket=currentBucket->next;
    }
  }
  return(NULL);
}

// the next 4 functions return a value from a hash table.
void *getHashValue(struct hashTable *hash, int key) {
  struct hashRecord *currentBucket=getHashRecord(hash,key);
  if (currentBucket==NULL) {
    return(NULL);
  } else {
    return (currentBucket->value.ptr);
  }
}

int getHashValueI(struct hashTable *hash, int key) {
  struct hashRecord *currentBucket=getHashRecord(hash,key);
  
  if (currentBucket==NULL) {
    return -1;
  } else {
    return (currentBucket->value.num);
  }

}

double getHashValueD(struct hashTable *hash, int key) {
  struct hashRecord *currentBucket=getHashRecord(hash,key);
  
  if (currentBucket==NULL) {
    return -1;
  } else {
    return (currentBucket->value.dnum);
  }

}

// getHashTotalI() returns the sum of all of the elements in a hash
// table storing integers.
int getHashTotalI(struct hashTable *hash) {

  struct hashRecord *currentRecord;
  int i,total;

  total=0;
  for (i=0; i<hash->numBuckets; i++) {
    currentRecord=hash->data[i];
    while (currentRecord != NULL) {
      total+=currentRecord->value.num;
      currentRecord=currentRecord->next;
    }
  }
  return total;
}

// removeHashValue() removes a record from a hash table.
void removeHashValue(struct hashTable *hash, int key) {

  struct hashRecord **prevRecord, *currentRecord;
  int hashVal = key % hash->numBuckets;


  prevRecord = &hash->data[hashVal];
  currentRecord = hash->data[hashVal];
  while (currentRecord !=NULL) {
    if (currentRecord->key==key) {
      *prevRecord = currentRecord->next;
      free(currentRecord);
      return;
    }
    prevRecord = &currentRecord->next;
    currentRecord=currentRecord->next;
  }
}

// putHashValue() puts a value in a hashtable.  It grows the hashtable
// if necessary.
void putHashValue(struct hashTable *hash, int key, void *value)
{
  struct hashRecord *newRecord;
  int hashVal = key % hash->numBuckets;

  if ((newRecord = getHashRecord(hash, key)) == NULL) {
    newRecord = newHashRecord(key, value, hash->data[hashVal]);
    hash->data[hashVal]=newRecord;
    hash->numRecords++;
    if (hash->numRecords>hash->numBuckets) {
      growHashTable(hash);
    }
  } else {
    newRecord->value.ptr=value;
  }
}

// putHashValueI() puts an integer into a hash table.
void putHashValueI(struct hashTable *hash, int key, int value)
{
  struct hashRecord *newRecord;
  int hashVal = key % hash->numBuckets;

  if ((newRecord = getHashRecord(hash, key)) == NULL) {
    newRecord = newHashRecordI(key, value, hash->data[hashVal]);
    hash->data[hashVal]=newRecord;
    hash->numRecords++;
    if (hash->numRecords>hash->numBuckets) {
      growHashTable(hash);
    }
  } else {
    newRecord->value.num=value;
  }
}

// putHashValueD() puts an double into a hash table.
void putHashValueD(struct hashTable *hash, int key, double value)
{
  struct hashRecord *newRecord;
  int hashVal = key % hash->numBuckets;

  if ((newRecord = getHashRecord(hash, key)) == NULL) {
    newRecord = newHashRecordD(key, value, hash->data[hashVal]);
    hash->data[hashVal]=newRecord;
    hash->numRecords++;
    if (hash->numRecords>hash->numBuckets) {
      growHashTable(hash);
    }
  } else {
    newRecord->value.dnum=value;
  }
}

// incHashValueI() increases the integer value corresponding to a key
// by 1.
void incHashValueI(struct hashTable *hash, int key)
{
  int currentValue=getHashValueI(hash,key);
  if (currentValue==-1) //If it is unseen, the current value is 0
    currentValue=0;
  putHashValueI(hash,key,currentValue+1);
}

// incHashValueD() increases the double value corresponding to a key
// by the amount.
void incHashValueD(struct hashTable *hash, int key, double amount)
{
  double currentValue=getHashValueD(hash,key);
  if (currentValue==-1) //If it is unseen, the current value is 0
    currentValue=0;
  putHashValueD(hash,key,currentValue+amount);
}
  

char *getHashValueC(struct hashTable *hash, int key) {
  void *hashValue=getHashValue(hash,key);
  char *retVal;

  if (hashValue==NULL) {
    return(NULL);
  } else {
    retVal=(char *)hashValue;
    return retVal;
  }
}

void putHashValueC(struct hashTable *hash, int key, char *value)
{
  void *hashValue;

  hashValue=(void *)value;

  putHashValue(hash,key,hashValue);
}


// initHashTable() initializes a hash table.
void initHashTable(struct hashTable *hash, int initBuckets, int TYPE) 
{

  struct hashRecord **data;
    
  hash->numBuckets=initBuckets;
  hash->numRecords=0;
  hash->type=TYPE;

  data=calloc(hash->numBuckets,sizeof(struct hashRecord *));
  
  hash->data=data;
}

// newHashTable() allocates memory for a hash table and initializes
// it.
struct hashTable *newHashTable(int initBuckets, int TYPE)
{


  struct hashTable *hash;
  

  hash=(struct hashTable *) malloc(sizeof(struct hashTable));

  initHashTable(hash, initBuckets, TYPE);
  return hash;
}

// copyHashTableSameData copies a hashtable keeping the same data.
void copyHashTableSameData(struct hashTable *hash, struct hashTable *newHash) {


  newHash->numBuckets=hash->numBuckets;
  newHash->numRecords=hash->numRecords;
  newHash->type=hash->type;
  newHash->data=hash->data;

  
}

// the next 6 functions are for initializing or creating specific
// types of hash tables.

void initHashTableI(struct hashTable *hash) {
  initHashTable(hash,INIT_BUCKETS,INT_HASH);
}

void initHashTableD(struct hashTable *hash) {
  initHashTable(hash,INIT_BUCKETS,DOUBLE_HASH);
}


// growHashTable() grows the hash table and puts all of the elements
// into the right buckets.
void growHashTable(struct hashTable *hash) {

  int i;
  struct hashRecord *currentRecord;
  struct hashRecord **oldData=hash->data;  // This stores the old data
  int oldNumBuckets=hash->numBuckets;      // The number of buckets in the old data.
  struct hashRecord **newData=(struct hashRecord **)calloc(hash->numBuckets*2,sizeof(struct hashRecord *));
  hash->data=newData;

  
  hash->numBuckets=hash->numBuckets*2;

  hash->numRecords=0;       // We set this to 0 because we are going 
                            // to reinsert all of the data.

  for (i=0; i<oldNumBuckets; i++) {
    currentRecord=oldData[i];
    while (currentRecord !=NULL) {
      if (hash->type==INT_HASH)
	putHashValueI(hash, currentRecord->key, currentRecord->value.num);
      else if (hash->type==DOUBLE_HASH)
	putHashValueD(hash, currentRecord->key, currentRecord->value.dnum);
      else if (hash->type==CHAR_HASH)
	putHashValueC(hash, currentRecord->key, (char *)currentRecord->value.ptr);
      currentRecord=currentRecord->next;
    }
  }
  freeData(oldData,oldNumBuckets,hash->type,FALSE);   // We free up the old data.


}

// freeData() frees the memory used by the data of a hash table.
void freeData(struct hashRecord **data, int numBuckets,int type, int freeAll) {

  int i;
  struct hashRecord *currentRecord,*nextRecord;
  
  for (i=0; i<numBuckets; i++) {
    currentRecord=data[i];
    while (currentRecord !=NULL) {
      nextRecord=currentRecord->next;
      if (type==INT_HASH) {
	//We do not need to do anything since it is not a pointer.
      }
      if (freeAll) {
	// Here we are freeing the data.
	if (type==CHAR_HASH) {
	  free(currentRecord->value.ptr);
	}
      }
      free(currentRecord);
      currentRecord=nextRecord;
    }
  }
  free(data);


}

// printHashTable() prints the values in a hashtable.  (Used only for
// debugging.)
void printHashTable(struct hashTable *hash) {

  int i;
  struct hashRecord *currentRecord;


  for(i=0; i<hash->numBuckets; i++) {

    printf("%d: ",i);
    currentRecord = hash->data[i];
    while (currentRecord!=NULL) {
      if (hash->type==INT_HASH)
	printf("%d=%d, ",currentRecord->key, currentRecord->value.num);
      if (hash->type==DOUBLE_HASH)
	printf("%d=%f, ",currentRecord->key, currentRecord->value.dnum);
      currentRecord=currentRecord->next;
    }

    printf("\n");
  }
}


// the next 6 functions deal with reading and writing hash tables to
// disk.
void writeHashTableI(struct hashTable *hash, FILE *fp) {

  int i;
  struct hashRecord *currentRec;


  fwrite(&hash->numBuckets, sizeof(int),1, fp);

  for (i=0; i<hash->numBuckets; i++) {
    
    currentRec=hash->data[i];
    while (currentRec!=NULL) {
      putc(254,fp); //code that there are more records.
      fwrite(currentRec, sizeof(struct hashRecord), 1, fp);
      fwrite(&currentRec->value.num, sizeof(int), 1, fp);
      currentRec=currentRec->next;
    }
    putc(255,fp);
  }
}

void readHashTableI(struct hashTable *newHash, FILE *fp) {

  int i, signal, *numBuckets;
  struct hashRecord *newRec,**currentRecP;


  numBuckets=malloc(sizeof(int));

  fread(numBuckets, sizeof(int), 1, fp);
  

  newHash->numBuckets=*numBuckets;
  free(numBuckets);

  newHash->data=(struct hashRecord **) calloc(newHash->numBuckets,sizeof(struct hashRecord *));

  for (i=0; i<newHash->numBuckets; i++) {

    currentRecP=&newHash->data[i];
    signal=getc(fp);

    while (signal!=255) {
      newRec=(struct hashRecord *)malloc(sizeof(struct hashRecord));      
      fread(newRec,sizeof(struct hashRecord),1,fp);

      fread(&newRec->value.num,sizeof(int),1,fp);
      *currentRecP=newRec;
      currentRecP=&newRec->next;
      *currentRecP=NULL;
      signal=getc(fp);

    }
    
  }

}   

void writeHashTableD(struct hashTable *hash, FILE *fp) {

  int i;
  struct hashRecord *currentRec;


  putc(hash->numBuckets,fp);

  for (i=0; i<hash->numBuckets; i++) {
    
    currentRec=hash->data[i];
    while (currentRec!=NULL) {
      putc(254,fp); //code that there are more records.
      fwrite(currentRec, sizeof(struct hashRecord), 1, fp);
      fwrite(&currentRec->value.dnum, sizeof(double), 1, fp);
      currentRec=currentRec->next;
    }
    putc(255,fp);
  }
}

void readHashTableD(struct hashTable *newHash, FILE *fp) {

  int i, signal, numBuckets;
  struct hashRecord *newRec,**currentRecP;

  numBuckets=getc(fp);

  newHash->numBuckets=numBuckets;

  newHash->data=(struct hashRecord **) calloc(newHash->numBuckets,sizeof(struct hashRecord *));

  for (i=0; i<newHash->numBuckets; i++) {

    currentRecP=&newHash->data[i];
    signal=getc(fp);

    while (signal!=255) {
      newRec=(struct hashRecord *)malloc(sizeof(struct hashRecord));      
      fread(newRec,sizeof(struct hashRecord),1,fp);

      fread(&newRec->value.dnum,sizeof(double),1,fp);
      *currentRecP=newRec;
      currentRecP=&newRec->next;
      *currentRecP=NULL;
      signal=getc(fp);

    }
    
  }

}   
 
void writeHashTableC(struct hashTable *hash, FILE *fp) {

  int i;
  struct hashRecord *currentRec;


  putc(hash->numBuckets,fp);

  for (i=0; i<hash->numBuckets; i++) {
    
    currentRec=hash->data[i];
    while (currentRec!=NULL) {
      putc(254,fp); //code that there are more records.
      fwrite(currentRec, sizeof(struct hashRecord), 1, fp);
      fwrite(currentRec->value.ptr, sizeof(char), 101, fp);
      currentRec=currentRec->next;
    }
    putc(255,fp);
  }
}

void readHashTableC(struct hashTable *newHash, FILE *fp) {

  int i, signal, numBuckets;
  struct hashRecord *newRec,**currentRecP;
  char inputBuffer[101];

  numBuckets=getc(fp);

  newHash->numBuckets=numBuckets;

  newHash->data=(struct hashRecord **) calloc(newHash->numBuckets,sizeof(struct hashRecord *));

  for (i=0; i<newHash->numBuckets; i++) {

    currentRecP=&newHash->data[i];
    signal=getc(fp);


    while (signal!=255) {
      newRec=(struct hashRecord *)malloc(sizeof(struct hashRecord));      
      fread(newRec,sizeof(struct hashRecord),1,fp);

      fread(inputBuffer,sizeof(char),101,fp);
      newRec->value.ptr=calloc(strlen(inputBuffer+1),sizeof(char));
      strcpy(newRec->value.ptr,inputBuffer);

      *currentRecP=newRec;
      currentRecP=&newRec->next;
      *currentRecP=NULL;
      signal=getc(fp);

    }
    
  }

}   

     
