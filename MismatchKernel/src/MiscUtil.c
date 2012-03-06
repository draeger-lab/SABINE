/******************************************************************************
 *
 * MiscUtil.h
 *
 * Utilities that should have implementations somewhere.
 *****************************************************************************/


#include <stdlib.h>

// this file provides implementations of auxillary functions.

  
int min(int a, int b) {

  if (a<b)
    return a;
  else
    return b;
}

int max(int a, int b) {
  if (a<b)
    return b;
  else
    return a;
}
