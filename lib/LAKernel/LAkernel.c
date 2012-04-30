/*********************************************************************
    LAkernel.c (version 0.3)                                          
                                                                    
    Compute the local alignment kernel

    Reference:
    H. Saigo, J.-P. Vert, T. Akutsu and N. Ueda, "Protein homology 
    detection using string alignment kernels", Bioinformatics, 
    vol.20, p.1682-1689, 2004.
                                                                    
                                                                    
    Copyright 2003 Jean-Philippe Vert                                 
    Copyright 2005 Jean-Philippe Vert, Hiroto Saigo                                 

    This file is part of LAkernel.

    LAkernel is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    Foobar is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Foobar; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

*********************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <ctype.h>
#include <string.h>

/* Usefule constants */
#define LOG0 -100000000         /*  log(0) */
#define INTSCALE 1000.0      /* critical for speed and precise computation*/
#define LOGSUM_TBL 10000      /* span of the logsum table */ 
#define LOGDELTA(x,y,i,j) ((x==i && y ==j) || (x==j && y==i)) ? 0 : LOG0

#define NAA 20                                  /* Number of amino-acids */
#define NLET 26                        /* Number of letters in the alphabet */
static char *aaList= "ARNDCQEGHILKMFPSTWYV";    /* The list of amino acids */
/* Index corresponding to the (i,j) entry (i,j=0..19) in the blosum matrix */
#define BINDEX(i,j) ( ((i)>(j)) ? (j)+(((i)*(i+1))/2) : (i)+(((j)*(j+1))/2))


/* You can switch here between SW and LA */

#define LOGP(x,y) (((x)>(y))?(x):(y)) /*Comment out if you need SW score */
/*#define LOGP(x,y) LogSum(x,y)  Comment if you need SW score */

/*********************/
/* Static variables  */
/*********************/


static int *isAA;                /* Indicates whether a char is an amino-acid */

static int *aaIndex;             /* The correspondance between amino-acid letter and index */

int opening,extension; /* Gap penalties */

/*********************/
/* Private Functions */
/*********************/



void terror(text)
  char *text;
  /* Print an error message and quit */
{
  fprintf(stderr,"\n\nERROR: %s.\n\n",text);
  exit(-1);
}


int *ArrayNew(size)
     /* Allocate an array of 'int' of size 'size' */
     int size;
{
  int *d;  
  if ((d=(int *)malloc(size*sizeof(int))) == NULL)
    terror("run out of memory");
  return d;
}

/* LogSum - default log funciotion. fast, but not exact */
/* LogSum2 - precise, but slow. Note that these two functions need different figure types  */

static int logsum_lookup[LOGSUM_TBL];
static void init_logsum(void){
  int i;
  for (i = 0; i < LOGSUM_TBL; i++) 
    logsum_lookup[i] = (int) (INTSCALE*
			       (log(1.+exp( (float) -i/INTSCALE))));
}
int LogSum(int p1, int p2){
  int diff;
  static int firsttime = 1;
  if (firsttime) {init_logsum(); firsttime = 0;}
  diff = p1 - p2;
  if      (diff >=  LOGSUM_TBL) return p1;
  else if (diff <= -LOGSUM_TBL) return p2;
  else if (diff > 0)            return p1+logsum_lookup[diff];
  else                          return p2+logsum_lookup[-diff];
}


float LogSum2(float p1, float p2)
{
  if (p1 > p2)
    return (p1-p2 > 50.) ? p1 : p1 + log(1. + exp(p2-p1));
  else
    return (p2-p1 > 50.) ? p2 : p2 + log(1. + exp(p1-p2));
}


void initialize(double *params)
     /* Initialize all static variables. This function should be called once before computing the first pair HMM score */
{
  register int i;
  double beta, gopn, gext;

  beta = params[210];
  gopn = params[211];
  gext = params[212];

  /* Initialization of the array which gives the position of each amino-acid in the set of amino-acid */
  if ((aaIndex=(int *)calloc(NLET,sizeof(int))) == NULL)
    terror("run out o memory");
  for (i=0;i<NAA;i++) 
    aaIndex[aaList[i]-'A']=i;
  
  /* Initialization of the array which indicates whether a char is an amino-acid */
  if ((isAA=(int *)calloc(256,sizeof(int))) == NULL)
    terror("run out of memory");
  for (i=0;i<NAA;i++) 
    isAA[(int)aaList[i]]=1;

  for (i=0 ; i<NAA*(NAA+1)/2; params[i++] *= beta *INTSCALE);

  /* Scale of gap penalties */
  opening = gopn * beta *INTSCALE;
  extension = gext * beta *INTSCALE;
}



void LAkernelcompute(aaX,aaY,nX,nY,params,kernel)
     /* Implementation of the convolution kernel which generalizes the Smith-Waterman algorithm */
     int *aaX , *aaY;   /* the two amino-acid sequences (as sequences of indexes in [0..NAA-1] indicating the position of the amino-acid in the variable 'aaList') */
     int nX , nY;       /* the lengths of both sequences */
     double *params, *kernel;
{
  double beta = params[210], gopn = params[211], gext = params[212];
   register int
    i,j,                /* loop indexes */
    cur, old,           /* to indicate the array to use (0 or 1) */
    curpos, frompos;    /* position in an array */
      /* reminder aa order  A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X */
   int           /* arrays to store the log-values of each state */
     *logX,    *logY,    *logM,    *logX2,    *logY2,
     aux , aux2;
     int
    cl;                /* length of a column for the dynamic programming */

  /*
  printf("now computing pairHMM between %d and %d:\n",nX,nY);
  for (i=0;i<nX;printf("%d ",aaX[i++]));
  printf("\n and \n");
  for (i=0;i<nY;printf("%d ",aaY[i++]));
  printf("\n");
  */

  /* Initialization of the arrays */
  /* Each array stores two successive columns of the (nX+1)x(nY+1) table used in dynamic programming */
  cl = nY+1;           /* each column stores the positions in the aaY sequence, plus a position at zero */

  logM = ArrayNew(2*cl);  logX = ArrayNew(2*cl);  logY = ArrayNew(2*cl); logX2 = ArrayNew(2*cl);  logY2 = ArrayNew(2*cl);

  /************************************************/
  /* First iteration : initialization of column 0 */
  /************************************************/
  /* The log=proabilities of each state are initialized for the first column (x=0,y=0..nY) */

  for (j=0;j<cl;j++) {
    logM[j] = logX[j] = logY[j] = logX2[j] = logY2[j] = LOG0;
  }

  /* Update column order */
  cur = 1;      /* Indexes [0..cl-1] are used to process the next column */
  old = 0;      /* Indexes [cl..2*cl-1] were used for column 0 */


  /************************************************/
  /* Next iterations : processing columns 1 .. nX */
  /************************************************/

  /* Main loop to vary the position in aaX : i=1..nX */
  for (i=1;i<=nX;i++) {

    /* Special update for positions (i=1..nX,j=0) */
    curpos = cur*cl;                  /* index of the state (i,0) */
    logM[curpos] = logX[curpos] = logY[curpos] = LOG0; logX2[curpos] = logY2[curpos] = LOG0;

    /* Secondary loop to vary the position in aaY : j=1..nY */
    for (j=1;j<=nY;j++) {

      curpos = cur*cl + j;            /* index of the state (i,j) */

      /* Update for states which emit X only */
      /***************************************/

      frompos = old*cl + j;            /* index of the state (i-1,j) */
      
      /* State RX */
      logX[curpos] = LOGP( - opening + logM[frompos] , - extension + logX[frompos] );

      /* State RX2 */
      logX2[curpos] = LOGP( logM[frompos] , logX2[frompos] );

      /* Update for states which emit Y only */
      /***************************************/

      frompos = cur*cl + j-1;          /* index of the state (i,j-1) */

      /* State RY */
      aux = LOGP( - opening + logM[frompos] , - extension + logY[frompos] );
      logY[curpos] = LOGP( aux , - opening + logX[frompos] );
      
      /* State RY2 */
      aux = LOGP( logM[frompos] , logY2[frompos] );
      logY2[curpos] = LOGP( aux , logX2[frompos] );
	
      /* Update for states which emit X and Y */
      /****************************************/

      frompos = old*cl + j-1;          /* index of the state (i-1,j-1) */

      aux = LOGP( logX[frompos] , logY[frompos] );
      aux2 = LOGP( 0 , logM[frompos] );
      logM[curpos] = LOGP( aux , aux2 ) + params[BINDEX(aaX[i-1], aaY[j-1])];

    }  /* end of j=1:nY loop */

    /* Update the culumn order */
    cur = 1-cur;
    old = 1-old;

  }  /* end of j=1:nX loop */


  /* Termination */
  /***************/

  curpos = old*cl + nY;                /* index of the state (nX,nY) */
  aux = LOGP( logX2[curpos] , logY2[curpos] );
  aux2 = LOGP( 0 , logM[curpos] );

  *kernel = LOGP( aux, aux2 )/(INTSCALE); /* Return the logarithm of the kernel */

  /* Memory release */
  free( logM );  free( logX );  free( logY );  free( logX2 );  free( logY2 );
}


/********************/
/* Public functions */
/********************/


void LAkernel(x,y,params,kernel)
     /* Return the log-probability of two sequences x and y under a pair HMM model */
     /* x and y are strings of aminoacid letters, e.g., "AABRS" */
     char *x , *y;
     double *params, *kernel;
{
  int *aax,*aay;  /* to convert x and y into sequences of amino-acid indexes */
  int lx,ly;       /* lengths of x and y */
  int i,j;

  /* If necessary, initialize static variables */
  if (isAA == NULL)
    initialize(params);

  lx = strlen(x);
  ly = strlen(y);

  if ((lx<1) || (ly<1))
    terror("empty chain");

  /* Create aax and aay */

  if ((aax=(int *)calloc(lx,sizeof(int))) == NULL)
    terror("run out of memory");
  if ((aay=(int *)calloc(ly,sizeof(int))) == NULL)
    terror("run out of memory");

  /* Extract the characters corresponding to aminoacids and keep their indexes */

  j=0;
  for (i=0 ; i<lx ; i++) 
    if (isAA[toupper(x[i])])
      aax[j++] = aaIndex[toupper(x[i])-'A'];
  lx = j;
  j=0;
  for (i=0 ; i<ly ; i++)
    if (isAA[toupper(y[i])])
      aay[j++] = aaIndex[toupper(y[i])-'A'];
  ly = j;

  /* Compute the pair HMM score */
  /*  return(LAkernelcompute(aax,aay,lx,ly));*/
  LAkernelcompute(aax,aay,lx,ly,params,kernel);

  /* Release memory */
  free(aax);
  free(aay);
}
