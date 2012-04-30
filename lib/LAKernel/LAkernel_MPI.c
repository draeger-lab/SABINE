/*********************************************************************
    LAkernel_MPI.c (version 0.3)                                          
                                                                    
    Compute the local alignment kernel (header file)

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

#include <stdio.h>
#include <stdlib.h>
#include "LAkernel.h"
#include "seqio.h"
#include "mpi.h"
#include "math.h"
#define LABEL_LENGTH 8 /* max length of label */
#define MASTER 0               /* taskid of first task */
#define FROM_MASTER 1          /* setting a message type */
#define FROM_WORKER 2          /* setting a message type */

int main(int argc,char *argv[]){
  SEQFILE *sfp;
  SEQINFO *sip;
  char *seq;
  int numtasks,              /* number of tasks in partition */
    taskid,                /* a task identifier */
    numworkers,            /* number of worker tasks */
    source,                /* task id of message source */
    dest,                  /* task id of message destination */
    mtype,                 /* message type */
    jobs,                  /* number of LA calculations sent to each worker */
    avejob, extra, offset, /* used to determine jobs sent to each worker */
    i, j, k, rc,           /* misc */
    num,
    data_size,             /* number of sequences in a file */
    len, maxlen=0;         /* length of sequence */ 
  char** M;                /* matrix to store all the protein sequence (memory consuming) */
  char **labels, **p, *q;
  int *X,*Y;
  double *scores;
  double **matrix;
  double startwtime,endwtime;
  MPI_Status status;
  double c, params[213];
  FILE *inp, *outp;

  rc = MPI_Init(&argc,&argv);
  rc|= MPI_Comm_size(MPI_COMM_WORLD,&numtasks);
  rc|= MPI_Comm_rank(MPI_COMM_WORLD,&taskid);
  if (rc != 0)
    printf ("error initializing MPI and obtaining task ID information\n");
   else
     /*     printf ("task ID = %d\n", taskid);*/
  numworkers = numtasks-1;

  if(argc == 3){
    if((inp=fopen(argv[2],"r")) == NULL){
      fprintf(stderr, "cannot open amino acid substitution matrix\n");
      exit(-2);
    }
  }else if(argc == 2){
    if((inp=fopen("parameters.dat","r")) == NULL){
      fprintf(stderr, "cannot open amino acid substitution matrix\n");
      exit(-2);
    }
  }else{
    fprintf(stderr, "usage: mpirun -np #cpu ./LAkernel_MPI database (parameter file)\n");
    exit(1);
  }

  /* Read the parameter file */
  i=0;
  while (fscanf(inp,"%lf", &c)!=EOF){
    params[i++]=c;
  }
  fclose(inp);

  /* examine the data size, then allocate memory */

  if ((sfp = seqfopen2(argv[1])) == NULL) {
    fprintf(stderr,"Unable to open database file\n");
    exit(1);
  }
  data_size=0;
  while ((seq = seqfgetseq(sfp, &len, 1)) != NULL){
    free(seq);
    if(len > maxlen){
      maxlen = len;
    }
    data_size++;
  }
  seqfclose(sfp);

  num = data_size*(data_size+1)/2;

  if((scores = (double*)malloc(sizeof(double)*num)) == NULL){
    fprintf(stderr,"Unable to allocate memory for matrix0 !\n");
    exit(1);
  }

  /* prepare index X and Y for M */

  if((X = (int *)malloc(sizeof(int)*num)) == NULL){
    fprintf(stderr,"Unable to allocate memory for X !\n");
    exit(1);
  }
  if((Y = (int *)malloc(sizeof(int)*num)) == NULL){
    fprintf(stderr,"Unable to allocate memory for Y !\n");
    exit(1);
  }

  k=0;
  for (i=0; i<data_size; i++){
    for (j=i; j<data_size; j++){
      X[k]=i;
      Y[k]=j;
      k++;
    }
  }


  /**************************** master task ************************************/
  if (taskid == MASTER)
    {
      
      if((labels = (char**)malloc(sizeof(char*)*data_size)) == NULL){
	fprintf(stderr,"Unable to allocate memory for labels0 !\n");
	exit(1);
      }
      for(i=0;i<data_size;i++){
	labels[i] = (char*)malloc(sizeof(char)*LABEL_LENGTH); 
	if(labels[i] == NULL){
	  fprintf(stderr,"Unable to allocate memory for labels %d!\n",i);
	  exit(1);
	}
      }
      
      if ((sfp = seqfopen2(argv[1])) == NULL) {
	fprintf(stderr,"Unable to open a file\n");
	exit(1);
      }
      
      k=0;
      while((seq = seqfgetseq(sfp, &len, 1)) != NULL){
	sip = seqfinfo(sfp,1);
	p=&sip->description;
	q=labels[k];
	for(i=0;i<LABEL_LENGTH;i++){
	  *(q+i)=*(*p+i);
	}
	*(q+LABEL_LENGTH)='\0';
	free(sip);
	free(seq);
	k++;
      }
      seqfclose(sfp);
      
      if((matrix = (double**)malloc(sizeof(double*)*data_size)) == NULL){
	fprintf(stderr,"Unable to allocate memory for matrix0 !\n");
	exit(1);
      }
      for(i=0;i<data_size;i++){
	matrix[i] = (double*)malloc(sizeof(double)*(data_size));
	if(matrix[i] == NULL){
	  fprintf(stderr,"Unable to allocate memory for matrix %d!\n",i);
	  exit(1);
	}
      }
      
      /*      printf("Number of worker tasks = %d\n",numworkers);*/
      startwtime=MPI_Wtime();
      /* send matrix data to the worker tasks */

      numworkers=numtasks-1;
      avejob = num/numworkers;
      extra = num%numworkers;
      offset = 0;
      mtype = FROM_MASTER;
      /*      printf("num=%d num_worker=%d avejob=%d extra=%d\n",num,numworkers,avejob,extra);*/
      for (dest=1; dest<=numworkers; dest++)
	{
	  jobs = (dest <= extra) ? avejob+1 : avejob;   
	  /*	  printf("   sending %d jobs to task %d\n",jobs,dest);*/
	  MPI_Send(&offset, 1, MPI_INT, dest, mtype, MPI_COMM_WORLD);
	  MPI_Send(&jobs, 1, MPI_INT, dest, mtype, MPI_COMM_WORLD);
	  offset = offset + jobs;
	}

      /* wait for results from all worker tasks */
      mtype = FROM_WORKER;
      k=0; /* index for scores */
      for (i=1; i<=numworkers; i++)
	{
	  source = i;
	  MPI_Recv(&offset, 1, MPI_INT, source, mtype, MPI_COMM_WORLD, &status);
	  MPI_Recv(&jobs, 1, MPI_INT, source, mtype, MPI_COMM_WORLD, &status);
	  MPI_Recv(&scores[offset], jobs, MPI_DOUBLE, source, mtype, MPI_COMM_WORLD, &status);
	  for(j=0;j<jobs;j++){
	    /*	 scores[k]=matrix[X[k]Y[k]]   */
	    if(X[k]<Y[k]){
              matrix[X[k]][Y[k]]=matrix[Y[k]][X[k]]=scores[k];
	    }else{
	      matrix[X[k]][Y[k]]=scores[k];
	    }
	    k++;
	  }
	}

      /*      printf("total %d cells\n",k);
      endwtime=MPI_Wtime();
      printf("mxm time = %f  %f Mflops\n",endwtime-startwtime, num/(endwtime-startwtime)*1.e-6); */

      /* print results */
      /*      printf("Here is the result matrix\n");*/
      printf("matrix");
      for(i=0;i<data_size;i++){
	printf("\t%s",labels[i]);
      }
      printf("\n");
      for(i=0;i<data_size;i++){
	printf("%s",labels[i]);
	for(j=0;j<data_size;j++){
	  printf("\t%.2f",matrix[i][j]);
	}
	printf("\n");
      }

      /*      for(i=0;i<num;i++)printf("\t%.2f",scores[i]);*/

      /* release memory */
      free(X);
      free(Y);
      free(scores);
      for(i=0;i<data_size;i++){
	free(matrix[i]);
      }
      free(matrix);
      for(i=0;i<data_size;i++){
	free(labels[i]);
      }
      free(labels);
    }

  /**************************** worker task ************************************/
  if (taskid > MASTER)
    {

      if((M = (char**)malloc(sizeof(char*)*data_size)) == NULL){
	fprintf(stderr,"Unable to allocate memory for M0 !\n");
	exit(1);
      }
      for(i=0;i<data_size;i++){
	M[i] = (char*)malloc(sizeof(char)*(maxlen+1)); 
	if(M[i] == NULL){
	  fprintf(stderr,"Unable to allocate memory for M %d!\n",i);
	  exit(1);
	}
      }
      
      if ((sfp = seqfopen2(argv[1])) == NULL) {
	fprintf(stderr,"Unable to open a file\n");
	exit(1);
      }
      
      k=0;
      while((seq = seqfgetseq(sfp, &len, 1)) != NULL){
	strcpy(M[k],seq);
	free(seq);
	k++;
      }
      seqfclose(sfp);

      mtype = FROM_MASTER;
      MPI_Recv(&offset, 1, MPI_INT, MASTER, mtype, MPI_COMM_WORLD, &status);
      MPI_Recv(&jobs, 1, MPI_INT, MASTER, mtype, MPI_COMM_WORLD, &status);
      /*      printf("offset=%d jobs=%d",offset,jobs);*/

      for(i=offset;i<offset+jobs;i++){
	LAkernel(M[X[i]],M[Y[i]],params,&kernel);
	scores[i] = kernel;
      }

      mtype = FROM_WORKER;
      MPI_Send(&offset, 1, MPI_INT, MASTER, mtype, MPI_COMM_WORLD);
      MPI_Send(&jobs, 1, MPI_INT, MASTER, mtype, MPI_COMM_WORLD);
      MPI_Send(&scores[offset], jobs, MPI_DOUBLE, MASTER, mtype, MPI_COMM_WORLD);
      free(X);
      free(Y);
      free(scores);
      for(i=0;i<data_size;i++){
	free(M[i]);
      }
      free(M);
    }
  MPI_Finalize();
}



