/*********************************************************************
    LAkernel_all_vs_all.c (version 0.3)                                          
                                                                    
    Compute the local alignment kernel between all pairs of proteins
    in a FASTA format file.

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
#include "math.h"

#define LABEL_LENGTH 8 /* max length of label */

int main(int argc, char *argv[])
{
  SEQFILE *sfp;
  SEQINFO *sip;
  char *seq1 , *seq2;
  int len, processed=0, cont=1, process, i, x, y, data_size;
  double **matrix;
  char **labels, **p, *q;
  FILE *inp, *outp;
  double c;
  double params[213];

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
    fprintf(stderr, "usage: ./LAkernel_direct database (parameter file)\n");
      exit(1);
  }

  /* Read the parameter file */
  i=0;
  while (fscanf(inp,"%lf", &c)!=EOF){
    params[i++]=c;
  }
  fclose(inp);

  /* count the number of sequences */
  if ((sfp = seqfopen2(argv[1])) == NULL) {
    fprintf(stderr,"Unable to open %s\n",argv[1]);
    exit(1);
  }
  data_size=0;
  while ((seq1 = seqfgetseq(sfp, &len, 1)) != NULL){
    free(seq1);
    data_size++;
  }

  if((matrix = (double**)malloc(sizeof(double*)*data_size)) == NULL){
    fprintf(stderr,"Unable to allocate memory for matrix0 !\n");
    exit(1);
  }
  for(i=0;i<data_size;i++){
    matrix[i] = (double*)malloc(sizeof(double)*data_size);
    if(matrix[i] == NULL){
      fprintf(stderr,"Unable to allocate memory for matrix %d!\n",i);
      exit(1);
    }
  }

  if((labels = (char**)malloc(sizeof(char*)*data_size)) == NULL){
    fprintf(stderr,"Unable to allocate memory for labels0 !\n");
    exit(1);
  }
  for(i=0;i<data_size;i++){
    labels[i] = (char*)malloc(sizeof(char)*(LABEL_LENGTH+1)); 
    if(labels[i] == NULL){
      fprintf(stderr,"Unable to allocate memory for labels %d!\n",i);
      exit(1);
    }
  }

  /* memory allocation done, main loop */

  x=0;
  while (cont) {
    /* open the database file and search the first query sequenced not processed yet */
    if ((sfp = seqfopen2(argv[1])) == NULL) {
      fprintf(stderr,"Unable to open %s\n",argv[2]);
      exit(1);
    }

    process=0;
    while (((seq1 = seqfgetseq(sfp, &len, 1)) != NULL) && (process<processed)){
      process++;
      free(seq1);
    }

    if (seq1 == NULL)
      cont = 0;
    else {
      processed++;
      /* seq1 contains the query protein. Retrieve and print its name */
      sip = seqfinfo(sfp,1);
      p=&sip->description;
      q=labels[x];
      for(i=0;i<LABEL_LENGTH;i++){
	*(q+i)=*(*p+i);
      }
      *(q+LABEL_LENGTH)='\0';
      free(sip);
      seqfclose(sfp);

      /* next open the database file */
      if ((sfp = seqfopen2(argv[1])) == NULL) {
	fprintf(stderr,"Unable to open %s\n",argv[1]);
	exit(1);
      }

      y=0;
      /* Read sequences one by one in the database */
      while ((seq2 = seqfgetseq(sfp, &len, 1)) != NULL) {
	/* Compute the pair HMM log-probability and print it */
	if(x<=y){
	  LAkernel(seq1,seq2,params,&kernel);
	  matrix[x][y] = kernel;
	}
	free(seq2);
	y++;
      }
      seqfclose(sfp);
      free(seq1);      
      x++;
    }
  }

  /*print the results*/
  printf("matrix");
  for(x=0;x<data_size;x++){
    printf("\t%s",labels[x]);
  }
  printf("\n");
  for(x=0;x<data_size;x++){
    printf("%s",labels[x]);
    for(y=0;y<data_size;y++){
      if(x>y){
	printf("\t%.2f",matrix[y][x]);
      }else{
      printf("\t%.2f",matrix[x][y]);
      }
    }
    printf("\n");
  }

  for(i=0;i<data_size;i++){
    free(matrix[i]);
  }
  free(matrix);

  for(i=0;i<data_size;i++){
    free(labels[i]);
  }
  free(labels);
  return 0;
}
