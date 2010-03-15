/*********************************************************************
    LAkernel_1_vs_all.c (version 0.3)                                          
n                                                                    
    Compute the local alignment kernel between one protein and all
    proteins in a database. The protein and the database are two
    FASTA format files.

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

FILE *inp, *outp;
double c;
static double params[213];

int main(int argc, char *argv[])
{
  SEQFILE *sfp;
  char *seq1 , *seq2;
  int i, len, data_size = 1;

  if(argc == 4){
    if((inp=fopen(argv[3],"r")) == NULL){
      fprintf(stderr, "cannot open amino acid substitution matrix\n");
      exit(-2);
    }
  }else if(argc == 3){
    if((inp=fopen("parameters.dat","r")) == NULL){
      fprintf(stderr, "cannot open amino acid substitution matrix\n");
      exit(-2);
    }
  }else{
    fprintf(stderr, "usage: ./LAkernel_1_vs_all sequence database (parameter file)\n");
      exit(1);
  }

  /* Read the parameter file */
  i=0;
  while (fscanf(inp,"%lf", &c)!=EOF){
    params[i++]=c;
  }
  fclose(inp);

  /* Read the protein sequence */
  if ((sfp = seqfopen2(argv[1])) == NULL) {
    fprintf(stderr,"Unable to open %s\n",argv[1]);
    exit(1);
  }
  if ((seq1 = seqfgetseq(sfp, &len, 1)) == NULL) {
    fprintf(stderr,"Unable to read the sequence in %s",argv[1]);
    exit(1);
  }
  seqfclose(sfp);

  /* Open the database */  
  if ((sfp = seqfopen2(argv[2])) == NULL) {
    fprintf(stderr,"Unable to open %s\n",argv[2]);
    exit(1);
  }

  /* Read sequences one by one in the database */

  while ((seq2 = seqfgetseq(sfp, &len, 1)) != NULL) {

    /* Compute the the score and print it */
    LAkernel(seq1,seq2,params,&kernel);
    printf("%.2f\n",kernel);

    free(seq2);
  }
  seqfclose(sfp);
  free(seq1);
  
  return 0;
}


