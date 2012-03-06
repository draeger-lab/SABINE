/*********************************************************************
    LAkernel_direct.c (version 0.3)                                          
                                                                    
    Compute the local alignment kernel between two sequences, input
    in the command line.

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

int main(int argc, char *argv[])
{
  FILE *inp;
  int i;
  double c;
  double params[213];

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
    fprintf(stderr, "usage: ./LAkernel_direct sequence1 sequence2 (parameter file)\n");
      exit(1);
  }
 
  i=0;
  while (fscanf(inp,"%lf", &c)!=EOF){
    params[i++]=c;
  }
  fclose(inp);

  LAkernel(argv[1],argv[2],params,&kernel);

  printf("%.5lf\n",kernel);
  return 0;
}


