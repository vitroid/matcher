#include <stdlib.h>
#include <string.h>
#include "common.h"


int
LoadGRO(FILE* file, double** Oatoms, double* cell){
  int nAtoms, nOatoms;
  char line[1024];
  fgets( line, sizeof(line), file);
  fgets( line, sizeof(line), file);
  nAtoms = atoi(line);
  nOatoms = 0;
  *Oatoms = (double*) malloc( sizeof(double) * nAtoms*3 );
  for(int i=0;i<nAtoms;i++){
    fgets(line, sizeof(line), file);
    double x,y,z;
    int   id;
    char  label[99],label2[99];
    sscanf(line, "%10s%5s%5d%8lf%8lf%8lf\n", label, label2, &id, &x, &y, &z);
    if ( strncasecmp( label2, "O", 1 ) == 0 ){
      (*Oatoms)[nOatoms*3+0] = x;
      (*Oatoms)[nOatoms*3+1] = y;
      (*Oatoms)[nOatoms*3+2] = z;
      nOatoms ++;
    }
  }
  //assert(nOatoms==nAtoms/3);
  fgets(line, sizeof(line), file);
  sscanf(line, "%lf %lf %lf\n", &cell[0], &cell[1], &cell[2]);
  return nOatoms;
}


double dot(double* x, double* y)
{
  double sum =0;
  for(int d=0;d<3;d++){
    sum += x[d]*y[d];
  }
  return sum;
}

  


void
sub(double*x, double* y, double* z)
{
  for(int d=0;d<3;d++){
    z[d] = x[d] -y[d];
  }
}




void
MakeNeighborList(int natoms, int npairs, int* pairs, 
		 //return values
		 bnode** nei)
{
  //make neighborlist
  for(int i=0;i<natoms; i++){
    nei[i] = NULL;
  }
  for(int i=0;i<npairs;i++){
    int r0 = pairs[i*2+0];
    int r1 = pairs[i*2+1];
    nei[r0] = insert(nei[r0], r1);
    nei[r1] = insert(nei[r1], r0);
  }
}
