#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <math.h>
#include "common.h"



double
det(double* A, double* B, double* C)
{
  // Ax Bx Cx   a b c
  // Ay By Cy = d e f
  // Az Bz Cz   g h k
  double a = A[0];
  double b = B[0];
  double c = C[0];
  double d = A[1];
  double e = B[1];
  double f = C[1];
  double g = A[2];
  double h = B[2];
  double k = C[2];
  return a * (e * k - f * h) - b * (d * k - f * g) + c * (d * h - e * g);
}





void
inv(double* A, double* B, double* C)
{
  // Ax Bx Cx   a b c
  // Ay By Cy = d e f
  // Az Bz Cz   g h k
  double D = det(A,B,C);
  double a = A[0];
  double b = B[0];
  double c = C[0];
  double d = A[1];
  double e = B[1];
  double f = C[1];
  double g = A[2];
  double h = B[2];
  double k = C[2];
  A[0] =  (e * k - f * h) / D;
  B[0] = -(b * k - c * h) / D;
  C[0] =  (b * f - c * e) / D;
  A[1] = -(d * k - f * g) / D;
  B[1] =  (a * k - c * g) / D;
  C[1] = -(a * f - c * d) / D;
  A[2] =  (d * h - e * g) / D;
  B[2] = -(a * h - b * g) / D;
  C[2] =  (a * e - b * d) / D;
}

int match_len(match* s)
{
  int n=0;
  while ( s != NULL ){
    n += 1;
    s = s->next;
  }
  return n;
}

int
LoadGRO0(FILE* file, double** Oatoms, double* cell){
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
    sscanf(line, "%10s%5s%5d", label, label2, &id);
    sscanf(&line[20], "%lf %lf %lf\n",  &x, &y, &z);
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

  
double vector_length(double x[3])
{
  return sqrt(dot(x,x));
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


//for hetero pairs
void
MakeNeighborList_hetero(int natoms, int npairs, int* pairs, 
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
    //nei[r1] = insert(nei[r1], r0);
  }
}





int
LoadGRO(FILE* file, double** Oatoms, double* cell, int rel){
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
    sscanf(line, "%10s%5s%5d", label, label2, &id);
    sscanf(&line[20], "%lf %lf %lf\n",  &x, &y, &z);
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
  //absolute to relative
  if (rel){
    for(int i=0; i<nOatoms; i++){
      for(int d=0; d<3; d++){
        double x = (*Oatoms)[i*3+d]/cell[d];
        (*Oatoms)[i*3+d] = x - floor(x+0.5);
      }
    }
  }
  return nOatoms;
}
