#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#include "pairlist.h"
#include "bst.h"

//read Gromacs-type data and find the OH covalent bonds.

//benchmark on vitroid-black 2017-3-31
//original matcher (no bst):             user	1m28.847s
//matcher with bst:                      user	1m48.723s
//matcher with bst (size() is improved): user	1m48.942s

#define max(A,B) ((A)>(B)?(A):(B))
int
isClose(float x, float y)
{
  float e = (x-y)/max(x,y);
  return (-1e-3 < e) && (e < +1e-3);
}




void
yaplot_linerel(float* a, float* d)
{
  printf("l %f %f %f %f %f %f\n", a[0],a[1],a[2], a[0]+d[0],a[1]+d[1], a[2]+d[2]);
}

int
LoadGRO(FILE* file, float** Oatoms, float* cell){
  int nAtoms, nOatoms;
  char line[1024];
  fgets( line, sizeof(line), file);
  fgets( line, sizeof(line), file);
  nAtoms = atoi(line);
  nOatoms = 0;
  *Oatoms = (float*) malloc( sizeof(float) * nAtoms*3 );
  for(int i=0;i<nAtoms;i++){
    fgets(line, sizeof(line), file);
    float x,y,z;
    int   id;
    char  label[99],label2[99];
    sscanf(line, "%10s%5s%5d%8f%8f%8f\n", label, label2, &id, &x, &y, &z);
    if ( strncasecmp( label2, "Ow", 2 ) == 0 ){
      (*Oatoms)[nOatoms*3+0] = x;
      (*Oatoms)[nOatoms*3+1] = y;
      (*Oatoms)[nOatoms*3+2] = z;
      nOatoms ++;
    }
  }
  //assert(nOatoms==nAtoms/3);
  fgets(line, sizeof(line), file);
  sscanf(line, "%f %f %f\n", &cell[0], &cell[1], &cell[2]);
  return nOatoms;
}


float dot(float* x, float* y)
{
  float sum =0;
  for(int d=0;d<3;d++){
    sum += x[d]*y[d];
  }
  return sum;
}

  
float vector_length(float* v)
{
  return sqrt(dot(v,v));
}

void
cross(float*x, float* y, float* z)
{
  z[0] = x[1]*y[2] - x[2]*y[1];
  z[1] = x[2]*y[0] - x[0]*y[2];
  z[2] = x[0]*y[1] - x[1]*y[0];
}


void
sub(float*x, float* y, float* z)
{
  for(int d=0;d<3;d++){
    z[d] = x[d] -y[d];
  }
}


void
mul(float*x, float a, float* z)
{
  for(int d=0;d<3;d++){
    z[d] = x[d] * a;
  }
}


void
regularize(float* x, float* y, float* z)
{
  float z0[3];
  for(int d=0;d<3;d++)
    z0[d] = z[d];
  float error = dot(x,y);
  float ye[3];
  mul(y,error/2,ye);
  float x_ort[3];
  sub(x,ye,x_ort);
  float xe[3];
  mul(x,error/2,xe);
  float y_ort[3];
  sub(y,xe,y_ort);
  float z_ort[3];
  cross(x_ort, y_ort, z_ort);
  mul(x_ort, 0.5*(3.0-dot(x_ort,x_ort)), x);
  mul(y_ort, 0.5*(3.0-dot(y_ort,y_ort)), y);
  mul(z_ort, 0.5*(3.0-dot(z_ort,z_ort)), z);
  //allow mirroring
  if (dot(z0, z) < 0){
    mul(z, -1, z);
  }
}


float
det(float* A, float* B, float* C)
{
  // Ax Bx Cx   a b c
  // Ay By Cy = d e f
  // Az Bz Cz   g h k
  float a = A[0];
  float b = B[0];
  float c = C[0];
  float d = A[1];
  float e = B[1];
  float f = C[1];
  float g = A[2];
  float h = B[2];
  float k = C[2];
  return a * (e * k - f * h) - b * (d * k - f * g) + c * (d * h - e * g);
}





void
inv(float* A, float* B, float* C)
{
  // Ax Bx Cx   a b c
  // Ay By Cy = d e f
  // Az Bz Cz   g h k
  float D = det(A,B,C);
  float a = A[0];
  float b = B[0];
  float c = C[0];
  float d = A[1];
  float e = B[1];
  float f = C[1];
  float g = A[2];
  float h = B[2];
  float k = C[2];
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

    





int
LoadAR3A(FILE* file, float** atoms){
  int nAtoms;
  char line[1024];
  while ( NULL != fgets( line, sizeof(line), file) ){
    if ( 0 == strncmp(line, "@AR3A", 5) ){
      fgets( line, sizeof(line), file);
      nAtoms = atoi(line);
      *atoms = (float*) malloc( sizeof(float) * nAtoms*3 );
      float dmin = 1e99;
      int   imin = -1;
      for(int i=0;i<nAtoms;i++){
	fgets(line, sizeof(line), file);
	float x,y,z;
	sscanf(line, "%f %f %f", &x, &y, &z);
	//atoms distribute around the origin.
	float rr = x*x + y*y + z*z;
	if ( rr < dmin ){
	  dmin = rr;
	  imin = i;
	}
	(*atoms)[i*3+0] = x / 10;
	(*atoms)[i*3+1] = y / 10;
	(*atoms)[i*3+2] = z / 10; //nm
      }
      //there must be a molecule at the origin
      for(int i=0;i<nAtoms;i++){
	(*atoms)[i*3+0] -= (*atoms)[imin*3+0];
	(*atoms)[i*3+1] -= (*atoms)[imin*3+1];
	(*atoms)[i*3+2] -= (*atoms)[imin*3+2];
      }
    }
  }
  return nAtoms;
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





int main(int argc, char* argv[])
{
  //usage: matcher grofile error msdmax
  //typical values: error = 0.03, msdmax == 0.35 for TPPI
  if ( argc != 5 ){
    fprintf(stderr, "usage: %s grofile error msdmax\n", argv[0]);
    exit(1);
  }
  float cell[3];
  float *Oatoms;
  FILE *file = fopen(argv[1], "r");
  float err, msdmax;
  sscanf(argv[3], "%f", &err);
  sscanf(argv[4], "%f", &msdmax);
  int nOatoms = LoadGRO(file, &Oatoms, cell);
  fclose(file);

  float *unitatoms; //relative
  file = fopen(argv[2], "r");
  int nunitatoms = LoadAR3A(file, &unitatoms);
  fclose(file);
  float rmax = 0;
  for(int i=0;i<nunitatoms;i++){
    float r = vector_length(&unitatoms[i*3]);
    if ( rmax < r ){
      rmax = r;
    }
  }
  float radius = rmax;
  fprintf(stderr,"System   %d %f %f %f\n", nOatoms, cell[0], cell[1], cell[2]);
  fprintf(stderr,"Template %d\n", nunitatoms);
  //Find the tetrahedra of a, b, and c cell vectors
  int* prox;
  //atoms of the proximity
  fprintf(stderr, "Preparing the neighbor lists.");
  int nProx   = pairlist(nOatoms, Oatoms, nOatoms, Oatoms, 0.0, radius*(1+err), cell, &prox);
  bnode* nei[nOatoms];
  MakeNeighborList(nOatoms, nProx, prox, nei);
  fprintf(stderr, "Done.\n");
  free(prox);
  fprintf(stderr,"%d nProx A\n", nProx);
  //find triangle PQR that matches the shape
  int ntet=0;
  for(int p=0; p<nOatoms; p++){
    fprintf(stderr, "\r%.1f %%", p*100./nOatoms);
    //this does not include p itself!
    int* neighbors = get_array(nei[p]);
    neighbors[size(nei[p])] = p; // add itself in place of the useless terminator
    float msd = 0.0;
    int partners[nunitatoms];
    for(int l=0; l<nunitatoms; l++){
      float u[3];
      u[0] = unitatoms[l*3+0] + Oatoms[p*3+0];
      u[1] = unitatoms[l*3+1] + Oatoms[p*3+1];
      u[2] = unitatoms[l*3+2] + Oatoms[p*3+2];
      float dmin = 1e99;
      for(int m=0; m<size(nei[p])+1; m++){
	float dd[3];
	int ne = neighbors[m];
	sub(u, &Oatoms[ne*3], dd);
	float L = dot(dd,dd);
	if ( L < dmin ){
	  dmin = L;
	  partners[l] = ne;
	}
      }
      msd += dmin;
    }
    if ( msd < msdmax * nunitatoms / 100 ){  //msd in AA
      printf("%d %f %f %f %f\n", p+1, Oatoms[p*3+0]*10, Oatoms[p*3+1]*10, Oatoms[p*3+2]*10, msd);
    }
    free(neighbors);
  }
  //dispose memory as soon as possible
  for(int i=0;i<nOatoms; i++){
    dispose(nei[i]);
  }
  free(Oatoms);
}




