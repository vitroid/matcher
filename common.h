#ifndef SMATCH_COMMON_H
#define SMATCH_COMMON_H
#include <stdio.h>
#include "bst.h"


typedef struct _matchtype
{
  double rmsd;
  int atom_gro;
  int atom_unitcell;
  double mat[9];
  int natom;
  int* atoms;
  struct _matchtype* next;
}
  match;

int LoadGRO(FILE* file, double** Oatoms, double* cell);
double dot(double* x, double* y);
void sub(double*x, double* y, double* z);
double vector_length(double x[3]);
void MakeNeighborList(int natoms, int npairs, int* pairs, 
		      //return values
		      bnode** nei);
void MakeNeighborList_hetero(int natoms, int npairs, int* pairs, 
		      //return values
		      bnode** nei);
#endif
