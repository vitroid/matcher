#ifndef SMATCH_COMMON_H
#define SMATCH_COMMON_H
#include <stdio.h>
#include "bst.h"

int LoadGRO(FILE* file, double** Oatoms, double* cell);
double dot(double* x, double* y);
void sub(double*x, double* y, double* z);
void MakeNeighborList(int natoms, int npairs, int* pairs, 
		      //return values
		      bnode** nei);
void MakeNeighborList_hetero(int natoms, int npairs, int* pairs, 
		      //return values
		      bnode** nei);
#endif
