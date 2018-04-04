#ifndef __PAIRLIST_H__
#define __PAIRLIST_H__
#define PL_FLOAT float
int pairlist(int nAtoms0, PL_FLOAT *atoms0, int nAtoms1, PL_FLOAT *atoms1, PL_FLOAT lower, PL_FLOAT higher, PL_FLOAT cell[3], int **pairs);
int Pairs(int npos, PL_FLOAT* rpos, int ngrid[3], int **pairs);
int Pairs2(int npos0, PL_FLOAT *rpos0, int npos1, PL_FLOAT *rpos1, int ngrid[3], int **pairs);


#endif
