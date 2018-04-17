#ifndef MATCHER_MATCHER_H
#define MATCHER_MATCHER_H

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
  matchtype;

int match_len(matchtype* s);
matchtype* matcher_core2(int nOatoms, double* Oatoms,
			 double cell[3],
			 int nunitatoms, double* unitatoms,
			 double unitcell[3],
			 double err,
			 double rprox,
			 int adjdens,
			 int nostore);

#endif
  
