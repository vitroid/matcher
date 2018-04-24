#ifndef SMATCH_SMATCHER_H
#define SMATCH_SMATCHER_H

typedef struct _smatch {
  int i,j;
  double radius;
  double d[3];
  double rmsd;
  struct _smatch *next;
}
  smatchtype;

smatchtype* smatcher_core(int nOatoms, double* Oatoms, double* cell, double radius, double rmsdmax, int every);
int smatch_len(smatchtype* s);

#endif
