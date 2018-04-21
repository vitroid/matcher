#ifndef MATCHER_MATCHER_H
#define MATCHER_MATCHER_H

#include "common.h"

int match_len(match* s);
match* matcher_core2(int nOatoms, double* Oatoms,
                     double cell[3],
                     int nunitatoms, double* unitatoms,
                     double unitcell[3],
                     double err,
                     double rprox,
                     int adjdens,
                     int nostore);

#endif
  
