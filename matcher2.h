#ifndef MATCHER_MATCHER2_H
#define MATCHER_MATCHER2_H

#include "common.h"

int match_len(match* s);
match* matcher2_core(int ngatoms,
                     double* gatoms,
                     double gcell[3],
                     int nuatoms,
                     double* uatoms,
                     double ucell[3],
                     int adjdens,
                     int nostore);

#endif
  
