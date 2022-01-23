#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <math.h>
#include "neighborlist.h"
#include <assert.h>



sAddressBook*
init_addressbook(int grid[3])
{
  sAddressBook* abook = (sAddressBook*) malloc (sizeof(sAddressBook));
  for(int d=0;d<3;d++){
    abook->grid[d] = grid[d];
  }
  int ngrid = grid[0]*grid[1]*grid[2];
  abook->nresidents = (int*) calloc(ngrid, sizeof(int));
  abook->residents  = (int**) malloc(sizeof(int*)*ngrid);
  return abook;
}

void
dispose_addressbook(sAddressBook* abook)
{
  int L = abook->grid[0]*abook->grid[1]*abook->grid[2];
  for(int i=0;i<L;i++){
    free(abook->residents[i]);
  }
  free(abook->residents);
  free(abook->nresidents);
  free(abook);
}


int
neighborhoods(int district[3], sAddressBook* abook, int** residents)
/*
returns the list of residents in the addressbook near addr.
*/
{
  int nnei = 0;
  int GX = abook->grid[0];
  int GY = abook->grid[1];
  int GZ = abook->grid[2];
  for(int i=-1; i<2; i++){
    int ii = (i + district[0] + GX) % GX;
    for(int j=-1; j<2; j++){
      int jj = (j + district[1] + GY) % GY;
      for(int k=-1; k<2; k++){
        int kk = (k + district[2] + GZ) % GZ;
        nnei += abook->nresidents[ADDRESS(ii,jj,kk)];
      }
    }
  }
  *residents = (int*) malloc(sizeof(int)*nnei);
  int nnei2 = 0;
  for(int i=-1; i<2; i++){
    int ii = (i + district[0] + GX) % GX;
    for(int j=-1; j<2; j++){
      int jj = (j + district[1] + GY) % GY;
      for(int k=-1; k<2; k++){
        int kk = (k + district[2] + GZ) % GZ;
        for(int l=0;l<abook->nresidents[ADDRESS(ii,jj,kk)]; l++){
          (*residents)[nnei2] = abook->residents[ADDRESS(ii,jj,kk)][l];
          nnei2++;
        }
      }
    }
  }
  assert (nnei == nnei2);
  return nnei;
}



int
near_origin(sAddressBook* abook, int** residents)
/*
returns the list of residents in the addressbook near the origin.
*/
{
  int nnei = 0;
  int GX = abook->grid[0];
  int GY = abook->grid[1];
  int GZ = abook->grid[2];
  for(int i=-1; i<1; i++){
    int ii = (i + GX) % GX;
    for(int j=-1; j<1; j++){
      int jj = (j + GY) % GY;
      for(int k=-1; k<1; k++){
        int kk = (k + GZ) % GZ;
        nnei += abook->nresidents[ADDRESS(ii,jj,kk)];
      }
    }
  }
  *residents = (int*) malloc(sizeof(int)*nnei);
  nnei = 0;
  for(int i=-1; i<1; i++){
    int ii = (i + GX) % GX;
    for(int j=-1; j<1; j++){
      int jj = (j + GY) % GY;
      for(int k=-1; k<1; k++){
        int kk = (k + GZ) % GZ;
        for(int l=0;l<abook->nresidents[ADDRESS(ii,jj,kk)]; l++){
          (*residents)[nnei] = abook->residents[ADDRESS(ii,jj,kk)][l];
          nnei++;
        }
      }
    }
  }
  return nnei;
}


sAddressBook* AddressBook(int grid[3], int npos, double* rpos)
/*
  make an address book from a list of atomic positions.
*/
{
  sAddressBook* abook = init_addressbook(grid);
  int GX = grid[0];
  int GY = grid[1];
  int GZ = grid[2];
  for(int i=0; i<npos; i++){
    int addr[3];
    for(int d=0; d<3; d++){
      double x = rpos[i*3+d];
      addr[d] = (int)(floor(x*grid[d]));
      if ( addr[d] < 0 ){
        addr[d] += grid[d];
      }
    }
    abook->nresidents[ADDRESS(addr[0], addr[1], addr[2])] ++;
  }
  int L = GX*GY*GZ;
  for(int i=0;i<L;i++){
    abook->residents[i] = (int*) malloc(sizeof(int)*(abook->nresidents[i]+1));
    abook->residents[i][abook->nresidents[i]] = -1; //terminator
    abook->nresidents[i] = 0;
  }
  for(int i=0; i<npos; i++){
    int addr[3];
    for(int d=0; d<3; d++){
      double x = rpos[i*3+d];
      addr[d] = (int)(floor(x*grid[d]));
      if ( addr[d] < 0 ){
        addr[d] += grid[d];
      }
    }
    int n = abook->nresidents[ADDRESS(addr[0], addr[1], addr[2])];
    abook->nresidents[ADDRESS(addr[0], addr[1], addr[2])] ++;
    abook->residents[ADDRESS(addr[0], addr[1], addr[2])][n] = i;
  }
  return abook;
}



int find_nearest(double rloc[3], sAddressBook* abook, double cell[3], double rpos[])
/*
  find the nearest atom near rloc.
  rloc: location (relative)
  abook: addressbook
  cell: cell dimension
  rpos: positions of atoms (relative)
*/
{
  int GX = abook->grid[0];
  int GY = abook->grid[1];
  int GZ = abook->grid[2];
  int addr[3];
  for(int d=0; d<3; d++){
    double x = rloc[d];
    addr[d] = (int)(floor(x*abook->grid[d]));
  }
  int* neis;
  int nneis = neighborhoods(addr, abook, &neis);
  double dmin = 1e99;
  int    jmin = -1;
  for(int i=0; i<nneis; i++){
    int j = neis[i];
    double dd = 0.0;
    for(int d=0;d<3;d++){
      double delta = rloc[d] - rpos[j*3+d];
      delta -= floor(delta+0.5);
      delta *= cell[d];
      dd += delta*delta;
    }
    if ( dd < dmin ){
      dmin = dd;
      jmin = j;
    }
  }
  free(neis);
  //fprintf(stderr, "jmin %d\n", jmin);
  return jmin;
}
