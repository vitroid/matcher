#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#include "pairlist.h"

#define ADDRESS(x,y,z) (((z)*GY + (y))*GX + (x))
#define True 1
#define False 0
int pairlist1(int nAtoms, float *atoms, float lower, float higher, float cell[3], int **pairs);
int pairlist2(int nAtoms0, float *atoms0, int nAtoms1, float *atoms1, float lower, float higher, float cell[3], int **pairs);

int
pairlist(int nAtoms0, float *atoms0, int nAtoms1, float *atoms1, float lower, float higher, float cell[3], int **pairs)
{
  if ( ( nAtoms1 == 0 ) || (atoms1 == NULL) || ( atoms0 == atoms1) ){
    if ( atoms0 == atoms1)
      assert (nAtoms0 == nAtoms1);
    return pairlist1(nAtoms0, atoms0, lower, higher, cell, pairs);
  }
  return pairlist2(nAtoms0, atoms0, nAtoms1, atoms1, lower, higher, cell, pairs);
}
   


int
gridpairlist(int ngrid[3], int single, int **pairs)
/* 
given:
    ngrid: number of grids in three dimensions
    single: if 1, make pair list for a single component system
returns:
    return value: number of grid pairs
    pairs: newly allocated list of grid pairs
*/
{
  const int GX = ngrid[0];
  const int GY = ngrid[1];
  const int GZ = ngrid[2];
  const int nTotalGrids = GX*GY*GZ;
  //assert ( GX > 2 );
  //assert ( GY > 2 );
  //assert ( GZ > 2 );
  int npair = 0;
  int neighborcells = 26;
  if ( single ) neighborcells = 13;
  (*pairs) = (int*) malloc(sizeof(int)*neighborcells*nTotalGrids*2);
  //make neighboring grid list
  for(int ix=0;ix<GX;ix++){
    for(int iy=0;iy<GY;iy++){
      for(int iz=0;iz<GZ;iz++){
        int a = ADDRESS(ix,iy,iz); //center
	//for all 26 neighbor cells
        int kx0 = (GX==2) ? ix : ix - 1;
        for(int kx=kx0;kx<=ix+1;kx++){
          int jx = (kx+GX) % GX;
          int ky0 = (GY==2) ? iy : iy - 1;
          for(int ky=ky0;ky<=iy+1;ky++){
            int jy = (ky+GY) % GY;
            int kz0 = (GZ==2) ? iz : iz - 1;
            for(int kz=kz0;kz<=iz+1;kz++){
              int jz = (kz+GZ) % GZ;
              int b = ADDRESS(jx,jy,jz);
              if ( (a != b) && ( ( ! single ) || ( a<b ) ) ){
                (*pairs)[npair*2+0] = a;
                (*pairs)[npair*2+1] = b;
                npair ++;
              }
            }
          }
        }
      }
    }
  }
  return npair;
}


int
pairlist1(int nAtoms, float *atoms, float lower, float higher, float cell[3], int **pairs)
/* 
given:
    nAtoms: number of atoms
    atoms: positions of atoms
    lower, higher: two atoms are considered to be paired when their distance is between them.
    cell: rectangular cell shape
returns:
    return value: number of pairs
    pairs: newly allocated list of pairs
*/
{
  int test = 0;
  //determine the grid size
  int ngrid[3];
  for(int d=0;d<3;d++){
    ngrid[d] = (int) floor(cell[d] / higher);
  }
  const int GX = ngrid[0];
  const int GY = ngrid[1];
  const int GZ = ngrid[2];
  const int nTotalGrids = GX*GY*GZ;
  
  //count the number of residents in each grid cells.
  int nResidents[nTotalGrids];
  for(int i=0;i<nTotalGrids;i++) nResidents[i] = 0;
  for(int i=0;i<nAtoms;i++){
    float r[3];  //relative location in the cell
    for(int d=0;d<3;d++){
      r[d] = atoms[i*3+d] / cell[d];
      r[d] -= floor(r[d]); //Periodic boundary condition
    }
    int grid[3];
    for(int d=0;d<3;d++){
      grid[d] = floor(r[d] * ngrid[d]);
    }
    nResidents[ADDRESS(grid[0],grid[1],grid[2])] ++;
  }

  //make the resident list
  //resident list is serialized to reduce the memory usage
  //head pointer to the list of residents in a grid cell is in headOfList.
  //atom ids are put in the residents[] and is terminated with -1.
  int residents[nAtoms + nTotalGrids];
  int headOfList[nTotalGrids];
  int pointer[nTotalGrids];
  int head = 0;
  for(int g=0;g<nTotalGrids;g++){
    headOfList[g] = head;
    pointer[g] = head;
    head += nResidents[g]+1;
    residents[head-1] = -1;
  }
  for(int i=0;i<nAtoms;i++){
    float r[3];  //relative location in the cell
    for(int d=0;d<3;d++){
      r[d] = atoms[i*3+d] / cell[d];
      r[d] -= floor(r[d]);
    }
    int grid[3];
    for(int d=0;d<3;d++){
      grid[d] = floor(r[d] * ngrid[d]);
    }
    int a = ADDRESS(grid[0],grid[1],grid[2]);
    residents[pointer[a]] = i;
    pointer[a] ++;
  }
  if(test){
    //just print
    for(int g=0;g<nTotalGrids;g++){
      printf("Grid %d\n", g);
      printf("Number of residents: %d\n", nResidents[g]);
      printf("First resident: %d\n", residents[headOfList[g]]);
      printf("Last resident: %d\n", residents[headOfList[g] + nResidents[g] - 1]);
      printf("Terminator (must be -1): %d\n", residents[headOfList[g] + nResidents[g]]);
    }
  }
  //make the neighboring grid pair list
  int* gridPairs;
  int nGridPairs = gridpairlist(ngrid, True, &gridPairs);
  //rough estimate of the number of atom pairs
  int rough = 0;
  for(int i=0;i<nGridPairs;i++){
    int g0 = gridPairs[i*2+0];
    int g1 = gridPairs[i*2+1];
    rough += nResidents[g0] * nResidents[g1];
  }
  for(int g=0;g<nTotalGrids;g++){
    rough += nResidents[g]*nResidents[g] / 2;
  }

  *pairs = (int*) malloc(sizeof(int) * rough*2);
  int nPairs=0;
  for(int i=0;i<nGridPairs;i++){
    int g0 = gridPairs[i*2+0];
    for(int j=0; j<nResidents[g0]; j++){
      int r0 = residents[headOfList[g0] + j];
      int g1 = gridPairs[i*2+1];
      for(int k=0; k<nResidents[g1]; k++){
        int r1 = residents[headOfList[g1] + k];
        float sum2 = 0.0;
        for(int d=0;d<3;d++){
          float delta = atoms[r0*3+d] - atoms[r1*3+d];
          delta -= floor( delta / cell[d] + 0.5 ) * cell[d];
          sum2 += delta*delta;
        }
        if ((lower*lower < sum2) && (sum2 < higher*higher)){
          (*pairs)[nPairs*2+0] = r0;
          (*pairs)[nPairs*2+1] = r1;
          nPairs ++;
        }
      }
    }
  }
  for(int g=0;g<nTotalGrids;g++){
    for(int j=0; j<nResidents[g]; j++){
      int r0 = residents[headOfList[g] + j];
      for(int k=j+1; k<nResidents[g]; k++){
        int r1 = residents[headOfList[g] + k];
        float sum2 = 0.0;
        for(int d=0;d<3;d++){
          float delta = atoms[r0*3+d] - atoms[r1*3+d];
          delta -= floor( delta / cell[d] + 0.5 ) * cell[d];
          sum2 += delta*delta;
        }
        if ((lower*lower < sum2) && (sum2 < higher*higher)){
          (*pairs)[nPairs*2+0] = r0;
          (*pairs)[nPairs*2+1] = r1;
          nPairs ++;
        }
      }
    }
  }
  //deallocate allocated memories
  free(gridPairs);
  return nPairs;
}
  

int
pairlist2(int nAtoms0, float *atoms0, int nAtoms1, float *atoms1, float lower, float higher, float cell[3], int **pairs)
/* 
given:
    nAtoms0: number of atoms of component 0
    atoms0: positions of atoms
    nAtoms1: number of atoms of component 1
    atoms1: positions of atoms
    higher: two atoms are considered to be paired when their distance is below it.
    cell: rectangular cell shape
returns:
    return value: number of pairs
    pairs: newly allocated list of pairs
*/
{
  int test = 0;
  //determine the grid size
  int ngrid[3];
  for(int d=0;d<3;d++){
    ngrid[d] = (int) floor(cell[d] / higher);
  }
  const int GX = ngrid[0];
  const int GY = ngrid[1];
  const int GZ = ngrid[2];
  const int nTotalGrids = GX*GY*GZ;
  printf("%dx%dx%d =%d nTotalGrids\n", GX,GY,GZ,nTotalGrids);
  //count the number of residents in each grid cells.
  int nResidents0[nTotalGrids];
  int nResidents1[nTotalGrids];
  for(int i=0;i<nTotalGrids;i++){
    nResidents0[i] = nResidents1[i] = 0;
  }
  for(int i=0;i<nAtoms0;i++){
    float r[3];  //relative location in the cell
    for(int d=0;d<3;d++){
      r[d] = atoms0[i*3+d] / cell[d];
      r[d] -= floor(r[d]); //Periodic boundary condition
    }
    int grid[3];
    for(int d=0;d<3;d++){
      grid[d] = floor(r[d] * ngrid[d]);
    }
    nResidents0[ADDRESS(grid[0],grid[1],grid[2])] ++;
  }
  for(int i=0;i<nAtoms1;i++){
    float r[3];  //relative location in the cell
    for(int d=0;d<3;d++){
      r[d] = atoms1[i*3+d] / cell[d];
      r[d] -= floor(r[d]); //Periodic boundary condition
    }
    int grid[3];
    for(int d=0;d<3;d++){
      grid[d] = floor(r[d] * ngrid[d]);
    }
    nResidents1[ADDRESS(grid[0],grid[1],grid[2])] ++;
  }

  //make the resident list
  //resident list is serialized to reduce the memory usage
  //head pointer to the list of residents in a grid cell is in headOfList.
  //atom ids are put in the residents[] and is terminated with -1.
  int pointer[nTotalGrids];
  int head;

  int residents0[nAtoms0 + nTotalGrids];
  int headOfList0[nTotalGrids];
  head = 0;
  for(int g=0;g<nTotalGrids;g++){
    headOfList0[g] = head;
    pointer[g] = head;
    head += nResidents0[g]+1;
    residents0[head-1] = -1;
  }
  for(int i=0;i<nAtoms0;i++){
    float r[3];  //relative location in the cell
    for(int d=0;d<3;d++){
      r[d] = atoms0[i*3+d] / cell[d];
      r[d] -= floor(r[d]);
    }
    int grid[3];
    for(int d=0;d<3;d++){
      grid[d] = floor(r[d] * ngrid[d]);
    }
    int a = ADDRESS(grid[0],grid[1],grid[2]);
    residents0[pointer[a]] = i;
    pointer[a] ++;
  }
  int residents1[nAtoms1 + nTotalGrids];
  int headOfList1[nTotalGrids];
  head = 0;
  for(int g=0;g<nTotalGrids;g++){
    headOfList1[g] = head;
    pointer[g] = head;
    head += nResidents1[g]+1;
    residents1[head-1] = -1;
  }
  for(int i=0;i<nAtoms1;i++){
    float r[3];  //relative location in the cell
    for(int d=0;d<3;d++){
      r[d] = atoms1[i*3+d] / cell[d];
      r[d] -= floor(r[d]);
    }
    int grid[3];
    for(int d=0;d<3;d++){
      grid[d] = floor(r[d] * ngrid[d]);
    }
    int a = ADDRESS(grid[0],grid[1],grid[2]);
    residents1[pointer[a]] = i;
    pointer[a] ++;
  }
  //make the neighboring grid pair list
  int* gridPairs;
  int nGridPairs = gridpairlist(ngrid, False, &gridPairs);
  //rough estimate of the number of atom pairs
  int rough = 0;
  for(int i=0;i<nGridPairs;i++){
    int g0 = gridPairs[i*2+0];
    int g1 = gridPairs[i*2+1];
    rough += nResidents0[g0] * nResidents1[g1];
  }
  for(int g=0;g<nTotalGrids;g++){
    rough += nResidents0[g]*nResidents1[g];
  }

  *pairs = (int*) malloc(sizeof(int) * rough*2);
  printf("rough: %d\n", rough);
  int nPairs=0;
  for(int i=0;i<nGridPairs;i++){
    int g0 = gridPairs[i*2+0];
    int g1 = gridPairs[i*2+1];
    for(int j=0; j<nResidents0[g0]; j++){
      int r0 = residents0[headOfList0[g0] + j];
      for(int k=0; k<nResidents1[g1]; k++){
        int r1 = residents1[headOfList1[g1] + k];
        float sum2 = 0.0;
        for(int d=0;d<3;d++){
          float delta = atoms0[r0*3+d] - atoms1[r1*3+d];
          delta -= floor( delta / cell[d] + 0.5 ) * cell[d];
          sum2 += delta*delta;
        }
        if ((lower*lower < sum2) && (sum2 < higher*higher)){
          (*pairs)[nPairs*2+0] = r0;
          (*pairs)[nPairs*2+1] = r1;
          nPairs ++;
        }
      }
    }
  }
  for(int g=0;g<nTotalGrids;g++){
    for(int j=0; j<nResidents0[g]; j++){
      int r0 = residents0[headOfList0[g] + j];
      for(int k=0; k<nResidents1[g]; k++){
        int r1 = residents1[headOfList1[g] + k];
        float sum2 = 0.0;
        for(int d=0;d<3;d++){
          float delta = atoms0[r0*3+d] - atoms1[r1*3+d];
          delta -= floor( delta / cell[d] + 0.5 ) * cell[d];
          sum2 += delta*delta;
        }
        if ((lower*lower < sum2) && (sum2 < higher*higher)){
          (*pairs)[nPairs*2+0] = r0;
          (*pairs)[nPairs*2+1] = r1;
          nPairs ++;
        }
      }
    }
  }
  //deallocate allocated memories
  free(gridPairs);
  return nPairs;
}
  
