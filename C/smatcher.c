#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#include "pairlist.h"
#include "bst.h"
#include "common.h"
#include "smatcher.h"

//Find similar coordinations.

int smatch_len(smatchtype* s)
{
  int n=0;
  while ( s != NULL ){
    n += 1;
    s = s->next;
  }
  return n;
}




smatchtype* smatcher_core(int nOatoms, double* Oatoms, double* cell, double radius, double rmsdmax, int every)
{
  smatchtype* smatch=NULL;
  
  int* prox;
  //atoms of the proximity
  fprintf(stderr, "Preparing the neighbor lists...\n");
  int nProx   = pairlist(nOatoms, Oatoms, nOatoms, Oatoms, 0.0, radius, cell, &prox);
  bnode* nei[nOatoms];
  MakeNeighborList(nOatoms, nProx, prox, nei);
  fprintf(stderr, "Done.\n");
  free(prox);
  fprintf(stderr,"%d nProx A\n", nProx);
  for(int p=0; p<nOatoms; p+=every){
    fprintf(stderr, "\n%.1f %%", p*100./nOatoms);
    //this does not include p itself!
    int* pnei = get_array(nei[p]);
    //for(int jq=0; jq<size(nei[p]); jq++){//q is in proximity of p
    //  int q = pnei[jq];
    for(int q=p+1; q<nOatoms; q++){
      int* qnei = get_array(nei[q]);

      double d[3];
      //offset between centers
      sub(&Oatoms[q*3], &Oatoms[p*3], d);
      double sumsqdev = 0.0;
      //fprintf(stderr, "%d\n", size(nei[q]));
      for(int ip=0; ip<size(nei[p]); ip++){
	double dmin = 1e99;
	int pn = pnei[ip];
	for(int iq=0; iq<size(nei[q]); iq++){
	  int qn = qnei[iq];
	  double dd[3];
	  sub(&Oatoms[qn*3], &Oatoms[pn*3], dd);
	  sub(dd, d, dd);
          for(int d=0;d<3;d++){
            dd[d] -= floor(dd[d]/cell[d]+0.5)*cell[d];
          }
	  double L = dot(dd,dd);
	  if ( L < dmin ){
	    dmin = L;
	  }
	}
	sumsqdev += dmin;
      }
      free(qnei);
      double msd = sumsqdev / size(nei[p]); //msd in nm**2
      double rmsd = sqrt(msd);
      if ( rmsd < rmsdmax ){
	//fprintf(stderr, "%d %d\n", p,q);
	smatchtype* s = (smatchtype*) malloc(sizeof(smatchtype));
	s->next = smatch;
	s->i = p;
	s->j = q;
	s->radius = radius;
	s->d[0] = d[0];
	s->d[1] = d[1];
	s->d[2] = d[2];
	s->rmsd = rmsd;
	smatch = s;
      }
      //return smatch;
    }
    free(pnei);
  }
  //dispose memory as soon as possible
  for(int i=0;i<nOatoms; i++){
    dispose(nei[i]);
  }
  free(Oatoms);
  return smatch;// return a chain.
}






int main(int argc, char* argv[])
{
  //usage: matcher grofile radius(nm) rmsd(nm)
  //typical values: rmsd == 0.06 nm for TPPI
  if ( argc != 4 ){
    fprintf(stderr, "usage: %s grofile radius rmsdmax\n", argv[0]);
    exit(1);
  }
  double cell[3];
  double *Oatoms;
  FILE *file = fopen(argv[1], "r");
  double rmsdmax;
  sscanf(argv[3], "%lf", &rmsdmax);
  int rel=0;
  int nOatoms = LoadGRO(file, &Oatoms, cell, rel);
  fclose(file);

  double radius = atof(argv[2]);
  fprintf(stderr,"System   %d %f %f %f\n", nOatoms, cell[0], cell[1], cell[2]);

  int N = nOatoms * nOatoms;
  int every = 1;
  if ( N > 100000 ){
    //evaluate only 10 M pairs.
    every = nOatoms*nOatoms/(2*100000)+1;
    fprintf(stderr, "Too many combinations! Will skip every %d to reduce time.\n",every);
  }

  smatchtype* smatch = smatcher_core(nOatoms, Oatoms, cell, radius, rmsdmax, every);
  while(smatch != NULL){
    int p = smatch->i;
    int q = smatch->j;
    double radius = smatch->radius;
    double d[3];
    d[0] = smatch->d[0];
    d[1] = smatch->d[1];
    d[2] = smatch->d[2];
    double rmsd = smatch->rmsd;
    printf("%d %d %f %f %f %f %f\n", p, q, radius, d[0], d[1], d[2], rmsd);
    smatchtype* s = smatch;
    smatch = smatch->next;
    free(s);
  }
}




