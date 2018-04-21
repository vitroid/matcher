#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#include "pairlist.h"
#include "bst.h"
#include "common.h"
#include "matcher.h"

//read Gromacs-type data and find the OH covalent bonds.

//benchmark on vitroid-black 2017-3-31
//original matcher (no bst):             user    1m28.847s
//matcher with bst:                      user    1m48.723s
//matcher with bst (size() is improved): user    1m48.942s

#define max(A,B) ((A)>(B)?(A):(B))
int
isClose(double x, double y)
{
  if ( ( x == 0.0 ) && (y == 0.0)){
    return TRUE;
  }
  double e = (x-y)/max(x,y);
  return (-1e-3 < e) && (e < +1e-3);
}

  

double cosine(double *a, double*b){
  return dot(a,b) / (vector_length(a)*vector_length(b));
}


void
cross(double*x, double* y, double* z)
{
  z[0] = x[1]*y[2] - x[2]*y[1];
  z[1] = x[2]*y[0] - x[0]*y[2];
  z[2] = x[0]*y[1] - x[1]*y[0];
}



double distance(double*a, double* b)
{
  double c[3];
  sub(a,b,c);
  return vector_length(c);
}

void
mul(double*x, double a, double* z)
{
  for(int d=0;d<3;d++){
    z[d] = x[d] * a;
  }
}


void
regularize(double* x, double* y, double* z)
{
  double z0[3];
  for(int d=0;d<3;d++)
    z0[d] = z[d];
  double error = dot(x,y);
  double ye[3];
  mul(y,error/2,ye);
  double x_ort[3];
  sub(x,ye,x_ort);
  double xe[3];
  mul(x,error/2,xe);
  double y_ort[3];
  sub(y,xe,y_ort);
  double z_ort[3];
  cross(x_ort, y_ort, z_ort);
  mul(x_ort, 0.5*(3.0-dot(x_ort,x_ort)), x);
  mul(y_ort, 0.5*(3.0-dot(y_ort,y_ort)), y);
  mul(z_ort, 0.5*(3.0-dot(z_ort,z_ort)), z);
  //allow mirroring
  if (dot(z0, z) < 0){
    mul(z, -1, z);
  }
}




double
det(double* A, double* B, double* C)
{
  // Ax Bx Cx   a b c
  // Ay By Cy = d e f
  // Az Bz Cz   g h k
  double a = A[0];
  double b = B[0];
  double c = C[0];
  double d = A[1];
  double e = B[1];
  double f = C[1];
  double g = A[2];
  double h = B[2];
  double k = C[2];
  return a * (e * k - f * h) - b * (d * k - f * g) + c * (d * h - e * g);
}





void
inv(double* A, double* B, double* C)
{
  // Ax Bx Cx   a b c
  // Ay By Cy = d e f
  // Az Bz Cz   g h k
  double D = det(A,B,C);
  double a = A[0];
  double b = B[0];
  double c = C[0];
  double d = A[1];
  double e = B[1];
  double f = C[1];
  double g = A[2];
  double h = B[2];
  double k = C[2];
  A[0] =  (e * k - f * h) / D;
  B[0] = -(b * k - c * h) / D;
  C[0] =  (b * f - c * e) / D;
  A[1] = -(d * k - f * g) / D;
  B[1] =  (a * k - c * g) / D;
  C[1] = -(a * f - c * d) / D;
  A[2] =  (d * h - e * g) / D;
  B[2] = -(a * h - b * g) / D;
  C[2] =  (a * e - b * d) / D;
}

    





int
LoadAR3R(FILE* file, double** atoms, double* cell){
  int nAtoms;
  char line[1024];
  while ( NULL != fgets( line, sizeof(line), file) ){
    if ( 0 == strncmp(line, "@BOX3", 5) ){
      fgets( line, sizeof(line), file);
      sscanf(line, "%lf %lf %lf\n", &cell[0], &cell[1], &cell[2]);
      cell[0] /= 10; //in nm, for gromacs
      cell[1] /= 10;
      cell[2] /= 10;
    }
    else if ( 0 == strncmp(line, "@AR3R", 5) ){
      fgets( line, sizeof(line), file);
      nAtoms = atoi(line);
      *atoms = (double*) malloc( sizeof(double) * nAtoms*3 );
      for(int i=0;i<nAtoms;i++){
        fgets(line, sizeof(line), file);
        double x,y,z;
        sscanf(line, "%lf %lf %lf", &x, &y, &z);
        //atoms distribute around the origin.
        x -= floor(x+0.5);
        y -= floor(y+0.5);
        z -= floor(z+0.5);
        (*atoms)[i*3+0] = x;
        (*atoms)[i*3+1] = y;
        (*atoms)[i*3+2] = z;
      }
    }
  }
  return nAtoms;
}





bnode**
NeighborAtoms(int nAtoms, double* Atoms, double lower, double upper, double* cell)
{
  int* prox;
  //memory leak test
  /*
  for(int i=0;i<10000;i++){
    int nProx   = pairlist(nAtoms, Atoms, nAtoms, Atoms, lower, upper, cell, &prox);
    free(prox);
    }*/
  
  int nProx   = pairlist(nAtoms, Atoms, nAtoms, Atoms, lower, upper, cell, &prox);
  //memory leak test
  /*
  for(int i=0;i<100000;i++){
    bnode** nei = (bnode**) malloc(sizeof(bnode*) * nAtoms);
    MakeNeighborList(nAtoms, nProx, prox, nei);
    for(int j=0;j<nAtoms;j++){
      dispose(nei[j]);
    }
    free(nei);
    }*/
  bnode** nei = (bnode**) malloc(sizeof(bnode*) * nAtoms);
  MakeNeighborList(nAtoms, nProx, prox, nei);
  free(prox);
  return nei;
}


bnode**
NeighborAtoms2(int nAtoms, double* Atoms, int nBtoms, double* Btoms, double lower, double upper, double* cell)
{
  int* prox;
  int nProx   = pairlist(nAtoms, Atoms, nBtoms, Btoms, lower, upper, cell, &prox);
  bnode** nei = (bnode**) malloc(sizeof(bnode*) * nAtoms);
  MakeNeighborList_hetero(nAtoms, nProx, prox, nei);
  free(prox);
  return nei;
}



int match_len(match* s)
{
  int n=0;
  while ( s != NULL ){
    n += 1;
    s = s->next;
  }
  return n;
}



match* matcher_core2(int nOatoms, double* Oatoms,
			 double cell[3],
			 int nunitatoms, double* unitatoms,
			 double unitcell[3],
			 double err,
			 double rprox,
			 int adjdens,
			 int nostore)
{
  double dens0 = nOatoms / (cell[0]*cell[1]*cell[2]);
  double dens1 = nunitatoms / (unitcell[0]*unitcell[1]*unitcell[2]);
  if (adjdens){
    double ratio = pow(dens1 / dens0, 1./3.);
    unitcell[0] *= ratio;
    unitcell[1] *= ratio;
    unitcell[2] *= ratio;
  }
  double unitcelli[3];
  //直方体セル以外に拡張することも容易だが、そこはまた必要が生じた時に。
  double radius = vector_length(unitcell)/2;
  for(int d=0;d<3;d++){
    unitcelli[d] = 1.0/unitcell[d];
  }
  fprintf(stderr,"System   %d %f %f %f\n", nOatoms, cell[0], cell[1], cell[2]);
  fprintf(stderr,"Template %d %f %f %f\n", nunitatoms, unitcell[0], unitcell[1], unitcell[2]);
  fprintf(stderr, "%f %f %d\n", err, rprox, adjdens);
  fprintf(stderr, "\n");
  for(int i=0;i<nunitatoms;i++){
    fprintf(stderr, "%f %f %f\n", Oatoms[i*3+0], Oatoms[i*3+1], Oatoms[i*3+2]);
  }
  fprintf(stderr, "\n");
  for(int i=0;i<nunitatoms;i++){
    fprintf(stderr, "%f %f %f\n", unitatoms[i*3+0], unitatoms[i*3+1], unitatoms[i*3+2]);
  }
  fprintf(stderr, "\n");
  //Find the tetrahedra of a, b, and c cell vectors
  double a = unitcell[0];
  double b = unitcell[1];
  double c = unitcell[2];
  double ab = sqrt(a*a + b*b);
  double ac = sqrt(a*a + c*c);
  double bc = sqrt(b*b + c*c);
  //atoms of the proximity
  fprintf(stderr, "Preparing the neighbor lists.");
  fprintf(stderr,".");
  bnode** nei = NeighborAtoms(nOatoms, Oatoms, 0.0, radius*(1+err), cell);
  //prepare pair lists
  fprintf(stderr, ".");
  bnode** neiA = NeighborAtoms(nOatoms, Oatoms, a*(1-err), a*(1+err), cell);
  fprintf(stderr, ".");
  bnode** neiB = NeighborAtoms(nOatoms, Oatoms, b*(1-err), b*(1+err), cell);
  fprintf(stderr, ".");
  bnode** neiC = NeighborAtoms(nOatoms, Oatoms, c*(1-err), c*(1+err), cell);
  fprintf(stderr, "Done.\n");
  //find triangle PQR that matches the shape
  //There might not be an atom at the origin of the unit cell...2018-4-7
  int ntet=0;
  match* matches = NULL;
  for(int p=0; p<nOatoms; p++){
    fprintf(stderr, "\r%.1f %%", p*100./nOatoms);
    int nnA = size(neiA[p]);
    int* nA = get_array(neiA[p]);
    int nnB = size(neiB[p]);
    int* nB = get_array(neiB[p]);
    int nnC = size(neiC[p]);
    int* nC = get_array(neiC[p]);
    for(int i=0; i<nnA; i++){
      int q = nA[i];
      double pq[3];
      for(int d=0;d<3;d++){
        pq[d] = Oatoms[q*3+d] - Oatoms[p*3+d];
        pq[d] -= floor( pq[d] / cell[d] + 0.5 ) * cell[d];
      }
      for(int j=0; j<nnB; j++){
        int r = nB[j];
        double pr[3];
        for(int d=0;d<3;d++){
          pr[d] = Oatoms[r*3+d] - Oatoms[p*3+d];
          pr[d] -= floor( pr[d] / cell[d] + 0.5 ) * cell[d];
        }
        double qrL = distance(pq, pr);
        if ( (ab*(1-err)<qrL) && (qrL<ab*(1+err)) ){
          for(int k=0; k<nnC; k++){
            int s = nC[k];
            double ps[3];
            for(int d=0;d<3;d++){
              ps[d] = Oatoms[s*3+d] - Oatoms[p*3+d];
              ps[d] -= floor( ps[d] / cell[d] + 0.5 ) * cell[d];
            }
            double qsL = distance(pq, ps);
            double rsL = distance(pr, ps);
            if ( (ac*(1-err)<qsL) && (qsL<ac*(1+err)) && (bc*(1-err)<rsL) && (rsL<bc*(1+err)) ){
                //mathcing without alignment.
                //find the rotation matrix from 
                //A is the unit cell
                //B is the unit cell in the atomic configuration
                //R is the rotation for A
                //Assume A and B has common origin. i.e. no translation
                //Then, if both A and B are not distorted,
                //B = RA, i.e. R = B A^-1
                //if A and B are stacks of row vectors,
                //B = AR, i.e. R = A^-1 B
                //In reality, R wont be a regular rotation matrix, so we have to regularize it.
                //Here is a simple idea.
                //http://stackoverflow.com/questions/23080791/eigen-re-orthogonalization-of-rotation-matrix
                //In the present case, unit cell is orthogonal, so inverse is trivial.
                //Assume A and B are row vectors
                double Rx[3], Ry[3], Rz[3];
                mul(pq, unitcelli[0], Rx);
                mul(pr, unitcelli[1], Ry);
                mul(ps, unitcelli[2], Rz);
                regularize(Rx,Ry,Rz);
                // check if it is really a unitary.
                //double t1[3];
                //cross(Rx,Ry,t1);
                //double vol = dot(t1,Rz);
                //fprintf(stderr, "%f %f %f\n", vector_length(pq)/a, vector_length(pr)/a, vector_length(qr)/ab);
                //fprintf(stderr, "%f %f %f\n", cosine(pq,pr), cosine(pr,ps), cosine(ps,pq));
                //fprintf(stderr, "%f %f %f %f\n", cosine(Rx,Ry), cosine(Ry,Rz), cosine(Rz,Rx), vol);
                //transposition 転置
                double Tx[3],Ty[3],Tz[3];
                Tx[0] = Rx[0];
                Tx[1] = Ry[0];
                Tx[2] = Rz[0];
                Ty[0] = Rx[1];
                Ty[1] = Ry[1];
                Ty[2] = Rz[1];
                Tz[0] = Rx[2];
                Tz[1] = Ry[2];
                Tz[2] = Rz[2];
                //We get R' by regularization, and B' = A R'
                //double ARx[3], ARy[3], ARz[3];
                //mul(Rx, unitcell[0], ARx);
                //mul(Ry, unitcell[1], ARy);
                //mul(Rz, unitcell[2], ARz);
                //this does not include p itself!
                int* neighbors = get_array(nei[p]);
                neighbors[size(nei[p])] = p; // add itself in place of the useless terminator
                // loop with center atoms
                // 単位胞の原子centerをpの位置に平行移動して重ねる。
                // 計算量が増える分、errを小さくしても見落しが少なくなるはず。
                for(int center=0; center<nunitatoms; center++){
                  //slide unit cell to the position of p
                  //with rotation
                  double slidunit[nunitatoms*3];
                  for(int l=0; l<nunitatoms; l++){
                    double rel[3];
                    for(int d=0; d<3; d++){
                      rel[d] = unitatoms[l*3+d]-unitatoms[center*3+d];
                      rel[d] -= floor(rel[d]+0.5);
                      rel[d] *= unitcell[d];
                    }
                    //affine transformation.
                    //Move and rotate center to p
                    slidunit[l*3+0] = dot(rel, Tx) + Oatoms[p*3+0];
                    slidunit[l*3+1] = dot(rel, Ty) + Oatoms[p*3+1];
                    slidunit[l*3+2] = dot(rel, Tz) + Oatoms[p*3+2];
                  }                       
                  double msd = 0.0;
                  int partners[nunitatoms];
		  if ( nunitatoms > 100 ){
		    //この部分が全分子の組みあわせになっていて遅い?
		    //そこで、nei[p]とunitatomsの間のpairlistをさらに生成する。
		    //まずneiの座標を別の配列にコピーする。
		    int nnei = size(nei[p])+1;
		    double pnei[3*nnei];
		    int   inei[nnei];
		    for(int m=0; m<nnei; m++){
		      int ne = neighbors[m];
		      for(int d=0;d<3;d++){
			pnei[3*m+d] = Oatoms[3*ne+d];
		      }
		      inei[m] = ne;
		    }
		    //近い分子はどれぐらいの距離にあるか。
		    //わからない。それをパラメータrproxで与えることにする。
		    bnode** neiX = NeighborAtoms2(nunitatoms, slidunit, nnei, pnei, 0.0, rprox, cell);
		    for(int l=0; l<nunitatoms; l++){
		      //近くにいる分子はneiXに入っている
		      int* nX = get_array(neiX[l]);
		      double dmin = 1e99;
		      //assert (size(neiX[l]) > 0); // causes error if rprox is too short
		      //fprintf(stderr, "NNEI: %d\n", size(neiX[l]));
		      for(int m=0; m<size(neiX[l]); m++){
			//molecular label in remaked list
			int re = nX[m];
			double dd[3];
			sub(&slidunit[l*3], &pnei[re*3], dd);
			for(int d=0;d<3;d++){
			  dd[d] -= floor(dd[d] / cell[d]+0.5) *cell[d];
			}
			//PBC should be here
			double L = dot(dd,dd);
			if ( L < dmin ){
			  dmin = L;
			  partners[l] = inei[re];
			}
			
		      }
		      free(nX);
		      msd += dmin;
		      if ( nunitatoms*rprox*rprox/25 < msd )
			break; //for l loop 
		    }
		    for(int l=0;l<nunitatoms; l++){
		      dispose(neiX[l]);
		    }
		    free(neiX);
		  }
		  else{
		    for(int l=0; l<nunitatoms; l++){
		      double dmin = 1e99;
		      for(int m=0; m<size(nei[p])+1; m++){
			int ne = neighbors[m];
			double dd[3];
			sub(&slidunit[l*3], &Oatoms[ne*3], dd);
			for(int d=0;d<3;d++){
			  dd[d] -= floor(dd[d] / cell[d]+0.5) *cell[d];
			}
			//PBC should be here
			double L = dot(dd,dd);
			if ( L < dmin ){
			  dmin = L;
			  partners[l] = ne;
			}
		      }
		      msd += dmin;
		      if ( nunitatoms*rprox*rprox/25 < msd )
			break;
		    }
		  }
		  double rmsd = sqrt(msd/nunitatoms);
		  if ( rmsd < rprox/5 ){
		    if ( nostore ){
		      printf("%f %d %d ", rmsd, p, center);
		      for(int o=0;o<3; o++){
			printf("%f ", Rx[o]);
		      }
		      for(int o=0;o<3; o++){
			printf("%f ", Ry[o]);
		      }
		      for(int o=0;o<3; o++){
			printf("%f ", Rz[o]);
		      }
		      printf("%d ", nunitatoms);
		      for(int o=0; o<nunitatoms;o++){
			printf("%d ", partners[o]);
		      }
		      printf("\n");
		    }
		    else{
		      match* m = (match*) malloc (sizeof(matches));
		      m->rmsd = rmsd;
		      fprintf(stderr, "%d %d %d %d %d %f\n", p,q,r,s,center,rmsd);
		      m->atom_gro = p;
		      m->atom_unitcell = center;
		      for(int o=0;o<3; o++){
			m->mat[o+0] = Rx[o];
			m->mat[o+3] = Ry[o];
			m->mat[o+6] = Rz[o];
		      }
		      m->natom= nunitatoms;
		      m->atoms = (int*) malloc (sizeof(int)*nunitatoms);
		      for(int l=0;l<nunitatoms;l++){
			m->atoms[l] = partners[l];
		      }
		      m->next = matches;
		      matches = m;
		    }//nostore
		  } //rmsd
		} // if ( !error )          
		free(neighbors);
		ntet += 1;
            } //if(range)
          } //for(k)
        } //if(range)
      } //for(j)
    } //for(i)
    free(nA);
    free(nB);
    free(nC);
  } // for(p)
  //dispose memory as soon as possible
  for(int i=0;i<nOatoms; i++){
    dispose(nei[i]);
    dispose(neiA[i]);
    dispose(neiB[i]);
    dispose(neiC[i]);
  }
  free(nei);
  free(neiA);
  free(neiB);
  free(neiC);
  fprintf(stderr,"%d ntet\n", ntet);
  fprintf(stderr, "%d nmatch\n", match_len(matches));

  free(Oatoms);
  return matches;
}





void usage(char *cmd)
{
  fprintf(stderr, "usage: %s [-e error][-r rprox][-a] grofile template.ar3r error rmsdmax\n", cmd);
  fprintf(stderr, " -e error   Allowance for p,q,r vectors (0.03)\n");
  fprintf(stderr, " -r value   Radius of proximity (0.4 nm)\n");
  fprintf(stderr, " -a         Automatic densityadjustment\n");
  exit(1);
}


int main(int argc, char* argv[])
{
  /*
    usage: matcher grofile error rmsdmax

    groファイルで構造データが与えられる。また、単位胞が別に与えられる。
    まず、ある原子pから、a軸の距離にある点q、b軸の距離にある点r、c軸の距離にある点sをさがす。(この時の距離の遊びがerror値)
    qr距離が|a-b|にほぼ一致し、qs距離が|a-c|にほぼ一致し、rs距離が|b-c|にほぼ一致するなら、そこには単位胞とおなじ周期構造がある可能性がある。
    そこで単位胞のなかから1原子centerを選び、それがpになるように、そして単位胞の基本
ベクトルがpq, pr, psに平行になるように単位胞を空間配置する。
単位胞の各原子に最も近い、構造データ内の原子をさがしだし、その対応表を出力する。
    これを、p,q,r,s,centerに関して繰り返す。
  */

  int arg = 1;
  double err = 0.03;
  double rprox = 0.4;
  int adjdens = 0;
  while( argv[arg][0] == '-' ){
    if ( 0 == strncmp(argv[arg], "-e", 2) ){
      sscanf(argv[arg+1], "%lf", &err);
      arg += 2;
    }
    else if ( 0 == strncmp(argv[arg], "-r", 2) ){
      sscanf(argv[arg+1], "%lf", &rprox);
      arg += 2;
    }
    else if ( 0 == strncmp(argv[arg], "-a", 2) ){
      adjdens = 1;
      arg += 1;
    }
    else {
      usage(argv[0]);
    }
  }      
  if ( arg +2 != argc ){
    usage(argv[0]);
  }
  double cell[3];
  double *Oatoms;
  FILE *file = fopen(argv[arg], "r");
  int nOatoms = LoadGRO(file, &Oatoms, cell);
  fclose(file);

  double unitcell[3];
  double *unitatoms; //relative
  file = fopen(argv[arg+1], "r");
  int nunitatoms = LoadAR3R(file, &unitatoms, unitcell);
  fclose(file);

  match* matches = matcher_core2(nOatoms,    Oatoms,    cell,
				   nunitatoms, unitatoms, unitcell,
				   err,
				   rprox,
				   adjdens,
				   1 //nostore
				   );
  //nostore then useless
  while ( matches != NULL ){
    match* m = matches;
    printf("%f %d %d ", m->rmsd, m->atom_gro, m->atom_unitcell);
    for(int i=0;i<9; i++){
      printf("%f ", m->mat[i]);
    }
    printf("%d ", m->natom);
    for(int i=0; i<m->natom;i++){
      printf("%d ", m->atoms[i]);
    }
    printf("\n");
    matches = m->next;
    free(m->atoms);
    free(m);
  }
}



