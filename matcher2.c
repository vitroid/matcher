#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "svd.h"
#include "neighborlist.h"
#include "common.h"



void multiply_t2(int m, int n, int p, double a[m*n], double b[p*n], double c[m*p])
//second arg b is transposed.
{
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < p; j++) {
            c[i*p+j] = 0;
            for (int k = 0; k < n; k++) {
                c[i*p+j] += a[i*n+k] * b[j*n+k];
            }
        }
    }
}


void multiply(int m, int n, int p, double a[m*n], double b[n*p], double c[m*p])
{
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < p; j++) {
            c[i*p+j] = 0;
            for (int k = 0; k < n; k++) {
                c[i*p+j] += a[i*n+k] * b[k*p+j];
            }
        }
    }
}




double
rot_trans(int nv, double uv[], double gv[], double R[9], double t[3])
/*
  determine affine projection from uv to gv

  returns:
  R: rotation matrix
  t: translation vector
  return: rmsd
*/
{
  double ucom[3], gcom[3];
  for(int d=0; d<3; d++){
    ucom[d] = 0.0;
    for(int i=0;i<nv; i++){
      ucom[d] += uv[i*3+d];
    }
    ucom[d] /= nv;
    gcom[d] = 0.0;
    for(int i=0;i<nv; i++){
      gcom[d] += gv[i*3+d];
    }
    gcom[d] /= nv;
  }
  double H[3*3];
  for(int j=0; j<3; j++){
    for(int k=0; k<3; k++){
      H[j*3+k] = 0.0;
      for(int i=0; i<nv; i++){
	H[j*3+k] += (uv[i*3+j]-ucom[j])*(gv[i*3+k]-gcom[k]);
      }
    }
  }
  double w[3];
  double v[3*3];
  dsvd(3, 3, H, w, v);
  //u is returned in H
  multiply_t2(3,3,3,H,v,R);
  multiply(1,3,3,ucom,R,t);
  for(int d=0;d<3;d++){
    t[d] = gcom[d] - t[d];
  }
  double I = 0.0;
  for(int i=0; i<nv; i++){
    double delta[3];
    multiply(1,3,3,&uv[i*3],R,delta);
    for(int d=0; d<3; d++){
      delta[d] += t[d] - gv[i*3+d];
      I += delta[d]*delta[d];
    }
  }
  double rmsd = sqrt(I/nv);
  return rmsd;
}


void quick_proximity(int n, int* order, double* distance, int ntop)
//orderにははじめ0ないしn-1の数字が並んでいる。
//これをdistance順に並べ、上位ntop個だけを返す。(あとは放置)
{
  int pivot = n/2;
  double pd = distance[order[pivot]];
  int head = 0;
  int tail = n-1;
  int count=4;
  if ( head >= tail )
    return;
  while ( 1 ){
    while ( distance[order[head]] < pd ){
      head += 1;
    }
    while (distance[order[tail]] > pd ){
      tail -= 1;
    }
    /*
    fprintf(stderr, "H%d T%d %d\n", head,tail,ntop);
    for(int i=0;i<n;i++){
      if ( i == head )
	fprintf(stderr, "H");
      if ( i == tail )
	fprintf(stderr, "T");
      if ( i == pivot )
	fprintf(stderr, "P");
      fprintf(stderr, "%f ", distance[order[i]]);
    }
    fprintf(stderr, "\n");
    */
    if ( head >= tail )
      break;
    int tmp = order[head];
    order[head] = order[tail];
    order[tail] = tmp;
    head += 1;
    tail -= 1;
    /*
    fprintf(stderr, "H%d T%d %d\n", head,tail,ntop);
    for(int i=0;i<n;i++){
      if ( i == head )
	fprintf(stderr, "H");
      if ( i == tail )
	fprintf(stderr, "T");
      if ( i == pivot )
	fprintf(stderr, "P");
      fprintf(stderr, "%f ", distance[order[i]]);
    }
    fprintf(stderr, "\n\n");
    */
  }
  if ( head > ntop ){
    quick_proximity(head, order, distance, ntop);
    return;
  }
  if ( tail < ntop ){
    quick_proximity(n-tail-1, &order[tail+1], distance, ntop-tail-1);
    return;
  }
  //quick_sort test
  /*quick_proximity(head, order, distance, 0);
  quick_proximity(n-tail-1, &order[tail+1], distance, 0);
  */
}
  



void quick_sort(int n, int* order, double* distance)
{
  int head = 0;
  int tail = n-1;
  if ( head >= tail )
    return;
  int pivot = n/2;
  double pd = distance[order[pivot]];
  while ( 1 ){
    while ( distance[order[head]] < pd ){
      head += 1;
    }
    while (distance[order[tail]] > pd ){
      tail -= 1;
    }
    if ( head >= tail )
      break;
    int tmp = order[head];
    order[head] = order[tail];
    order[tail] = tmp;
    head += 1;
    tail -= 1;
  }
  quick_sort(head, order, distance);
  quick_sort(n-tail-1, &order[tail+1], distance);
}
  




match*
matcher2_core(int nuatoms,
              double* uatoms,
              double ucell[3],
              int ngatoms,
              double* gatoms,
              double gcell[3],
              int adjdens,
              int nostore)
{
  //密度の調整
  if ( adjdens ){
    double gdens = ngatoms / (gcell[0]*gcell[1]*gcell[2]);
    double udens = nuatoms / (ucell[0]*ucell[1]*ucell[2]);
    double ratio = pow(udens/gdens, 1./3.);
    fprintf(stderr, "Density correction: %f\n", ratio);
    ucell[0] *= ratio;
    ucell[1] *= ratio;
    ucell[2] *= ratio;
  }
  match* matches = NULL;
  // 単位胞の外接球の半径
  double uR = (ucell[0]+ucell[1]+ucell[2])/2.;
  // 近接距離。グリッドのサイズ
  double rprox = 0.4;
  // まず重ねる点を決める。
  int ucenter = 0;
  for(int gcenter=0; gcenter<ngatoms; gcenter++){
    if (gcenter*100/ngatoms != (gcenter-1)*100/ngatoms)
      fprintf(stderr, "Progress %d%%\n", gcenter*100/ngatoms);
    // それらの点を原点に平行移動する。
    // 絶対座標に変換する
    double sgatoms[ngatoms*3];
    double sgratoms[ngatoms*3];
    for(int i=0;i<ngatoms; i++){
      for(int d=0;d<3;d++){
	double x = gatoms[i*3+d] - gatoms[gcenter*3+d];
	x -= floor(x+0.5);
	sgratoms[i*3+d] = x;
	sgatoms[i*3+d] = x*gcell[d];
      }
    }
    double suatoms[nuatoms*3];
    for(int i=0;i<nuatoms; i++){
      for(int d=0;d<3;d++){
	double x = uatoms[i*3+d] - uatoms[ucenter*3+d];
	x -= floor(x+0.5);
	suatoms[i*3+d] = x*ucell[d];
      }
    }
    // 住所録を作る
    int ggrid[3];
    for(int d=0;d<3;d++){
      ggrid[d] = (int)floor(gcell[d]/rprox);
    }
    //fprintf(stderr, "%d %d %d\n", ggrid[0], ggrid[1], ggrid[2]);
    sAddressBook* abook = AddressBook(ggrid, ngatoms, sgratoms);
    // 原点からの距離をリストにする。
    // gについては、全部は要らない。
    int* gorig;
    int ngorig = near_origin( abook, &gorig );
    double distance[ngorig];
    for(int i=0; i<ngorig; i++){
      int j = gorig[i];
      double sum=0.0;
      for(int d=0; d<3; d++){
	double x = sgatoms[j*3+d];
	sum += x*x;
      }
      distance[i] = sum;
    }
    int gorder[ngorig];
    for(int i=0;i<ngorig;i++)
      gorder[i] = i;
    quick_sort(ngorig, gorder, distance);
    /*
    fprintf(stderr, "%d ngorig\n", ngorig);
    for(int i=0;i<ngorig;i++){
      fprintf(stderr, "%f ", distance[gorder[i]]);
    }
    fprintf(stderr, "\n");
    exit(0);
    */
    for(int i=0;i<ngorig;i++)
      gorder[i] = gorig[gorder[i]];
    free(gorig);
    double udistance[nuatoms];
    for(int i=0; i<nuatoms; i++){
      double sum=0.0;
      for(int d=0; d<3; d++){
	double x = suatoms[i*3+d];
	sum += x*x;
      }
      udistance[i] = sum;
    }
    int uorder[nuatoms];
    for(int i=0;i<nuatoms;i++)
      uorder[i] = i;
    quick_sort(nuatoms, uorder, udistance);
    /* it works
    fprintf(stderr, "%d nuatoms\n", nuatoms);
    for(int i=0;i<nuatoms;i++){
      fprintf(stderr, "%d %f ", uorder[i], udistance[uorder[i]]);
    }
    fprintf(stderr, "\n");
    exit(0);
    exit(0);
    */
    for(int mi=1; mi<5; mi++){
      int ni = gorder[mi];
      for(int mj=1; mj<5; mj++){
	int nj = gorder[mj];
	if ( ni != nj ){
	  for(int mk=1; mk<5; mk++){
	    int nk = gorder[mk];
	    if ( ( ni != nk ) && ( nj != nk ) ){
	      int glist[nuatoms];
	      glist[0] = gcenter;
	      glist[1] = ni;
	      glist[2] = nj;
	      glist[3] = nk;
	      // その状態で並進と回転を最適化する。
	      double uv[nuatoms*3];
	      double gv[nuatoms*3];
	      int nv = 4;
	      for(int i=0;i<nv;i++){
		for(int d=0;d<3;d++){
		  uv[i*3+d] = suatoms[uorder[i]*3+d];
		  gv[i*3+d] = sgatoms[glist[i]*3+d];
		}
	      }
	      double R[9], t[3];
	      // uvをRだけ回転しtだけ並進するとgvに重なる。
	      double rmsd = rot_trans(nv, uv, gv, R, t);
	      //printf("%d %f\n", nv,rmsd);
	      // 粒子を1つ増やす。
	      for(; nv<nuatoms; nv++){
		if ( rmsd > 0.1 )
		  break;
		// 次のuはソートされたリストから選ぶ。
		//abs
		for(int d=0;d<3;d++){
		  uv[nv*3+d] = suatoms[uorder[nv]*3+d];
		}
		// gは近接点のなかからさがす。重複してもよい。
		double uvRt[3];
		multiply(1,3,3,&uv[nv*3], R, uvRt);
		for(int d=0;d<3;d++){
		  uvRt[d] += t[d];
		  uvRt[d] /= gcell[d]; //abs to rel
		}
		int nearest = find_nearest(uvRt, abook, gcell, sgratoms);
		/*
		  printf("%f %f %f\n", uvRt[0], uvRt[1], uvRt[2]);
		  printf("%f %f %f\n",
		  sgratoms[nearest*3+0],
		  sgratoms[nearest*3+1],
		  sgratoms[nearest*3+2]);
		*/
		glist[nv] = nearest;
		for(int d=0;d<3;d++){
		  gv[nv*3+d] = sgatoms[nearest*3+d];
		}
		// uvをRだけ回転しtだけ並進するとgvに重なる。
		rmsd = rot_trans(nv+1, uv, gv, R, t);
                if (rmsd > 0.1 ){
                  break;
                }
		//printf("%d %d %f\n", nv+1, nearest, rmsd);
	      } //for nv
	      if ( rmsd < 0.1 ){
                if ( nostore ){
                  printf("%f %d %d ", rmsd, gcenter, ucenter);
                  for(int i=0; i<9; i++){
                    printf("%f ", R[i]);
                  }
                  printf("%d ", nuatoms);
                  int corr[nuatoms];
                  for(int i=0; i<nuatoms; i++){
                    //printf("%d %d\n", uorder[i], glist[i]);
                    corr[uorder[i]] = glist[i];
                  }
                  for(int i=0; i<nuatoms; i++){
                    printf("%d ", corr[i]);
                  }
                  printf("\n");
                }
                else{
                  match* ma = (match*)malloc(sizeof(match));
                  ma->rmsd = rmsd;
                  ma->atom_gro = gcenter;
                  ma->atom_unitcell = ucenter;
                  for(int i=0;i<9;i++){
                    ma->mat[i] = R[i];
                  }
                  ma->natom = nuatoms;
                  ma->atoms = (int*) malloc(sizeof(int)*nuatoms);
                  int corr[nuatoms];
                  for(int i=0; i<nuatoms; i++){
                    //printf("%d %d\n", uorder[i], glist[i]);
                    corr[uorder[i]] = glist[i];
                  }
                  for(int i=0;i<nuatoms;i++){
                    ma->atoms[i] = corr[i];
                  }
                  ma->next = matches;
                  matches = ma;
                }
	      } //if rmsd
	    } // if
	  } //for mk
	}// if
      }//for mj
    }//for mi
    
    dispose_addressbook(abook);
  }//for gcenter;
  //no return value
  return matches;
}



int main(int argc, char *argv[]){
  int arg = 1;

  double gcell[3];
  double *gatoms;
  FILE *file = fopen(argv[arg], "r");
  int ngatoms = LoadGRO(file, &gatoms, gcell);
  fclose(file);
  //absolute to relative
  for(int i=0; i<ngatoms; i++){
    for(int d=0; d<3; d++){
      double x = gatoms[i*3+d]/gcell[d];
      gatoms[i*3+d] = x - floor(x+0.5);
    }
  }
  /* test of abook and find_nearest. it just works.
  int grid[3] = {10,10,10};
  sAddressBook* abook = AddressBook(grid, ngatoms, gatoms);
  double rloc[3] = {0.3, 0.5, 0.5};
  for(int d=0;d<3;d++){
    rloc[d] /= gcell[d];
  }
  int nearest = find_nearest(rloc, abook, gcell, gatoms);
  printf("%d %f %f %f\n",nearest,
	 gatoms[nearest*3+0]*gcell[0],
	 gatoms[nearest*3+1]*gcell[1],
	 gatoms[nearest*3+2]*gcell[2]);
  exit(0);
  */
  
  double ucell[3];
  double *uatoms;
  file = fopen(argv[arg+1], "r");
  int nuatoms = LoadGRO(file, &uatoms, ucell);
  fclose(file);
  //absolute to relative
  for(int i=0; i<nuatoms; i++){
    for(int d=0; d<3; d++){
      double x = uatoms[i*3+d]/ucell[d];
      uatoms[i*3+d] = x - floor(x+0.5);
    }
  }
  /*
  double distance[]={ 3.0, 9.0, 4.0, 9.0, 5.0, 9.0, 2.0, 6.0, 5.0, 3.0, 5.0, 8.0 };
  int order[] = {0,1,2,3,4,5,6,7,8,9,10,11};
  int norder = 12;
  quick_proximity(norder, order, distance, 4);
  for(int i=0;i<norder;i++){
    fprintf(stderr, "%f ", distance[order[i]]);
  }
  fprintf(stderr, "\n");
  exit(0);
  */
  int adjdens = 1;
  int nostore = 1;
  matcher2_core(nuatoms, uatoms, ucell, ngatoms, gatoms, gcell, adjdens, nostore);
}

  


void
test(int argc, char *argv[]){
  //compare with the result of test() in matcher2.py
  double H[3*3] = { 2.0, 3.0, 5.0, 1.0, 4.0, 2.0, 3.0, 5.0, 1.0 };
  double w[3];
  double v[3*3];
  double R[3*3];
  dsvd(3, 3, H, w, v);
  //u is returned in H
  multiply_t2(3,3,3,H,v,R);
  printf("H\n");
  for(int i=0;i<3;i++){
    for(int j=0;j<3;j++){
      printf("%f ", H[i*3+j]);
    }
    printf("\n");
  }
  printf("V\n");
  for(int i=0;i<3;i++){
    for(int j=0;j<3;j++){
      printf("%f ", v[i*3+j]);
    }
    printf("\n");
  }
  printf("R\n");
  for(int i=0;i<3;i++){
    for(int j=0;j<3;j++){
      printf("%f ", R[i*3+j]);
    }
    printf("\n");
  }
}
  
  
      
  
