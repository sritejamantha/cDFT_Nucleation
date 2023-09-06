#include <complex>
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iostream>
#include "omp.h"
using namespace std ;

#define PI   3.141592653589793238462643383

#define pow2(x) ((x)*(x))
#define pow3(x) ((x)*(x)*(x))
#define pow4(x) ((x)*(x)*(x)*(x))
#define min(A,B) ((A)<(B) ? (A) : (B) )
#define Kdelta(i,j) ((i==j) ? 1 : 0 )
#define Dim 1
#define TNsp 4
#define crsTNsp 6
#define MAXClen 200

#ifndef MAIN
extern
#endif
double *Grandr,amuP,amuSrf,Xa,Xb,*Xar,*Xbr,lamb_mu,**wwall,**q_pgt,**qs_fwd_pgt,**qs_bkd_pgt,**nwrhok,**wk,**extwk,amu[TNsp],**wn0,**wn1,**wn2,**wn3,***wnv1,***wnv2,**knltz1,**knlt1,**kntz1,rf_len,Adisp[7][3],Bdisp[7][3],Lh[Dim],L[Dim],dx[Dim],Lknl[Dim],dxknl[Dim],Tmpr,kijpara[crsTNsp][2],epsln[TNsp],crsepsln[8],sigma[TNsp],rsigma3[TNsp],rsigma[TNsp],Nsig[TNsp],N2sig[TNsp], blkrho[TNsp][2],**rhoK, **dfloc,**extdflocdn0,**extdflocdn3,**extrhok2,*extz,**extrhok,**qk,**tmpN2,**tmpNhf,**extqk,*qtmp,*extqtmp,*tmpar1,*exttmpar1,lamb,**dfchdn3,**dfhsdn3,**dfttdrho,**dfdispndr,***dfdispndr_pir,**dfdispl2dn0,**dfdispl2dn3,**dfdispl1dn0,**dfdispl1dn3,*ttF,*Fhs,*Fch, *Fdispl1,*Fdispl2, *Fdispnl,rhoLA,rhoLB;


#ifndef MAIN
extern
#endif
int ir_intrf,cntpad,wall_para,Nbar,padN2,padN,ext,extM,extM2,M,nthreads, nstot, Nx[Dim],extNx2[Dim],extNx[Dim],ext2Nx[Dim],Nknl[Dim],Nknl2[Dim],NK[TNsp],nK[TNsp],Nbr,sigNhf[TNsp],sigN[TNsp],otp_freq2,otp_freq,ttstep;


#ifndef MAIN
extern
#endif
long idum;

void write_extgrid_data( const char * , double * );
void write_grid_data( const char * , double * );
void write_ext2grid_data( const char *, double * );
double integrate( double * );
double get_k(int , double *);
double integ_trap( int , double*  );
void read_one_resume_file(FILE *, double* );
void read_resume_files() ;
void interp_values(double*  , double ,int );
