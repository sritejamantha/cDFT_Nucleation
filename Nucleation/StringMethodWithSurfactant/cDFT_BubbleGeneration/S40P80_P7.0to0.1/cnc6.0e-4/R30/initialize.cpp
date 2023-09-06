#include "globals.h"
void allocate( void ) ;
void read_input( void ) ;
void unstack_input(int , int* , int* ) ;
void read_matdisp(void);

void initialize() {
  int i,j,nn[Dim] ;

  idum =   -long( time(0) ) ; // 9 ; //

  read_input() ;

  //cout<<"here"<<endl;
  read_matdisp();

  lamb_mu = -0.5;
 
  rf_len = rsigma[0];

  M = 1;
 
  for(i=0;i<Dim;i++){
	L[i] /= rf_len;
	dx[i] /= rf_len;
	Nx[i] = floor(L[i]/dx[i])+1;
  	M *=  Nx[i];
  	Lh[i] = L[i] *0.5;
  }

  padN = 0; 
  padN2 = 0; 

  for(i=0;i<TNsp;i++){
 	rsigma[i] /= rf_len;
	sigma[i] = rsigma[i]*(1. - 0.12*exp(-3.0*epsln[i]/Tmpr));
	
	rsigma3[i] = pow(rsigma[i],3);	
	Nsig[i] = round(sigma[i]/dx[0]/2.0)*dx[0]*2.0;
	N2sig[i] = 2.0*round(sigma[i]/dx[0])*dx[0];
	sigNhf[i]  = round(sigma[i]/2./dx[0]);
	sigN[i] = round(sigma[i]/dx[0]);

	if(sigN[i] > padN){
		
		padN2 =5*sigN[i];
		padN = sigN[i] ; 
	}

        cout<<"type i: sigma "<<sigma[i] <<" Nsig "<<Nsig[i]<<" N2sig "<<N2sig[i]<<" sigNhf "<<sigNhf[i]<<" sigN "<<sigN[i]<<endl;
  }

  for(i=0;i<Dim;i++){
	extNx[i] = Nx[i];
  	extNx2[i] = Nx[i] ;
  }

  extNx[Dim-1] += 2*padN;
  extNx2[Dim-1] += 2*padN2;

  cout<<"padN "<<padN<<"; padN2 "<<padN2<<endl;

  extM = M + 2*padN;
  extM2 = M + 2*padN2;

  cout<<"M "<<M<< " extM "<<extM<<" extM2 "<<extM2<<endl;
	

  int tmp_count=0;
  for(i=0;i<TNsp-1;i++)
  {
   for(j=i+1;j<TNsp;j++)
   {
        crsepsln[tmp_count] = (kijpara[tmp_count][0]*Tmpr - kijpara[tmp_count][1]);

        crsepsln[tmp_count] = (1.0-  crsepsln[tmp_count] )*sqrt(epsln[i]*epsln[j]);

        cout<<"crsepsln type "<<tmp_count<<": "<<crsepsln[tmp_count]<<endl;
        tmp_count=tmp_count+1;
   }
  }

  allocate( ) ;


  double tmp_dist,xInp,cofac;
  //int mid_ind = int(2.0*M/3.);
  int mid_ind =150;

  /*FILE *f0,*f1,*f2,*f3;
  f0 = fopen("rho0_Inp.dat","r");
  f1 = fopen("rho1_Inp.dat","r");
  f2 = fopen("rho2_Inp.dat","r");
  f3 = fopen("rho3_Inp.dat","r");*/
//  printf("%lf \t %d \t %lf \t %d\n",RgA,NRgA,RgB,NRgB);

  /*cofac = blkrho[1][1]/blkrho[3][1];
  for(j=0;j<M;j++)
  {
   fscanf(f0,"%lf %lf",&xInp,&rhoK[0][j]);
   fscanf(f1,"%lf %lf",&xInp,&rhoK[1][j]);
   fscanf(f2,"%lf %lf",&xInp,&rhoK[2][j]);
   fscanf(f3,"%lf %lf",&xInp,&rhoK[3][j]);
  }
  fclose(f0);fclose(f1);fclose(f2);fclose(f3);*/

    for(i=0;i<TNsp;i++){

        for(j=0;j<M;j++){


                unstack_input(j,nn,Nx);
                if(wall_para == 0){
                        //tmp_dist = j;//nn[Dim-1]*dx[Dim-1] - L[Dim-1]/2.0;//+ (i != 1 ?Nsig[i]/2.0 : 0.0);
                        //tmp_dist = nn[Dim-1]*dx[Dim-1] -  L[Dim-1]/2.0;
                        //rhoK[i][j] = (tmp_dist > 0 ? blkrho[i][1]: blkrho[i][0]);
                        if(j == mid_ind)
                                rhoK[i][j] = blkrho[i][1]*0.5 + blkrho[i][0]*0.5;
                        else if(j <mid_ind)
                                rhoK[i][j] = blkrho[i][0];
                        else
                                rhoK[i][j] = blkrho[i][1];
                        //(blkrho[i][0]-blkrho[i][1])*erfc( tmp_dist/sigma[i])/2.0 + blkrho[i][1];

                }
                else{


                        tmp_dist = nn[Dim-1]*dx[Dim-1] - Nsig[i];
                         rhoK[i][j] = (blkrho[i][0]-blkrho[i][1])*erfc( 5*tmp_dist/sigma[i])/2.0 + blkrho[i][1];
                }

                q_pgt[0][j] = 1;

                if((nn[Dim-1]<(sigN[i]-1)) and (wall_para > 0))
                        wwall[i][j] = 1.e10;
                else
                        wwall[i][j] = 0.0;
        }

        for(j=0;j<(2*sigNhf[i]+1);j++){

          knlt1[i][j] = -sigNhf[i]*dx[Dim-1]+ j*dx[Dim-1];

          knltz1[i][j] = Nsig[i]*Nsig[i]/4.0 -pow(knlt1[i][j],2);
//        cout<<"t1 "<<knlt1[i][j]<<" t2  "<<knltz1[i][j]<<endl;
        }


    }
  


   for(j=0;j<M;j++)
   {
     Xar[j] = rhoK[TNsp-2][j]/(rhoK[TNsp-2][j] + rhoK[TNsp-1][j]);
     Grandr[j] = 0.0;
   }

   read_resume_files();

   if(Dim ==1)
   {
    for(i=0 ;i<extM2; i++)
    {	
       extz[i] = -1.*dx[Dim-1]*padN2 + i*dx[Dim-1]; 

    }
     
   }

   for(j=0;j<M;j++)
    {
     qs_fwd_pgt[0][j]=1;
     qs_bkd_pgt[0][j]=1;
    }

//   exit(-1);


}


void allocate(){


// int tmp_chl=NK[1],i,j; 
 int tmp_chl,i,j;
  
// if (NK[1]> MAXClen)
 	tmp_chl = MAXClen;
 
 Xar =  (double*)  calloc(M, sizeof( double ) ) ;
 Grandr = (double*)  calloc(M, sizeof( double ) ) ;
 q_pgt = (double**)  calloc(tmp_chl,sizeof( double* ) ) ;
 qs_fwd_pgt = (double**)  calloc(tmp_chl,sizeof( double* ) ) ;
 qs_bkd_pgt = (double**)  calloc(tmp_chl,sizeof( double* ) ) ;
 nwrhok =  (double**)  calloc(TNsp,sizeof( double* ) ) ;
 wk = (double**)  calloc(TNsp,sizeof( double* ) ) ;
 wwall = (double**)  calloc(TNsp,sizeof( double* ) ) ;
 rhoK = (double**)  calloc(TNsp+1,sizeof( double* ) ) ; 
 extrhok = (double**) calloc(TNsp,sizeof( double* ) ) ;
 extrhok2 = (double**) calloc(TNsp,sizeof( double* ) ) ;
 extwk = (double**) calloc(TNsp,sizeof( double* ) ) ;

 knlt1 = (double**) calloc(TNsp,sizeof( double* ) ) ;
 knltz1 = (double**) calloc(TNsp,sizeof( double* ) ) ;
 wn0 = (double**)  calloc(TNsp+1,sizeof( double* ) ) ; 
 wn1 = (double**)  calloc(TNsp+1,sizeof( double* ) ) ; 
 wn2 = (double**)  calloc(TNsp+1,sizeof( double* ) ) ; 
 wn3 = (double**)  calloc(TNsp+1,sizeof( double* ) ) ; 
 
 wnv1 = (double***)  calloc(TNsp+1,sizeof( double** ) ) ; 
 wnv2 = (double***)  calloc(TNsp+1,sizeof( double** ) ) ; 

 extdflocdn0 = (double**) calloc(TNsp,sizeof( double* ) ) ;
 extdflocdn3 = (double**) calloc(TNsp,sizeof( double* ) ) ;

 dfloc = (double**) calloc(TNsp,sizeof( double* ) ) ;
 
 tmpN2 = (double**) calloc(nthreads,sizeof( double* ) ) ;


 tmpNhf = (double**) calloc(nthreads,sizeof( double* ) ) ;

 dfdispndr_pir = (double***)  calloc(TNsp+1,sizeof( double** ) ) ;
 dfdispndr = (double**)  calloc(TNsp+1,sizeof( double* ) ) ;

 dfttdrho = (double**)  calloc(TNsp,sizeof( double* ) ) ;

 dfdispl1dn0 =  (double**)  calloc(TNsp+1,sizeof( double* ) ) ;
 dfdispl1dn3 =  (double**)  calloc(TNsp+1,sizeof( double* ) ) ;

 dfdispl2dn0 =  (double**)  calloc(TNsp+1,sizeof( double* ) ) ;
 dfdispl2dn3 =  (double**)  calloc(TNsp+1,sizeof( double* ) ) ;
 


 dfhsdn3 = (double**)  calloc(TNsp+1,sizeof( double* ) ) ;

 dfchdn3 = (double**)  calloc(TNsp+1,sizeof( double* ) ) ;

 extz =  (double*)  calloc(extM2, sizeof( double )); 

 ttF =  (double*)  calloc(M, sizeof( double ));

for(i=0;i< tmp_chl;i++)
{
	q_pgt[i] =  (double*)  calloc(M, sizeof( double ) ) ;
        qs_fwd_pgt[i] = (double*)  calloc(M, sizeof( double ) ) ;
        qs_bkd_pgt[i] = (double*)  calloc(M, sizeof( double ) ) ;
}

for(i=0;i< TNsp;i++){
	wwall[i] = (double*)  calloc(M, sizeof( double ) ) ;
	wk[i] =  (double*)  calloc(M, sizeof( double ) ) ;
	nwrhok[i] = (double*)  calloc(M, sizeof( double ) ) ;
	dfttdrho[i] = (double*)  calloc(M, sizeof( double ) ) ;

	knlt1[i] = (double*)  calloc(2*sigNhf[i]+1, sizeof( double ) ) ;
 	knltz1[i] = (double*)  calloc(2*sigNhf[i]+1, sizeof( double ) ) ;
  	tmpar1 = (double*)  calloc(M, sizeof( double ) ) ;
	extrhok[i] = (double*)  calloc(extM, sizeof( double ) ) ;
	extrhok2[i] = (double*)  calloc(extM2, sizeof( double ) ) ;
	extwk[i] =  (double*)  calloc(extM2, sizeof( double ) ) ;

//	tmpN2[i] = (double* )calloc(padN2, sizeof( double ) ) ;
	//tmpNhf[i] = (double* )calloc(2*sigNhf[i]+1, sizeof( double ) ) ;	
	extdflocdn0[i] = (double*) calloc(extM,sizeof( double ) ) ; 
 	extdflocdn3[i] = (double*) calloc(extM,sizeof( double ) ) ;
 	dfloc[i] =  (double*) calloc(M,sizeof( double )) ;
 }

 for(i=0; i< nthreads; i++){
    tmpNhf[i] = (double* )calloc(padN+3, sizeof( double ) );
    tmpN2[i] = (double* )calloc(padN2, sizeof( double ) );
 }

 Fhs = (double*)  calloc(M, sizeof( double ) ) ;
 Fch = (double*)  calloc(M, sizeof( double ) ) ;
 Fdispl1 = (double*)  calloc(M, sizeof( double ) ) ;
 Fdispl2 = (double*)  calloc(M, sizeof( double ) ) ;
 Fdispnl =(double*)  calloc(M, sizeof( double ) ) ;

 for(i=0;i< TNsp+1;i++){
	rhoK[i] = (double*)  calloc(M, sizeof( double ) ) ;

	wn0[i] =  (double*)  calloc(M, sizeof( double ) ) ;
	wn1[i] =  (double*)  calloc(M, sizeof( double ) ) ;

	wn2[i] =  (double*)  calloc(M, sizeof( double ) ) ;

	wn3[i] =  (double*)  calloc(M, sizeof( double ) ) ;
 	wnv1[i] = (double**)  calloc(Dim, sizeof( double *) ) ;
 	wnv2[i] = (double**)  calloc(Dim, sizeof( double *) ) ;

	dfhsdn3[i] = (double*)  calloc(M, sizeof( double ) ) ;

 	dfchdn3[i] = (double*)  calloc(M, sizeof( double ) ) ;
	
	dfdispndr[i] = (double*)  calloc(M, sizeof( double ) ) ;

	dfdispndr_pir[i] = (double**)  calloc(TNsp+1, sizeof( double* ) ) ;
	dfdispl1dn0[i] = (double*)  calloc(M, sizeof( double ) ) ;
	dfdispl1dn3[i] = (double*)  calloc(M, sizeof( double ) ) ;
	dfdispl2dn0[i] = (double*)  calloc(M, sizeof( double ) ) ;
	dfdispl2dn3[i] = (double*)  calloc(M, sizeof( double ) ) ;


	for(j=0;j<Dim ; j++){
		wnv1[i][j] =  (double*)  calloc(M, sizeof( double ) ) ;

		wnv2[i][j] =  (double*)  calloc(M, sizeof( double ) ) ;
	}//j == Dim


	for(j=0;j<TNsp+1 ; j++){
		 dfdispndr_pir[i][j] = (double*)  calloc(M, sizeof( double ) ) ;
	}

	//exttmpar1[i] = (double*)  calloc(M, sizeof( double ) ) ;


 }//i == Tnsp 

}
