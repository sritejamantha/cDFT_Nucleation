#include "globals.h"
void prof_shift(int, int);
void read_rst() ;
double find_blkmu();
void strg_dfree();
void dfree();
void allocate( void ) ;
void read_termianl() ;
void read_input( void ) ;
void unstack_input(int , int* , int* ) ;
void read_matdisp(void);

void initialize() 
{
  int i,j,k,nn[Dim] ;

  idum =   -long( time(0) ) ; // 9 ; //

  read_input() ;

  read_matdisp();

  lamb_mu = -0.0;
 
  rf_len = rsigma[0];

  M = 1;

 
  for(i=0;i<Dim;i++)
  {
	L[i] /= rf_len;
	dx[i] /= rf_len;
	Nx[i] = floor(L[i]/dx[i])+1;
  	M *=  Nx[i];
  	Lh[i] = L[i] *0.5;
  	cut_r /= rf_len;
  }

  padN = 0; 
  padN2 = 0; 

  cut_ind = int(cut_r/dx[Dim-1]+0.5);

  cout<<"strg cut ind at "<<cut_ind<<endl;
  for(i=0;i<TNsp;i++)
  {
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

  for(i=0;i<Dim;i++)
  {
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


  read_termianl() ;

   // exit(1);
  if(rst_ind > 0)
  {

	read_rst();

  }
  else
  {

	for(i=0;i<TNsp;i++)
        {
		for(j=0;j<M;j++)
                {
			tmil_rho0[i][j] = blkrho[i][1];
			unstack_input(j,nn,Nx);
			for(k=0; k<N_strg; k++)
                        {
			 	if(k<(N_strg/3))
				strg_rhoK[i][k][j] = tmil_rho0[i][j]*(double(N_strg/3-k-1)/double(N_strg/3-1)) + tmil_rho1[i][j]*(double(k)/double(N_strg/3-1));
				else
				    strg_rhoK[i][k][j] = tmil_rho1[i][j];
			}
		}

  	}

        /*for(j=0;j<M;j++)
        {
         printf("%d \t %.15e \t %.15e \t %.15e\n",j,strg_rhoK[0][28][j],strg_rhoK[1][28][j],strg_rhoK[2][28][j]+strg_rhoK[3][28][j]);
        }

        exit(-1);*/

        
	for(k=1;k<N_strg-1;k++)
        {
		prof_shift(k,cut_ind);
	}

        /*for(j=0;j<M;j++)
        {
         printf("%d \t %.15e \t %.15e \t %.15e\n",j,strg_rhoK[0][18][j],strg_rhoK[1][18][j],strg_rhoK[2][18][j]+strg_rhoK[3][18][j]);
        }

        exit(-1);*/


}


  double tmp_dist;
  int mid_ind = int(2*M/4.);
  for(i=0;i<TNsp;i++)
    {
        for(j=0;j<M;j++)
        {

                tmil_rho0[i][j] = blkrho[i][1];

                unstack_input(j,nn,Nx);

                for(k=0; k<N_strg; k++)
                {

                 if(k>0 and k < (N_strg-1))
                        drhodstrg[i][k][j] = (strg_rhoK[i][k+1][j] - strg_rhoK[i][k][j]);
                 else
                        drhodstrg[i][k][j] = 0;

                }//k
                q_pgt[0][j] = 1;

                if((nn[Dim-1]<(sigN[i]-1)) and (wall_para > 0))
                        wwall[i][j] = 1.e10;
                else
                        wwall[i][j] = 0.0;
        }//j

        /*for(k=0; k < N_strg ; k++)
        {
                nmdrds[i][k] = integ_inpd(M,drhodstrg[i][k],drhodstrg[i][k])*dx[Dim-1];
                strg_s[k]= double(k)/double(N_strg-1);
        }*///k = N-trg


        for(j=0;j<(2*sigNhf[i]+1);j++)
        {

          knlt1[i][j] = -sigNhf[i]*dx[Dim-1]+ j*dx[Dim-1];

          knltz1[i][j] = Nsig[i]*Nsig[i]/4.0 -pow(knlt1[i][j],2);
        }
    

    }//i

//    double tmpdrhodstrg[1][60][750];
    for(k=0;k<N_strg;k++)
    {
     for(j=0;j<M;j++)
     {
     tmpdrhodstrg[0][k][j]=0.0;
     }
    }

    for(i=0;i<TNsp;i++)
    {
     for(k=0; k < N_strg ; k++)
     {
      if(i<TNsp-2)
      {
       nmdrds[i][k] = integ_inpd(M,drhodstrg[i][k],drhodstrg[i][k])*dx[Dim-1];
      }
      else if(i==TNsp-2)
      {
       for(j=0;j<M;j++)
       {
        tmpdrhodstrg[0][k][j] = drhodstrg[TNsp-2][k][j]+drhodstrg[TNsp-1][k][j];
       }
       nmdrds[TNsp-2][k] = integ_inpd(M,tmpdrhodstrg[0][k],tmpdrhodstrg[0][k])*dx[Dim-1];
      }
      else
      {
       nmdrds[TNsp-1][k]=nmdrds[TNsp-2][k];
      }
      strg_s[k]= double(k)/double(N_strg-1);
     }///k = N-trg
    }

    /*for(j=0;j<M;j++)
    {
     printf("%d \t %.16f \t %.16f \t %.16f\n",j,drhodstrg[0][18][j],drhodstrg[1][18][j],drhodstrg[2][18][j]+drhodstrg[3][18][j]);
    }

    exit(-1);*/

    /*for(k=0;k<N_strg;k++)
    {
     printf("%d \t %.16f \t %.16f \t %.16f \t %.16f\n",k,nmdrds[0][k],nmdrds[1][k],nmdrds[2][k],nmdrds[3][k]);
    }

    exit(-1);*/
   /*for(j=0;j<M;j++)
   {
     Xar[j] = rhoK[TNsp-2][j]/(rhoK[TNsp-2][j] + rhoK[TNsp-1][j]);
   }*/


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

    //exit(1);

    double blkpres,rdc;
    rdc = (6.022e23)*pow(rf_len,3)*(1e-30);
    blkpres = find_blkmu();
    
    printf("%lf\n",blkpres*8.314*Tmpr/rdc);
//    exit(-1);    
    strg_dfree();

    //exit(1);

/*    for(i=0; i<TNsp; i++)
    	for(k=0; k<N_strg; k++)
        {
		inpdmudrds[i][k] = integ_inpd(M,drhodstrg[i][k],strg_dwdrho[i][k])*dx[Dim-1];
		//cout<<"inpdmu "<<inpdmudrds[i][k]<<" nmds "<<nmdrds[i][k]<<endl;
	}
*/
     //exit(1);
}


void allocate()
{

 int tmp_chl,i,j; 
  
// if (NK[1]> MAXClen)
 	tmp_chl = MAXClen;

 Xar =  (double*)  calloc(M, sizeof( double ) ) ;
 strg_dwdrho =  (double***)  calloc(TNsp,sizeof( double** ) ) ;
 tmpstrg_dwdrho =  (double***)  calloc(TNsp,sizeof( double** ) ) ;
 strg_dfttdrho =  (double***)  calloc(TNsp,sizeof( double** ) ) ; 
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
 tmil_rho0 = (double**) calloc(TNsp,sizeof( double* ) ) ;
 tmil_rho1 =  (double**) calloc(TNsp,sizeof( double* ) ) ;
 strg_newrhoK = (double***) calloc(TNsp,sizeof( double** ) ) ;
 strg_rhoK =  (double***) calloc(TNsp,sizeof( double** ) ) ;
 drhodstrg = (double***) calloc(TNsp,sizeof( double** ) ) ;
 tmpdrhodstrg = (double***) calloc(1,sizeof( double** ) ) ;
 nmdrds = (double**) calloc(TNsp,sizeof( double* ) ) ;
 inpdmudrds = (double**) calloc(TNsp,sizeof( double* ) ) ;

 strg_s = (double*) calloc(N_strg,sizeof(double));
 strg_s_ind = (int**) calloc(N_strg,sizeof(int*));  
 strg_s_w = (double**) calloc(N_strg,sizeof(double*));

 mxc_v = (double*) calloc(N_strg,sizeof(double));

 for(i=0;i<N_strg;i++){
	strg_s_ind[i] = (int*) calloc(2,sizeof(int));
	strg_s_w[i] = (double*) calloc(2,sizeof(double));
 }


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

 strg_itgF = (double*)  calloc(N_strg, sizeof( double ));
 strg_ttF = (double**)  calloc(N_strg, sizeof( double* ));
 
 strg_lgrgk = (double**)  calloc(N_strg, sizeof( double* ));
  
for(i=0;i< N_strg;i++){
	strg_ttF[i] =  (double*)  calloc(M, sizeof( double ) ) ;
	strg_lgrgk[i] = (double*)  calloc(M, sizeof( double ) ) ;
}

for(i=0;i< tmp_chl;i++)
{
  q_pgt[i] =  (double*)  calloc(M, sizeof( double ) ) ;
  qs_fwd_pgt[i] = (double*)  calloc(M, sizeof( double ) ) ;
  qs_bkd_pgt[i] = (double*)  calloc(M, sizeof( double ) ) ;
}	
qtmp = (double*)  calloc(M, sizeof( double ) ) ;

for(i=0;i< TNsp;i++){
	wwall[i] = (double*)  calloc(M, sizeof( double ) ) ;
	wk[i] =  (double*)  calloc(M, sizeof( double ) ) ;
	nwrhok[i] = (double*)  calloc(M, sizeof( double ) ) ;
	dfttdrho[i] = (double*)  calloc(M, sizeof( double ) ) ;

	knlt1[i] = (double*)  calloc(2*sigNhf[i]+1, sizeof( double ) ) ;
 	knltz1[i] = (double*)  calloc(2*sigNhf[i]+1, sizeof( double ) ) ;
  	tmpar1 = (double*)  calloc(M, sizeof( double ) ) ;		/*SM: Why here? Does't have any i loop*/
	extrhok[i] = (double*)  calloc(extM, sizeof( double ) ) ;
	extrhok2[i] = (double*)  calloc(extM2, sizeof( double ) ) ;
	extwk[i] =  (double*)  calloc(extM2, sizeof( double ) ) ;

//	tmpN2[i] = (double* )calloc(padN2, sizeof( double ) ) ;
	//tmpNhf[i] = (double* )calloc(2*sigNhf[i]+1, sizeof( double ) ) ;	
	extdflocdn0[i] = (double*) calloc(extM,sizeof( double ) ) ; 
 	extdflocdn3[i] = (double*) calloc(extM,sizeof( double ) ) ;
 	dfloc[i] =  (double*) calloc(M,sizeof( double )) ;
 
 	tmil_rho0[i] = (double*) calloc(M,sizeof( double )) ;
	tmil_rho1[i] =(double*) calloc(M,sizeof( double )) ;

	strg_newrhoK[i] = (double**) calloc(N_strg,sizeof( double* )) ;
	
	strg_dfttdrho[i] = (double**) calloc(N_strg,sizeof( double* )) ;

	strg_rhoK[i] = (double**) calloc(N_strg,sizeof( double* )) ;

	drhodstrg[i] = (double**) calloc(N_strg,sizeof( double* )) ;

	nmdrds[i] = (double*) calloc(N_strg,sizeof( double )) ;
	inpdmudrds[i] = (double*) calloc(N_strg,sizeof( double )) ;

	strg_dwdrho[i] = (double**) calloc(N_strg,sizeof( double* )) ;
	
	
	for(j=0;j<N_strg; j++){

		strg_newrhoK[i][j] = (double*) calloc(M,sizeof( double )) ;
		strg_rhoK[i][j] = (double*) calloc(M,sizeof( double )) ;
		drhodstrg[i][j] = (double*) calloc(M,sizeof( double )) ;

		strg_dfttdrho[i][j] = (double*) calloc(M,sizeof( double )) ;
		strg_dwdrho[i][j] = (double*) calloc(M,sizeof( double )) ;
	}
 }

 tmpdrhodstrg[0] = (double**) calloc(N_strg,sizeof( double* )) ;
 tmpstrg_dwdrho[0] = (double**) calloc(N_strg,sizeof( double* )) ;
 for(j=0;j<N_strg;j++)
 {
  tmpdrhodstrg[0][j] = (double*) calloc(M,sizeof( double )) ;
  tmpstrg_dwdrho[0][j] = (double*) calloc(M,sizeof( double )) ;
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

void prof_shift(int str, int tt_sftind)
{

 int tmp_ttsftind = int(double(N_strg - str +1)/double(N_strg)*tt_sftind) ,i, j, k;
// printf("%d \t %d \t %d\n",str,tt_sftind,tmp_ttsftind);

  for(j=0;j<TNsp;j++)
  {
   	for(i=0;i<M;i++)
        {
		qtmp[i] = strg_rhoK[j][str][i];
   	}

	for(i=0;i<(M-tmp_ttsftind);i++)
		strg_rhoK[j][str][i] = qtmp[i+tmp_ttsftind];
  
  	for(i=(M-tmp_ttsftind);i<M;i++)
		strg_rhoK[j][str][i] = qtmp[M-1];
  }

}
