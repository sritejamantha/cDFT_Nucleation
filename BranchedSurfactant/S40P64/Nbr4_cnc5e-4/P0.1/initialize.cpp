#include "globals.h"
void read_resume_files() ;
void allocate( void ) ;
void read_input( void ) ;
void unstack_input(int , int* , int* ) ;
void read_matdisp(void);
void read_blkrho(double);
void gen_rnnum(int,int,int,int*);

void initialize(double inp_pres) {
  int i,j,k,nn[Dim] ;

  idum =   -long( time(0) ) ; // 9 ; //

  read_input() ;

  read_blkrho(inp_pres);


  read_matdisp();

 

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


  double tmp_dist;
  
  int mid_ind = int(M/2.);

    for(i=0;i<TNsp;i++){
	
	for(j=0;j<M;j++){


		unstack_input(j,nn,Nx);
		if(wall_para == 0){
			tmp_dist = nn[Dim-1]*dx[Dim-1] -  L[Dim-1]/2.0;//+ (i != 1 ?Nsig[i]/2.0 : 0.0);
			
			rhoK[i][j] = (tmp_dist > 0 ? blkrho[i][1]: blkrho[i][0]);
	
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
//	  cout<<"t1 "<<knlt1[i][j]<<" t2  "<<knltz1[i][j]<<endl;
	}
	

    }

    for(j=0;j<M;j++)
    {
     Xar[j] = rhoK[TNsp-2][j]/(rhoK[TNsp-2][j] + rhoK[TNsp-1][j]);
    }

    read_resume_files() ;

    //write_grid_data("rho2_test.dat", rhoK[2]);

    if(Dim ==1){
     for(i=0 ;i<extM2; i++){
	
	 extz[i] = -1.*dx[Dim-1]*padN2 + i*dx[Dim-1]; 

     }
     
    }


    int tmp_nl=1;
    int tmp_nu=NK[TNsp-2];
    int tmp_ct = Nbr;
    for(i=0;i<tmp_nu;i++)
    {
     brflag[i]=0;
    }
    tmp_nu=tmp_nu-2;
    srand(time(0));
    gen_rnnum(tmp_nl,tmp_nu,tmp_ct,brflag);
    for(i=0;i<tmp_nu+2;i++)
    {
     printf("i:%d \t flag:%d\n",i,brflag[i]);
    }
//    exit(-1);

    int tmp_b=NK[TNsp-2];
    int tmp_s=NK[TNsp-1]/Nbr;
    for(j=0;j<M;j++)
    {
     qp_bb_pgt[0][j]=1;
     if(brflag[tmp_b-1]==0)      /*No branch at the end backbone bead*/
     {
      qc_bb_pgt[tmp_b-1][j]=1;
     }
     //printf("%d \t %.15e \t %.15e\n",j,qp_bb_pgt[0][j],qc_bb_pgt[tmp_b-1][j]);
     for(k=0;k<Nbr;k++)
     {
      //qp_br_pgt[k][0][j]=1;
      qc_br_pgt[k][tmp_s-1][j]=1;
     }
    }

//   exit(-1);

   prev_dif = 1.0e-2;
}


void allocate(){


 int tmp_chl=NK[1],i,j,tmp_brl=Nbr,tmp_bbl=NK[TNsp-2]; 
  
 if (NK[1]> MAXClen)
 	tmp_chl = MAXClen;

 Xar =  (double*)  calloc(M, sizeof( double ) ) ;
 brflag = (int*) calloc(tmp_bbl,sizeof( int* ) ) ;
 q_pgt = (double**)  calloc(tmp_chl,sizeof( double* ) ) ;
 qp_bb_pgt = (double**)  calloc(tmp_chl,sizeof( double* ) ) ;
 qc_bb_pgt = (double**)  calloc(tmp_chl,sizeof( double* ) ) ;
 qc_bbr_pgt = (double**)  calloc(tmp_brl,sizeof( double* ) ) ;
 qp_br_pgt = (double***)  calloc(tmp_brl,sizeof( double** ) ) ;
 qc_br_pgt = (double***)  calloc(tmp_brl,sizeof( double** ) ) ;
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

int tmp_brsl = NK[TNsp-1]/Nbr;
for(i=0;i<tmp_brl;i++)
{
 qp_br_pgt[i] = (double**)  calloc(tmp_brsl, sizeof( double* ) ) ;
 qc_br_pgt[i] = (double**)  calloc(tmp_brsl, sizeof( double* ) ) ;
 qc_bbr_pgt[i] = (double*)  calloc(M, sizeof( double ) ) ; 
 for(j=0;j<tmp_brsl;j++)
 {
  qp_br_pgt[i][j] = (double*)  calloc(M, sizeof( double ) ) ;
  qc_br_pgt[i][j] = (double*)  calloc(M, sizeof( double ) ) ;
 }
}

for(i=0;i< tmp_chl;i++)
{
	q_pgt[i] =  (double*)  calloc(M, sizeof( double ) ) ;
        qp_bb_pgt[i] = (double*)  calloc(M, sizeof( double ) ) ;
        qc_bb_pgt[i] = (double*)  calloc(M, sizeof( double ) ) ;
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

	dfdispl1dn0[i] = (double*)  calloc(M, sizeof( double ) ) ;
	dfdispl1dn3[i] = (double*)  calloc(M, sizeof( double ) ) ;
	dfdispl2dn0[i] = (double*)  calloc(M, sizeof( double ) ) ;
	dfdispl2dn3[i] = (double*)  calloc(M, sizeof( double ) ) ;


	for(j=0;j<Dim ; j++){
		wnv1[i][j] =  (double*)  calloc(M, sizeof( double ) ) ;

		wnv2[i][j] =  (double*)  calloc(M, sizeof( double ) ) ;
	}//j == Dim


	//exttmpar1[i] = (double*)  calloc(M, sizeof( double ) ) ;


 }//i == Tnsp 

}

void gen_rnnum(int nl,int nu,int ct,int *flg)
{
  int i,j,fl,tmp,tmpnum,rnum[ct];

  for(i=0;i<ct;)
  {
   fl=1;
   tmpnum = rand()%(nu-nl+1)+nl;
   if(i==0)
   {
    rnum[i]=tmpnum;
    i=i+1;
   }
   else
   {
    for(j=0;j<i;j++)
    {
     if(rnum[j]==tmpnum)
     {
      fl=-1; 
     }
    }
    if(fl>0)
    {
     rnum[i]=tmpnum;
     i=i+1;
    }
   }   
  }

  /*for(i=0;i<ct;i++)
  {
   printf("rn-%d:%d\n",i,rnum[i]);
  }*/

  for(i=0;i<ct-1;i++)
  {
   for(j=0;j<ct-i-1;j++)
   {
    if(rnum[j]>rnum[j+1])
    {
     tmp = rnum[j];
     rnum[j] = rnum[j+1];
     rnum[j+1] = tmp;
    }
   }
  }

  /*for(i=0;i<ct;i++)
  {
   printf("rn-%d:%d\n",i,rnum[i]);
  }*/

  j=0;
  for(i=0;i<nu+2;i++)
  {
   if(i==rnum[j])
   {
    flg[i]=1;
    j=j+1;
   }
  }


  // For block copolymer
  /*for(i=0;i<nu+1;i++)
  {
   if(i==nu)
   {
    flg[i]=1;
   }
   else
   {
    flg[i]=0;
   }
  }*/

  /*for(i=0;i<nu+1;i++)
  {
   printf("i:%d \t flag:%d\n",i,flg[i]);
  }*/



}
