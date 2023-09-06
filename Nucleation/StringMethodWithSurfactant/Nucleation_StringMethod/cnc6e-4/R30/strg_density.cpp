#include "globals.h"
int stack(int*);
double integ_simpson(int , double* );
void unstack_ext2(int, int*  );
void unstack(int, int*  );

double strg_density(int sti, double final_rmax)
{
    int adj_len,i,j,k,l,ind1,ind2,indz,nn[Dim];
    int tmpTNsp,tmpNs,tmpsigNbs,Nbp,tmpsigNhfbs;    
    double *q_tmp,k2, vect[Dim],cent_v,tmpN2sigbs,tmp_rsigab=0.0,tmp_sigab=0.0,mubp=0.0;
   
    tmpNs = NK[TNsp-1]/Nbr;

    tmp_rsigab = (rsigma[TNsp-1]+rsigma[TNsp-2])/2.0;
    tmp_sigab = tmp_rsigab*(1. - 0.12*exp(-3.0*crsepsln[crsTNsp-1]/Tmpr));

    tmpsigNbs = round(tmp_sigab/dx[0]);
    tmpN2sigbs = 2.0*round(tmp_sigab/dx[0])*dx[0];
    tmpsigNhfbs = round(tmp_sigab/2./dx[0]);

    tmpTNsp = TNsp-2;
    for (i=0;i<tmpTNsp;i++) /*Loop over Molecules other than the surfactant (here Polyol)*/
    {
#pragma omp parallel for private(k,nn)
	for(j=0;j<extM2;j++)
	{
		unstack_ext2(j,nn);
		if((nn[Dim-1] >= padN2 ) and  (nn[Dim-1]< (Nx[Dim-1]+padN2)) )
		{
			nn[Dim-1] -= padN2;
			k = stack(nn);
			extwk[i][j] = exp(amu[i]-strg_dfttdrho[i][sti][k] );
			wk[i][k] = extwk[i][j] ;
		}
		else if (nn[Dim-1]>= (Nx[Dim-1] +padN2))
		{
			extwk[i][j] = exp(amu[i]-strg_dfttdrho[i][sti][M-1] );
		}
		else
		{
			extwk[i][j] = exp(amu[i]-strg_dfttdrho[i][sti][0] );
		}
	
	}//j extM2
      
    }//i tmpTNsp

    /*if(sti==20)
    {
     for(j=0;j<extM2;j++)
     {
      printf("%d \t %.15e \t %.15e\n",j,extwk[0][j],extwk[1][j]);
     }
     exit(-1);
    }*/


    Nbp = NK[TNsp-1]+NK[TNsp-2];
    mubp = (double)Nbp*(amu[TNsp-2]+amu[TNsp-1]);
  /*For surfactant*/
    #pragma omp parallel for private(k,nn)
    for(j=0;j<extM2;j++)
    {
        unstack_ext2(j,nn);
        if((nn[Dim-1] >= padN2 ) and  (nn[Dim-1]< (Nx[Dim-1]+padN2)) )
        {
           nn[Dim-1] -= padN2;
           k = stack(nn);
           
           extwk[TNsp-2][j] = exp(amu[TNsp-2]+amu[TNsp-1]-strg_dfttdrho[TNsp-2][sti][k]);
           extwk[TNsp-1][j] = exp(amu[TNsp-2]+amu[TNsp-1]-strg_dfttdrho[TNsp-1][sti][k]);
           wk[TNsp-2][k] = extwk[TNsp-2][j] ;
           wk[TNsp-1][k] = extwk[TNsp-1][j] ;
        }
        else if (nn[Dim-1]>= (Nx[Dim-1] +padN2))
        {
           extwk[TNsp-2][j] = exp(amu[TNsp-2]+amu[TNsp-1]-strg_dfttdrho[TNsp-2][sti][M-1]);
           extwk[TNsp-1][j] = exp(amu[TNsp-2]+amu[TNsp-1]-strg_dfttdrho[TNsp-1][sti][M-1]);
        }
        else
        {
          extwk[TNsp-2][j] = exp(amu[TNsp-2]+amu[TNsp-1]-strg_dfttdrho[TNsp-2][sti][0]);
          extwk[TNsp-1][j] = exp(amu[TNsp-2]+amu[TNsp-1]-strg_dfttdrho[TNsp-1][sti][0]);
        }
    }//j extM2   


    /*if(sti==20)
    {
     for(j=0;j<extM2;j++)
     {
      printf("%d \t %.15e \t %.15e\n",j,extwk[2][j],extwk[3][j]);
     }
     exit(-1);
    }*/
      
    q_tmp = (double*) malloc(extM2*sizeof(double));
  

	/*calculate chain propagator*/
   
    int tmp_double; 

    for(i=0;i<tmpTNsp;i++)    /*Loop over molecules other than the surfactant (Polyol and CO2 here)*/
    {
    	for(l=1;l<NK[i];l++)
        {		
		if( l >= MAXClen )
                {

		    /*for(j=0;j<M;j++){
			  q_pgt[l][j]  = q_pgt[l-1][j] ;
		    }*/
		    continue;
		}

		// calc q+1 = q* wk
		for(j=0;j<extM2;j++)
                {
		    unstack_ext2(j,nn);
		    if((nn[Dim-1] >= padN2 ) and  (nn[Dim-1]< (Nx[Dim-1]+padN2)) )
                    {
		    	 nn[Dim-1] -= padN2;
			 k = stack(nn);
			 q_tmp[j] = extwk[i][j]*q_pgt[l-1][k];
		    }
		    else if (nn[Dim-1]>= (Nx[Dim-1] +padN2))
                    {
		    	q_tmp[j] = extwk[i][j]*q_pgt[l-1][M-1];
		    }
		    else
                    {
			//q_tmp[j] = extwk[i][j]*q_pgt[l-1][0];
		    }
			
		}//j extM2

		//convolv q+1 with delta function
#pragma omp parallel for private(nn,k,ind1,ind2,adj_len)
		for(j=0;j<M;j++)
		{
		    unstack(j,nn);

		    if(nn[Dim-1] < sigNhf[i])
		    	continue;



		    int tid = omp_get_thread_num() ;
		    if(j < sigN[i])
		    {
		    	ind1 =  padN2 + sigN[i] - j;
		    	ind2 = j+padN2 + sigN[i];
		    	adj_len = 2*j +1;
		    }
		    else
		    {
		    	ind1 = j+padN2 - sigN[i];
		    	ind2 = j+padN2 + sigN[i];
		    	adj_len = 2*sigN[i]+1;
		    }

		    
		    for(k=0;k<adj_len;k++)
                    {
			tmpN2[tid][k] = q_tmp[ind1+k]*extz[ind1+k];  
		    }//k
		    q_pgt[l][j] = dx[Dim-1]*integ_simpson(adj_len,tmpN2[tid])/N2sig[i]/extz[j+padN2];
		}//j  M

		cent_v = q_tmp[padN2 +sigN[i]] ;

		interp_values(q_pgt[l],cent_v,sigNhf[i]); 

		/*if(i==0)
                {
                  printf("\n");
                  for(j=0;j<M;j++)
                  {
                   printf("%d \t %d \t %d \t %.15e\n",i,l,j,q_pgt[l][j]);
                  }
                }*/	

	}//l NK

        /*if(sti==20 && i==0)
        {
          printf("\n");
          for(j=0;j<M;j++)
          {
            printf("%d \t %.15e\n",j,q_pgt[60][j]);
          }
          exit(-1);
         }*/


	int lft_n;
	double tmp_rmax = final_rmax; 

	//accum particle density 
//#pragma omp parallel for reduction(max:tmp_rmax)
	for(j=0;j<M;j++)
	{
		//strg_newrhoK[i][sti][j]  = 0;
		double	tmp_rhoh = 0;; 
		int lt ,ind1t, ind2t; 
		
		for( lt=0;lt<NK[i];lt++)
		{
			
			ind1t = (lt<MAXClen ? lt : MAXClen-1);
			ind2t = ((NK[i]-lt-1) < MAXClen ? (NK[i]-lt-1) : MAXClen-1);
			
			 tmp_rhoh += wk[i][j]*q_pgt[ind1t][j]*q_pgt[ind2t][j]; 	

		}//l NK[i]
	

		if( (tmp_rhoh > (0.01/lamb))  and (tmp_rhoh > tmp_rmax))
                {
		     tmp_rmax  =  tmp_rhoh ;  
//                     printf("%d \t %d \t %d \t %.14f\n",sti,i,j,tmp_rhoh);
                }
		strg_newrhoK[i][sti][j] = tmp_rhoh;
	}//j M

	final_rmax = tmp_rmax;

   }//i tmpTNsp

   /*if(sti==20)
   {
    for(j=0;j<M;j++)
    {
     printf("%d \t %.15e \t %.15e\n",j,strg_newrhoK[0][sti][j],strg_newrhoK[1][sti][j]);
    }
    exit(-1);
   }*/

    /*For surfactant*/

   for(l=1;l<NK[TNsp-2];l++)
   {
    for(j=0;j<extM2;j++)
    {
     unstack_ext2(j,nn);
     if((nn[Dim-1] >= padN2 ) and  (nn[Dim-1]< (Nx[Dim-1]+padN2)) )
     {
      nn[Dim-1] -= padN2;
      k = stack(nn);
      q_tmp[j] = extwk[TNsp-2][j]*qs_fwd_pgt[l-1][k];
     }
     else if (nn[Dim-1]>= (Nx[Dim-1] +padN2))
     {
      q_tmp[j] = extwk[TNsp-2][j]*qs_fwd_pgt[l-1][M-1];
     }
     else
     {
      q_tmp[j] = extwk[TNsp-2][j]*qs_fwd_pgt[l-1][0];
     }

    }//j extM2

    //convolv q+1 with delta function
#pragma omp parallel for private(k,ind1,ind2)
    for(j=0;j<M;j++)
    {
     unstack(j,nn);

     if(nn[Dim-1] < sigNhf[TNsp-2])
        continue;

     int tid = omp_get_thread_num() ;
     if(j < sigN[TNsp-2])
     {
       ind1 =  padN2 + sigN[TNsp-2] - j;
       ind2 = j+padN2 + sigN[TNsp-2];
       adj_len = 2*j +1;
     }
     else
     {
       ind1 = j+padN2 - sigN[TNsp-2];
       ind2 = j+padN2 + sigN[TNsp-2];
       adj_len = 2*sigN[TNsp-2]+1;
     }

     for(k=0;k<adj_len;k++)
     {
       tmpN2[tid][k] = q_tmp[ind1+k]*extz[ind1+k];
     }//k
     qs_fwd_pgt[l][j] = dx[Dim-1]*integ_simpson(adj_len,tmpN2[tid])/N2sig[TNsp-2]/extz[j+padN2];

    }//j  M

    cent_v = q_tmp[padN2 +sigN[TNsp-2]] ;

    interp_values(qs_fwd_pgt[l],cent_v,sigNhf[TNsp-2]);

   }// l NK


   l=NK[TNsp-2];
   for(j=0;j<extM2;j++)
   {
    unstack_ext2(j,nn);
    if((nn[Dim-1] >= padN2 ) and  (nn[Dim-1]< (Nx[Dim-1]+padN2)) )
    {
     nn[Dim-1] -= padN2;
     k = stack(nn);
     q_tmp[j] = extwk[TNsp-2][j]*qs_fwd_pgt[l-1][k];
    }
    else if (nn[Dim-1]>= (Nx[Dim-1] +padN2))
    {
     q_tmp[j] = extwk[TNsp-2][j]*qs_fwd_pgt[l-1][M-1];
    }
    else
    {
     q_tmp[j] = extwk[TNsp-2][j]*qs_fwd_pgt[l-1][0];
    }

   }//j extM2

   //convolv q+1 with delta function
#pragma omp parallel for private(k,ind1,ind2)
   for(j=0;j<M;j++)
   {

     unstack(j,nn);

     if(nn[Dim-1] < tmpsigNhfbs)
       continue;

    int tid = omp_get_thread_num() ;
    if(j < tmpsigNbs)
    {
      ind1 =  padN2 + tmpsigNbs - j;
      ind2 = j+padN2 + tmpsigNbs;
      adj_len = 2*j +1;
    }
    else
    {
      ind1 = j+padN2 - tmpsigNbs;
      ind2 = j+padN2 + tmpsigNbs;
      adj_len = 2*tmpsigNbs+1;
    }
    for(k=0;k<adj_len;k++)
    {
      tmpN2[tid][k] = q_tmp[ind1+k]*extz[ind1+k];
    }//k
    qs_fwd_pgt[l][j] = dx[Dim-1]*integ_simpson(adj_len,tmpN2[tid])/tmpN2sigbs/extz[j+padN2];

   }//j  M

   cent_v = q_tmp[padN2 +tmpsigNbs] ;

   interp_values(qs_fwd_pgt[l],cent_v,tmpsigNhfbs);

   
   for(l=NK[TNsp-2]+1;l<NK[TNsp-2]+NK[TNsp-1];l++)
   {
    for(j=0;j<extM2;j++)
    {
     unstack_ext2(j,nn);
     if((nn[Dim-1] >= padN2 ) and  (nn[Dim-1]< (Nx[Dim-1]+padN2)) )
     {
      nn[Dim-1] -= padN2;
      k = stack(nn);
      q_tmp[j] = extwk[TNsp-1][j]*qs_fwd_pgt[l-1][k];
     }
     else if (nn[Dim-1]>= (Nx[Dim-1] +padN2))
     {
      q_tmp[j] = extwk[TNsp-1][j]*qs_fwd_pgt[l-1][M-1];
     }
     else
     {
      q_tmp[j] = extwk[TNsp-1][j]*qs_fwd_pgt[l-1][0];
     }

    }//j extM2

    //convolv q+1 with delta function
#pragma omp parallel for private(k,ind1,ind2)
    for(j=0;j<M;j++)
    {
     unstack(j,nn);

     if(nn[Dim-1] < sigNhf[TNsp-1])
       continue;

     int tid = omp_get_thread_num() ;
     if(j < sigN[TNsp-1])
     {
       ind1 =  padN2 + sigN[TNsp-1] - j;
       ind2 = j+padN2 + sigN[TNsp-1];
       adj_len = 2*j +1;
     }
     else
     {
       ind1 = j+padN2 - sigN[TNsp-1];
       ind2 = j+padN2 + sigN[TNsp-1];
       adj_len = 2*sigN[TNsp-1]+1;
     }

     for(k=0;k<adj_len;k++)
     {
       tmpN2[tid][k] = q_tmp[ind1+k]*extz[ind1+k];
     }//k
     qs_fwd_pgt[l][j] = dx[Dim-1]*integ_simpson(adj_len,tmpN2[tid])/N2sig[TNsp-1]/extz[j+padN2];

     cent_v = q_tmp[padN2 +sigN[TNsp-1]] ;

     interp_values(qs_fwd_pgt[l],cent_v,sigNhf[TNsp-1]);

    }//j  M

   }// l NK


   for(l=1;l<NK[TNsp-1];l++)
   {
    for(j=0;j<extM2;j++)
    {
     unstack_ext2(j,nn);
     if((nn[Dim-1] >= padN2 ) and  (nn[Dim-1]< (Nx[Dim-1]+padN2)) )
     {
      nn[Dim-1] -= padN2;
      k = stack(nn);
      q_tmp[j] = extwk[TNsp-1][j]*qs_bkd_pgt[l-1][k];
     }
     else if (nn[Dim-1]>= (Nx[Dim-1] +padN2))
     {
      q_tmp[j] = extwk[TNsp-1][j]*qs_bkd_pgt[l-1][M-1];
     }
     else
     {
      q_tmp[j] = extwk[TNsp-1][j]*qs_bkd_pgt[l-1][0];
     }

    }//j extM2

    //convolv q+1 with delta function
#pragma omp parallel for private(k,ind1,ind2)
    for(j=0;j<M;j++)
    {
     unstack(j,nn);

     if(nn[Dim-1] < sigNhf[TNsp-1])
        continue;
     int tid = omp_get_thread_num() ;
     if(j < sigN[TNsp-1])
     {
       ind1 =  padN2 + sigN[TNsp-1] - j;
       ind2 = j+padN2 + sigN[TNsp-1];
       adj_len = 2*j +1;
     }
     else
     {
       ind1 = j+padN2 - sigN[TNsp-1];
       ind2 = j+padN2 + sigN[TNsp-1];
       adj_len = 2*sigN[TNsp-1]+1;
     }
     for(k=0;k<adj_len;k++)
     {
       tmpN2[tid][k] = q_tmp[ind1+k]*extz[ind1+k];
     }//k
     qs_bkd_pgt[l][j] = dx[Dim-1]*integ_simpson(adj_len,tmpN2[tid])/N2sig[TNsp-1]/extz[j+padN2];

    }//j  M

    cent_v = q_tmp[padN2 +sigN[TNsp-1]] ;

    interp_values(qs_bkd_pgt[l],cent_v,sigNhf[TNsp-1]);

   }// l NK


   l=NK[TNsp-1];
   for(j=0;j<extM2;j++)
   {
    unstack_ext2(j,nn);
    if((nn[Dim-1] >= padN2 ) and  (nn[Dim-1]< (Nx[Dim-1]+padN2)) )
    {
     nn[Dim-1] -= padN2;
     k = stack(nn);
     q_tmp[j] = extwk[TNsp-1][j]*qs_bkd_pgt[l-1][k];
    }
    else if (nn[Dim-1]>= (Nx[Dim-1] +padN2))
    {
     q_tmp[j] = extwk[TNsp-1][j]*qs_bkd_pgt[l-1][M-1];
    }
    else
    {
     q_tmp[j] = extwk[TNsp-1][j]*qs_bkd_pgt[l-1][0];
    }

   }//j extM2

   //convolv q+1 with delta function
#pragma omp parallel for private(k,ind1,ind2)
   for(j=0;j<M;j++)
   {
    unstack(j,nn);

    if(nn[Dim-1] < tmpsigNhfbs)
      continue;

    int tid = omp_get_thread_num() ;
    if(j < tmpsigNbs)
    {
      ind1 =  padN2 + tmpsigNbs - j;
      ind2 = j+padN2 + tmpsigNbs;
      adj_len = 2*j +1;
    }
    else
    {
      ind1 = j+padN2 - tmpsigNbs;
      ind2 = j+padN2 + tmpsigNbs;
      adj_len = 2*tmpsigNbs+1;
    }
    for(k=0;k<adj_len;k++)
    {
      tmpN2[tid][k] = q_tmp[ind1+k]*extz[ind1+k];
    }//k
    qs_bkd_pgt[l][j] = dx[Dim-1]*integ_simpson(adj_len,tmpN2[tid])/tmpN2sigbs/extz[j+padN2];

   }//j  M

   cent_v = q_tmp[padN2 +tmpsigNbs] ;

   interp_values(qs_bkd_pgt[l],cent_v,tmpsigNhfbs);

   for(l=NK[TNsp-1]+1;l<NK[TNsp-1]+NK[TNsp-2];l++)
   {
    for(j=0;j<extM2;j++)
    {
     unstack_ext2(j,nn);
     if((nn[Dim-1] >= padN2 ) and  (nn[Dim-1]< (Nx[Dim-1]+padN2)) )
     {
      nn[Dim-1] -= padN2;
      k = stack(nn);
      q_tmp[j] = extwk[TNsp-2][j]*qs_bkd_pgt[l-1][k];
     }
     else if (nn[Dim-1]>= (Nx[Dim-1] +padN2))
     {
      q_tmp[j] = extwk[TNsp-2][j]*qs_bkd_pgt[l-1][M-1];
     }
     else
     {
      q_tmp[j] = extwk[TNsp-2][j]*qs_bkd_pgt[l-1][0];
     }

    }//j extM2

    //convolv q+1 with delta function
#pragma omp parallel for private(k,ind1,ind2)
    for(j=0;j<M;j++)
    {
     unstack(j,nn);

     if(nn[Dim-1] < sigNhf[TNsp-2])
       continue;
     int tid = omp_get_thread_num() ;
     if(j < sigN[TNsp-2])
     {
       ind1 =  padN2 + sigN[TNsp-2] - j;
       ind2 = j+padN2 + sigN[TNsp-2];
       adj_len = 2*j +1;
     }
     else
     {
       ind1 = j+padN2 - sigN[TNsp-2];
       ind2 = j+padN2 + sigN[TNsp-2];
       adj_len = 2*sigN[i]+1;
     }
     for(k=0;k<adj_len;k++)
     {
       tmpN2[tid][k] = q_tmp[ind1+k]*extz[ind1+k];
     }//k
     qs_bkd_pgt[l][j] = dx[Dim-1]*integ_simpson(adj_len,tmpN2[tid])/N2sig[TNsp-2]/extz[j+padN2];

    }//j  M

    cent_v = q_tmp[padN2 +sigN[TNsp-2]] ;

    interp_values(qs_bkd_pgt[l],cent_v,sigNhf[TNsp-2]);

   }// l NK


   double tmp_rmax = final_rmax;
   //accum particle density 
//#pragma omp parallel for reduction(max:tmp_rmax)
//   Nbp = NK[TNsp-1]+NK[TNsp-2];
   for(j=0;j<M;j++)
   {
    // strg_newrhok[TNsp-2][sti][j] = 0.;
    // strg_newrhok[TNsp-1][sti][j] = 0.;
     double tmp_rhoh1 = 0.0;
     double tmp_rhoh2 = 0.0;
     for(l=0;l<NK[TNsp-2];l++)
     {
      tmp_rhoh1 += wk[TNsp-2][j]*qs_fwd_pgt[l][j]*qs_bkd_pgt[Nbp-l-1][j];
     }//l NK[i]

     for(l=NK[TNsp-2];l<Nbp;l++)
     {
      tmp_rhoh2 += wk[TNsp-1][j]*qs_fwd_pgt[l][j]*qs_bkd_pgt[Nbp-l-1][j];
     }//l NK[i]

     if( (tmp_rhoh1 > (0.01/lamb))  and (tmp_rhoh1 > tmp_rmax))
     {
        tmp_rmax  =  tmp_rhoh1;
  //      printf("%d \t %d \t %.14f\n",sti,j,tmp_rhoh1);
     }
     if( (tmp_rhoh2 > (0.01/lamb))  and (tmp_rhoh2 > tmp_rmax))
     {
        tmp_rmax  =  tmp_rhoh2;
//        printf("%d \t %d \t %.14f\n",sti,j,tmp_rhoh1);
     }
     strg_newrhoK[TNsp-2][sti][j] = tmp_rhoh1;
     strg_newrhoK[TNsp-1][sti][j] = tmp_rhoh2;

     final_rmax = tmp_rmax;
//     printf("%d \t %.15e \t %.15e \t %.15e\n",j,strg_newrhok[TNsp-2][j],strg_newrhok[TNsp-1][j],strg_newrhok[TNsp-2][j]+strg_newrhok[TNsp-1][j]);
   }//j M


   /*if(sti==55)
   {
    for(j=0;j<M;j++)
    {
     printf("%d \t %.12e \t %.12e \t %.12e\n",j,strg_newrhoK[0][sti][j],strg_newrhoK[1][sti][j],strg_newrhoK[2][sti][j]+strg_newrhoK[3][sti][j]);
    }
    exit(-1);
   }*/

   /*if(final_rmax>0)
   {
    printf("%d \t %.14f\n",sti,final_rmax);
   } */
  
   for(i=0;i<TNsp;i++)
   {
      for(j=0;j<M; j++)
      {
	strg_newrhoK[i][sti][j] = (1-(final_rmax > 0 ?  0.01/final_rmax: lamb))*strg_rhoK[i][sti][j] + strg_newrhoK[i][sti][j]*( final_rmax > 0 ?  0.01/final_rmax : lamb );
//        strg_newrhoK[i][sti][j] = (1.0-lamb)*strg_rhoK[i][sti][j] + strg_newrhoK[i][sti][j]*( lamb );
      }//j M
   }//i TNsp

   /*if(sti==48)
   {
    for(j=0;j<M;j++)
    {
     printf("%d \t %.12e \t %.12e \t %.12e\n",j,strg_newrhoK[0][sti][j],strg_newrhoK[1][sti][j],strg_newrhoK[2][sti][j]+strg_newrhoK[3][sti][j]);
    }
    exit(-1);
   }*/



  
   free(q_tmp);
   return final_rmax;
}


