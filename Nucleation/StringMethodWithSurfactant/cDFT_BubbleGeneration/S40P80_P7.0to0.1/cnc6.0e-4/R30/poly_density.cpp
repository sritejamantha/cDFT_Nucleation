#include "globals.h"
int stack(int*);
double integ_simpson(int , double* );
void unstack_ext2(int, int*  );
void unstack(int, int*  );

void poly_density()
{
    int adj_len,i,j,k,l,ind1,ind2,indz,nn[Dim];
    int tmpTNsp,tmpNs,tmpsigNbs,Nbp,tmpsigNhfbs;    
    double *q_tmp,k2, vect[Dim],cent_v,tmpN2sigbs,tmp_rsigab=0.0,tmp_sigab=0.0,mubp=0.0;
   
    
    tmpNs = NK[TNsp-1]/Nbr;

  //  tmp_rsigab = 2.0*(sigma[TNsp-1]*sigma[TNsp-2])/(sigma[TNsp-1]+sigma[TNsp-2]);
    tmp_rsigab = (rsigma[TNsp-1]+rsigma[TNsp-2])/2.0;
    tmp_sigab = tmp_rsigab*(1. - 0.12*exp(-3.0*crsepsln[crsTNsp-1]/Tmpr));

    //tmpsigNbs = round((sigma[TNsp-1]+sigma[TNsp-2])/dx[0]);
    //tmpN2sigbs = 2.0*round((sigma[TNsp-1]+sigma[TNsp-2])/dx[0])*dx[0];

   tmpsigNbs = round(tmp_sigab/dx[0]);
   tmpN2sigbs = 2.0*round(tmp_sigab/dx[0])*dx[0];
   tmpsigNhfbs = round(tmp_sigab/2./dx[0]);

//    printf("%d \t %lf\n",tmpsigNbs,tmpN2sigbs);
//    exit(-1);


    tmpTNsp = TNsp-2;   
    for (i=0;i<tmpTNsp;i++) /*Loop over Molecules other than Surfactant (Polyol here)*/
    {
#pragma omp parallel for private(k,nn)
	for(j=0;j<extM2;j++)
	{
		unstack_ext2(j,nn);
		if((nn[Dim-1] >= padN2 ) and  (nn[Dim-1]< (Nx[Dim-1]+padN2)) )
		{
			nn[Dim-1] -= padN2;
			k = stack(nn);
			extwk[i][j] = exp(amu[i]-dfttdrho[i][k]-wwall[i][k]);
			wk[i][k] = extwk[i][j] ;
		}
		else if (nn[Dim-1]>= (Nx[Dim-1] +padN2))
		{
				extwk[i][j] = exp(amu[i]-dfttdrho[i][M-1]-wwall[i][M-1]);
		}
		else
		{
				extwk[i][j] = exp(amu[i]-dfttdrho[i][0]-wwall[i][0]);
		}
	
	}//j extM2

       /*if(i==0)
       {
        for(j=0;j<M;j++)
        {
         printf("%d \t %d \t %.15e \t %.15e\n",i,j,dfttdrho[i][j],wk[i][j]);
        }
       }*/


    }//i tmpTNsp

    /*for(j=0;j<M;j++)
    {
      printf("%d \t %.15e \t %.15e\n",j,wk[0][j],wk[1][j]);
    }*/
 
   // exit(-1);
    
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
           //extwk[TNsp-2][j] = exp(amu[TNsp-2]+amu[TNsp-1]-dfttdrho[TNsp-2][k]-dfttdrho[TNsp-1][k]-wwall[TNsp-1][k]-wwall[TNsp-2][k]);
           extwk[TNsp-2][j] = exp(amu[TNsp-2]+amu[TNsp-1]-dfttdrho[TNsp-2][k]);
           extwk[TNsp-1][j] = exp(amu[TNsp-2]+amu[TNsp-1]-dfttdrho[TNsp-1][k]);
           wk[TNsp-2][k] = extwk[TNsp-2][j] ;
           wk[TNsp-1][k] = extwk[TNsp-1][j] ;
        }
        else if (nn[Dim-1]>= (Nx[Dim-1] +padN2))
        {
         //  extwk[TNsp-2][j] = exp(amu[TNsp-2]+amu[TNsp-1]-dfttdrho[TNsp-2][M-1]-dfttdrho[TNsp-1][M-1]-wwall[TNsp-2][M-1]-wwall[TNsp-1][M-1]);
           extwk[TNsp-2][j] = exp(amu[TNsp-2]+amu[TNsp-1]-dfttdrho[TNsp-2][M-1]);
           extwk[TNsp-1][j] = exp(amu[TNsp-2]+amu[TNsp-1]-dfttdrho[TNsp-1][M-1]);
        }
        else
        {
         //  extwk[TNsp-2][j] = exp(amu[TNsp-2]+amu[TNsp-1]-dfttdrho[TNsp-2][0]-dfttdrho[TNsp-1][0]-wwall[TNsp-2][0]-wwall[TNsp-1][0]);
          extwk[TNsp-2][j] = exp(amu[TNsp-2]+amu[TNsp-1]-dfttdrho[TNsp-2][0]);
          extwk[TNsp-1][j] = exp(amu[TNsp-2]+amu[TNsp-1]-dfttdrho[TNsp-1][0]);
        }
    }//j extM2  

    /*for(j=0;j<M;j++)
    {
     printf("%d \t %.15e \t %.15e\n",j, amu[TNsp-2]+amu[TNsp-1]-dfttdrho[TNsp-2][j]-dfttdrho[TNsp-1][j],wk[TNsp-2][j]);
    }*/

    /*for(j=0;j<M;j++)
    {
     printf("%d \t %.15e \t %.15e\n",j, dfttdrho[TNsp-2][j],dfttdrho[TNsp-1][j]);
    }*/

   // exit(-1);

    q_tmp = (double*) malloc(extM2*sizeof(double));
  

	/*calculate chain propagator*/
   
    int tmp_double; 

    
    for (i=0;i<tmpTNsp;i++)  /*Loop over molecules other than surfactant (Polyol here)*/
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
		    else if (nn[Dim-1]>= (Nx[Dim-1] +padN2)){
		    	q_tmp[j] = extwk[i][j]*q_pgt[l-1][M-1];
		    }
		    else
		    {
			q_tmp[j] = extwk[i][j]*q_pgt[l-1][0];
		    }

		    /*if(i==0)
                    {
                     printf("%d \t %d \t %d \t %.15e\n",i,l,j,q_tmp[j]);
                    }*/
			
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

		    
		    for(k=0;k<adj_len;k++){
			tmpN2[tid][k] = q_tmp[ind1+k]*extz[ind1+k];  
		    }//k
		    q_pgt[l][j] = dx[Dim-1]*integ_simpson(adj_len,tmpN2[tid])/N2sig[i]/extz[j+padN2];
                    /*if(i==0 && l<2)
                    {
                     printf("%d \t %d \t %d \t %.15e\n",i,l,j,q_pgt[l][j]);
                    }*/
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
        
	int lft_n; 
	//accum particle density 
#pragma omp parallel for private(l,lft_n,ind1,ind2)
	for(j=0;j<M;j++)
	{
		nwrhok[i][j] = 0;
		for(l=0;l<NK[i];l++)
		{
			
			ind1 = (l<MAXClen ? l : MAXClen-1);
			ind2 = ((NK[i]-l-1) < MAXClen ? (NK[i]-l-1) : MAXClen-1);
			
			nwrhok[i][j] += wk[i][j]*q_pgt[ind1][j]*q_pgt[ind2][j]; 	

		}//l NK[i]
//		if(i==1) 
//                printf("%d \t %d \t %.15e\n",i,j,nwrhok[i][j]);
	}//j M

    }//i tmpTNsp

    /*for(j=0;j<M;j++)
    {
     printf("%d \t %.15e \t %.15e\n",j,nwrhok[0][j],nwrhok[1][j]);
    }*/
   // exit(-1);
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


   //accum particle density 
//   Nbp = NK[TNsp-1]+NK[TNsp-2];
   for(j=0;j<M;j++)
   {
     nwrhok[TNsp-2][j] = 0.;
     nwrhok[TNsp-1][j] = 0.;
     for(l=0;l<NK[TNsp-2];l++)
     {
      nwrhok[TNsp-2][j] += wk[TNsp-2][j]*qs_fwd_pgt[l][j]*qs_bkd_pgt[Nbp-l-1][j];
     }//l NK[i]

     for(l=NK[TNsp-2];l<Nbp;l++)
     {
      nwrhok[TNsp-1][j] += wk[TNsp-1][j]*qs_fwd_pgt[l][j]*qs_bkd_pgt[Nbp-l-1][j];
     }//l NK[i]

  //   printf("%d \t %.15e \t %.15e \t %.15e\n",j,nwrhok[TNsp-2][j],nwrhok[TNsp-1][j],nwrhok[TNsp-2][j]+nwrhok[TNsp-1][j]);
   }//j M

  /*for(i=0;i<M;i++)
    {
     printf("%d \t %.15e \t %.15e\n",i,nwrhok[0][i],nwrhok[1][i]);
    }*/
   free(q_tmp);

  // exit(-1);
//   free(q_tmp);
}

 
