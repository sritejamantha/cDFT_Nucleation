#include "globals.h"
int stack(int*);
double integ_simpson(int , double* );
void unstack_ext2(int, int*  );

void poly_density()
{
    int i,j,k,l,fl,ind1,ind2,nn[Dim];
    int tmpTNsp,tmpNs,tmpsigNbs,Nbp;
    double *q_tmp,k2, vect[Dim],tmpN2sigbs,tmp_rsigab=0.0,tmp_sigab=0.0,mubp=0.0;
    double tmp_rsigbc=0.0,tmp_sigbc=0.0,tmpsigNbc,tmpN2sigbc;
//    tmpNs = NK[TNsp-1]/Nbr;

    tmp_rsigab = (rsigma[TNsp-2]+rsigma[TNsp-3])/2.0;
    tmp_sigab = tmp_rsigab*(1. - 0.12*exp(-3.0*crsepsln[crsTNsp-3]/Tmpr));

    tmp_rsigab = (rsigma[TNsp-2]+rsigma[TNsp-1])/2.0;
    tmp_sigbc = tmp_rsigab*(1. - 0.12*exp(-3.0*crsepsln[crsTNsp-1]/Tmpr));

    tmpsigNbs = round(tmp_sigab/dx[0]);
    tmpN2sigbs = 2.0*round(tmp_sigab/dx[0])*dx[0];

    tmpsigNbc = round(tmp_sigbc/dx[0]);
    tmpN2sigbc = 2.0*round(tmp_sigbc/dx[0])*dx[0];
//    printf("%d \t %lf\n",tmpsigNbs,tmpN2sigbs);
//    exit(-1);


    tmpTNsp = TNsp-3;
    for (i=0;i<tmpTNsp;i++) /*Loop over Molecules other than surfactant*/
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


    }//i tmpTNsp

    //exit(-1);

   Nbp = NK[TNsp-1]+NK[TNsp-2]+NK[TNsp-3];
   mubp = (double)Nbp*(amu[TNsp-3]+amu[TNsp-2]+amu[TNsp-1]);

  /*For surfactant*/
    #pragma omp parallel for private(k,nn)
    for(j=0;j<extM2;j++)
    {
        unstack_ext2(j,nn);
        if((nn[Dim-1] >= padN2 ) and  (nn[Dim-1]< (Nx[Dim-1]+padN2)) )
        {
           nn[Dim-1] -= padN2;
           k = stack(nn);
           extwk[TNsp-3][j] = exp(amu[TNsp-3]+amu[TNsp-2]+amu[TNsp-1]-dfttdrho[TNsp-3][k]); 
           extwk[TNsp-2][j] = exp(amu[TNsp-3]+amu[TNsp-2]+amu[TNsp-1]-dfttdrho[TNsp-2][k]);
           extwk[TNsp-1][j] = exp(amu[TNsp-3]+amu[TNsp-2]+amu[TNsp-1]-dfttdrho[TNsp-1][k]);
           wk[TNsp-3][k] = extwk[TNsp-3][j] ;
           wk[TNsp-2][k] = extwk[TNsp-2][j] ;
           wk[TNsp-1][k] = extwk[TNsp-1][j] ;
        }
        else if (nn[Dim-1]>= (Nx[Dim-1] +padN2))
        {
           extwk[TNsp-3][j] = exp(amu[TNsp-3]+amu[TNsp-2]+amu[TNsp-1]-dfttdrho[TNsp-3][M-1]);
           extwk[TNsp-2][j] = exp(amu[TNsp-3]+amu[TNsp-2]+amu[TNsp-1]-dfttdrho[TNsp-2][M-1]);
           extwk[TNsp-1][j] = exp(amu[TNsp-3]+amu[TNsp-2]+amu[TNsp-1]-dfttdrho[TNsp-1][M-1]);
        }
        else
        {
          extwk[TNsp-3][j] = exp(amu[TNsp-3]+amu[TNsp-2]+amu[TNsp-1]-dfttdrho[TNsp-3][0]);
          extwk[TNsp-2][j] = exp(amu[TNsp-3]+amu[TNsp-2]+amu[TNsp-1]-dfttdrho[TNsp-2][0]);
          extwk[TNsp-1][j] = exp(amu[TNsp-3]+amu[TNsp-2]+amu[TNsp-1]-dfttdrho[TNsp-1][0]);
        }
    }//j extM2
  
    q_tmp = (double*) malloc(extM2*sizeof(double));


	/*calculate chain propagator*/
   
    int tmp_double; 

    for (i=0;i<tmpTNsp;i++) /*Loop over molecules other than surfactant*/
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
			q_tmp[j] = extwk[i][j]*q_pgt[l-1][0];
		    }

		}//j extM2

		//convolv q+1 with delta function
#pragma omp parallel for private(k,ind1,ind2)
		for(j=0;j<M;j++)
                {
		    int tid = omp_get_thread_num() ;
		    ind1 = j+padN2 - sigN[i];
		    ind2 = j+padN2 + sigN[i];
		    for(k=0;k<(2*sigN[i]+1);k++){
			tmpN2[tid][k] = q_tmp[ind1+k];  
		    }//k
		    q_pgt[l][j] = dx[Dim-1]*integ_simpson(2*sigN[i]+1,tmpN2[tid])/N2sig[i];
		}//j  M
		
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
//                if(i==0) 
//                printf("%d \t %d \t %.15e\n",i,j,nwrhok[i][j]);
        }//j M

        //printf("\n");
    }//i tmpTNsp

//    exit(-1);
   /*For surfactant*/
    
   for(l=1;l<NK[TNsp-3];l++)
   {
    for(j=0;j<extM2;j++)
    {
     unstack_ext2(j,nn);
     if((nn[Dim-1] >= padN2 ) and  (nn[Dim-1]< (Nx[Dim-1]+padN2)) )
     {
      nn[Dim-1] -= padN2;
      k = stack(nn);
      q_tmp[j] = extwk[TNsp-3][j]*qs_fwd_pgt[l-1][k];
     }
     else if (nn[Dim-1]>= (Nx[Dim-1] +padN2))
     {
      q_tmp[j] = extwk[TNsp-3][j]*qs_fwd_pgt[l-1][M-1];
     }
     else
     {
      q_tmp[j] = extwk[TNsp-3][j]*qs_fwd_pgt[l-1][0];
     }

    }//j extM2

    //convolv q+1 with delta function
#pragma omp parallel for private(k,ind1,ind2)
    for(j=0;j<M;j++)
    {
     int tid = omp_get_thread_num() ;
     ind1 = j+padN2 - sigN[TNsp-3];
     ind2 = j+padN2 + sigN[TNsp-3];
     for(k=0;k<(2*sigN[TNsp-3]+1);k++)
     {
       tmpN2[tid][k] = q_tmp[ind1+k];
     }//k
     qs_fwd_pgt[l][j] = dx[Dim-1]*integ_simpson(2*sigN[TNsp-3]+1,tmpN2[tid])/N2sig[TNsp-3];
 
    }//j  M
 
   }// l NK


   l=NK[TNsp-3];
   for(j=0;j<extM2;j++)
   {
    unstack_ext2(j,nn);
    if((nn[Dim-1] >= padN2 ) and  (nn[Dim-1]< (Nx[Dim-1]+padN2)) )
    {
     nn[Dim-1] -= padN2;
     k = stack(nn);
     q_tmp[j] = extwk[TNsp-3][j]*qs_fwd_pgt[l-1][k];
    }
    else if (nn[Dim-1]>= (Nx[Dim-1] +padN2))
    {
     q_tmp[j] = extwk[TNsp-3][j]*qs_fwd_pgt[l-1][M-1];
    }
    else
    {
     q_tmp[j] = extwk[TNsp-3][j]*qs_fwd_pgt[l-1][0];
    }

   }//j extM2

   //convolv q+1 with delta function
#pragma omp parallel for private(k,ind1,ind2)
   for(j=0;j<M;j++)
   {
    int tid = omp_get_thread_num() ;
    ind1 = j+padN2 - tmpsigNbs;
    ind2 = j+padN2 + tmpsigNbs;
    for(k=0;k<(2*tmpsigNbs+1);k++)
    {
      tmpN2[tid][k] = q_tmp[ind1+k];
    }//k
    qs_fwd_pgt[l][j] = dx[Dim-1]*integ_simpson(2*tmpsigNbs+1,tmpN2[tid])/tmpN2sigbs;

   }//j  M

   

   for(l=NK[TNsp-3]+1;l<NK[TNsp-3]+NK[TNsp-2];l++)
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
     int tid = omp_get_thread_num() ;
     ind1 = j+padN2 - sigN[TNsp-2];
     ind2 = j+padN2 + sigN[TNsp-2];
     for(k=0;k<(2*sigN[TNsp-2]+1);k++)
     {
       tmpN2[tid][k] = q_tmp[ind1+k];
     }//k
     qs_fwd_pgt[l][j] = dx[Dim-1]*integ_simpson(2*sigN[TNsp-2]+1,tmpN2[tid])/N2sig[TNsp-2];

    }//j  M

   }// l NK


   l=NK[TNsp-3]+NK[TNsp-2];
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
    int tid = omp_get_thread_num() ;
    ind1 = j+padN2 - tmpsigNbc;
    ind2 = j+padN2 + tmpsigNbc;
    for(k=0;k<(2*tmpsigNbc+1);k++)
    {
      tmpN2[tid][k] = q_tmp[ind1+k];
    }//k
    qs_fwd_pgt[l][j] = dx[Dim-1]*integ_simpson(2*tmpsigNbc+1,tmpN2[tid])/tmpN2sigbc;

   }//j  M

   for(l=NK[TNsp-3]+NK[TNsp-2]+1;l<NK[TNsp-3]+NK[TNsp-2]+NK[TNsp-1];l++)
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
     int tid = omp_get_thread_num() ;
     ind1 = j+padN2 - sigN[TNsp-1];
     ind2 = j+padN2 + sigN[TNsp-1];
     for(k=0;k<(2*sigN[TNsp-1]+1);k++)
     {
       tmpN2[tid][k] = q_tmp[ind1+k];
     }//k
     qs_fwd_pgt[l][j] = dx[Dim-1]*integ_simpson(2*sigN[TNsp-1]+1,tmpN2[tid])/N2sig[TNsp-1];

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
     int tid = omp_get_thread_num() ;
     ind1 = j+padN2 - sigN[TNsp-1];
     ind2 = j+padN2 + sigN[TNsp-1];
     for(k=0;k<(2*sigN[TNsp-1]+1);k++)
     {
       tmpN2[tid][k] = q_tmp[ind1+k];
     }//k
     qs_bkd_pgt[l][j] = dx[Dim-1]*integ_simpson(2*sigN[TNsp-1]+1,tmpN2[tid])/N2sig[TNsp-1];

    }//j  M

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
    int tid = omp_get_thread_num() ;
    ind1 = j+padN2 - tmpsigNbc;
    ind2 = j+padN2 + tmpsigNbc;
    for(k=0;k<(2*tmpsigNbc+1);k++)
    {
      tmpN2[tid][k] = q_tmp[ind1+k];
    }//k
    qs_bkd_pgt[l][j] = dx[Dim-1]*integ_simpson(2*tmpsigNbc+1,tmpN2[tid])/tmpN2sigbc;

   }//j  M

   
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
     int tid = omp_get_thread_num() ;
     ind1 = j+padN2 - sigN[TNsp-2];
     ind2 = j+padN2 + sigN[TNsp-2];
     for(k=0;k<(2*sigN[TNsp-2]+1);k++)
     {
       tmpN2[tid][k] = q_tmp[ind1+k];
     }//k
     qs_bkd_pgt[l][j] = dx[Dim-1]*integ_simpson(2*sigN[TNsp-2]+1,tmpN2[tid])/N2sig[TNsp-2];

    }//j  M

   }// l NK

   l=NK[TNsp-1]+NK[TNsp-2];
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
    int tid = omp_get_thread_num() ;
    ind1 = j+padN2 - tmpsigNbs;
    ind2 = j+padN2 + tmpsigNbs;
    for(k=0;k<(2*tmpsigNbs+1);k++)
    {
      tmpN2[tid][k] = q_tmp[ind1+k];
    }//k
    qs_bkd_pgt[l][j] = dx[Dim-1]*integ_simpson(2*tmpsigNbs+1,tmpN2[tid])/tmpN2sigbs;

   }//j  M

   for(l=NK[TNsp-1]+NK[TNsp-2]+1;l<NK[TNsp-1]+NK[TNsp-2]+NK[TNsp-3];l++)
   {
    for(j=0;j<extM2;j++)
    {
     unstack_ext2(j,nn);
     if((nn[Dim-1] >= padN2 ) and  (nn[Dim-1]< (Nx[Dim-1]+padN2)) )
     {
      nn[Dim-1] -= padN2;
      k = stack(nn);
      q_tmp[j] = extwk[TNsp-3][j]*qs_bkd_pgt[l-1][k];
     }
     else if (nn[Dim-1]>= (Nx[Dim-1] +padN2))
     {
      q_tmp[j] = extwk[TNsp-3][j]*qs_bkd_pgt[l-1][M-1];
     }
     else
     {
      q_tmp[j] = extwk[TNsp-3][j]*qs_bkd_pgt[l-1][0];
     }

    }//j extM2

    //convolv q+1 with delta function
#pragma omp parallel for private(k,ind1,ind2)
    for(j=0;j<M;j++)
    {
     int tid = omp_get_thread_num() ;
     ind1 = j+padN2 - sigN[TNsp-3];
     ind2 = j+padN2 + sigN[TNsp-3];
     for(k=0;k<(2*sigN[TNsp-3]+1);k++)
     {
       tmpN2[tid][k] = q_tmp[ind1+k];
     }//k
     qs_bkd_pgt[l][j] = dx[Dim-1]*integ_simpson(2*sigN[TNsp-3]+1,tmpN2[tid])/N2sig[TNsp-3];

    }//j  M

   }// l NK
   
   //accum particle density 
   for(j=0;j<M;j++)
   {
     nwrhok[TNsp-3][j] = 0.;
     nwrhok[TNsp-2][j] = 0.;
     nwrhok[TNsp-1][j] = 0.;
     for(l=0;l<NK[TNsp-3];l++)
     {
      nwrhok[TNsp-3][j] += wk[TNsp-3][j]*qs_fwd_pgt[l][j]*qs_bkd_pgt[Nbp-l-1][j];
     }//l NK[i]

     for(l=NK[TNsp-3];l<NK[TNsp-3]+NK[TNsp-2];l++)
     {
      nwrhok[TNsp-2][j] += wk[TNsp-2][j]*qs_fwd_pgt[l][j]*qs_bkd_pgt[Nbp-l-1][j];
     }//l NK[i]

     for(l=NK[TNsp-3]+NK[TNsp-2];l<Nbp;l++)
     {
      nwrhok[TNsp-1][j] += wk[TNsp-1][j]*qs_fwd_pgt[l][j]*qs_bkd_pgt[Nbp-l-1][j];
     }//l NK[i]

//     printf("%d \t %.15e \t %.15e \t %.15e \t %.15e\n",j,nwrhok[TNsp-3][j],nwrhok[TNsp-2][j],nwrhok[TNsp-1][j],nwrhok[TNsp-3][j]+nwrhok[TNsp-2][j]+nwrhok[TNsp-1][j]);
   }//j M



//   exit(-1);
	
   free(q_tmp);
}

 
