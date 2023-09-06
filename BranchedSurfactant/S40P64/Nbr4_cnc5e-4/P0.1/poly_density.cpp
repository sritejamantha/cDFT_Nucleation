#include "globals.h"
int stack(int*);
double integ_simpson(int , double* );
void unstack_ext2(int, int*  );

void poly_density()
{
    int i,j,k,l,fl,ind1,ind2,nn[Dim];
    int tmpTNsp,tmpNs,tmpsigNbs;
    double *q_tmp,k2, vect[Dim],tmpN2sigbs;
   
    tmpNs = NK[TNsp-1]/Nbr;

    tmpsigNbs = round((sigma[TNsp-1]+sigma[TNsp-2])/dx[0]);
    tmpN2sigbs = 2.0*round((sigma[TNsp-1]+sigma[TNsp-2])/dx[0])*dx[0];



    tmpTNsp = TNsp-2;
    for (i=0;i<tmpTNsp;i++) /*Molecules other than the surfactant*/
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

  /*For surfactant*/
    #pragma omp parallel for private(k,nn)
    for(j=0;j<extM2;j++)
    {
        unstack_ext2(j,nn);
        if((nn[Dim-1] >= padN2 ) and  (nn[Dim-1]< (Nx[Dim-1]+padN2)) )
        {
           nn[Dim-1] -= padN2;
           k = stack(nn);
           extwk[TNsp-2][j] = exp(amu[TNsp-2]+amu[TNsp-1]-dfttdrho[TNsp-2][k]-dfttdrho[TNsp-1][k]-wwall[TNsp-1][k]-wwall[TNsp-2][k]);
           wk[TNsp-2][k] = extwk[TNsp-2][j] ;
        }
        else if (nn[Dim-1]>= (Nx[Dim-1] +padN2))
        {
           extwk[TNsp-2][j] = exp(amu[TNsp-2]+amu[TNsp-1]-dfttdrho[TNsp-2][M-1]-dfttdrho[TNsp-1][M-1]-wwall[TNsp-2][M-1]-wwall[TNsp-1][M-1]);
        }
        else
        {
           extwk[TNsp-2][j] = exp(amu[TNsp-2]+amu[TNsp-1]-dfttdrho[TNsp-2][0]-dfttdrho[TNsp-1][0]-wwall[TNsp-2][0]-wwall[TNsp-1][0]);
        }
    }//j extM2
  
    q_tmp = (double*) malloc(extM2*sizeof(double));



    /*Calculate Propagators for Molecules other than the surfactant */   
    int tmp_double; 

    for (i=0;i<tmpTNsp;i++) /*Loop over Molecules other than the Surfactant*/
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

                    /*if(i==1)
                    {
                     printf("%d \t %d \t %d \t %.15e\n",i,l,j,q_tmp[j]);
                    }*/
                    
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
        }//j M

        //printf("\n");
    }//i tmpTNsp
    //exit(-1);
 

  /*For Surfactant*/

   //Weigths of child propagators for the branches
   for(i=0;i<Nbr;i++) /*Loop over number of branches*/
   {
    for(l=tmpNs-2;l>=0;l--) /*Loop over number of segments in each branch*/
    {

     // calc q^{C\alpha}_l = q^{C\beta}_{l+1}* wk
     for(j=0;j<extM2;j++)
     {
       unstack_ext2(j,nn);
       if((nn[Dim-1] >= padN2 ) and  (nn[Dim-1]< (Nx[Dim-1]+padN2)) )
       {
         nn[Dim-1] -= padN2;
         k = stack(nn);
         q_tmp[j] = extwk[TNsp-2][j]*qc_br_pgt[i][l+1][k];
       }
       else if (nn[Dim-1]>= (Nx[Dim-1] +padN2))
       {
         q_tmp[j] = extwk[TNsp-2][j]*qc_br_pgt[i][l+1][M-1];
       }
       else
       {
         q_tmp[j] = extwk[TNsp-2][j]*qc_br_pgt[i][l+1][0];
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
        qc_br_pgt[i][l][j] = dx[Dim-1]*integ_simpson(2*sigN[TNsp-1]+1,tmpN2[tid])/N2sig[TNsp-1];
      }//j  M


     } //l segments of each branch

    /*Weight for backbone-branch connecting point*/
     for(j=0;j<extM2;j++)
     {
       unstack_ext2(j,nn);
       if((nn[Dim-1] >= padN2 ) and  (nn[Dim-1]< (Nx[Dim-1]+padN2)) )
       {
         nn[Dim-1] -= padN2;
         k = stack(nn);
         q_tmp[j] = extwk[TNsp-2][j]*qc_br_pgt[i][0][k];
       }
       else if (nn[Dim-1]>= (Nx[Dim-1] +padN2))
       {
         q_tmp[j] = extwk[TNsp-2][j]*qc_br_pgt[i][0][M-1];
       }
       else
       {
         q_tmp[j] = extwk[TNsp-2][j]*qc_br_pgt[i][0][0];
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
        qc_bbr_pgt[i][j] = dx[Dim-1]*integ_simpson(2*tmpsigNbs+1,tmpN2[tid])/tmpN2sigbs;
      }


   }//Nbr number of branches

   

   //Weights of child propagators for the backbone beads
   fl=Nbr-1;
   for(l=NK[TNsp-2]-2;l>=0;l--) /*Loop over backbone segments*/
   {   
     
     // calc q^{C\alpha}_l = q^{C\beta}_{l+1}* wk
     if(brflag[l+1]>0)
     {
      for(j=0;j<extM2;j++)
      {
        unstack_ext2(j,nn);
        if((nn[Dim-1] >= padN2 ) and  (nn[Dim-1]< (Nx[Dim-1]+padN2)) )
        {
          nn[Dim-1] -= padN2;
          k = stack(nn);
          if(l==NK[TNsp-2]-2)
          {
           q_tmp[j] = extwk[TNsp-2][j]*qc_bbr_pgt[fl][k];
          }
          else
          {
           q_tmp[j] = extwk[TNsp-2][j]*qc_bb_pgt[l+1][k]*qc_bbr_pgt[fl][k];
          }
        }
        else if (nn[Dim-1]>= (Nx[Dim-1] +padN2))
        {
          if(l==NK[TNsp-2]-2)
          {
           q_tmp[j] = extwk[TNsp-2][j]*qc_bbr_pgt[fl][M-1];
          }
          else
          {
           q_tmp[j] = extwk[TNsp-2][j]*qc_bb_pgt[l+1][M-1]*qc_bbr_pgt[fl][M-1];
          }
        }
        else
        {
         if(l==NK[TNsp-2]-2)
         {
          q_tmp[j] = extwk[TNsp-2][j]*qc_bbr_pgt[fl][0];
         }
         else
         {
          q_tmp[j] = extwk[TNsp-2][j]*qc_bb_pgt[l+1][0]*qc_bbr_pgt[fl][0];
         }
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
        /*sigN and N2sig corresponds to backbone bead*/
        qc_bb_pgt[l][j] = dx[Dim-1]*integ_simpson(2*sigN[TNsp-2]+1,tmpN2[tid])/N2sig[TNsp-2];
      }//j  M

      fl=fl-1;
     } /*if over brflag*/

     else
     {
      for(j=0;j<extM2;j++)
      { 
        unstack_ext2(j,nn);
        if((nn[Dim-1] >= padN2 ) and  (nn[Dim-1]< (Nx[Dim-1]+padN2)) )
        { 
          nn[Dim-1] -= padN2;
          k = stack(nn);
          q_tmp[j] = extwk[TNsp-2][j]*qc_bb_pgt[l+1][k];
        } 
        else if (nn[Dim-1]>= (Nx[Dim-1] +padN2))
        {
          q_tmp[j] = extwk[TNsp-2][j]*qc_bb_pgt[l+1][M-1];
        }
        else
        {
         q_tmp[j] = extwk[TNsp-2][j]*qc_bb_pgt[l+1][0];
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
        qc_bb_pgt[l][j] = dx[Dim-1]*integ_simpson(2*sigN[TNsp-2]+1,tmpN2[tid])/N2sig[TNsp-2];
      }

     } /*else over brflag*/


    } //l backbone segments

    //Weights of parent propagators for the backbone
   fl=0;
   for(l=1;l<NK[TNsp-2];l++) /*Loop over backbone segments*/
   {

     // calc q^{P}_l = q^{P}_{l-1}* wk
     if(brflag[l-1]>0)
     {
      for(j=0;j<extM2;j++)
      {
       unstack_ext2(j,nn);
       if((nn[Dim-1] >= padN2 ) and  (nn[Dim-1]< (Nx[Dim-1]+padN2)) )
       {
         nn[Dim-1] -= padN2;
         k = stack(nn);
         q_tmp[j] = extwk[TNsp-2][j]*qp_bb_pgt[l-1][k]*qc_bbr_pgt[fl][k];
       }
       else if (nn[Dim-1]>= (Nx[Dim-1] +padN2))
       {
         q_tmp[j] = extwk[TNsp-2][j]*qp_bb_pgt[l-1][M-1]*qc_bbr_pgt[fl][M-1];
       }
       else
       {
         q_tmp[j] = extwk[TNsp-2][j]*qp_bb_pgt[l-1][0]*qc_bbr_pgt[fl][0];
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
  
        /*sigN and N2sig corresponds to backbone beads*/
        qp_bb_pgt[l][j] = dx[Dim-1]*integ_simpson(2*sigN[TNsp-2]+1,tmpN2[tid])/N2sig[TNsp-2];
      }//j  M

      fl = fl+1;
     } /*if over brflag*/

     else
     {
      for(j=0;j<extM2;j++)
      {
       unstack_ext2(j,nn);
       if((nn[Dim-1] >= padN2 ) and  (nn[Dim-1]< (Nx[Dim-1]+padN2)) )
       {
         nn[Dim-1] -= padN2;
         k = stack(nn);
         q_tmp[j] = extwk[TNsp-2][j]*qp_bb_pgt[l-1][k];
       }
       else if (nn[Dim-1]>= (Nx[Dim-1] +padN2))
       {
         q_tmp[j] = extwk[TNsp-2][j]*qp_bb_pgt[l-1][M-1];
       }
       else
       {
         q_tmp[j] = extwk[TNsp-2][j]*qp_bb_pgt[l-1][0];
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
        qp_bb_pgt[l][j] = dx[Dim-1]*integ_simpson(2*sigN[TNsp-2]+1,tmpN2[tid])/N2sig[TNsp-2];
      }//j  M

     } /*Else for brflag*/

    } //l backbone segments

     //Weights of parent propagators for the branch segments
   fl =0;
   for(i=0;i<Nbr;i++) /*Loop over number of branches*/
   {
    for(l=fl;l<NK[TNsp-2];l++)
    {
     if(brflag[l]>0)
     {
      fl = l;
      l=9999;
     }
    }
    for(l=0;l<tmpNs;l++) /*Loop over branch segments*/
    {

     // calc q^{C\alpha}_l = q^{C\beta}_{l+1}* wk
     if(l>0)
     {
      for(j=0;j<extM2;j++)
      {
       unstack_ext2(j,nn);
       if((nn[Dim-1] >= padN2 ) and  (nn[Dim-1]< (Nx[Dim-1]+padN2)) )
       {
         nn[Dim-1] -= padN2;
         k = stack(nn);
         q_tmp[j] = extwk[TNsp-2][j]*qp_br_pgt[i][l-1][k];
       }
       else if (nn[Dim-1]>= (Nx[Dim-1] +padN2))
       {
         q_tmp[j] = extwk[TNsp-2][j]*qp_br_pgt[i][l-1][M-1];
       }
       else
       {
         q_tmp[j] = extwk[TNsp-2][j]*qp_br_pgt[i][l-1][0];
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
        qp_br_pgt[i][l][j] = dx[Dim-1]*integ_simpson(2*sigN[TNsp-1]+1,tmpN2[tid])/N2sig[TNsp-1];
      }//j  M

     }/*if over l>0*/

     else
     {
      for(j=0;j<extM2;j++)
      {
       unstack_ext2(j,nn);
       if((nn[Dim-1] >= padN2 ) and  (nn[Dim-1]< (Nx[Dim-1]+padN2)) )
       {
         nn[Dim-1] -= padN2;
         k = stack(nn);
         if(fl==(NK[TNsp-2]-1))
         {
          q_tmp[j] = extwk[TNsp-2][j]*qp_bb_pgt[fl][k];
         }
         else
         {
          q_tmp[j] = extwk[TNsp-2][j]*qp_bb_pgt[fl][k]*qc_bb_pgt[fl][k];
         }
       }
       else if (nn[Dim-1]>= (Nx[Dim-1] +padN2))
       {
         if(fl==(NK[TNsp-2]-1))
         {
          q_tmp[j] = extwk[TNsp-2][j]*qp_bb_pgt[fl][M-1];
         }
         else
         {
          q_tmp[j] = extwk[TNsp-2][j]*qp_bb_pgt[fl][M-1]*qc_bb_pgt[fl][M-1];
         }
       }
       else
       {
         if(fl==(NK[TNsp-2]-1))
         {
          q_tmp[j] = extwk[TNsp-2][j]*qp_bb_pgt[fl][0];
         }
         else
         {
          q_tmp[j] = extwk[TNsp-2][j]*qp_bb_pgt[fl][0]*qc_bb_pgt[fl][0];
         }
       }

      }//j extM2

      //convolv q+1 with delta function
#pragma omp parallel for private(k,ind1,ind2)
      for(j=0;j<M;j++)
      {
        int tid = omp_get_thread_num() ;
        ind1 = j+padN2 - tmpsigNbs;  /*Modified sigN*/
        ind2 = j+padN2 + tmpsigNbs;  /*Modified sigN*/
        for(k=0;k<(2*tmpsigNbs+1);k++)
        {
          tmpN2[tid][k] = q_tmp[ind1+k];
        }//k
        /*Modify sigN adn N2sig accordingly*/
        qp_br_pgt[i][l][j] = dx[Dim-1]*integ_simpson(2*tmpsigNbs+1,tmpN2[tid])/tmpN2sigbs;
      }

     }

    } //l backbone segments
    fl=fl+1;
   } //Nbr number of brnaches

//Accumulating surfactant segment densities
   for(j=0;j<M;j++)
   {
      nwrhok[TNsp-2][j] = 0;
      nwrhok[TNsp-1][j] = 0;
      fl=0;
      for(l=0;l<NK[TNsp-2];l++)
      {
        ind1 = (l<MAXClen ? l : MAXClen-1);
        if(brflag[l]>0)
        {
         if(l==(NK[TNsp-2]-1))
         {
          nwrhok[TNsp-2][j] += wk[TNsp-2][j]*qp_bb_pgt[ind1][j]*qc_bbr_pgt[fl][j];
         }
         else
         {
          nwrhok[TNsp-2][j] += wk[TNsp-2][j]*qp_bb_pgt[ind1][j]*qc_bb_pgt[ind1][j]*qc_bbr_pgt[fl][j];
         }
         fl=fl+1;
        }
        else
        {
         nwrhok[TNsp-2][j] += wk[TNsp-2][j]*qp_bb_pgt[ind1][j]*qc_bb_pgt[ind1][j];
        }
//        printf("%d \t %d \t %d \t %.15e \t %.15e\n",i,ind1,j,qp_bb_pgt[ind1][j],qc_bb_pgt[ind1][j]);
      }//l NK[i]


      for(i=0;i<Nbr;i++)
      {
       for(l=0;l<tmpNs;l++)
       {
        ind1 = (l<MAXClen ? l : MAXClen-1);
        nwrhok[TNsp-1][j] += wk[TNsp-2][j]*qp_br_pgt[i][ind1][j]*qc_br_pgt[i][ind1][j];
//        printf("%d \t %d \t %d \t %.15e\n",i,ind1,j,qp_br_pgt[i][ind1][j]);
       }
      }//i Nbr
//      printf("%d \t %.15e \t %.15e \t %.15e\n",j,nwrhok[TNsp-2][j],nwrhok[TNsp-1][j],nwrhok[TNsp-2][j]+nwrhok[TNsp-1][j]);
    }//j M

//    exit(-1);
	
   free(q_tmp);
}

 
