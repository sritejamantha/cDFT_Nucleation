#include "globals.h"

double integ_simpson(int , double* );
void unstack_ext2(int, int* ) ;
int stack(int*);
void unstack(int, int*  );
void unstack_ext(int, int*  );


void get_dfnloc()
{


	int i,j,k,l,m, ind1,ind2,nn[Dim];

	double tmp_epsln,tmp_sig,tmp_disp,sigN_12,sigN_23,sigN_13,sigma13,sigma12,sigma23;
        double sigN_ij[TNsp][TNsp],sigmaij[TNsp][TNsp],epsij[TNsp][TNsp];
       

        k=0;
        for(i=0;i<TNsp;i++)
        {
         for(j=i;j<TNsp;j++)
         {
          sigN_ij[i][j] = round((sigma[i]+sigma[j])/2./dx[Dim-1]);
          sigmaij[i][j] = round((sigma[i]+sigma[j])/4./dx[Dim-1])*dx[Dim-1]*2;
          if(i==j)
          {
           epsij[i][j] = epsln[i];
          }
          else
          {
           epsij[i][j] = crsepsln[k];
           epsij[j][i] = epsij[i][j]; 
           sigN_ij[j][i] = sigN_ij[i][j];
           sigmaij[j][i] = sigmaij[i][j];
           k=k+1;
          }
         }
        }

 
  	for (i=0;i<TNsp;i++)
	{
#pragma omp parallel for private(k,nn)
		for(j=0;j<extM2;j++)
		{
			
			unstack_ext2(j,nn);
			if((nn[Dim-1] >= padN2 ) and  (nn[Dim-1]< (Nx[Dim-1]+padN2)) )
			{
				nn[Dim-1] -= padN2;
				k = stack(nn);
				extrhok2[i][j] = rhoK[i][k];

			}
			else if (nn[Dim-1]>= (Nx[Dim-1] +padN2))
			{
				extrhok2[i][j] = rhoK[i][M-1];

			}
			else
			{
				extrhok2[i][j] = rhoK[i][0];

			}

		}//j, extM

 	}//i TNsp



	int tmp_length,tmp_flag;
	double tmp_sigr;

	for (i=0;i<TNsp;i++)
	{
#pragma omp parallel for private(tmp_sigr,tmp_length,tmp_epsln,tmp_sig,tmp_flag,k,l,ind1,ind2)
		for(j=0;j<M;j++)
		{
		   
		   dfdispndr[i][j] = 0;

		   for (k=0;k<TNsp;k++)
		   {
                        tmp_length  = sigN_ij[i][k];
                        tmp_flag = 0;
                        tmp_sig = sigmaij[i][k];
                        tmp_sigr = (sigma[i]+sigma[k])/2.;
                        tmp_epsln = epsij[i][k];
	
			int tid = omp_get_thread_num() ;
			ind1 = j + padN2 - tmp_length;
			ind2 = j + padN2 + tmp_length;

			for(l=0; l<(2*tmp_length+1);l++)
			{
				
				tmpN2[tid][l] = (extrhok2[k][ind1+l]-rhoK[k][j])*tmp_epsln*tmp_sigr*tmp_sigr;	

			}

			dfdispndr[i][j] += integ_simpson(2*tmp_length+1,tmpN2[tid])*tmp_sigr/tmp_sig;
	
			ind1 = j + padN2 + tmp_length;
			ind2 = j +padN2 + 5*tmp_length;

			for(l=0; l<(4*tmp_length+1);l++)
			{
				
				tmpN2[tid][l] = (extrhok2[k][ind1+l]-rhoK[k][j])/pow(extz[ind1+l]-extz[j+padN2],4)*tmp_epsln*pow(tmp_sigr,6);	

			}

			 dfdispndr[i][j] += integ_simpson(4*tmp_length+1,tmpN2[tid])*pow(tmp_sig/tmp_sigr,3);
                        
			ind1 = j +padN2 - 5*tmp_length;
			ind2 = j +padN2 - tmp_length;

			for(l=0; l<(4*tmp_length+1);l++)
			{
				
				tmpN2[tid][l] = (extrhok2[k][ind1+l]-rhoK[k][j])/pow(extz[ind1+l]-extz[j+padN2],4)*tmp_epsln*pow(tmp_sigr,6);	

			}

			 dfdispndr[i][j] += integ_simpson(4*tmp_length+1,tmpN2[tid])*pow(tmp_sig/tmp_sigr,3);
			
                        

		   }//k TNsp
			dfdispndr[i][j] *= -PI*dx[Dim-1]/Tmpr/4.;
		}//j M

	
	}//i TNsp

      //write_grid_data("dfdispndr1.dat",dfdispndr[0]);




}
