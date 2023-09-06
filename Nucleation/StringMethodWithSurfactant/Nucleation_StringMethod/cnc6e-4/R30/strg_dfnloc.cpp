#include "globals.h"

double integ_simpson(int , double* );
void unstack_ext2(int, int* ) ;
int stack(int*);
void unstack(int, int*  );
void unstack_ext(int, int*  );


void get_strdfnloc(int sti){


	int i,sigN_12,sigN_23, sigN_13, j,k,l,m, ind1,ind2,nn[Dim];

	double cent_v,cent_v_tmp,tmp_epsln,tmp_sig,tmp_disp,sigma13,sigma12,sigma23;
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
				extrhok2[i][j] = strg_rhoK[i][sti][k];

			}
			else if (nn[Dim-1]>= (Nx[Dim-1] +padN2))
			{
				extrhok2[i][j] = strg_rhoK[i][sti][M-1];

			}
			else
			{
				extrhok2[i][j] = strg_rhoK[i][sti][0];

			}

		}//j, extM

 	}//i TNsp


	int tmp_length,tmp_flag,adj_len;

	for (i=0;i<TNsp;i++)
	{
	   for(j=0;j<M;j++)
		dfdispndr[i][j] = 0;

           for (k=0;k<TNsp;k++)
           {
		
		tmp_length  = sigN_ij[i][k];
                tmp_flag = 0;
                tmp_sig = sigmaij[i][k];
                tmp_epsln = epsij[i][k];	
	
#pragma omp parallel for private(nn,adj_len,l,ind1,ind2)
		for(j=0;j<M;j++)
		{
			unstack(j,nn);
			dfdispndr_pir[i][k][j] = 0;
			int tid = omp_get_thread_num() ;
			
			if(nn[Dim-1] <tmp_length )
			{
				continue;


			}
			
			if(j<= tmp_length)
			{
				ind1 = padN2 + tmp_length - j; 
				ind2 = padN2 + j+ tmp_length;
				adj_len = 2*j+1;
			}
			else
			{
				ind1 = j + padN2 - tmp_length;
				ind2 = j + padN2 + tmp_length;
				adj_len = 2*tmp_length+1;
			}

			
			for(l=0; l<adj_len;l++)
			{
				

				tmpN2[tid][l] = (extrhok2[k][ind1+l]-strg_rhoK[k][sti][j])*tmp_epsln*(tmp_sig*tmp_sig-pow(tmp_sig,6.)/pow(extz[j+padN2]+extz[ind1+l],4.))*extz[ind1+l];	

			}

			dfdispndr_pir[i][k][j] += integ_trap(adj_len,tmpN2[tid])/extz[j + padN2];
			
			ind1 = j + padN2 + tmp_length;
			ind2 = j +padN2 + 3*tmp_length;

			for(l=0; l<(2*tmp_length+1);l++)
			{
				
				tmpN2[tid][l] = (extrhok2[k][ind1+l]-strg_rhoK[k][sti][j])*(1./pow(extz[ind1+l]-extz[j+padN2],4)- 1./pow(extz[ind1+l]+extz[j+padN2],4))*tmp_epsln*pow(tmp_sig,6)*extz[ind1+l];	

			}

			dfdispndr_pir[i][k][j] += integ_simpson(2*tmp_length+1,tmpN2[tid])/extz[j + padN2];

		

			if(j<= 3*tmp_length)
			{
				ind1 = padN2 ; 
				ind2 = padN2 + ( (j - tmp_length) >0 ? j-tmp_length : -1);
				adj_len = ind2-ind1+1;
			}
			else
			{
				ind1 = j +padN2 - 3*tmp_length;
				ind2 = j +padN2 - tmp_length;
				adj_len = 2*tmp_length+1;
			}

			for(l=0; l<adj_len;l++)
			{
				
				tmpN2[tid][l] = (extrhok2[k][ind1+l]-strg_rhoK[k][sti][j])*(1./pow(extz[ind1+l]-extz[j+padN2],4)- 1./pow(extz[ind1+l]+extz[j+padN2],4))*tmp_epsln*pow(tmp_sig,6)*extz[ind1+l];	

			}

			if(adj_len>0)
			 	dfdispndr_pir[i][k][j] += integ_trap(adj_len,tmpN2[tid])/extz[j + padN2];

			
			dfdispndr_pir[i][k][j] *= -PI*dx[Dim-1]/Tmpr/4.;


		}//j M

		/* calcu lim r -0 dfdispndr*/
		cent_v = 0; 
		/*for (k=0;k<TNsp;k++){

			if(k== i ){
				tmp_length = sigN[i];
				tmp_sig = Nsig[i];
				tmp_epsln = epsln[i];
			}
			else if((k==0 and i==1) or (k==1 and i == 0)){
				tmp_length = sigN_12;
				tmp_sig = sigma12;
				tmp_epsln = crsepsln[0];
			}
			else if((k==0 and i==2) or (k==2 and i == 0)){
				tmp_length = sigN_13;
				tmp_sig = sigma13;
				tmp_epsln = crsepsln[1];
			}
			else{
				tmp_length = sigN_23;
				tmp_sig = sigma23;
				tmp_epsln = crsepsln[2];	
			}

			ind1 = padN2 + tmp_length;
			ind2 = padN2 + 5*tmp_length;

			for(l=0; l<(4*tmp_length+1);l++){
				
				tmpN2[0][l] = (extrhok2[k][ind1+l]-strg_rhoK[k][sti][0])*tmp_epsln*pow(tmp_sig,6)/pow(extz[ind1+l],4);	

			}

			cent_v_tmp = -2*PI*dx[Dim-1]*integ_simpson(4*tmp_length+1,tmpN2[0]);
			cent_v += cent_v_tmp; 
		}// k TNsp
		*/	
		interp_values(dfdispndr_pir[i][k],cent_v,tmp_length);

		for(j=0;j<M; j++) 
			dfdispndr[i][j] += dfdispndr_pir[i][k][j];	
	      }//k TNsp	
	}//i TNsp

      //write_grid_data("dfdispndr2.dat",dfdispndr[1]);




}
