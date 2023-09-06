#include "globals.h"

double integ_simpson(int , double* );
void unstack_ext2(int, int* ) ;
int stack(int*);
void unstack(int, int*  );
void unstack_ext(int, int*  );


void get_dfnloc(){


	int i,sigN_12,sigN_23, sigN_13, j,k,l,m, ind1,ind2,nn[Dim];

	double cent_v,cent_v_tmp,tmp_epsln,tmp_sig,tmp_disp,sigma13,sigma12,sigma23;

	sigN_12 = round((sigma[0]+sigma[1])/2./dx[Dim-1]);
	sigma12 = round((sigma[0]+sigma[1])/4./dx[Dim-1])*dx[Dim-1]*2;
	
	if(TNsp >2){

	sigN_23 = round((sigma[2]+sigma[1])/2./dx[Dim-1]);
	sigma23 = round((sigma[2]+sigma[1])/4./dx[Dim-1])*dx[Dim-1]*2;
	
	sigN_13 = round((sigma[2]+sigma[0])/2./dx[Dim-1]);
	sigma13 = round((sigma[2]+sigma[0])/4./dx[Dim-1])*dx[Dim-1]*2;
	


	}
	
  	for (i=0;i<TNsp;i++){
#pragma omp parallel for private(k,nn)
		for(j=0;j<extM2;j++){
			
			unstack_ext2(j,nn);
			if((nn[Dim-1] >= padN2 ) and  (nn[Dim-1]< (Nx[Dim-1]+padN2)) ){
				nn[Dim-1] -= padN2;
				k = stack(nn);
				extrhok2[i][j] = rhoK[i][k];

			}
			else if (nn[Dim-1]>= (Nx[Dim-1] +padN2)){
				extrhok2[i][j] = rhoK[i][M-1];

			}
			else{
				extrhok2[i][j] = rhoK[i][0];

			}

		}//j, extM

 	}//i TNsp


	int tmp_length,tmp_flag,adj_len;

	for (i=0;i<TNsp;i++){
#pragma omp parallel for private(nn,adj_len,tmp_length,tmp_epsln,tmp_sig,tmp_flag,k,l,ind1,ind2)
		for(j=0;j<M;j++){

			unstack(j,nn);

		/*	if(nn[Dim-1] <sigNhf[i] ){
		  
		   		continue;
		   	}
		  */ 
		   dfdispndr[i][j] = 0;

		   for (k=0;k<TNsp;k++){
			
			if(k== i ){
			
				tmp_length  = sigN[i];
				tmp_flag = 0; 
				tmp_sig = Nsig[i];
				tmp_epsln = epsln[i];
			}
			else if((k==0 and i==1) or (k==1 and i == 0)){
				tmp_length = sigN_12;
				tmp_flag =1;
				tmp_sig = sigma12;
				tmp_epsln = crsepsln[0];
			}
			else if((k==0 and i==2) or (k==2 and i == 0)){
				tmp_length = sigN_13;
				tmp_flag =2;
				tmp_sig = sigma13;
				tmp_epsln = crsepsln[1];
			}
			else{
				tmp_length = sigN_23;
				tmp_flag =3; 
				tmp_sig = sigma23;
				tmp_epsln = crsepsln[2];	
			}
			int tid = omp_get_thread_num() ;
			
			if(nn[Dim-1] <sigNhf[i] ){
		  
		   		continue;
		   	}
	

			if(j<= tmp_length){
				ind1 = padN2 + tmp_length - j; 
				ind2 = padN2 + j+ tmp_length;
				adj_len = 2*j+1;
			}
			else{
				ind1 = j + padN2 - tmp_length;
				ind2 = j + padN2 + tmp_length;
				adj_len = 2*tmp_length+1;
			}

			
			for(l=0; l<adj_len;l++){
				

				tmpN2[tid][l] = (extrhok2[k][ind1+l]-rhoK[k][j])*tmp_epsln*(tmp_sig*tmp_sig-pow(tmp_sig,6.)/pow(extz[j+padN2]+extz[ind1+l],4.))*extz[ind1+l];	

			}

			dfdispndr[i][j] += integ_trap(adj_len,tmpN2[tid])/extz[j + padN2];
			
			ind1 = j + padN2 + tmp_length;
			ind2 = j +padN2 + 5*tmp_length;

			for(l=0; l<(4*tmp_length+1);l++){
				
				tmpN2[tid][l] = (extrhok2[k][ind1+l]-rhoK[k][j])*(1./pow(extz[ind1+l]-extz[j+padN2],4)- 1./pow(extz[ind1+l]+extz[j+padN2],4))*tmp_epsln*pow(tmp_sig,6)*extz[ind1+l];	

			}

			dfdispndr[i][j] += integ_simpson(4*tmp_length+1,tmpN2[tid])/extz[j + padN2];

		

			if(j<= 5*tmp_length){
				ind1 = padN2 ; 
				ind2 = padN2 + ( (j - tmp_length) >0 ? j-tmp_length : -1);
				adj_len = ind2-ind1+1;
			}
			else{
				ind1 = j +padN2 - 5*tmp_length;
				ind2 = j +padN2 - tmp_length;
				adj_len = 4*tmp_length+1;
			}

			for(l=0; l<adj_len;l++){
				
				tmpN2[tid][l] = (extrhok2[k][ind1+l]-rhoK[k][j])*(1./pow(extz[ind1+l]-extz[j+padN2],4)- 1./pow(extz[ind1+l]+extz[j+padN2],4))*tmp_epsln*pow(tmp_sig,6)*extz[ind1+l];	

			}

			if(adj_len>0)
			 	dfdispndr[i][j] += integ_trap(adj_len,tmpN2[tid])/extz[j + padN2];

			


		   }//k TNsp
			dfdispndr[i][j] *= -PI*dx[Dim-1]/Tmpr/4.;
		}//j M

		/* calcu lim r -0 dfdispndr*/
		cent_v = 0; 
		for (k=0;k<TNsp;k++){

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
				
				tmpN2[0][l] = (extrhok2[k][ind1+l]-rhoK[k][0])*tmp_epsln*pow(tmp_sig,6)/pow(extz[ind1+l],4);	

			}

			cent_v_tmp = -2*PI*dx[Dim-1]*integ_simpson(4*tmp_length+1,tmpN2[0]);
			cent_v += cent_v_tmp; 
		}// k TNsp

		interp_values(dfdispndr[i],cent_v,sigNhf[i]);

	
	}//i TNsp

      //write_grid_data("dfdispndr2.dat",dfdispndr[1]);




}
