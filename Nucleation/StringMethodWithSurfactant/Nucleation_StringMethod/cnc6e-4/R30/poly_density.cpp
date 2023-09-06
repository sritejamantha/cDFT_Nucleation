#include "globals.h"
int stack(int*);
double integ_simpson(int , double* );
void unstack_ext2(int, int*  );
void unstack(int, int*  );

void poly_density(){
    int adj_len,i,j,k,l,ind1,ind2,indz,nn[Dim];
    
    double *q_tmp,k2, vect[Dim],cent_v;
   
    for (i=0;i<TNsp;i++){
#pragma omp parallel for private(k,nn)
	for(j=0;j<extM2;j++){
		unstack_ext2(j,nn);
		if((nn[Dim-1] >= padN2 ) and  (nn[Dim-1]< (Nx[Dim-1]+padN2)) ){
			nn[Dim-1] -= padN2;
			k = stack(nn);
			extwk[i][j] = exp(amu[i]-dfttdrho[i][k]-wwall[i][k]);
			wk[i][k] = extwk[i][j] ;
		}
		else if (nn[Dim-1]>= (Nx[Dim-1] +padN2)){
				extwk[i][j] = exp(amu[i]-dfttdrho[i][M-1]-wwall[i][M-1]);
		}
		else{
				extwk[i][j] = exp(amu[i]-dfttdrho[i][0]-wwall[i][0]);
		}
	
	}//j extM2


    }//i TNsp
      

    //write_grid_data("wk1.dat", wk[1]);

    //exit(1);
    q_tmp = (double*) malloc(extM2*sizeof(double));
  

	/*calculate chain propagator*/
   
    int tmp_double; 

    for (i=0;i<TNsp;i++){

    	for(l=1;l<NK[i];l++){
		
		if( l >= MAXClen ){

		    /*for(j=0;j<M;j++){
			  q_pgt[l][j]  = q_pgt[l-1][j] ;
		    }*/
		    continue;
		}

		// calc q+1 = q* wk
		for(j=0;j<extM2;j++){
		    unstack_ext2(j,nn);
		    if((nn[Dim-1] >= padN2 ) and  (nn[Dim-1]< (Nx[Dim-1]+padN2)) ){
		    	 nn[Dim-1] -= padN2;
			 k = stack(nn);
			 q_tmp[j] = extwk[i][j]*q_pgt[l-1][k];
		    }
		    else if (nn[Dim-1]>= (Nx[Dim-1] +padN2)){
		    	q_tmp[j] = extwk[i][j]*q_pgt[l-1][M-1];
		    }
		    else{
			//q_tmp[j] = extwk[i][j]*q_pgt[l-1][0];
		    }
			
		}//j extM2

		//convolv q+1 with delta function
#pragma omp parallel for private(nn,k,ind1,ind2,adj_len)
		for(j=0;j<M;j++){
		    unstack(j,nn);

		    if(nn[Dim-1] < sigNhf[i])
		    	continue;



		    int tid = omp_get_thread_num() ;
		    if(j < sigN[i]){
		    	ind1 =  padN2 + sigN[i] - j;
		    	ind2 = j+padN2 + sigN[i];
		    	adj_len = 2*j +1;
		    }
		    else{
		    	ind1 = j+padN2 - sigN[i];
		    	ind2 = j+padN2 + sigN[i];
		    	adj_len = 2*sigN[i]+1;
		    }

		    
		    for(k=0;k<adj_len;k++){
			tmpN2[tid][k] = q_tmp[ind1+k]*extz[ind1+k];  
		    }//k
		    q_pgt[l][j] = dx[Dim-1]*integ_simpson(adj_len,tmpN2[tid])/N2sig[i]/extz[j+padN2];
		}//j  M

		cent_v = q_tmp[padN2 +sigN[i]] ;

		interp_values(q_pgt[l],cent_v,sigNhf[i]); 

		
	}//l NK

	int lft_n; 
	//accum particle density 
#pragma omp parallel for private(l,lft_n,ind1,ind2)
	for(j=0;j<M;j++){
		nwrhok[i][j] = 0;
		for(l=0;l<NK[i];l++){
			
			ind1 = (l<MAXClen ? l : MAXClen-1);
			ind2 = ((NK[i]-l-1) < MAXClen ? (NK[i]-l-1) : MAXClen-1);
			
		/*	if(wk[i][j]>1.){
				lft_n = -1 ;
			}
			else
				lft_n = l - ind1;
				lft_n += (NK[i]-l-1) - ind2;
*/
			nwrhok[i][j] += wk[i][j]*q_pgt[ind1][j]*q_pgt[ind2][j]; 	

		}//l NK[i]


	}//j M

/*	if(i==1){
	  	k2 = nwrhok[i][int(M-1)];
	 	for(j=0;j<M;j++) 
			nwrhok[i][j] *=blkrho[i][1]/k2;// (blkrho[i][1]+blkrho[i][0])/k2/2.;

		//	nwrhok[i][j] *= (blkrho[i][1]+blkrho[i][0])/k2/2.;


	}*/

    }//i TNsp
 
	
   free(q_tmp);
}

 
