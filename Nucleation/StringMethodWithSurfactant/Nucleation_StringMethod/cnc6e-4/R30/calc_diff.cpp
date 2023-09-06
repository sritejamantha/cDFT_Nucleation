#include "globals.h"
int stack(int*);
void unstack_ext2(int, int*  );
void unstack_ext(int, int*  );

double calc_diff(){
    int mid_ind = int(M/2.),i,j,k,l,ind1,ind2,nn[Dim];

    double *tmp_mxc_r ,ttdiff =0,tmp_end[2]; 

    tmp_mxc_r = (double*) calloc(M,sizeof(double));

    for(i=0; i<N_strg; i++)
    {

	strg_itgF[i] = integ_simpson_sph(M,strg_ttF[i])*dx[Dim-1]*dx[Dim-1]*dx[Dim-1];

	for(j =0 ; j<M ;j++)
		tmp_mxc_r[j] = (strg_rhoK[0][0][M-1] - strg_rhoK[0][i][j])/ strg_rhoK[0][0][M-1] ; 
	
	mxc_v[i] = integ_simpson_sph(M,tmp_mxc_r)*dx[Dim-1]*dx[Dim-1]*dx[Dim-1];
//        printf("%d \t %.14f \t %.9f\n",i,mxc_v[i],strg_itgF[i]);
 	//cout<<strg_itgF[i]<<endl;
    }

//    exit(-1);


    free(tmp_mxc_r);
    

   
    return integ_trap(N_strg, strg_itgF)/double(N_strg) ;
 


}
