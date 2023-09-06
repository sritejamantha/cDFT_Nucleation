#include "globals.h"

void calc_lgrg(){
int i , j, k;
double *tmp_lgrgk ; 

tmp_lgrgk = (double*) calloc(M, sizeof(double));
	for(i=0; i<TNsp;i++){
		for(j=1;j<(N_strg-1);j++){
		 for(k=0; k<M;k++){
			tmp_lgrgk[k] = strg_rhoK[i][j+1][k] - strg_rhoK[i][j][k];
		 }
		strg_lgrgk[i][j] = sqrt(integ_inpd(M,tmp_lgrgk,tmp_lgrgk)*dx[Dim-1]);
		
		for(k=0; k<M;k++){
		      tmp_lgrgk[k] = strg_rhoK[i][j][k] - strg_rhoK[i][j-1][k];
		}

		strg_lgrgk[i][j] -= sqrt(integ_inpd(M,tmp_lgrgk,tmp_lgrgk)*dx[Dim-1]);
		
	       	 strg_lgrgk[i][j] /= NK[i];
	       }
	}

free(tmp_lgrgk);
}
