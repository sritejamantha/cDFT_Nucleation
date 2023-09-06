#include "globals.h" 
void calc_lgrg();
void strg_reparamrz();
void update_strg_para();
double strg_density(int,double);
void strg_update_eular2();
void strg_dfree();

void strg_update()
{

  int i, j, k;

  double final_rmax = 0.0;

  for(j=1; j<(N_strg-1); j++)
  {
    strg_density(j,final_rmax);
  }//j N_strg
  for(i=0;i<TNsp;i++)
  {
    for(k=0 ; k<M; k++)
    {
	strg_newrhoK[i][0][k] =  strg_rhoK[i][0][k] ;
	strg_newrhoK[i][N_strg-1][k] =  strg_rhoK[i][N_strg-1][k] ;
    }
  }//i TNsp

//  exit(-1);
  
  strg_reparamrz();

  strg_dfree();
}


void strg_update_eular2(){

  int i, j, k;

  //double ;

  calc_lgrg();

  for(i=0; i<TNsp; i++){
	for(j=1; j<(N_strg-1); j++){
#pragma omp parallel for 
	    for(k=0 ; k<M; k++){

		 strg_newrhoK[i][j][k] =  strg_rhoK[i][j][k]  -lamb*(strg_dwdrho[i][j][k]  - drhodstrg[i][j][k]*inpdmudrds[i][j]/nmdrds[i][j]); //+lamb*(drhodstrg[i][j][k]*strg_lgrgk[i][j]/sqrt(nmdrds[i][j]));
		//strg_newrhoK[i][j][k] =  strg_rhoK[i][j][k]  -lamb*(strg_dwdrho[i][j][k]) ;
	    
	    
	    	if(strg_newrhoK[i][j][k] < 0 )
			strg_newrhoK[i][j][k] *=-1;
	    }
	//cout<<i<<" "<<j<<endl;
  
 	} 
 
 	for(k=0 ; k<M; k++){
		strg_newrhoK[i][0][k] =  strg_rhoK[i][0][k] ;
		strg_newrhoK[i][N_strg-1][k] =  strg_rhoK[i][N_strg-1][k] ;
	}
  
  }//i TNsp


  //strg_reparamrz();

  //update_strg_para();

  //write_grid_data("test_strg12_rho1.dat",strg_rhoK[1][29]);
  //write_grid_data("test_strg12_rho0.dat",strg_rhoK[0][29]);
}





void strg_update_eular()
{

  int i, j, k;

  //double ;

  calc_lgrg();

  for(i=0; i<TNsp; i++){
	for(j=1; j<(N_strg-1); j++){
#pragma omp parallel for 
	    for(k=0 ; k<M; k++){

		 strg_newrhoK[i][j][k] =  strg_rhoK[i][j][k]  -lamb*(strg_dwdrho[i][j][k]  - drhodstrg[i][j][k]*inpdmudrds[i][j]/nmdrds[i][j]) ;//+lamb*(drhodstrg[i][j][k]*strg_lgrgk[i][j]/sqrt(nmdrds[i][j]));
		//strg_newrhoK[i][j][k] =  strg_rhoK[i][j][k]  -lamb*(strg_dwdrho[i][j][k]) ;
	    
	    
	    	if(strg_newrhoK[i][j][k] < 0 )
			strg_newrhoK[i][j][k] = 1e-60;
	    }
	//cout<<i<<" "<<j<<endl;
  
 	} 
 
 	for(k=0 ; k<M; k++){
		strg_newrhoK[i][0][k] =  strg_rhoK[i][0][k] ;
		strg_newrhoK[i][N_strg-1][k] =  strg_rhoK[i][N_strg-1][k] ;
	}
  
  }//i TNsp


  strg_reparamrz();

  update_strg_para();

  //write_grid_data("test_strg12_rho1.dat",strg_rhoK[1][29]);
  //write_grid_data("test_strg12_rho0.dat",strg_rhoK[0][29]);
}





void update_strg_para(){

    int i,j,k;

     for(i=0;i<TNsp;i++){ 
	for(j=0;j<M;j++){ 
	 for(k=0; k<N_strg; k++){


		if(k>0 and k < (N_strg-1))
			drhodstrg[i][k][j] = (strg_rhoK[i][k+1][j] - strg_rhoK[i][k][j]);
		else if(k==0)
			drhodstrg[i][k][j] = (strg_rhoK[i][1][j] - strg_rhoK[i][0][j]); 
	
		else
			drhodstrg[i][k][j] = (strg_rhoK[i][N_strg-1][j] - strg_rhoK[i][N_strg-2][j]);
		

	 }//k
	
	
	
	}//j

	for(k=0; k < N_strg ; k++){
	    
		nmdrds[i][k] = integ_inpd(M,drhodstrg[i][k],drhodstrg[i][k])*dx[Dim-1];  
	}//k = N-trg
	

     }//i


}
