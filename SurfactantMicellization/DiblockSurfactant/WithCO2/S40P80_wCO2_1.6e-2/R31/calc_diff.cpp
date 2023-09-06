#include "globals.h"
int stack(int*);
void unstack_ext2(int, int*  );
void unstack_ext(int, int*  );
double integ_simpson(int , double* );

double calc_diff(){
    int mid_ind = 450,i,j,k,l,ind1,ind2,nn[Dim];

    double ttdiff =0,tmp_end[2]; 


    for(i=0;i<M;i++)
    {
        for(j=0;j<TNsp;j++)
        {
		nwrhok[j][i] = (1-lamb)*rhoK[j][i] + lamb*nwrhok[j][i];
	}
    }

   
    for(j=0;j<TNsp;j++)
    {
    	//tmp_end[0] = rhoK[j][0];
	//tmp_end[1] = rhoK[j][M-1];

	/*for(i=0;i<M;i++){
		rhoK[j][i] /=tmp_end[1] / blkrho[j][1]; //(rhoK[j][i]- tmp_end[0])/(tmp_end[1]-tmp_end[0])*(blkrho[j][1]-blkrho[j][0])  + blkrho[j][0];
    	}*/	
       
        for(i=0;i<M;i++)
        {
	  
	 // if(abs(rhoK[j][i] - nwrhok[j][i])> ttdiff)
	  	//ttdiff = abs(rhoK[j][i] - nwrhok[j][i]);
	  if(i>= sigNhf[j])
	  	ttdiff += pow(rhoK[j][i] - nwrhok[j][i],2);
	 
	  rhoK[j][i] = nwrhok[j][i];
	}
	
	/*tmp_end[1] = rhoK[j][M-1];
	tmp_end[0] = rhoK[j][0];

	for(i=0;i<M;i++)
		rhoK[j][i] = (rhoK[j][i]- tmp_end[0])/(tmp_end[1]-tmp_end[0])*(blkrho[j][1]-blkrho[j][0])  + blkrho[j][0];
*/
    }

    for(i=0;i<M;i++)
    {
    	rhoK[TNsp][i] = 0;
    	for(j=0;j<TNsp;j++)
        {

		rhoK[TNsp][i] += rhoK[j][i];
    	}
        Xar[j] = rhoK[TNsp-2][j]/(rhoK[TNsp-2][j] + rhoK[TNsp-1][j]);
    }

   //for(i=0;i<M;i++){
   //	tmpar1[i] = extz[i+padN2]*extz[i+padN2]*(blkrho[1][1]-rhoK[1][i])*3;
   //}  
 
   //double intg1 = integ_simpson(M,tmpar1)*dx[Dim-1];
   //cout<<"lg mtlp "<<intg1<<" "<<pow(extz[padN2 + 438],3)*blkrho[1][1]<<endl;//pow(extz[padN2 + 432],3)*blkrho[1][1]<<endl;
  // lamb_mu += (1e-8*lamb)*(integ_simpson(M,tmpar1)*dx[Dim-1] - pow(extz[padN2 + 438],3)*blkrho[1][1]);

   
   //10*lamb*(rhoK[1][int(2*M/4)] - (blkrho[1][1]+blkrho[1][0])/2.);  

    return ttdiff;
}

void MaxShift(double*datr)
{
 int i,j,adj_v2,Max_ri=0;
 double tmpMax=0.0;
 
 tmpMax=0.0;
 for(i=0;i<M;i++)
 {
  if(datr[i]>tmpMax)
  {
   tmpMax=datr[i];
   Max_ri = i;
  }
 }

    
    for(i=0;i<TNsp;i++){
#pragma omp parallel for 
        for(j=0;j<Max_ri;j++){
                rhoK[i][j] = rhoK[i][Max_ri];
        }
    }

// exit(-1);

}


int shift_curv(double* datr,double cent_v,int aim_ind){
	int adj_v2,nn[Dim],i, j,k,mid_ind,mid_ind2;

	double excess , adj_v;
//        printf("%lf \t %d\n",cent_v,aim_ind);
	for(i=0; i< M-1 ;i++){
		adj_v = (datr[i]- cent_v)*(datr[i+1]-cent_v);
                	    

	    	if(adj_v == 0){
			mid_ind = ( (datr[i]- cent_v) == 0 ? i : i+1);
			excess = 0 ;
			break;
		
		}
		else if( adj_v < 0){

			mid_ind =  ( abs(datr[i]- cent_v) > abs(datr[i+1]- cent_v) ? i+1 : i    );
			//mid_ind2 = 1+i; ( abs(datr[i]- cent_v) > abs(datr[i+1]- cent_v) ? i+1 : i    );
			
			excess = abs(datr[i+1]- cent_v)/abs(datr[i+1]-datr[i]);//(abs(datr[i]- cent_v) > abs(datr[i+1]- cent_v) ? i+1 : i   )	
			break;
			

		}
		else{
			mid_ind = i;
		}

	}

  //      printf("mid_ind:%d\n",mid_ind);
    //    exit(-1);
	if(mid_ind == M-2){

		cout<<"can not find interf "<<endl;
		exit(1);
	}


  for(i=0;i<TNsp;i++){

#pragma omp parallel for private(k,nn)
   for(j=0;j<extM;j++){
	
	unstack_ext(j,nn);
  	if((nn[Dim-1] >= padN ) and  (nn[Dim-1]< (Nx[Dim-1]+padN)) ){
		nn[Dim-1] -= padN;
		k = stack(nn);
		extrhok[i][j] = rhoK[i][k];
	}
	else if (nn[Dim-1]>= (Nx[Dim-1] +padN)){
		extrhok[i][j] = rhoK[i][M-1];
	}
	else 
	{
	         extrhok[i][j] = rhoK[i][0];  
	}
   }//j

  }//i



    adj_v2 = aim_ind - mid_ind; 
//    printf("adj_v2:%d \n",adj_v2);			
//    exit(-1);
//    if(adj_v2 >= 0 )
//    	return 0;

    for(i=0;i<TNsp;i++){
#pragma omp parallel for 
	for(j=0;j<M;j++){ 
		rhoK[i][j] = extrhok[i][j+padN-adj_v2];
                /*if(i==0)
		{
		  printf("%d \t %.15e\n",j,rhoK[i][j]);
                }*/
	}
    }

//    exit(-1);
//    printf("%d \t %d \t %d\n",aim_ind,mid_ind,adj_v2);
    return adj_v2;

}
