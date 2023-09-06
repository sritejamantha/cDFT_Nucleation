#include "globals.h"
int stack(int*);
void unstack_ext2(int, int*  );

double calc_diff()
{
    int i,j,k,l,ind1,ind2,nn[Dim];

    double ttdiff =0,tmp_end[2]; 

    for(i=0;i<M;i++)
    {
	for(j=0;j<TNsp;j++)
        {
		nwrhok[j][i] = (1-tmp_lamb)*rhoK[j][i] + tmp_lamb*nwrhok[j][i];
	}
    }

    for(j=0;j<TNsp;j++)
    {
        for(i=0;i<M;i++)
        {
	  
	  if(abs(rhoK[j][i] - nwrhok[j][i])> ttdiff)
	  	ttdiff = abs(rhoK[j][i] - nwrhok[j][i]);
	  rhoK[j][i] = nwrhok[j][i];
	}
	
    }

    for(i=0;i<M;i++)
    {
    	rhoK[TNsp][i] = 0;
    	for(j=0;j<TNsp;j++)
        {

		rhoK[TNsp][i] += rhoK[j][i];
    	}
        Xar[j] = rhoK[TNsp-3][j]/(rhoK[TNsp-3][j] + rhoK[TNsp-2][j] + rhoK[TNsp-1][j]);
        Xbr[j] = rhoK[TNsp-2][j]/(rhoK[TNsp-3][j] + rhoK[TNsp-2][j] + rhoK[TNsp-1][j]);
    }

    if(prev_dif < ttdiff)
    {
     //tmp_lamb = 0.5*lamb;
     tmp_lamb = 0.99*tmp_lamb;
    }
    /*else
    {
     tmp_lamb=lamb;
    }*/

   prev_dif =ttdiff;
   return ttdiff;
}




int shift_curv(double* datr,double cent_v,int aim_ind)
{
	int adj_v2,nn[Dim],i, j,k,mid_ind,mid_ind2;

	double excess , adj_v;
	for(i=0; i< M-1 ;i++)
        {
		adj_v = (datr[i]- cent_v)*(datr[i+1]-cent_v);
	    
	    	if(adj_v == 0){
			mid_ind = ( (datr[i]- cent_v) == 0 ? i : i+1);
			excess = 0 ;
			break;
		
		}
		else if( adj_v < 0)
                {

			mid_ind =  ( abs(datr[i]- cent_v) > abs(datr[i+1]- cent_v) ? i+1 : i    );
			//mid_ind2 = 1+i; ( abs(datr[i]- cent_v) > abs(datr[i+1]- cent_v) ? i+1 : i    );
			
			excess = abs(datr[i+1]- cent_v)/abs(datr[i+1]-datr[i]);//(abs(datr[i]- cent_v) > abs(datr[i+1]- cent_v) ? i+1 : i   )	
			break;
			

		}
		else
                {
			mid_ind = i;
		}

	}


	if(mid_ind == M-2)
        {

		cout<<"can not find interf "<<endl;
		exit(1);
	
	}



  for(i=0;i<TNsp;i++)
  {

#pragma omp parallel for private(k,nn)
   for(j=0;j<extM;j++)
   {
	
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
			
    if(adj_v2 >= 0 )
    	return 0;

    for(i=0;i<TNsp;i++)
    {
#pragma omp parallel for 
	for(j=0;j<M;j++)
        { 
		rhoK[i][j] = extrhok[i][j+padN-adj_v2];
	}
    }

    return adj_v2;

}

