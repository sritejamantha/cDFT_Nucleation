#include "globals.h"


double calc_press(double blk_pres){

	double rdc,*grand, Leq, surftsn, gbv ,val;
	int i,j,k;

	grand =  (double*) malloc(M*sizeof(double));
	
	rdc = (6.022e23)*pow(rf_len,3)*(1e-30);
	
        for(i=0;i<M;i++)
        {
         Grandr[i]=0.0;
        }
	for(i=0;i<M;i++)
	{
		grand[i] = ttF[i] ; 
		tmpar1[i] = 0 ;
		for(j=0;j<TNsp-2;j++)
		{
			if(rhoK[j][i] > 1.0e-9)
			{
				grand[i] += (wwall[j][i]-amu[j])*rhoK[j][i];
			}
		}
		if(rhoK[TNsp-2][i]+rhoK[TNsp-1][i] > 1.0e-9)
                {
                  grand[i] +=  (wwall[j][i]-amu[TNsp-2]-amu[TNsp-1])*(rhoK[TNsp-2][i]+rhoK[TNsp-1][i]);
                }
                Grandr[i]=grand[i];
		tmpar1[i] = extz[padN2 +i]*extz[padN2 +i]*grand[i];
//		printf("%d \t %.15f\n",i,grand[i]);
	}

        

	gbv = pow(extz[padN2+M-1],3)*grand[M-1]*4.*PI/3.;

	//for(i=0;i<Dim;i++)
	//	gbv *= dx[i] ;

	surftsn = (dx[Dim-1]*integ_trap( M,tmpar1 )*4.*PI- gbv)*8.314*Tmpr/rdc*rf_len*(1.e-10);

	write_grid_data("grand.dat",grand);

	free(grand);
	
	return surftsn;
}

