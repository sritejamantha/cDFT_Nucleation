#include "globals.h"


double calc_press(double blk_pres)
{

	double rdc,*grand, Leq, surftsn, gbv ,val;
	int i,j,k;

	grand =  (double*) malloc(M*sizeof(double));
	
	rdc = (6.022e23)*pow(rf_len,3)*(1e-30);
	

	for(i=0;i<M;i++)
        {
		grand[i] = ttF[i] ; 
		for(j=0;j<TNsp-2;j++)
                {
			if(rhoK[j][i] > 10e-9)
				grand[i] +=  (wwall[j][i]-amu[j])*rhoK[j][i];  

		}
                if(rhoK[TNsp-2][i]+rhoK[TNsp-1][i] > 10e-9)
                {
                  grand[i] +=  (wwall[j][i]-amu[TNsp-2]-amu[TNsp-1])*(rhoK[TNsp-2][i]+rhoK[TNsp-1][i]);
                }
	}

	gbv = M*grand[M-1];//blk_pres;//grand[M-1];

	for(i=0;i<Dim;i++)
		gbv *= dx[i] ;

	surftsn = (integrate( grand )- gbv)*8.314*Tmpr/rdc*rf_len*(1.e-10);

	write_grid_data("grand.dat",grand);

	free(grand);
	
	return surftsn;
}

