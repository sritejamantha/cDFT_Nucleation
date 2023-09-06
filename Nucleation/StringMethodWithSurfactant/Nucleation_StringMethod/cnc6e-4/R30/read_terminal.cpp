#include "globals.h"

void read_one_resume_file(FILE *inp, double* w ) {

        int i, j,nn[Dim];

	        double dr, di,dm[Dim];
		        for (i=0; i<M; i++) {

			                for (j=0; j<Dim; j++)
					      fscanf(inp,"%lf ", &dm[j]);						fscanf(inp, "%lf\n", &dr);

					w[i] = dr;						
			}	
}

void read_termianl() 
{
    char fname[100];
    FILE *inp; 
    int i;

    for(i=0;i<TNsp;i++)
    {
     sprintf(fname,"ter_rho%d.dat",i);
     inp = fopen(fname,"r");

     if (inp !=NULL)
     {

	read_one_resume_file(inp, tmil_rho1[i]) ;
    	 fclose(inp);
	 printf("Read ter_rho%d.dat!\n",i);
     }
     else
     {
	printf("can not read in ter_rho%d.dat!\n",i);
    	exit(1);
     }
    
    }

}


void read_rst() {

    FILE *inp; 
    char nm[80];

    int i,j,k;

    for(i=0;i<TNsp;i++)
	for(j=0;j<N_strg;j++){

    	sprintf(nm,"Nstrg%d_rho%d.dat",j,i);
	inp = fopen(nm,"r");

    	if (inp !=NULL){
		read_one_resume_file(inp,strg_rhoK[i][j]) ;
    	 	fclose(inp);
	 	cout << "Read restart rho "<<i<<" Nstrg "<<j<<" !" << endl;
    	}
    	else{
		cout<<"can not read restart rho "<<i<<" Nstrg "<<j<<" !" << endl;
    		exit(1);
    	}



		
	}
	


}
