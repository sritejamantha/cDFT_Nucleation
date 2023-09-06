#define MAIN
#include "globals.h"
void poly_density() ;
void initialize( double ) ;
void dfree();
double calc_diff();
double find_blkmu();
double calc_press( double );
int  shift_curv(double* , double ,int );

int main( int argc , char** argv ) 
{

  int i,j,k,track ;
  
  FILE *otp_tsn;
  
  double rdc,blkpres,chg_rho, surftsn;

  double inp_pres;

  if ( argc == 4 && !strcmp( "-nt" , argv[1] ) ) 
  {
    nthreads = atoi( argv[2] ) ;
    inp_pres = atof(argv[3]);

    omp_set_num_threads( nthreads ) ;
    printf("\nNumber of threads set to %d!\n" , nthreads ) ;
   
  }
  else 
  {
  	exit(1);
  }

  initialize(inp_pres) ;  


  char nm[40];

  snprintf (nm, sizeof nm, "%.0f_surftsn.dat", inp_pres);
  otp_tsn = fopen(nm,"w");

  rdc = (6.022e23)*pow(rf_len,3)*(1e-30);

  blkpres = find_blkmu();//*8.314*Tmpr/rdc;
  printf("Bulk Pressure:%.15f\n",blkpres*8.314*Tmpr/rdc);
//  exit(-1);
  track=0;
  for(i=0; i<ttstep; i++)
  {
        if(i==0)
        {
         tmp_lamb=1.0e-5;
        }

        if(i==1000 && chg_rho>1.0e-6)
        {
         tmp_lamb=1.0e-3;
        }

        if(chg_rho<1.0e-9 && track<3 && i>500)
        {
         tmp_lamb=1.0e-3;
         track=track+1;
        }
        
   	dfree();
	
	poly_density();

	chg_rho = calc_diff();


	if(i%otp_freq ==0)
        {

		for(j=0;j<TNsp;j++)
                {
		 snprintf (nm, sizeof nm, "rho%d_%.0f.dat", j,inp_pres);
		 write_grid_data(nm, rhoK[j]);

		}

		surftsn = calc_press(blkpres);
		
		fprintf(otp_tsn, "%d %lf %lf\n",i,surftsn,blkpres*8.314*Tmpr/rdc);fflush(stdout);
                //cout<<i<<" "<<chg_rho<<" "<<surftsn<<endl;
                printf("%d \t %d \t %.15e \t %.15e\n",i,track,chg_rho,surftsn);
	}

	if((chg_rho<1.e-10) and (i > 500)){
		cout<<"chg_rho "<<chg_rho<<endl;
	 	cout<<"sim finished "<<endl;
		break;
	}
  	

  }
}
