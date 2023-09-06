#define MAIN
#include "globals.h"
void write_strg_profile();
int shift_curv(double* , double );
void poly_density() ;
void initialize( void ) ;
void strg_update();
void strg_update_eular();
void strg_update_eular2();
void strg_dfree();
double calc_diff();
double find_blkmu();
double calc_press( double );


int main( int argc , char** argv ) {

  int tmp_sft,  i,j,k ;
  
  FILE *otp_tsn;
  
  double rdc,blkpres,chg_rho, surftsn;


  if ( argc == 3 && !strcmp( "-nt" , argv[1] ) ) {
    nthreads = atoi( argv[2] ) ;
    omp_set_num_threads( nthreads ) ;
    printf("\nNumber of threads set to %d!\n" , nthreads ) ;
   
  }
  else {
    nthreads = 1 ;
    omp_set_num_threads( nthreads ) ;
    printf("\nNumber of threads set to %d!\n" , nthreads ) ;
  }

  initialize() ;  


  rdc = (6.022e23)*pow(rf_len,3)*(1e-30);

   double W_0 , W_1, dW;
   W_1 = calc_diff();
   
//   printf("%.11e\n",W_1);  
//   exit(-1); 

  for(i=0; i<ttstep; i++)
  {
      	if(i > strg_initstp)
        {
         strg_update();
	}
	else
        {
         strg_update_eular();
	 strg_dfree();
	}
	
	if(i%otp_freq ==0)
        {
	 W_0 = W_1 ;
	 W_1 = calc_diff();

	 dW = W_1 - W_0;
	
	 cout<<"step "<<i<<" delt f "<<dW<<endl;


	 otp_tsn = fopen("strg_ttf.dat","w");
	 for(j=0; j<N_strg; j++)
         {
			
	   fprintf(otp_tsn, "%f %1.8e \n", mxc_v[j]*pow(rf_len/10.,3),strg_itgF[j]); 

	 }
	 fclose(otp_tsn);
	
	}

	if(i%otp_freq ==0)
        {
	 write_strg_profile();
	}

        if((abs(dW)<1.e-8) and (i > 5))
        {
                cout<<"dW "<<dW<<endl;
                cout<<"sim finished "<<endl;
                break;
        }

  }
}
