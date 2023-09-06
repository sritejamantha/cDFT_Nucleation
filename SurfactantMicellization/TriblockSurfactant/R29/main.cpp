#define MAIN
#include "globals.h"
void MaxShift(double*);
int shift_curv(double* , double , int);
void poly_density() ;
void initialize( void ) ;
void dfree();
double calc_diff();
double find_blkmu();
double calc_press( double );


int main( int argc , char** argv ) {

  int tmp_sft=0,  i,j,k ;

  char nm[80];
  FILE *otp_tsn;
  otp_tsn = fopen("surftsn.dat","w");
  
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

  blkpres = find_blkmu();//*8.314*Tmpr/rdc;
  printf("%lf\n",blkpres*8.314*Tmpr/rdc);
//  amu[0] = amuP;
//  exit(1);

  for(i=0; i<ttstep; i++){
   
   	dfree();
	
	poly_density();

	chg_rho = calc_diff(); 

	if((i>= 5) and (i%5 == 0)){
//		MaxShift(rhoK[1]);
		tmp_sft = shift_curv(rhoK[0], (blkrho[0][0]+blkrho[0][1])/2.,145 );
	}

	if( (i%otp_freq2 == 0) and (i>0)){
		sprintf(nm,"rho0_%d.dat",i);
		write_grid_data(nm, rhoK[0]);

		sprintf(nm,"rho1_%d.dat",i);
		write_grid_data(nm, rhoK[1]);
	
		sprintf(nm,"rho2_%d.dat",i);
                write_grid_data(nm, rhoK[2]);

		sprintf(nm,"rho3_%d.dat",i);
                write_grid_data(nm, rhoK[3]);

                sprintf(nm,"Grand_%d.dat",i);
                write_grid_data(nm, Grandr);
	}

	if(i%otp_freq ==0){

		write_grid_data("rho0.dat", rhoK[0]);
		write_grid_data("rho1.dat", rhoK[1]);
		write_grid_data("rho2.dat", rhoK[2]);
		write_grid_data("rho3.dat", rhoK[3]);
                write_grid_data("Grand.dat", Grandr);
// 		write_grid_data("rhott.dat", rhoK[TNsp]);

		surftsn = calc_press(blkpres);
		
		fprintf(otp_tsn, "%lf %lf\n",blkpres*8.314*Tmpr/rdc, surftsn);fflush(stdout);
		if(i < 0)	
			cout<<i<<" "<<chg_rho<<" "<<surftsn<<" "<<lamb_mu<<endl;
		else
			cout<<i<<" "<<chg_rho<<" "<<surftsn<<" shift "<<tmp_sft<<endl;
	}

	if((chg_rho<1.e-15) and (i > 5)){
		cout<<"chg_rho "<<chg_rho<<endl;
	 	cout<<"sim finished "<<endl;
		break;
	}

  	

  }
}
