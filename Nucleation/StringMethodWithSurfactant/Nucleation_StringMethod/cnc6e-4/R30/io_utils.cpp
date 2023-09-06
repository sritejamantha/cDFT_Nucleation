#include "globals.h"
void unstack(int, int*  );
void unstack_ext(int, int*);
void unstack_ext2(int, int*);


void write_kspace_data( const char *nm , complex<double> *kdt ) {
  int i, j , nn[Dim] ;
  FILE *otp ;
  double kv[Dim], k2 ;
  
  otp = fopen( nm , "w" ) ;

  for ( i=1 ; i<M ; i++ ) {
    unstack( i , nn ) ;

    k2 = get_k( i , kv ) ;

    for ( j=0 ; j<Dim ; j++ ) 
      fprintf( otp , "%lf " , kv[j] ) ;

    fprintf( otp , "%1.5e %1.5e %1.5e %1.5e\n" , abs(kdt[i]), sqrt(k2), 
        real(kdt[i]) , imag(kdt[i]) ) ;

    if ( Dim == 2 && nn[0] == Nx[0]-1 )
      fprintf( otp , "\n" ) ;
  }

  fclose( otp ) ;


}

void write_grid_data( const char *nm , double *dat ) {

  int i, j, nn[Dim] ;
  FILE *otp ;
  double r[Dim] ;
  otp = fopen( nm , "w" ) ;

  for ( i=0 ; i<M ; i++ ) {
    unstack( i , nn ) ;
    
    for ( j=0 ; j<Dim ; j++ )
      fprintf( otp , "%lf " , double(nn[j]) * dx[j]* rf_len ) ;
    
    fprintf( otp , "%1.16e \n" , dat[i] ) ;
    
    if ( Dim == 2 && nn[0] == Nx[0]-1 )
      fprintf( otp , "\n" ) ;
  }

  fclose( otp ) ;

}

void write_extgrid_data( const char *nm , double *dat ) {

  int i, j, nn[Dim] ;
  FILE *otp ;
  double r[Dim] ;
  otp = fopen( nm , "w" ) ;

  for ( i=0 ; i<extM ; i++ ) {
    unstack_ext( i , nn ) ;
    
    for ( j=0 ; j<Dim ; j++ )
      fprintf( otp , "%lf " , double(nn[j]) * dx[j]* rf_len ) ;
    
    fprintf( otp , "%1.16e \n" , dat[i] ) ;
    
    if ( Dim == 2 && nn[0] == extNx[0]-1 )
      fprintf( otp , "\n" ) ;
  }

  fclose( otp ) ;

}

void write_ext2grid_data( const char *nm , double *dat ) {

  int i, j, nn[Dim] ;
  FILE *otp ;
  double r[Dim] ;
  otp = fopen( nm , "w" ) ;

  for ( i=0 ; i<extM2 ; i++ ) {
    unstack_ext2( i , nn ) ;
    
    for ( j=0 ; j<Dim ; j++ )
      fprintf( otp , "%lf " , double(nn[j]) * dx[j]* rf_len ) ;
    
    fprintf( otp , "%1.16e \n" , dat[i] ) ;
    
    if ( Dim == 2 && nn[0] == ext2Nx[0]-1 )
      fprintf( otp , "\n" ) ;
  }

  fclose( otp ) ;

}

void write_strg_profile(){

   int i, j, k;

   char nm[80];

   for(i=0;i<N_strg; i++){
	for(j=0; j<TNsp; j++){
		sprintf(nm, "Nstrg%d_rho%d.dat",i,j);
		 write_grid_data(nm,strg_rhoK[j][i]);
	}
   }

}
