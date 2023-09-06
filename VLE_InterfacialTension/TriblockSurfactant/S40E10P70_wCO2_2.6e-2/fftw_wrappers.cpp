#include "globals.h"
//int fftw_init_threads( void ) ;
//void fftw_plan_with_nthreads( int );


void fftw_fwd( double* in , complex<double>* out ) {

  int i;
  double norm = 1.0 / double( M );

#pragma omp parallel for
  for ( i=0 ; i<M ; i++ ) {
    fin[i][0] = in[i];
    fin[i][1] = 0.0 ;
  }

  fftw_execute( ft_fwd ) ;


#pragma omp parallel for
  for ( i=0 ; i<M ; i++ ) 
    out[i] = ( fout[i][0] + I * fout[i][1] ) * norm ;

}



void fftw_back( complex<double>* in , double* out ) {

  int i;

#pragma omp parallel for
  for ( i=0 ; i<M ; i++ ) {
    fin[i][0] = real( in[i] ) ;
    fin[i][1] = imag( in[i] ) ;
  }

  fftw_execute( ft_bck ) ;

#pragma omp parallel for
  for ( i=0 ; i<M ; i++ )
    out[i] = fout[i][0] ;

}





void fft_init() {

  int rtn, Nfp[Dim], j ;

  fin = ( fftw_complex* ) fftw_malloc( M * sizeof( fftw_complex ) ) ;
  fout = ( fftw_complex* ) fftw_malloc( M * sizeof( fftw_complex ) ) ;

  rtn = fftw_init_threads() ;

  if ( !rtn ) die("Error in fftw initialization!\n") ;

  fftw_plan_with_nthreads( nthreads ) ;

  for ( j=0 ; j<Dim ; j++ )
    Nfp[j] = Nx[Dim-j-1] ;

  ft_fwd = fftw_plan_dft( Dim , Nfp , fin, fout, FFTW_FORWARD , FFTW_MEASURE ) ;

  ft_bck = fftw_plan_dft( Dim , Nfp , fin, fout, FFTW_BACKWARD , FFTW_MEASURE ) ;

  mem_use += M * sizeof( fftw_complex ) * 2 ;

}
