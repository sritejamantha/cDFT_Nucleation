#include "globals.h"

double integrate( double *inp ) {
  double sum = 0.0 ;
  int i ;

#pragma omp parallel for reduction(+:sum)
  for ( i=0 ; i<M ; i++ )
    sum += inp[i] ;

  for ( i=0 ; i<Dim ; i++ )
    sum *= dx[i] ;

  return sum ;
}

double integ_simpson(int np, double* dat){

  int i, j; 
  double sum=0;

  if(np<7){
	cout<<"error! pnt in sig is less than 7!"<<endl;
        exit(1);
  }

  sum = 3.0*dat[0]/8.0 + 7.0*dat[1]/6.0 + 23.0*dat[2]/24.0;

  for (i=3; i<np-3; i++)
	sum += dat[i];

  sum += 23.0*dat[np-3]/24.0 + 7.0*dat[np-2]/6.0 + 3.0*dat[np-1]/8.0;

  return sum;
}

