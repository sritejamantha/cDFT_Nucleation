#include "globals.h"

double pbc_mdr2( double *r1 , double *r2 , double *dr ) {

  int j ;
  double mdr2 = 0.0 ;

  for ( j=0 ; j<Dim ; j++ ) {
    dr[j] = r1[j] - r2[j] ;
    
    if ( dr[j] > Lh[j] )
      dr[j] -= L[j] ;

    else if ( dr[j] <= -Lh[j] )
      dr[j] += L[j] ;

    mdr2 += dr[j] * dr[j] ;
  }

  return mdr2 ;

}



int stack_input(int x[Dim], int Nxx[Dim]) {
  if (Dim==1)
    return x[0];
  else if (Dim==2)
    return (x[0] + x[1]*Nxx[0]);
  else
    return  (x[0] + (x[1] + x[2]*Nxx[1])*Nxx[0] );
}

// Stacks x using only local values
int stack_local(int x[Dim]) {
  if (Dim==1)
    return x[0];
  else if (Dim==2)
    return (x[0] + x[1]*Nx[0]);
  else
    return  (x[0] + (x[1] + x[2]*Nx[1])*Nx[0] );
}

// Stacks vector x into 1D array index in [ 0, M ]
int stack( int x[Dim] ) { 
  if (Dim==1)
    return x[0];
  else if (Dim==2)
    return (x[0] + x[1]*Nx[0]);
  else
    return  (x[0] + (x[1] + x[2]*Nx[1])*Nx[0] );
}

void unstack_input(int id, int nn[Dim], int Nxx[Dim]) {

  if (Dim==1) {
    nn[0] = id;
    return;
  }
  else if (Dim==2) {
    nn[1] = id/Nxx[0];
    nn[0] = (id - nn[1]*Nxx[0]);
    return;
  }
  else if (Dim==3) {
    nn[2] = id/Nxx[1]/Nxx[0];
    nn[1] = id/Nxx[0] - nn[2]*Nxx[1];
    nn[0] = id - (nn[1] + nn[2]*Nxx[1])*Nxx[0];
  }
  else {
    cout << "Dim is goofy!" << endl;
    return;
  }
}

// Receives index id in [0 , M ] and makes array
// nn[Dim] in [ 0 , Nx[Dim] ]
void unstack(int id, int nn[Dim] ) {

  if (Dim==1) {
    nn[0] = id;
    return;
  }
  else if (Dim==2) {
    nn[1] = id/Nx[0];
    nn[0] = (id - nn[1]*Nx[0]);
    return;
  }
  else if (Dim==3) {
    nn[2] = id/Nx[1]/Nx[0];
    nn[1] = id/Nx[0] - nn[2]*Nx[1];
    nn[0] = id - (nn[1] + nn[2]*Nx[1])*Nx[0];
  }
  else {
    cout << "Dim is goofy!" << endl;
    return;
  }
}


void unstack_ext2(int id, int nn[Dim] ) {

  if (Dim==1) {
    nn[0] = id;
    return;
  }
  else if (Dim==2) {
    nn[1] = id/extNx2[0];
    nn[0] = (id - nn[1]*extNx2[0]);
    return;
  }
  else if (Dim==3) {
    nn[2] = id/extNx2[1]/extNx2[0];
    nn[1] = id/extNx2[0] - nn[2]*extNx2[1];
    nn[0] = id - (nn[1] + nn[2]*extNx2[1])*extNx2[0];
  }
  else {
    cout << "Dim is goofy!" << endl;
    return;
  }
}







void unstack_ext(int id, int nn[Dim] ) {

  if (Dim==1) {
    nn[0] = id;
    return;
  }
  else if (Dim==2) {
    nn[1] = id/extNx[0];
    nn[0] = (id - nn[1]*extNx[0]);
    return;
  }
  else if (Dim==3) {
    nn[2] = id/extNx[1]/extNx[0];
    nn[1] = id/extNx[0] - nn[2]*extNx[1];
    nn[0] = id - (nn[1] + nn[2]*extNx[1])*extNx[0];
  }
  else {
    cout << "Dim is goofy!" << endl;
    return;
  }
}




void get_r( int id , double r[Dim] ) {
  int i, id2, n[Dim];

  unstack(id, n);


  for ( i=0; i<Dim; i++) {
    r[i] = dx[i] * double( n[i] );
    
    if ( r[i] > L[i]/2.0 )
      r[i] -= L[i];
    else if ( r[i] <= -L[i]/2.0 )
      r[i] += L[i];
 
  }
}

double get_k_alias( int id , double k[Dim] ) {

  double kmag = 0.0;
  int i, id2, n[Dim] , j , has_nyquist = 0;
  for ( i=0 ; i<Dim ; i++ )
    if ( Nx[i] % 2 == 0 )
      has_nyquist = 1;

  unstack(id, n);

  if ( Nx[0] % 2 == 0 && n[0] == Nx[0] / 2 )
    k[0] = 0.0 ;
  else if ( double(n[0]) < double(Nx[0]) / 2.)
   k[0] = 2*PI*double(n[0])/L[0];
  else
   k[0] = 2*PI*double(n[0]-Nx[0])/L[0];

  if (Dim>1) {
    if ( Nx[1] % 2 == 0 && n[1] == Nx[1] / 2 )
      k[1] = 0.0 ;
    else if ( double(n[1]) < double(Nx[1]) / 2.)
      k[1] = 2*PI*double(n[1])/L[1];
    else
      k[1] = 2*PI*double(n[1]-Nx[1])/L[1];
  }

  if (Dim==3) {
    if ( Nx[2] % 2 == 0 && n[2] == Nx[2] / 2 )
      k[2] = 0.0 ;
    else if ( double(n[2]) < double(Nx[2]) / 2.)
      k[2] = 2*PI*double(n[2])/L[2];
    else
      k[2] = 2*PI*double(n[2]-Nx[2])/L[2];
  }

  // Kills off the Nyquist modes
  if ( id2 != 0 && has_nyquist ) {
    for ( i=0 ; i<Dim ; i++ ) {
      if ( k[i] == 0.0 ) {
        for ( j=0 ; j<Dim ; j++ )
          k[j] = 0.0 ;
        kmag = 0.0;
        break;
      }
    }
  }
  
  for (i=0; i<Dim; i++)
    kmag += k[i]*k[i];

  return kmag;

}



// Receives index id in [ 0 , ML ] and returns 
// proper k-value, whether running in parallel or not
double get_k(int id, double k[Dim]) {

  double kmag = 0.0;
  int i, id2, n[Dim];

  unstack(id, n);

  for ( i=0 ; i<Dim ; i++ ) {
    if ( double( n[i] ) < double( Nx[i] ) / 2. )
      k[i] = 2 * PI * double( n[i] ) / L[i] ;
    else
      k[i] = 2 * PI * double( n[i] - Nx[i] ) / L[i] ;

    kmag += k[i] * k[i] ;
  }

  return kmag;

}

