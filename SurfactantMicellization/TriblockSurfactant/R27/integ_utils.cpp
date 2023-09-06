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


double integ_trap( int np, double* dat ) {
  double sum = 0.0 ;
  int i ;

  if(np <2){


    cout<<"too few point! np: "<<np<<endl;
    exit(1);


  }

  if(np ==2){
    sum = (dat[0] +dat[1])/2.0;
  }
  else if(np ==3){
    sum = (1.*dat[0]+4.* dat[1] + 1.*dat[2])/3.;

  }
  else if(np ==4){
    sum = 3.*(1.*dat[0] + 3.*dat[1] + 3.*dat[2]+1*dat[3])/8;
  }
  else if(np ==5){

    sum = 2.*(7*dat[0] +32.*dat[1] +12.*dat[2] +32.*dat[3]+7.*dat[4])/45.;
  }
  else if(np ==6){

    sum = 5.*(19*dat[0] +75.*dat[1] +50.*dat[2] +50.*dat[3]+75.*dat[4] + 19*dat[5])/288.;
  }
  else{

    sum = integ_simpson(np,dat);
  }


  return sum ;

}


void interp_values( double* out_v, double cent_v, int np ) {
  double xrp, x1,x2,x3,y1,y2,y3 ;
  int i ;

  x1 = 0;
  x2 = np*dx[Dim-1];
  x3 = (np+1)*dx[Dim-1];

  y1 = cent_v;
  y2 = out_v[np];
  y3 = out_v[np+1];

  out_v[0] = y1;


  for(i=1 ;i<np; i++){

    xrp = i*dx[Dim-1];

    out_v[i] = y1*(xrp - x2)*(xrp-x3)/(x1- x2)/(x1-x3)+y2*(xrp - x1)*(xrp-x3)/(x2- x1)/(x2-x3)+ y3*(xrp - x1)*(xrp-x2)/(x3- x1)/(x3-x2);

  }

  
}



