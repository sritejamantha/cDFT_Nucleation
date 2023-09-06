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

void read_one_resume_file(FILE *inp, double* w ) {
	
	int i, j,nn[Dim];

	double dr, di,dm[Dim];
	for (i=0; i<M; i++) {
	
		for (j=0; j<Dim; j++)
			fscanf(inp,"%lf ", &dm[j]);
		fscanf(inp, "%lf\n", &dr);
		
		w[i] = dr;
	}
}



void read_resume_files() {
	FILE *inp;
	char nm[20];
	int flag = 0; 
	inp = fopen("rst0.dat","r");
	if (inp!=NULL) {
		read_one_resume_file(inp, rhoK[0] ) ;
		   fclose(inp);
		   cout << "Read rst0.dat!" << endl;
		flag +=1;
	}
	
	inp = fopen("rst1.dat","r");
	if (inp!=NULL) {
		read_one_resume_file(inp, rhoK[1] ) ;
		   fclose(inp);
		   cout << "Read rst1.dat!" << endl;
		
		 flag +=1;
	
	}

	inp = fopen("rst2.dat","r");
	if (inp!=NULL) {
		read_one_resume_file(inp, rhoK[2] ) ;
		   fclose(inp);
		   cout << "Read rst2.dat!" << endl;
	
		 flag +=1;
	}


	int i, j; 
	double k2;

	if( flag == TNsp){
        for(i=0; i<TNsp; i++){
	    if(i==0)
	    	k2 = blkrho[0][1]/rhoK[i][0];
	    else if(i==1)
	    	k2 = blkrho[1][1]/rhoK[i][M-1];
	    else
	    	 k2 = blkrho[2][1]/rhoK[i][M-1];
	    for(j=0;j<M;j++){

	    		rhoK[i][j] *= k2;



	    }
	}
	}
	else{
	

	cout<<"can not readin full sets of rst files!"<<endl;
	}

}
