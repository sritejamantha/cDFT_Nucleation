#include "globals.h"

void read_matdisp( void ) {

  FILE *inp ;
  int i,j;
  double d1 ;

  char tt[80] ;

  inp = fopen( "dft_mat.dat" , "r" ) ;
  
  if ( inp == NULL ){
	printf("Can not find dft_mat.dat!\n");
	exit(1);
  }


  for(i=0; i<7; i++){
  	for(j=0;j<6;j++){
		fscanf( inp , "%lf" , &d1);
		if(j<3){
		  Adisp[i][j] = d1;
		}
		else{
  		  Bdisp[i][j-3] = d1;

		}
		//printf("readin mat val: %lf\n",d1);
  	}

	fgets( tt , 80 , inp ) ;
  }
  fclose( inp ) ;
  printf("dft_mat.dat loaded!\n");
}
