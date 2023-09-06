#include "globals.h"

void read_input( void ) {

  FILE *inp ;
  int i,j;
  double d1 ;

  char tt[80] ;

  inp = fopen( "dft.input" , "r" ) ;

  fscanf( inp , "%lf %d %lf" , &Tmpr, &N_strg, &cut_r  ) ;
  printf("current simulation temprt %lf\n",Tmpr);
  printf("with N string points: %d\n",N_strg);
  printf("with cut r at: %lf\n",cut_r);

  fgets( tt , 80 , inp ) ;

  for(i=0; i<TNsp; i++){
  	fscanf( inp , "%d " , &NK[i]  ) ;
  	printf("type %d mol length %d\n",i,NK[i]);
  }
  fgets( tt , 80 , inp ) ;

  Xa = double(NK[TNsp-2])/double(NK[TNsp-1]+NK[TNsp-2]);
  printf("Fraction of Backbone beads per chain=%lf\n",Xa);

  fscanf( inp,"%d ",&Nbr);
  printf("Number of PEO Branches per chain=%d\n",Nbr);
  fgets( tt , 80 , inp ) ;

  for(i=0; i<TNsp; i++){
  	fscanf( inp , "%lf " , &rsigma[i]  ) ;
	printf("type %d sigma  %lf\n",i,rsigma[i]);	
  
  }
  fgets( tt , 80 , inp ) ;

  for(i=0; i<TNsp; i++){
  	fscanf( inp , "%lf " , &epsln[i]  ) ;
 	printf("type %d epsln  %lf\n",i,epsln[i]);	
  
  }
  fgets( tt , 80 , inp ) ;

  for(i=0; i<crsTNsp; i++){
  	fscanf( inp , "%lf %lf " , &kijpara[i][0],&kijpara[i][1]  ) ;
  	printf("crskij  %lf %lf\n",kijpara[i][0],kijpara[i][1]);	
  

  }
  fgets( tt , 80 , inp) ;
  
  for(i=0; i<TNsp; i++){
  	fscanf( inp , "%lf %lf" ,&blkrho[i][0],&blkrho[i][1] ) ;
  	printf("bulk density read in %lf %lf\n", blkrho[i][0],blkrho[i][1]);	
        fgets( tt , 80 , inp) ;

  }
  fgets( tt , 80 , inp) ;

  for(i=0; i<Dim; i++){
  	fscanf( inp , "%lf " , &L[i]  ) ;
  	printf("dim %d w length %lf\n",i,L[i]);
  }
  fgets( tt , 80 , inp);

  for(i=0; i<Dim; i++){
  	fscanf( inp , "%lf " , &dx[i]  ) ;
  	printf("dim %d w dx %lf\n",i,dx[i]);
  }

  fgets( tt , 80 , inp) ;
  fgets( tt , 80 , inp) ;

  fscanf( inp , "%d %lf " , &ttstep, &lamb  ) ;
  fgets( tt , 80 , inp) ;
  
  fscanf( inp , "%d %d %d" , &otp_freq, &otp_freq2, &strg_initstp ) ;
  fgets( tt , 80 , inp) ;
  printf("The ttstep  = %d and otp freq = %d, otp freq2 = %d, strg_init = %d\n",ttstep,otp_freq, otp_freq2,strg_initstp);

  fgets( tt , 80 , inp) ;
  for(i=0; i<TNsp; i++){
           fscanf( inp , "%lf" , &amu[i]  ) ;

	   printf("type %d w amu: %lf\n",i,amu[i]);
   }
   
   fgets( tt , 80 , inp) ;
   fgets( tt , 80 , inp) ;
   fscanf( inp, "%d" , &wall_para  ) ;
   cout<<"wall para "<<wall_para<<endl;
   //printf("wall para:  %d  \n",wall_para);
   fgets( tt , 80 , inp) ;
   fgets( tt , 80 , inp) ;
   fscanf( inp, "%d" , &rst_ind ) ;
   cout<<"restart ind "<<rst_ind<<endl;

 fclose( inp ) ;
  

  

  
  



}
