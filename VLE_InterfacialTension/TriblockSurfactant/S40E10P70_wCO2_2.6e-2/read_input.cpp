#include "globals.h"

void read_input( void ) {

  FILE *inp ;
  int i,j;
  double d1 ;

  char tt[80] ;

  inp = fopen( "dft.input" , "r" ) ;

  fscanf( inp , "%lf " , &Tmpr  ) ;
  printf("current simulation temprt %lf\n",Tmpr);
  fgets( tt , 80 , inp ) ;

  for(i=0; i<TNsp; i++){
  	fscanf( inp , "%d " , &NK[i]  ) ;
  	printf("type %d mol length %d\n",i,NK[i]);
  }
  fgets( tt , 80 , inp ) ;
  
  Xa = double(NK[TNsp-3])/double(NK[TNsp-1]+NK[TNsp-2]+NK[TNsp-3]);
  Xb = double(NK[TNsp-2])/double(NK[TNsp-1]+NK[TNsp-2]+NK[TNsp-3]);
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
        fgets( tt , 80 , inp) ; 

  }
  
  for(i=0; i<TNsp; i++){
  	//fscanf( inp , "%lf %lf" ,&blkrho[i][0],&blkrho[i][1] ) ;
  	//printf("bulk density read in %lf %lf\n", blkrho[i][0],blkrho[i][1]);	
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

  tmp_lamb = lamb;
  
  fscanf( inp , "%d " , &otp_freq ) ;
  fgets( tt , 80 , inp) ;
  printf("The ttstep  = %d and otp freq = %d\n",ttstep,otp_freq);

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
	
 fclose( inp ) ;
  

  

  
  



}

void read_blkrho(double inp_pres){

int count = 0, j,k ; 


char c, tt[380];
FILE* fp;
 fp = fopen("blkrho.inp","r");

if(fp == NULL){

	cout<<"can not read blkrho.inp"<<endl;
	exit(1);

}
	for (c = getc(fp); c != EOF; c = getc(fp))
		if (c == '\n') 
			count = count + 1;

	fclose(fp);
   
	fp = fopen("blkrho.inp","r");

	double tmp_pres;
	int flag = 0; 

	for(k = 0; k<count ;k++){

	    fscanf(fp,"%lf ",&tmp_pres);
	    cout<<tmp_pres<<endl;
	    if( abs(tmp_pres - inp_pres)<1e-4){
	//	cout<<"test rho "<<tmp_pres<<" "<<inp_pres<<endl;
		fscanf(fp, "%lf ", &blkrho[0][0]);
		fscanf(fp, "%lf ", &blkrho[0][1]);
		fscanf(fp, "%lf ", &blkrho[1][0]);
		fscanf(fp, "%lf ", &blkrho[1][1]);
		fscanf(fp, "%lf ", &blkrho[2][0]);
		fscanf(fp, "%lf ", &blkrho[2][1]);
                fscanf(fp, "%lf ", &blkrho[3][0]);
                fscanf(fp, "%lf ", &blkrho[3][1]);
                fscanf(fp, "%lf ", &blkrho[4][0]);
                fscanf(fp, "%lf ", &blkrho[4][1]);
		//blkrho[1][0] = 1e-185;
	//	cout<<blkrho[0][1]<<endl;
		flag = 1; 
	    	break;
	    }
	    else{
		fgets(tt,380,fp);

	    }


	}


	fclose(fp);

        for(k=0;k<TNsp;k++)
        {
         printf("%.15e \t %.15e\n",blkrho[k][0],blkrho[k][1]);
        }

	if(flag != 1){
  		
		exit(1);

  	}

}
