#include "globals.h"

void strg_reparamrz()
{
   int i, j, k ;

   double *tmp_slic,tt_slen = 0;

   tmp_slic = (double*) calloc(M,sizeof(double));

    // calc distance    
    strg_s[0] = 0;

    for(i=1 ;i<N_strg; i++)
    {

	strg_s[i] = strg_s[i-1];

	for(j=0;j<TNsp-2; j++)
        {

	    for(k=0;k<M ;k++)
            {
		tmp_slic[k] = strg_newrhoK[j][i][k] - strg_newrhoK[j][i-1][k];

	    }//k

	   strg_s[i] +=  pow(integ_inpd(M,tmp_slic,tmp_slic)*dx[Dim-1],0.5);
	
	}//j

        for(k=0;k<M ;k++)
        {
          tmp_slic[k] = (strg_newrhoK[TNsp-2][i][k] + strg_newrhoK[TNsp-1][i][k]) - (strg_newrhoK[TNsp-2][i-1][k] + strg_newrhoK[TNsp-1][i-1][k]);

        }//k

        strg_s[i] +=  pow(integ_inpd(M,tmp_slic,tmp_slic)*dx[Dim-1],0.5);
	
//        printf("%d \t %.14f\n",i,strg_s[i]);
    }//i

//    exit(-1);
    //normlize string
//    printf("%d \t %.14f\n",0,strg_s[0]);
    for(i=1 ;i<N_strg; i++)
    {
	strg_s[i] /= strg_s[N_strg-1];
//        printf("%d \t %.14f\n",i,strg_s[i]);
    }

//   exit(-1);


    //intepolation 

    int tmp_indx=1;
    double cur_s;

    for(i=1 ;i<N_strg-1; i++)
    {
	
	strg_s_ind[i][0] = 0;
	strg_s_ind[i][1] = 0;
	cur_s = i*(1./(double(N_strg)-1));
	
	// find the weights in linear interp 

	for(j=tmp_indx;j< N_strg; j++ )
        {

		if( strg_s[j]>= cur_s)
                {

			strg_s_ind[i][0] = j-1;
			strg_s_ind[i][1] = j;

			strg_s_w[i][0] = strg_s[j] - cur_s;
			strg_s_w[i][1] = cur_s-strg_s[j-1] ;
			//tmp_indx = j;	 		
			//cout<<strg_s_ind[i][0]<<" "<<strg_s_ind[i][1]<<endl;
			break;
		}
	}

	if(strg_s_ind[i][1] == 0 )
        {
	  cout<<"can not find near point for str = "<<i<<" !"<<endl;
	  exit(1);
	}

	//linear interp 

	for(j=0; j<TNsp;j++)
        {
#pragma omp parallel for 
		for(k=0; k<M;k++)
                {

		    strg_rhoK[j][i][k] = strg_newrhoK[j][strg_s_ind[i][0]][k]*strg_s_w[i][0] + strg_newrhoK[j][strg_s_ind[i][1]][k]*strg_s_w[i][1];
		    strg_rhoK[j][i][k] /= ( strg_s_w[i][0] +  strg_s_w[i][1]);	
		}
	}//j=TNsp


     /*if(i==15)
     {
      for(k=0;k<M;k++)
      {
       printf("%d \t %.12e \t %.12e \t %.12e\n",k,strg_rhoK[0][i][k],strg_rhoK[1][i][k],strg_rhoK[2][i][k]+strg_rhoK[3][i][k]);
      }
     }*/
    }//i = sti
	    
 //  exit(-1);
   free(tmp_slic);

}
