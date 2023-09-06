#include "globals.h"
double integ_simpson(int , double* );
void get_weightn(void);
void get_dfhsdn(void);
void get_dfchdn();
void get_dfdispl1dn();
void get_dfdispl2dn();
void get_dfnloc();
int stack(int*);
void unstack(int, int*  );
void unstack_ext(int, int*  );


void dfree(){

  int i,j,l,nn[Dim],k ;

  double sigv12,sigv13,sigv23,tmprhoN3,tmpNN3; 

  get_weightn();
  
  get_dfhsdn();
  
  get_dfchdn();
 

  get_dfdispl1dn();

  get_dfdispl2dn();


  for(i=0;i<TNsp;i++)
  {
   
#pragma omp parallel for private(k,nn)
   for(j=0;j<extM;j++)
   {
	
	unstack_ext(j,nn);
  	if((nn[Dim-1] >= padN ) and  (nn[Dim-1]< (Nx[Dim-1]+padN)) )
	{
		nn[Dim-1] -= padN;
		k = stack(nn);
		extdflocdn0[i][j] = dfdispl1dn0[i][k]/Nsig[i]+dfdispl2dn0[i][k]/Nsig[i] ;
		extdflocdn3[i][j] =  PI*(dfhsdn3[i][k] + dfchdn3[i][k]+dfdispl1dn3[i][k]+dfdispl2dn3[i][k])*pow(sigma[i],3)/pow(Nsig[i],3);
	}
	else if (nn[Dim-1]>= (Nx[Dim-1] +padN))
	{
		extdflocdn0[i][j] = dfdispl1dn0[i][M-1]/Nsig[i]+dfdispl2dn0[i][M-1]/Nsig[i];
		extdflocdn3[i][j] =  PI*(dfhsdn3[i][M-1] + dfchdn3[i][M-1]+dfdispl1dn3[i][M-1]+dfdispl2dn3[i][M-1])*pow(sigma[i],3)/pow(Nsig[i],3);

	}
	else 
	{
	        extdflocdn0[i][j] = dfdispl1dn0[i][0]/Nsig[i] + dfdispl2dn0[i][0]/Nsig[i];  
		extdflocdn3[i][j] =  PI*(dfhsdn3[i][0] + dfchdn3[i][0]+dfdispl1dn3[i][0]+dfdispl2dn3[i][0])*pow(sigma[i],3)/pow(Nsig[i],3);

	}
   }//j


  }//i
  //exit(-1);

	int ind1, ind2;


  if(Dim ==1 )
  {
  	for(i=0;i<TNsp;i++)
	{
#pragma omp parallel for private(k,l,ind1,ind2)
     		for(j=0;j<M;j++)
		{
		   int tid = omp_get_thread_num() ; 
		   ind1 = j+padN - sigNhf[i];
		   ind2 = j+padN + sigNhf[i]; 
		   dfloc[i][j] = 0; 

		
		   for(l=0 ;l<(2*sigNhf[i]+1);l++)
		   {
			tmpNhf[tid][l] = extdflocdn0[i][ind1+l];
		   }//l == 2sigN+1    		   
		   dfloc[i][j] = integ_simpson(2*sigNhf[i]+1, tmpNhf[tid])*dx[Dim-1]; //dx[Dim-1]*PI*sigma[i]*sigma[i]/Nsig[i]*integ_simpson(2*sigNhf[i]+1, tmpNhf[i]);

		   for(l=0 ;l<(2*sigNhf[i]+1);l++)
		   {
			tmpNhf[tid][l] = extdflocdn3[i][ind1+l]*knltz1[i][l];
		   }//l == 2sigN+1    		   
		   dfloc[i][j] += dx[Dim-1]*integ_simpson(2*sigNhf[i]+1, tmpNhf[tid]);
		   //if(i==1)
		   	//cout<<"dfloc "<<dfloc[i][j]<<endl;
		}//j 

  	}//i
  }//Dim ==1 

  //exit(-1);

  get_dfnloc();


  tmpNN3 = NK[TNsp-2]+NK[TNsp-1];
  for(i=0;i<M;i++)
  {

	ttF[i] =0;
        tmprhoN3=rhoK[TNsp-2][i]+rhoK[TNsp-1][i];
	for(j=0;j<TNsp;j++)
        {
                if(j==TNsp-2)
                {
                 ttF[i] += Xar[i]*tmprhoN3/tmpNN3*(log(tmprhoN3/tmpNN3 + 1.e-30)-1.);
                 ttF[i] += dfdispndr[j][i]*rhoK[j][i];

                 dfttdrho[j][i] = Xar[i]*(dfloc[j][i] + dfdispndr[j][i]);
                }
  
                else if(j==TNsp-1)
                {
                 ttF[i] += (1.0 - Xar[i])*tmprhoN3/tmpNN3*(log(tmprhoN3/tmpNN3 + 1.e-30)-1.);
                 ttF[i] += dfdispndr[j][i]*rhoK[j][i];

                 dfttdrho[j][i] = (1.0 - Xar[i])*(dfloc[j][i] + dfdispndr[j][i]);
                }

                else
                {
		 ttF[i] += rhoK[j][i]/NK[j]*(log(rhoK[j][i]/NK[j]+1.e-30)-1.);
		 ttF[i] += dfdispndr[j][i]*rhoK[j][i];
	
		 dfttdrho[j][i] = dfloc[j][i] + dfdispndr[j][i];
                }

	}

	ttF[i] += Fhs[i];

	ttF[i] += Fch[i];

	ttF[i] += Fdispl1[i] + Fdispl2[i];

	

  }

//  write_grid_data("dfloc1.dat",dfttdrho[0]);

  //write_grid_data("dfloc2.dat",dfttdrho[1]);
 

//  exit(-1);
}


void get_dfhsdn(){

  int i,j,k,nn[Dim];
  double fhspart[6] ,tmp_lg, tmp_1m, tmp_dd[3],tmp_eta;


#pragma omp parallel for private(tmp_lg,tmp_1m,tmp_eta,tmp_dd,i,k,nn,fhspart)
    for(j=0;j<M;j++)
    {
	Fhs[j] = 0;

	unstack(j,nn);
	for(i =0 ;i<3; i++)
	    tmp_dd[i] = 0;

  	for(i =0 ;i<TNsp; i++)
        {

	    tmp_dd[0] += wn3[i][j]/sigma[i];
	    tmp_dd[1] += wn3[i][j]/sigma[i]/sigma[i]/PI;
	    tmp_dd[2] += wn3[i][j]/pow(sigma[i],3)/PI*6.;

	}
	tmp_eta = wn3[TNsp][j];
	if(!(tmp_eta>0))
        {
	    Fhs[j] = 0 ; 
	    for(i =0 ;i<TNsp; i++)
            {
		 dfhsdn3[i][j] = 0; 	
	    }
	    continue;
	}
	
	
	tmp_1m = 1.-wn3[TNsp][j];
	       
	tmp_lg = log(tmp_1m);

	Fhs[j] = -1*tmp_lg*tmp_dd[2] + 18*tmp_dd[0]*tmp_dd[1]/tmp_1m + 6.0*pow(tmp_dd[0],3)*(tmp_lg/tmp_eta/tmp_eta + 1./tmp_eta/tmp_1m/tmp_1m)/PI;

        //printf("%d \t %.15e\n",j,Fhs[j]);

	fhspart[0] = tmp_dd[2]/tmp_1m;
	fhspart[1] = 18*tmp_dd[0]/tmp_1m;
	fhspart[2] = 18*tmp_dd[1]/tmp_1m;
	fhspart[3] = 18*tmp_dd[1]*tmp_dd[0]/tmp_1m/tmp_1m;
   	fhspart[4] = 18*tmp_dd[0]*tmp_dd[0]*(tmp_lg/pow(tmp_eta,2.)+ 1./(tmp_eta*tmp_1m*tmp_1m))/PI;
	fhspart[5] = 6*pow(tmp_dd[0],3)*(2*tmp_lg/pow(tmp_eta,3.) + (2. -5*tmp_eta+ tmp_eta*tmp_eta)/pow(tmp_eta,2)/pow(tmp_1m,3) )/PI;
	
	for(i =0 ;i<TNsp; i++)
        {
	    dfhsdn3[i][j] = fhspart[0] - tmp_lg/pow(sigma[i],3)/PI*6.+fhspart[1]/sigma[i]/sigma[i]/PI+ fhspart[2]/sigma[i] + fhspart[3] + fhspart[4]/sigma[i] - fhspart[5];
        }
    }//j 

    //write_grid_data("dfhsdn31.dat",dfhsdn3[0]);
    //write_grid_data("Fhs.dat",Fhs);
    //exit(-1);
}

void get_dfchdn(){

  int i,j,k,nn[Dim];
  double fhspart[6] ,tmp_gg[TNsp],tmp_lg, tmp_1m, tmp_dd[4],tmp_eta;
  double tmp_gab=0.0,tmp_lgab,tmp_rsigab=0.0,tmp_sigab=0.0,sigv[crsTNsp],tmprhoN3,tmpNN3,tmpch_cf3b,tmpch_cf3ab;
  double tmp_alpha,tmp_beta,tmp_dd0,tmp_dd1,tmp_dd2;
 
  k=0;
  for(i=0;i<TNsp-1;i++)
  {
   for(j=i+1;j<TNsp;j++)
   {
    sigv[k]=(rsigma[i] + rsigma[j])/2.0;
    k=k+1;
   }
  }
  tmp_rsigab = (sigma[TNsp-1]*sigma[TNsp-2])/(sigma[TNsp-1]+sigma[TNsp-2]);
  tmp_sigab = tmp_rsigab*(1. - 0.12*exp(-3.0*crsepsln[crsTNsp-1]/Tmpr));
//  tmp_sigab = sigma[TNsp-1];
  tmpNN3 = (double)(NK[TNsp-2]+NK[TNsp-1]);
#pragma omp parallel for private(tmp_gg,tmp_lg,tmp_1m,tmp_eta,tmp_dd,i,k,nn,tmp_gab,tmprhoN3,tmpch_cf3b,tmpch_cf3ab,tmp_lgab,tmp_alpha,tmp_beta,tmp_dd0,tmp_dd1,tmp_dd2)
    for(j=0;j<M;j++)
    {
	Fch[j] = 0;
	
	tmp_eta = wn3[TNsp][j];
	tmp_1m = 1.-wn3[TNsp][j];

        tmprhoN3 = wn0[TNsp-2][j]+wn0[TNsp-1][j];
       

	for(i =0 ;i<4; i++)
	    tmp_dd[i] = 0;
	
	for(i =0 ;i<TNsp; i++)
        {
		tmp_dd[3] += wn3[i][j]/sigma[i];

	}

  	for(i =0 ;i<TNsp; i++)
        {
	    
             tmp_gg[i] = 1./tmp_1m + 3*tmp_dd[3]*sigma[i]/2.0/pow(tmp_1m,2) + pow(tmp_dd[3]*sigma[i],2)/2./pow(tmp_1m,3);
        }

        tmp_gab =  1./tmp_1m + 3*tmp_dd[3]*tmp_sigab/2.0/pow(tmp_1m,2) + pow(tmp_dd[3]*tmp_sigab,2)/2./pow(tmp_1m,3);


        tmp_alpha = wn3[TNsp-2][j]/pow(sigma[TNsp-2],3) + wn3[TNsp-1][j]/pow(sigma[TNsp-1],3);
        tmp_beta = (1. - double(NK[TNsp-2]))*log( tmp_gg[TNsp-2]) + (1. - (double)NK[TNsp-1])*log( tmp_gg[TNsp-1]) - (1. - (double)Nbr)*log( tmp_gg[TNsp-1]) - (double)Nbr*log(tmp_gab);
        tmp_beta = tmp_beta/tmpNN3;

  	for(i =0 ;i<TNsp-2; i++)
        {

	    tmp_dd[0] += (1.-double(NK[i]))/double(NK[i])*wn3[i][j]/sigma[i]/tmp_gg[i];
	    tmp_dd[1] +=(1.-double(NK[i]))/double(NK[i])*wn3[i][j]/sigma[i]/sigma[i]/PI/tmp_gg[i];
	    tmp_dd[2] += (1.-double(NK[i]))/double(NK[i])*wn3[i][j]/pow(sigma[i],3)/PI/tmp_gg[i]*6.;
	
	}

       tmp_dd2 = (1.0 - double(NK[TNsp-2]))/tmp_gg[TNsp-2] + (1.0 - double(NK[TNsp-1]))/tmp_gg[TNsp-1] - (1.0 - double(Nbr))/tmp_gg[TNsp-1] - double(Nbr)/tmp_gab;
   tmp_dd2 = 6.0*tmp_alpha*tmp_dd2/tmpNN3/PI;
   tmp_dd[2] = tmp_dd[2] + tmp_dd2;


       tmp_dd1 = sigma[TNsp-2]*(1.0 - double(NK[TNsp-2]))/tmp_gg[TNsp-2] +  sigma[TNsp-1]*(1.0 - double(NK[TNsp-1]))/tmp_gg[TNsp-1] - sigma[TNsp-1]*(1.0 - double(Nbr))/tmp_gg[TNsp-1] - tmp_sigab*double(Nbr)/tmp_gab;
   tmp_dd1 = tmp_alpha*tmp_dd1/tmpNN3/PI;
   tmp_dd[1] = tmp_dd[1] + tmp_dd1;

        tmp_dd0 = pow(sigma[TNsp-2],2.0)*(1.0 - double(NK[TNsp-2]))/tmp_gg[TNsp-2] +  pow(sigma[TNsp-1],2.0)*(1.0 - double(NK[TNsp-1]))/tmp_gg[TNsp-1] - pow(sigma[TNsp-1],2.0)*(1.0 - double(Nbr))/tmp_gg[TNsp-1] - pow(tmp_sigab,2.0)*double(Nbr)/tmp_gab;
   tmp_dd0 = tmp_alpha*tmp_dd0/tmpNN3;
   tmp_dd[0] = tmp_dd[0] + tmp_dd0;


	       
        Fch[j]=-(tmprhoN3/tmpNN3)*( (1. - double(Nbr))*log(tmp_gg[TNsp-1]) + double(Nbr)*log(tmp_gab));	
	for(i =0; i<TNsp; i++)
        {

	    tmp_lg = log(tmp_gg[i]);
 

	    Fch[j] += 6.*(1.-double(NK[i]))/double(NK[i])*wn3[i][j]*tmp_lg/PI/pow(sigma[i],3);
	    
            if(i<TNsp-2)
            {	    
	     dfchdn3[i][j] = 6*(1.-double(NK[i]))/double(NK[i])*tmp_lg/PI/pow(sigma[i],3);
	    }
            else
            {
             dfchdn3[i][j] = 6*tmp_beta/PI/pow(sigma[i],3);
            }

            dfchdn3[i][j] += 1./tmp_1m/tmp_1m*tmp_dd[2];
	    dfchdn3[i][j] +=(9./tmp_1m/tmp_1m/sigma[i]+ 18.*tmp_dd[3]/pow(tmp_1m,3) )*tmp_dd[1];
	    dfchdn3[i][j] +=(6.*tmp_dd[3]/pow(tmp_1m,3)/PI/sigma[i]  +  9.*pow(tmp_dd[3],2)/pow(tmp_1m,4)/PI)*tmp_dd[0];
	

	}

    //    printf("%d \t %.15e\n",j,Fch[j]);
    }//j 

    //exit(-1);
    //write_grid_data("Fch.dat",Fch);
   //write_grid_data("dfchdn31.dat",dfchdn3[1]);
}

void get_dfdispl1dn(){

  int i,j,l,k,nn[Dim];
  double loc_tt,loc_x[TNsp],tmp_I2,loc_nbar,dI1dn0[TNsp],tmp_I1,dmn0[TNsp],loc_da[7],loc_a[7],tmp_gg[TNsp],tmp_lg, tmp_pars, tmp_n0,tmp_dd[4],tmp_eta;

  double sigv[crsTNsp],tmpsig[TNsp][TNsp],tmpeps[TNsp][TNsp];
 
  k=0;
  for(i=0;i<TNsp-1;i++)
  {
   for(j=i+1;j<TNsp;j++)
   {
    sigv[k]=(rsigma[i] + rsigma[j])/2.0;
    tmpsig[i][j] = sigv[k];
    tmpsig[j][i]=sigv[k];
    tmpeps[i][j] = crsepsln[k];
    tmpeps[j][i] = crsepsln[k];
    k=k+1;
   }
  }


#pragma omp parallel for private(tmp_pars,tmp_gg,tmp_dd,tmp_I1,tmp_eta,dmn0,dI1dn0,loc_a,loc_da,loc_tt,loc_x,loc_nbar,k,l,i)
  for(j=0; j<M; j++)
  {
	loc_tt=0;

	loc_nbar = 0; 
	
	tmp_n0 = 0;
	
	tmp_eta = wn3[TNsp][j];

	for(i=0;i<TNsp;i++)
        {
		tmp_n0 += wn0[i][j];
                if(i==TNsp-2)
        	{
         	  loc_x[i] =Xar[j]*wn0[i][j]/NK[i];
        	}
        	else if(i==TNsp-1)
        	{
         	  loc_x[i] = (1.0 - Xar[j])*wn0[i][j]/NK[i];
        	}
        	else
        	{	
		  loc_x[i] =wn0[i][j]/NK[i];
                }
		loc_tt += loc_x[i];
		tmp_gg[i] =0;

	}//i
   
        if(loc_tt > 0 )
        {
		for(i=0;i<TNsp;i++)
                {
			loc_x[i] /= loc_tt;
                        if(i==TNsp-2)
                	{
                 	  loc_nbar += NK[i]*loc_x[i]/Xar[j];
               		}
                	else if(i==TNsp-1)
                	{
                 	  loc_nbar += NK[i]*loc_x[i]/(1.0-Xar[j]);
                	}
                	else	
                	{
			  loc_nbar += NK[i]*loc_x[i];
			}
		}//i
		
	}
	else{
	  Fdispl1[j] = 0;
	   
	  for(i=0;i<TNsp;i++)
          {
	 	dfdispl1dn0[i][j] = 0;
		dfdispl1dn3[i][j] = 0;
	  }
	  continue;
	 //cout<<"no density @ "<<j<<" !"<<endl;
	 //exit(1);
	}
	//get dnbar/dn0	
	for(i=0;i<TNsp;i++)
        {
                if(i==TNsp-2)
        	{
         	  dmn0[i] = -1.*Xar[j]*tmp_n0/loc_tt/loc_tt/NK[i] + 1./loc_tt;
        	}
        	else if(i==TNsp-1)
        	{
         	  dmn0[i] = -1.*(1.0-Xar[j])*tmp_n0/loc_tt/loc_tt/NK[i] + 1./loc_tt;
        	}
        	else
        	{
		  dmn0[i] = -1.*tmp_n0/loc_tt/loc_tt/NK[i]+ 1./loc_tt;
		}
	
	}//i,TNsp
	
	//get coeff a, I1 ,Fdsipl1 and da/dn0
	tmp_I1 = 0; 
	Fdispl1[j] = 0;
	tmp_dd[0] = 0; //dI1dnbar
	tmp_dd[1] = 0; //dI1dn3
	tmp_dd[2]= 0 ; //factors ass w  dI1dn0 & dI1dn3
	// tmp_gg: factors in dfdispldn0 

	for(i=0;i<7;i++)
        {
	   	loc_a[i] = Adisp[i][0] + +(loc_nbar-1.)/loc_nbar*Adisp[i][1]+(loc_nbar-1.)*(loc_nbar-2)/loc_nbar/loc_nbar*Adisp[i][2]; 
	   	
	   	loc_da[i] = 1./loc_nbar/loc_nbar*Adisp[i][1]+ (3./loc_nbar/loc_nbar- 4./pow(loc_nbar,3))*Adisp[i][2];
	   	tmp_I1 += loc_a[i]*pow(tmp_eta,i);
 		
		tmp_dd[0] += loc_da[i]*pow(tmp_eta,i);
	
		if(i>0)
			tmp_dd[1] += loc_a[i]*i*pow(tmp_eta,i-1);
	
		
	}

	for(i=0;i<TNsp;i++)
        {
             tmp_pars = wn0[i][j]*wn0[i][j]*epsln[i]*rsigma3[i]/Tmpr;
             Fdispl1[j] += tmp_pars;
             tmp_dd[2] += tmp_pars;

             dI1dn0[i] = tmp_dd[0]*dmn0[i];
                
        }//i TNsp


	for(i=0;i<TNsp-1;i++)
	{
		for(k=i+1;k<TNsp;k++)
		{
			
			tmp_pars = 2.*wn0[i][j]*wn0[k][j]*tmpeps[i][k]*pow(tmpsig[i][k],3)/Tmpr;
			Fdispl1[j] += tmp_pars;

			tmp_dd[2] += tmp_pars;
		}

	
	}//i TNsp
	
       

	for(i=0;i<TNsp;i++)
        {
		for(k=0;k<TNsp;k++)
		{
			 if(i== k)
			 	tmp_gg[i] += wn0[i][j]*epsln[i]*rsigma3[i]/Tmpr;
			 else
			 {
			 	tmp_gg[i] +=  wn0[k][j]*tmpeps[i][k]*pow(tmpsig[i][k],3)/Tmpr;
			 }
 		}
	}

	Fdispl1[j] *= -2.*PI*tmp_I1;
       // printf("%d \t %.15e\n",j,Fdispl1[j]);

	for(i=0;i<TNsp;i++)
        {
	  dfdispl1dn0[i][j] = -4.*PI*tmp_I1*tmp_gg[i]-2.*PI*dI1dn0[i]*tmp_dd[2]; 
	  dfdispl1dn3[i][j] = -2.*PI*tmp_dd[1]*tmp_dd[2];
	}
 

  }//j,M

  //exit(-1);

//  write_grid_data("dfdispl1dn01.dat",dfdispl1dn0[1]);
//  write_grid_data("dfdispl1dn31.dat",dfdispl1dn3[1]);
}


void get_dfdispl2dn()
{

  int i,j,l,k,nn[Dim];
  double dmn0[TNsp],loc_tt,loc_x[TNsp],tmp_I2,loc_nbar,dMdn3,dMdn0[TNsp],dI2dn0[TNsp],tmp_I1,loc_da[7],loc_a[7],tmp_1m ,tmp_gg[TNsp],tmp_lg, tmp_pars, tmp_n0,tmp_dd[4],tmp_eta;


  double tmp_M,sigv[crsTNsp],tmpsig[TNsp][TNsp],tmpeps[TNsp][TNsp];

  k=0;
  for(i=0;i<TNsp-1;i++)
  {
   for(j=i+1;j<TNsp;j++)
   {
    sigv[k]=(rsigma[i] + rsigma[j])/2.0;
    tmpsig[i][j] = sigv[k];
    tmpsig[j][i]=sigv[k];
    tmpeps[i][j] = crsepsln[k];
    tmpeps[j][i] = crsepsln[k];
    k=k+1;
   }
  }

#pragma omp parallel for private(tmp_M,tmp_pars,tmp_gg,tmp_dd,tmp_I2,tmp_eta,tmp_1m,dmn0,dI2dn0,dMdn3,dMdn0,loc_a,loc_da,loc_tt,loc_x,loc_nbar,i,k,l)
  for(j=0; j<M; j++)
  {

	loc_tt=0;

	loc_nbar = 0; 
	
	tmp_n0 = 0;
	
	tmp_eta = wn3[TNsp][j];
	tmp_1m = 1. - tmp_eta;

	for(i=0;i<TNsp;i++)
        {
		tmp_n0 += wn0[i][j];	
		if(i==TNsp-2)
                {
                  loc_x[i] =Xar[j]*wn0[i][j]/NK[i];
                }
                else if(i==TNsp-1)
                {
                  loc_x[i] = (1.0 - Xar[j])*wn0[i][j]/NK[i];
                }
                else
                {
		  loc_x[i] =wn0[i][j]/NK[i];
		}
		loc_tt += loc_x[i];
		tmp_gg[i] =0;

	}//i
   
        if(loc_tt > 0 )
	{
		for(i=0;i<TNsp;i++)
		{
		        	
			loc_x[i] /= loc_tt;
			if(i==TNsp-2)
                        {
                          loc_nbar += NK[i]*loc_x[i]/Xar[j];
                        }
                        else if(i==TNsp-1)
                        {
                          loc_nbar += NK[i]*loc_x[i]/(1.0-Xar[j]);
                        }
                        else
                        {
			  loc_nbar += NK[i]*loc_x[i];
			}

		}//i
		
	}
	else
	{
	 
	 Fdispl2[j] = 0;
	   
	 for(i=0;i<TNsp;i++)
	 {
	 	dfdispl2dn0[i][j] = 0;
		dfdispl2dn3[i][j] = 0;
	 }
	 continue;
	 
	 
	 
	 //cout<<"no density @ "<<j<<" !"<<endl;
	 //exit(1);
	}
	//get dnbar/dn0	
	for(i=0;i<TNsp;i++)
        {
		if(i==TNsp-2)
                {
                  dmn0[i] = -1.*Xar[j]*tmp_n0/loc_tt/loc_tt/NK[i] + 1./loc_tt;
                }
                else if(i==TNsp-1)
                {
                  dmn0[i] = -1.*(1.0-Xar[j])*tmp_n0/loc_tt/loc_tt/NK[i] + 1./loc_tt;
                }
                else
                {
		  dmn0[i] = -1.*tmp_n0/loc_tt/loc_tt/NK[i]+ 1./loc_tt;
		}
	
	}//i,TNsp
	
	//get coeff a, I1 ,Fdsipl1 and da/dn0
	tmp_I2 = 0; 
	Fdispl2[j] = 0;
	tmp_dd[0] = 0; //factors in dMdn0 (dMdnbar)
	tmp_dd[1] = 0; //factors in dI2dn0
	tmp_dd[2]= 0 ; //dI2dn3
	tmp_dd[3] = 0; //factors in dfdisp2dn0 and dfdisp2dn3
	// tmp_gg: factors(2) in dfdispld2n0 


	for(i=0;i<7;i++)
	{
	   	loc_a[i] = Bdisp[i][0] +(loc_nbar-1.)/loc_nbar*Bdisp[i][1]+(loc_nbar-1.)*(loc_nbar-2)/loc_nbar/loc_nbar*Bdisp[i][2]; 
	   	
	   	loc_da[i] = 1./loc_nbar/loc_nbar*Bdisp[i][1]+ (3./loc_nbar/loc_nbar- 4./pow(loc_nbar,3))*Bdisp[i][2];
	   	tmp_I2 += loc_a[i]*pow(tmp_eta,i);
 		
	
		
		tmp_dd[1] += loc_da[i]*pow(tmp_eta,i);
		
		if(i>0)
			tmp_dd[2] += i*loc_a[i]*pow(tmp_eta,i-1);
	}


	
	tmp_M = 1.+loc_nbar*(8.*tmp_eta - 2.*tmp_eta*tmp_eta)/pow(tmp_1m,4) + (1.- loc_nbar)*(20.*tmp_eta - 27*tmp_eta*tmp_eta + 12.*pow(tmp_eta,3)-2.*pow(tmp_eta,4.))/tmp_1m/tmp_1m/pow(2.-tmp_eta,2);

	tmp_dd[0] = (8.*tmp_eta - 2.*tmp_eta*tmp_eta)/pow(tmp_1m,4) - (20.*tmp_eta - 27.*tmp_eta*tmp_eta+12.*pow(tmp_eta,3)-2.*pow(tmp_eta,4))/tmp_1m/tmp_1m/pow(2.-tmp_eta,2);
	
	dMdn3 = loc_nbar*(8.+20.*tmp_eta-4.*tmp_eta*tmp_eta)/pow(tmp_1m,5) + (1.- loc_nbar)*(2*pow(tmp_eta,3)+12.*tmp_eta*tmp_eta-48.*tmp_eta+40)/pow(tmp_1m,3)/pow(2.-tmp_eta,3);
	
	for(i=0;i<TNsp;i++)
	{
		tmp_pars = wn0[i][j]*wn0[i][j]*pow(epsln[i]/Tmpr,2)*rsigma3[i];
		Fdispl2[j] += tmp_pars;
		tmp_dd[3] += tmp_pars;
		
		dMdn0[i] = tmp_dd[0] * dmn0[i]; 
		dI2dn0[i] = tmp_dd[1]*dmn0[i];
			
	 }//i TNsp

	 for(i=0;i<TNsp-1;i++)
         {
                for(k=i+1;k<TNsp;k++)
                {

                    tmp_pars = 2.*wn0[i][j]*wn0[k][j]*pow(tmpeps[i][k]/Tmpr,2)*pow(tmpsig[i][k],3);
                    Fdispl2[j] += tmp_pars;

                    tmp_dd[3] += tmp_pars;
                
                }


         }//i TNsp

         for(i=0;i<TNsp;i++)
	 {
		for(k=0;k<TNsp;k++)
		{
			 if(i== k)
			 	tmp_gg[i] += wn0[i][j]*pow(epsln[i]/Tmpr,2)*rsigma3[i];
			 else
                         {
			 	tmp_gg[i] +=  wn0[k][j]*pow(tmpeps[i][k]/Tmpr,2)*pow(tmpsig[i][k],3);
			 }
		}
	}//i /TNsp

	Fdispl2[j] *= -1.*PI*tmp_I2/tmp_M*loc_nbar;
  //      printf("%d \t %.15e\n",j,Fdispl2[j]);
       
	for(i=0;i<TNsp;i++)
	{

	  dfdispl2dn0[i][j]  = -2.*PI/tmp_M*tmp_I2*loc_nbar*tmp_gg[i] + PI/tmp_M/tmp_M*dMdn0[i]*tmp_I2*loc_nbar*tmp_dd[3] -PI/tmp_M*dI2dn0[i]*loc_nbar*tmp_dd[3] - PI/tmp_M*tmp_I2*dmn0[i]*tmp_dd[3];
	  
	  
	  dfdispl2dn3[i][j] = PI/tmp_M/tmp_M*dMdn3*tmp_I2*loc_nbar*tmp_dd[3]- PI/tmp_M*tmp_dd[2]*loc_nbar*tmp_dd[3];
	}

	
  }//j,M

//  exit(-1);
//  write_grid_data("dfdispl2dn01.dat",dfdispl2dn0[1]);
 // write_grid_data("dfdispl2dn31.dat",dfdispl2dn3[1]);

}







void get_weightn(){

  int i, j ,k, l,m,ind1,ind2,nn[Dim];

  for(i=0;i<TNsp;i++){

#pragma omp parallel for private(k,nn)
   for(j=0;j<extM;j++){
	
	unstack_ext(j,nn);
  	if((nn[Dim-1] >= padN ) and  (nn[Dim-1]< (Nx[Dim-1]+padN)) ){
		nn[Dim-1] -= padN;
		k = stack(nn);
		extrhok[i][j] = rhoK[i][k];
	}
	else if (nn[Dim-1]>= (Nx[Dim-1] +padN)){
		extrhok[i][j] = rhoK[i][M-1];
	}
	else 
	{
	         extrhok[i][j] = rhoK[i][0];  
	}
   }//j

  }//i

 
  
  if(Dim ==1 ){
  	for(i=0;i<TNsp;i++){
#pragma omp parallel for private(k,l,ind1,ind2)
     		for(j=0;j<M;j++){
		   int tid = omp_get_thread_num() ; 
		   ind1 = j+padN - sigNhf[i];
		   ind2 = j+padN + sigNhf[i]; 
		   wn2[i][j] = 0; 
		   wn3[i][j] = 0 ;
		   //wnv2[i][Dim-1][j] = 0; 
		
		   for(l=0 ;l<(2*sigNhf[i]+1);l++){
			tmpNhf[tid][l] = extrhok[i][ind1+l];
		   }//l == 2sigN+1    		   
		   wn2[i][j] = integ_simpson(2*sigNhf[i]+1, tmpNhf[tid])*sigma[i]*sigma[i]/Nsig[i]*PI *dx[Dim-1]; //dx[Dim-1]*PI*sigma[i]*sigma[i]/Nsig[i]*integ_simpson(2*sigNhf[i]+1, tmpNhf[i]);
		   wn0[i][j] =  wn2[i][j]/PI/sigma[i]/sigma[i];

		   for(l=0 ;l<(2*sigNhf[i]+1);l++){
			tmpNhf[tid][l] = extrhok[i][ind1+l]*knltz1[i][l];
		   }//l == 2sigN+1    		   
		   wn3[i][j] = dx[Dim-1]*PI*pow(sigma[i],3)/pow(Nsig[i],3)*integ_simpson(2*sigNhf[i]+1, tmpNhf[tid]);

		 /*  for(l=0 ;l<(2*sigNhf[i]+1);l++){
			tmpNhf[i][l] = extrhok[i][ind1+l]*knlt1[i][l];
		   }//l == 2sigN+1    		   
		   wnv2[i][Dim-1][j] = dx[Dim-1]*2.*PI*integ_simpson(2*sigNhf[i]+1, tmpNhf[i]);
		*/
                
		}//j 	
      }//i
  }//Dim ==1 

  for(j=0;j<M;j++){

	wn3[TNsp][j] = 0;

	for(i=0;i<TNsp;i++){

		wn3[TNsp][j] += wn3[i][j];
	}
  }
  
 
}


