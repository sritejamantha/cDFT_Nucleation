#include "globals.h"

double find_blkmu()
{

 int i,j,k;

 double sigv[crsTNsp],fblkr[TNsp],tmpsig[TNsp][TNsp],tmpeps[TNsp][TNsp];
 double tmpfid,tmprhoN3,tmpNN3,blkn2[TNsp],blkn3[TNsp+1],blkn0[TNsp],blkdfhsn3[TNsp],blkdfhs[TNsp],blkfhs,blkfid[TNsp], blkdfid[TNsp];


 for(i=0; i< TNsp ; i++)
 {
  fblkr[i] = blkrho[i][1]; //load blk density 
//  printf("%d \t %.15f\n",i,fblkr[i]);
 }

 tmprhoN3 = fblkr[TNsp-3]+fblkr[TNsp-2]+fblkr[TNsp-1];
 tmpNN3 = (double)(NK[TNsp-3]+NK[TNsp-2]+NK[TNsp-1]);
 
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
   cout<<"sigvs "<<sigv[k]<<endl;
   k=k+1;
  }
 }

/* for(i=0;i<TNsp;i++)
 {
  printf("%d \t %lf\n",i,epsln[i]);
 }*/
// exit(-1);
// epsln[0]=170.5;epsln[1]=226.5;epsln[2]=204.58;epsln[3]=257.73;epsln[4]=226.5;
 for(i=0;i<TNsp;i++)
 {
  tmpsig[i][i] = rsigma[i];
  tmpeps[i][i] = epsln[i];
 }

 //exit(-1);


 blkn3[TNsp] = 0;
 
  tmpfid=0;
  for(i=0; i< TNsp ; i++)
  {
  
   if(fblkr[i] >0 )
   {
    if(i<TNsp-3)
    {
     blkdfid[i] = log(fblkr[i]/NK[i])/NK[i];
     blkfid[i] = fblkr[i]*(log(fblkr[i]/NK[i])-1)/NK[i];
     tmpfid = tmpfid+blkfid[i];
//     printf("%d \t %.15f \t %.15f\n",i,blkdfid[i],blkfid[i]);
    }
    else
    {
     if(i==TNsp-3)
     {
      blkdfid[i] = Xa*log(tmprhoN3/tmpNN3)/tmpNN3;
      blkfid[i] = Xa*tmprhoN3*(log(tmprhoN3/tmpNN3)-1)/tmpNN3;
      tmpfid = tmpfid+blkfid[i];
//      printf("%d \t %.15f \t %.15f\n",i,blkdfid[i],blkfid[i]);
     }
     else if(i==TNsp-2)
     {
      blkdfid[i] = Xb*log(tmprhoN3/tmpNN3)/tmpNN3;
      blkfid[i] = Xb*tmprhoN3*(log(tmprhoN3/tmpNN3)-1)/tmpNN3;
      tmpfid = tmpfid+blkfid[i];
//      printf("%d \t %.15f \t %.15f\n",i,blkdfid[i],blkfid[i]);
     }
     else if(i==TNsp-1)
     {
      blkdfid[i] = (1.0-Xa-Xb)*log(tmprhoN3/tmpNN3)/tmpNN3;
      blkfid[i] = (1.0-Xa-Xb)*tmprhoN3*(log(tmprhoN3/tmpNN3)-1)/tmpNN3;
      tmpfid = tmpfid+blkfid[i];
//      printf("%d \t %.15f \t %.15f\n",i,blkdfid[i],blkfid[i]);
     }
    }
   }
   else
   {
    blkdfid[i] = 0;
    blkfid[i] =0;
 
   }
	
   blkn0[i] = fblkr[i];
   blkn3[i] = fblkr[i]*PI/6*pow(sigma[i],3); 
   blkn3[TNsp] += blkn3[i];
  } //i TNsp

//  printf("Fid=%.15f\n",tmpfid);
//  printf("dFiddrho3=%.15f\n",blkdfid[TNsp-3]+blkdfid[TNsp-2]+blkdfid[TNsp-1]);
  // calculate bulk values from fhs
  
//  exit(-1);
  double tmp_1m,tmp_lg, tmp_dd[4], loc_x[TNsp], loc_nbar, loc_tt, tmp_gg[TNsp],tmp_eta, fhspart[6];
  
  for(i =0 ;i<4; i++)
              tmp_dd[i] = 0;
    
   for(i=0; i< TNsp ; i++)
   {
    tmp_dd[0] += blkn3[i]/sigma[i];
    tmp_dd[1] += blkn3[i]/sigma[i]/sigma[i]/PI;
    tmp_dd[2] += blkn3[i]/pow(sigma[i],3)/PI*6.;
   }

   tmp_eta = blkn3[TNsp];
   tmp_1m = 1- tmp_eta;
   tmp_lg =  log(tmp_1m);
  
   blkfhs = -1*tmp_lg*tmp_dd[2] + 18*tmp_dd[0]*tmp_dd[1]/tmp_1m + 6.0*pow(tmp_dd[0],3)*(tmp_lg/tmp_eta/tmp_eta + 1./tmp_eta/tmp_1m/tmp_1m)/PI; 

//   printf("Fhs:%.15f\n",blkfhs);
   
   fhspart[0] = tmp_dd[2]/tmp_1m;
   fhspart[1] = 18*tmp_dd[0]/tmp_1m;
   fhspart[2] = 18*tmp_dd[1]/tmp_1m;
   fhspart[3] = 18*tmp_dd[1]*tmp_dd[0]/tmp_1m/tmp_1m;
   fhspart[4] = 18*tmp_dd[0]*tmp_dd[0]*(tmp_lg/pow(tmp_eta,2.)+ 1./(tmp_eta*tmp_1m*tmp_1m))/PI;
   fhspart[5] = 6*pow(tmp_dd[0],3)*(2*tmp_lg/pow(tmp_eta,3.) + (2. -5*tmp_eta+ tmp_eta*tmp_eta)/pow(tmp_eta,2)/pow(tmp_1m,3) )/PI;
       
   for(i =0 ;i<TNsp; i++)
   {
    blkdfhsn3[i] =  fhspart[0] - tmp_lg/pow(sigma[i],3)/PI*6.+fhspart[1]/sigma[i]/sigma[i]/PI+ fhspart[2]/sigma[i] + fhspart[3] + fhspart[4]/sigma[i] - fhspart[5];

    if(i==TNsp-3)
    {
     blkdfhs[i] = Xa*PI*blkdfhsn3[i]/6*pow(sigma[i],3);
//     printf("%d \t %.15f\n",i,blkdfhs[i]);
    }

    else if(i==TNsp-2)
    {
     blkdfhs[i] = Xb*PI*blkdfhsn3[i]/6*pow(sigma[i],3);
//     printf("%d \t %.15f\n",i,blkdfhs[i]);
    }

    else if(i==TNsp-1)
    {
     blkdfhs[i] = (1.0 - Xa - Xb)*PI*blkdfhsn3[i]/6*pow(sigma[i],3);
//     printf("%d \t %.15f\n",i,blkdfhs[i]);
    }

    else
    {
     blkdfhs[i] = PI*blkdfhsn3[i]/6*pow(sigma[i],3);
//     printf("%d \t %.15f\n",i,blkdfhs[i]);
    }
 	     
	       //cout<<blkfhs<<endl;
   }

//   printf("3a-3b-3c \t %.15f\n",blkdfhs[TNsp-3]+blkdfhs[TNsp-2]+blkdfhs[TNsp-1]);
//   printf("Fhs:%.15f\n",blkfhs);
//   exit(-1);

   //calculate values from fch
   double blkfch, blkdfchn3[TNsp],blkdfch[TNsp],tmp_rsigab,tmp_sigab,tmp_gab,tmp_alpha,tmp_beta,tmp_dd0,tmp_dd1,tmp_dd2;
   double tmp_rsigbc,tmp_sigbc,tmp_gbc;

   tmp_rsigab = 2.0*rsigma[TNsp-3]*rsigma[TNsp-2]/(rsigma[TNsp-3]+rsigma[TNsp-2]);
   tmp_sigab = tmp_rsigab*(1. - 0.12*exp(-3.0*crsepsln[crsTNsp-3]/Tmpr));
   
   tmp_rsigbc = 2.0*rsigma[TNsp-2]*rsigma[TNsp-1]/(rsigma[TNsp-2]+rsigma[TNsp-1]);
   tmp_sigbc = tmp_rsigbc*(1. - 0.12*exp(-3.0*crsepsln[crsTNsp-1]/Tmpr));
//   printf("%.15e \t %.15e\n",tmp_sigab,tmp_sigbc); 
 
  // tmp_sigab = sigma[TNsp-1];
   for(i =0 ;i<4; i++)
	tmp_dd[i] = 0;
	
   for(i =0 ;i<TNsp; i++)
	tmp_dd[3] +=  blkn3[i]/sigma[i];
	
   for(i =0 ;i<TNsp; i++)
	tmp_gg[i] =  1./tmp_1m + 3*tmp_dd[3]*sigma[i]/2.0/pow(tmp_1m,2) + pow(tmp_dd[3]*sigma[i],2)/2./pow(tmp_1m,3);

   tmp_gab =  1./tmp_1m + 3*tmp_dd[3]*tmp_sigab/2.0/pow(tmp_1m,2) + pow(tmp_dd[3]*tmp_sigab,2)/2./pow(tmp_1m,3);
   tmp_gbc =  1./tmp_1m + 3*tmp_dd[3]*tmp_sigbc/2.0/pow(tmp_1m,2) + pow(tmp_dd[3]*tmp_sigbc,2)/2./pow(tmp_1m,3);
//   printf("%.15e \t %.15e\n", tmp_gab,tmp_gbc); 	

   tmp_alpha = blkn3[TNsp-3]/pow(sigma[TNsp-3],3) + blkn3[TNsp-2]/pow(sigma[TNsp-2],3) + blkn3[TNsp-1]/pow(sigma[TNsp-1],3);
   tmp_beta = (1. - double(NK[TNsp-3]))*log( tmp_gg[TNsp-3]) + (1. - double(NK[TNsp-2]))*log( tmp_gg[TNsp-2]) + (1. - (double)NK[TNsp-1])*log( tmp_gg[TNsp-1]) - (1. - (double)Nbr)*( log( tmp_gg[TNsp-2])+log( tmp_gg[TNsp-1]) ) - (double)Nbr*( log(tmp_gab) + log(tmp_gbc) );
   tmp_beta = tmp_beta/tmpNN3;


   for(i =0 ;i<TNsp-3; i++)
   {

    tmp_dd[0] += (1.-double(NK[i]))/double(NK[i])*blkn3[i]/sigma[i]/tmp_gg[i];
    tmp_dd[1] += (1.-double(NK[i]))/double(NK[i])*blkn3[i]/sigma[i]/sigma[i]/PI/tmp_gg[i];
    tmp_dd[2] += (1.-double(NK[i]))/double(NK[i])*blkn3[i]/pow(sigma[i],3)/PI/tmp_gg[i]*6.;

   }//i TNsp


   tmp_dd2 = (1.0 - double(NK[TNsp-3]))/tmp_gg[TNsp-3] + (1.0 - double(NK[TNsp-2]))/tmp_gg[TNsp-2] + (1.0 - double(NK[TNsp-1]))/tmp_gg[TNsp-1] - (1.0 - double(Nbr))/tmp_gg[TNsp-2] - (1.0 - double(Nbr))/tmp_gg[TNsp-1] - double(Nbr)/tmp_gab - double(Nbr)/tmp_gbc;
   tmp_dd2 = 6.0*tmp_alpha*tmp_dd2/tmpNN3/PI;
   tmp_dd[2] = tmp_dd[2] + tmp_dd2;

  
   tmp_dd1 = sigma[TNsp-3]*(1.0 - double(NK[TNsp-3]))/tmp_gg[TNsp-3] + sigma[TNsp-2]*(1.0 - double(NK[TNsp-2]))/tmp_gg[TNsp-2] +  sigma[TNsp-1]*(1.0 - double(NK[TNsp-1]))/tmp_gg[TNsp-1] - sigma[TNsp-2]*(1.0 - double(Nbr))/tmp_gg[TNsp-2] - sigma[TNsp-1]*(1.0 - double(Nbr))/tmp_gg[TNsp-1] - tmp_sigab*double(Nbr)/tmp_gab - tmp_sigbc*double(Nbr)/tmp_gbc;
   tmp_dd1 = tmp_alpha*tmp_dd1/tmpNN3/PI;
   tmp_dd[1] = tmp_dd[1] + tmp_dd1;	

   tmp_dd0 = pow(sigma[TNsp-3],2.0)*(1.0 - double(NK[TNsp-3]))/tmp_gg[TNsp-3] + pow(sigma[TNsp-2],2.0)*(1.0 - double(NK[TNsp-2]))/tmp_gg[TNsp-2] +  pow(sigma[TNsp-1],2.0)*(1.0 - double(NK[TNsp-1]))/tmp_gg[TNsp-1] - pow(sigma[TNsp-2],2.0)*(1.0 - double(Nbr))/tmp_gg[TNsp-2] - pow(sigma[TNsp-1],2.0)*(1.0 - double(Nbr))/tmp_gg[TNsp-1] - pow(tmp_sigab,2.0)*double(Nbr)/tmp_gab - pow(tmp_sigbc,2.0)*double(Nbr)/tmp_gbc;
   tmp_dd0 = tmp_alpha*tmp_dd0/tmpNN3;
   tmp_dd[0] = tmp_dd[0] + tmp_dd0; 
 

   blkfch =0 ;
   blkfch = -(tmprhoN3/tmpNN3)*(1. - double(Nbr))*( log(tmp_gg[TNsp-2]) + log(tmp_gg[TNsp-1]) );
   blkfch = blkfch - (tmprhoN3/tmpNN3)*double(Nbr)*( log(tmp_gab) + log(tmp_gbc) );
   for(i =0; i<TNsp; i++)
   {
     tmp_lg = log(tmp_gg[i]);
     blkfch += 6.*(1.-double(NK[i]))/double(NK[i])*blkn3[i]*tmp_lg/PI/pow(sigma[i],3);

     if(i<TNsp-3)
     {
      blkdfchn3[i] = 6*(1.-double(NK[i]))/double(NK[i])*tmp_lg/PI/pow(sigma[i],3);
     }
     else
     {
      blkdfchn3[i] = 6*tmp_beta/PI/pow(sigma[i],3);
     }

     blkdfchn3[i] += 1./tmp_1m/tmp_1m*tmp_dd[2];
     blkdfchn3[i] += (9./tmp_1m/tmp_1m/sigma[i]+ 18.*tmp_dd[3]/pow(tmp_1m,3) )*tmp_dd[1];
     blkdfchn3[i] += (6.*tmp_dd[3]/pow(tmp_1m,3)/PI/sigma[i]  +  9.*pow(tmp_dd[3],2)/pow(tmp_1m,4)/PI)*tmp_dd[0];

     if(i==TNsp-3)
     {
      blkdfch[i] =  Xa*PI*pow(sigma[i],3)/6*blkdfchn3[i];
//      printf("%d \t %.15f\n",i,blkdfch[i]);
     }
     else if(i==TNsp-2)
     {
      blkdfch[i] =  Xb*PI*pow(sigma[i],3)/6*blkdfchn3[i];
//      printf("%d \t %.15f\n",i,blkdfch[i]); 
     }
     else if(i==TNsp-1)
     {
      blkdfch[i] =  (1.0 - Xa - Xb)*PI*pow(sigma[i],3)/6*blkdfchn3[i];
//      printf("%d \t %.15f\n",i,blkdfch[i]);
     }
     else
     {
      blkdfch[i] =  PI*pow(sigma[i],3)/6*blkdfchn3[i];
//      printf("%d \t %.15f\n",i,blkdfch[i]);
     }
	 
   }
 
//    printf("Fch:%.15f\n",blkfch);
//    printf("dFchdr3:%.15f\n",blkdfch[TNsp-1]+blkdfch[TNsp-2]+blkdfch[TNsp-3]);
  
//    exit(-1);
   //calculate values from fdisp1

   double blkfdisp1,blkdfdisp1n0[TNsp],blkdfdisp1n3,blkdfdisp1[TNsp],tmp_I1, dmn0[TNsp], dI1dn0[TNsp], loc_da[7],loc_a[7], tmp_pars,tmp_n0;

   loc_tt=0;

   loc_nbar = 0; 
	
   tmp_n0 = 0;
	
   for(i=0;i<TNsp;i++)
   {
	tmp_n0 += blkn0[i];	
        if(i==TNsp-3)
        {
	 loc_x[i] =Xa*blkn0[i]/NK[i];
        }
        else if(i==TNsp-2)
        {
         loc_x[i] = Xb*blkn0[i]/NK[i];
        }
        else if(i==TNsp-1)
        {
         loc_x[i] = (1.0 - Xa - Xb)*blkn0[i]/NK[i];
        }
        else
        {
         loc_x[i] = blkn0[i]/NK[i];
        }
	loc_tt += loc_x[i];
	tmp_gg[i] =0;

   }//i
   
   if(loc_tt >0)
   {
	for(i=0;i<TNsp;i++)
        { 
		loc_x[i] /= loc_tt;
                if(i==TNsp-3)
                {
		 loc_nbar += NK[i]*loc_x[i]/Xa;
                }
                else if(i==TNsp-2)
                {
                 loc_nbar += NK[i]*loc_x[i]/Xb;
                }
                else if(i==TNsp-1)
                {
                 loc_nbar += NK[i]*loc_x[i]/(1.0-Xa-Xb);
                }
                else
                {
                 loc_nbar += NK[i]*loc_x[i];
                }
	}
	
   }
   else
   {
	cout<<"wrong blk setup!"<<endl;
	exit(1);

   }
   
//   printf("Nbar:%.15f\n",loc_nbar);

   for(i=0;i<TNsp;i++)
   {
        if(i==TNsp-3)
        {
	 dmn0[i] = -1.*Xa*tmp_n0/loc_tt/loc_tt/NK[i] + 1./loc_tt;
        }
        else if(i==TNsp-2)
        {
         dmn0[i] = -1.*Xb*tmp_n0/loc_tt/loc_tt/NK[i] + 1./loc_tt;
        }
        else if(i==TNsp-1)
        {
         dmn0[i] = -1.*(1.0-Xa-Xb)*tmp_n0/loc_tt/loc_tt/NK[i] + 1./loc_tt;
        }
        else
        {
         dmn0[i] = -1.*tmp_n0/loc_tt/loc_tt/NK[i] + 1./loc_tt;
        }

   }

   tmp_I1 = 0; 
   blkfdisp1 = 0;
   tmp_dd[0] = 0; //dI1dnbar
   tmp_dd[1] = 0; //dI1dn3
   tmp_dd[2]= 0 ; //factors ass w  dI1dn0 & dI1dn3

   // tmp_gg: factors in dfdispldn0 
   for(i=0;i<7;i++)
   {
        loc_a[i] = Adisp[i][0] +(loc_nbar-1.)/loc_nbar*Adisp[i][1]+(loc_nbar-1.)*(loc_nbar-2)/loc_nbar/loc_nbar*Adisp[i][2]; 
	   	
	loc_da[i] = 1./loc_nbar/loc_nbar*Adisp[i][1]+ (3./loc_nbar/loc_nbar- 4./pow(loc_nbar,3))*Adisp[i][2];
	tmp_I1 += loc_a[i]*pow(tmp_eta,i);
 		
	tmp_dd[0] += loc_da[i]*pow(tmp_eta,i);
	
	if(i>0)
	  tmp_dd[1] += loc_a[i]*i*pow(tmp_eta,i-1);
	
		
   }

//   printf("J1:%.15f\n",tmp_I1);

   for(i=0;i<TNsp;i++)
   {
	tmp_pars = blkn0[i]*blkn0[i]*epsln[i]*rsigma3[i]/Tmpr;
	blkfdisp1 += tmp_pars;
	tmp_dd[2] += tmp_pars;
	
	dI1dn0[i] = tmp_dd[0]*dmn0[i];
	
   }//i TNsp
	
   int tmp_count=0;
   for(i=0;i<TNsp-1;i++)
   {
        for(k=i+1;k<TNsp;k++)
        {

                tmp_pars = 2.*blkn0[i]*blkn0[k]*tmpeps[i][k]*pow(tmpsig[i][k],3)/Tmpr;
                blkfdisp1 += tmp_pars;
                tmp_dd[2] += tmp_pars;
               

        }


    }//i TNsp



   tmp_count=0;
    for(i=0;i<TNsp;i++)
    {
        for(k=0;k<TNsp;k++)
        {
                if(i== k)
                  tmp_gg[i] += blkn0[k]*epsln[i]*rsigma3[i]/Tmpr;
                else
                {
                  tmp_gg[i] +=  blkn0[k]*tmpeps[i][k]*pow(tmpsig[i][k],3)/Tmpr;
                }
        }

    }


    blkfdisp1*= -2.*PI*tmp_I1;
//    printf("Fdisp1:%.15f\n",blkfdisp1);

    for(i=0;i<TNsp;i++)
    {
	blkdfdisp1n0[i] =  -4.*PI*tmp_I1*tmp_gg[i] -2.*PI*dI1dn0[i]*tmp_dd[2]; 
	
		//cout<<tmp_I1<<" "<<tmp_gg[i]<<" "<<dI1dn0[i]<<" "<<tmp_dd[2]<<endl;

    }

    blkdfdisp1n3 = -2.*PI*tmp_dd[1]*tmp_dd[2];

    for(i=0;i<TNsp;i++)
    {
        if(i==TNsp-3)
        {
         blkdfdisp1[i] =  Xa*(blkdfdisp1n0[i] + PI/6*pow(sigma[i],3)*blkdfdisp1n3);
        }
        else if(i==TNsp-2)
        {
         blkdfdisp1[i] =  Xb*(blkdfdisp1n0[i] + PI/6*pow(sigma[i],3)*blkdfdisp1n3);
        }
        else if(i==TNsp-1)
        {
         blkdfdisp1[i] =  (1.0-Xa-Xb)*(blkdfdisp1n0[i] + PI/6*pow(sigma[i],3)*blkdfdisp1n3);
        }
        else
        {
         blkdfdisp1[i] =  (blkdfdisp1n0[i] + PI/6*pow(sigma[i],3)*blkdfdisp1n3);
        }
    }
    //calculate values from fdisp2

    double tmp_M, blkdfdisp2[TNsp],blkfdisp2, blkdfdisp2n0[TNsp],blkdfdisp2n3,tmp_I2,dMdn3,dMdn0[TNsp],dI2dn0[TNsp];

    for(i=0; i<TNsp; i++)
	tmp_gg[i] = 0;

    //get coeff a, I1 ,Fdsipl1 and da/dn0
    tmp_I2 = 0; 
    blkfdisp2 = 0;
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

//    printf("J2:%.15f\n",tmp_I2);
 	
    tmp_M = 1.+loc_nbar*(8.*tmp_eta - 2.*tmp_eta*tmp_eta)/pow(tmp_1m,4) + (1.- loc_nbar)*(20.*tmp_eta - 27*tmp_eta*tmp_eta + 12.*pow(tmp_eta,3)-2.*pow(tmp_eta,4.))/tmp_1m/tmp_1m/pow(2.-tmp_eta,2);

//    printf("MInv:%.15f\n",1/tmp_M);

    tmp_dd[0] = (8.*tmp_eta - 2.*tmp_eta*tmp_eta)/pow(tmp_1m,4) - (20.*tmp_eta - 27.*tmp_eta*tmp_eta+12.*pow(tmp_eta,3)-2.*pow(tmp_eta,4))/tmp_1m/tmp_1m/pow(2.-tmp_eta,2);
	
    dMdn3 = loc_nbar*(8.+20.*tmp_eta-4.*tmp_eta*tmp_eta)/pow(tmp_1m,5) + (1.- loc_nbar)*(2*pow(tmp_eta,3)+12.*tmp_eta*tmp_eta-48.*tmp_eta+40)/pow(tmp_1m,3)/pow(2.-tmp_eta,3);
	

    for(i=0;i<TNsp;i++)
    {
	tmp_pars = blkn0[i]*blkn0[i]*pow(epsln[i]/Tmpr,2)*rsigma3[i];
	blkfdisp2 += tmp_pars;
	tmp_dd[3] += tmp_pars;

	dMdn0[i] = tmp_dd[0] * dmn0[i]; 
	dI2dn0[i] = tmp_dd[1]*dmn0[i];
			
    }//i TNsp

   tmp_count=0;
    for(i=0;i<TNsp-1;i++)
    {
        for(k=i+1;k<TNsp;k++)
        {
               tmp_pars = 2.*blkn0[i]*blkn0[k]*pow(tmpeps[i][k]/Tmpr,2)*pow(tmpsig[i][k],3);
               blkfdisp2 += tmp_pars;

               tmp_dd[3] += tmp_pars;
        }


    }//i TNsp

   for(i=0;i<TNsp;i++)
   {
        for(k=0;k<TNsp;k++)
        {
                if(i== k)
                  tmp_gg[i] += blkn0[i]*pow(epsln[i]/Tmpr,2)*rsigma3[i];
                         
                else
                {
                  tmp_gg[i] += blkn0[k]*pow(tmpeps[i][k]/Tmpr,2)*pow(tmpsig[i][k],3);
                }
        }
   }//i /TNsp

    blkfdisp2 *= -1.*PI*tmp_I2/tmp_M*loc_nbar;

//    printf("Fdisp2:%.15f\n",blkfdisp2);
//    printf("Fdisp:%.15f\n",blkfdisp2+blkfdisp1);

    for(i=0;i<TNsp;i++)
    {   
	blkdfdisp2n0[i] = -2.*PI/tmp_M*tmp_I2*loc_nbar*tmp_gg[i] + PI/tmp_M/tmp_M*dMdn0[i]*tmp_I2*loc_nbar*tmp_dd[3] -PI/tmp_M*dI2dn0[i]*loc_nbar*tmp_dd[3] - PI/tmp_M*tmp_I2*dmn0[i]*tmp_dd[3];
	  

    }
	
    blkdfdisp2n3 = PI/tmp_M/tmp_M*dMdn3*tmp_I2*loc_nbar*tmp_dd[3]- PI/tmp_M*tmp_dd[2]*loc_nbar*tmp_dd[3];
	
    for(i=0;i<TNsp;i++)
    {
        if(i==TNsp-3)
        {
		blkdfdisp2[i] = Xa*(blkdfdisp2n0[i] + PI/6*pow(sigma[i],3)*blkdfdisp2n3);
        }
        else if(i==TNsp-2)
        {
        	blkdfdisp2[i] = Xb*(blkdfdisp2n0[i] + PI/6*pow(sigma[i],3)*blkdfdisp2n3);
        }
        else if(i==TNsp-1)
        {
                blkdfdisp2[i] = (1.0 - Xa - Xb)*(blkdfdisp2n0[i] + PI/6*pow(sigma[i],3)*blkdfdisp2n3);
        }
        else
        {
         	blkdfdisp2[i] = (blkdfdisp2n0[i] + PI/6*pow(sigma[i],3)*blkdfdisp2n3);
        }

	amu[i] = blkdfid[i] + blkdfhs[i] + blkdfch[i] + blkdfdisp1[i] +blkdfdisp2[i];

	//cout<<blkdfid[i]<<" "<<blkdfhs[i]<<" "<<blkdfch[i]<<" "<<blkdfdisp1[i] <<" "<<blkdfdisp2[i]<<endl;
	//	cout<<i<<" w amu "<<amu[i]<<endl;
        printf("Mu_%d:%.15f\n",i,amu[i]);
    }

    printf("Mu3a-3b:%.15f\n",amu[TNsp-3]+amu[TNsp-2]+amu[TNsp-1]);
//    exit(-1);

    /*for(i=0;i<TNsp;i++)
    {
     printf("dfdisp-%d:%.15f\n",i,blkdfdisp1[i] +blkdfdisp2[i]);
    }*/

   // printf("dfdisp-3a3b3c:%.15f\n",blkdfdisp1[TNsp-1] +blkdfdisp2[TNsp-1]+blkdfdisp1[TNsp-2] +blkdfdisp2[TNsp-2]+blkdfdisp1[TNsp-3] +blkdfdisp2[TNsp-3]); 
    
   // exit(-1);
    double ttpres=0;
  
    for(i=0;i<TNsp;i++)
    {
        if(i==TNsp-3)
        {
	 ttpres += amu[i]*fblkr[i]/Xa;
        }
        else if(i==TNsp-2)
        {
         ttpres += amu[i]*fblkr[i]/Xb;
        }
        else if(i==TNsp-1)
        {
         ttpres += amu[i]*fblkr[i]/(1.0-Xa-Xb);
        }
        else
        {
         ttpres += amu[i]*fblkr[i];
        }        

	ttpres -= blkfid[i];
    }

    ttpres += -1* blkfhs- blkfch -blkfdisp1 - blkfdisp2;

    return ttpres;

}
