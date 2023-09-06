#include<stdio.h>
#include<stdlib.h>
#include<math.h>

int main()
{
 int i,k,it,nit=40,dt=5000,count=0,FnCount=0;
 double ExFE[40]={0.0},Exrho1[40]={0.0};
 double FnExFE[40]={0.0},FnExrho1[40]={0.0};
 double dx=0.2,x=0.0,gp[750]={0.0},gb;
 double xi=0.0,xb=0.0,Integ=0.0,rf_len=2.99;
 double diff=0.0,Integ_a=0.0,r1_b=0.0,rho1[750]={0.0};
 double PI=3.14159265359;
 double AvgExFE=0.0,AvgExrho1=0.0,VarFE=0.0,Varrho1=0.0;
 double FnAvgExFE=0.0,FnAvgExrho1=0.0;
 char fname0[100],fname1[100];

 FILE *f1,*f2;

 dx = dx/rf_len;
 xb = 749*dx;

 count=0;
 for(it=0;it<nit-1;it++)
 {
  ExFE[it]=0.0;Exrho1[it]=0.0;
  count=count+1;
 }

 AvgExFE=0.0;AvgExrho1=0.0;VarFE=0.0;Varrho1=0.0;
 for(it=0;it<nit-1;it++)
 {
  gb=0.0;Integ=0.0;
  r1_b=0.0;Integ_a=0.0;

  sprintf(fname0,"Grand_%d.dat",it*dt + 50000);
  f1 = fopen(fname0,"r");

  sprintf(fname1,"rho1_%d.dat",it*dt + 50000);
  f2 = fopen(fname1,"r");

//  printf("%s \t %s\n",fname0,fname1);
  for(i=0;i<750;i++)
  {
   gp[i]=0.0;
   rho1[i]=0.0;
  }

  for(i=0;i<750;i++)
  {
   fscanf(f1,"%lf %lf\n",&x,&gp[i]);
   fscanf(f2,"%lf %lf\n",&x,&rho1[i]);
  }
  fclose(f1);fclose(f2);


  gb = 4.0*PI*gp[749]*pow(xb,3.0)/3.0;
  r1_b = 4.0*PI*rho1[749]*pow(xb,3.0)/3.0;

  Integ=0.5*gp[749]*xb*xb;
  Integ_a=0.5*rho1[749]*xb*xb;
  for(i=1;i<749;i++)
  {
   xi=(double)i*dx;
   Integ = Integ + gp[i]*xi*xi;
   Integ_a = Integ_a + rho1[i]*xi*xi;
  }
 
  Integ = Integ*4.0*PI*dx;
  Integ_a = Integ_a*4.0*PI*dx;

  ExFE[it] = Integ-gb; 
  Exrho1[it] = Integ_a - r1_b;
 
  AvgExFE = AvgExFE + ExFE[it];
  AvgExrho1 = AvgExrho1 + Exrho1[it];
//  printf("Diff:%.15e\tSystem:%.15e\tBlk:%.15e\n",Integ-gb,Integ,gb);
  printf("%d \t %.15e \t %.15e\n",it,ExFE[it],Exrho1[it]);
 }

 AvgExFE = AvgExFE/((double)nit-1.0);
 AvgExrho1 = AvgExrho1/((double)nit-1.0);
 printf("%.15e \t %.15e\n",AvgExrho1,AvgExFE);
// exit(-1);
 k=0;diff=0.0;
 for(it=0;it<nit-1;it++)
 {
  diff = (ExFE[it]-AvgExFE)/AvgExFE;
  /*if(diff<0)
  {
   diff=-diff;
  }
  if(diff<0.25)*/
  if(ExFE[it]>0.0 && diff<0.25)
  {
   FnExFE[k] = ExFE[it];
   FnExrho1[k] = Exrho1[it];
//   printf("Final:%.15e\n",FnExFE[k]);
   printf("Fin:%d \t %.15e \t %.15e\n",k,FnExFE[k],FnExrho1[k]);
   k=k+1;
   FnCount = FnCount+1;
  }

 }

 printf("%d \t %d\n",nit-1,FnCount);
 
 FnAvgExFE=0.0;FnAvgExrho1=0.0;
 for(it=0;it<FnCount;it++)
 {
  FnAvgExFE = FnAvgExFE + FnExFE[it];
  FnAvgExrho1 = FnAvgExrho1 + FnExrho1[it];
 } 

 FnAvgExFE = FnAvgExFE/((double)FnCount);
 FnAvgExrho1 = FnAvgExrho1/((double)FnCount);

 for(it=0;it<FnCount;it++)
 {
  VarFE = VarFE + pow(FnExFE[it]-FnAvgExFE,2.0);
  Varrho1 = Varrho1 + pow(FnExrho1[it]-FnAvgExrho1,2.0);
 }

 VarFE = sqrt(VarFE/((double)FnCount));
 Varrho1 = sqrt(Varrho1/((double)FnCount));
 
 printf("%.15e \t %.15e \t %.15e \t %.15e\n",FnAvgExrho1,FnAvgExFE,Varrho1,VarFE);
// printf("%d\n",count);
}
