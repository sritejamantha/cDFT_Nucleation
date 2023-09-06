#include<stdio.h>
#include<stdlib.h>
#include<math.h>

int main()
{
 int i,j;
 double x[800]={0.0},y[800]={0.0};
 double Test,diff,dx=0.25,L=200.0,rhoA_blk=0.0,rhoB_blk=0.0,x_GDS=118.25,ExAds=0.0;
 int Lhlf,ntest,L_GDS;
 double PDMS=0.0,PPO=0.0;

 FILE *f1,*f2;

 f1 = fopen("rho2_0.dat","r");

 for(i=0;i<800;i++)
 {
  fscanf(f1,"%lf %lf",&x[i],&y[i]);
 }

 fclose(f1);


 rhoA_blk = y[0];
 rhoB_blk = y[799];

 printf("%.15e \t %.15e\n",rhoA_blk,rhoB_blk);

 L_GDS = (int)(x_GDS/dx);


 ExAds=0.0;
 ExAds = 0.5*( (y[0]-rhoA_blk) + (y[L_GDS]-rhoA_blk) );

 for(i=1;i<L_GDS;i++)
 {
  ExAds = ExAds + (y[i]-rhoA_blk);
 }

 ExAds += 0.5*( (y[L_GDS]-rhoB_blk) + (y[799]-rhoB_blk) );

 for(i=L_GDS+1;i<799;i++)
 {
  ExAds = ExAds + (y[i]-rhoB_blk);
 }

 ExAds = ExAds*dx;
 PDMS = ExAds;
 printf("PDMS:%.15f\n",ExAds);


 f1 = fopen("rho3_0.dat","r");

 for(i=0;i<800;i++)
 {
  x[i]=0.0;y[i]=0.0;
  fscanf(f1,"%lf %lf",&x[i],&y[i]);
 }

 fclose(f1);


 rhoA_blk = y[0];
 rhoB_blk = y[799];

 printf("%.15e \t %.15e\n",rhoA_blk,rhoB_blk);

 L_GDS = (int)(x_GDS/dx);


 ExAds=0.0;
 ExAds = 0.5*( (y[0]-rhoA_blk) + (y[L_GDS]-rhoA_blk) );

 for(i=1;i<L_GDS;i++)
 {
  ExAds = ExAds + (y[i]-rhoA_blk);
 }

 ExAds += 0.5*( (y[L_GDS]-rhoB_blk) + (y[799]-rhoB_blk) );

 for(i=L_GDS+1;i<799;i++)
 {
  ExAds = ExAds + (y[i]-rhoB_blk);
 }

 ExAds = ExAds*dx;
 PPO = ExAds;
 printf("PPO:%.15f\n",ExAds);
 
 printf("PDMS+PO:%.15f\n",PDMS+PPO);
}
