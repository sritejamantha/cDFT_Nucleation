#include<stdio.h>
#include<stdlib.h>
#include<math.h>

int main()
{
 int i;
 double dx=0.2,x=0.0,gp[750]={0.0},gb;
 double xi=0.0,xb=0.0,Integ=0.0,rf_len=2.99;
 double PI=3.14159265359;

 FILE *f1;

 f1 = fopen("Grand_5000.dat","r");
 for(i=0;i<750;i++)
 {
  fscanf(f1,"%lf %lf\n",&x,&gp[i]);
 }
 fclose(f1);


 dx = dx/rf_len;
 xb = 749*dx;
 gb = 4.0*PI*gp[749]*pow(xb,3.0)/3.0;

 Integ=0.5*gp[749]*xb*xb;
 for(i=1;i<749;i++)
 {
  xi=(double)i*dx;
  Integ = Integ + gp[i]*xi*xi;
 }
 
 Integ = Integ*4.0*PI*dx;
 printf("Diff:%.15e\tSystem:%.15e\tBlk:%.15e\n",Integ-gb,Integ,gb);

}
