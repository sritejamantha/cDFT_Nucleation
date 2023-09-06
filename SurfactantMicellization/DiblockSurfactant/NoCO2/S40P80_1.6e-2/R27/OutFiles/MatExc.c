#include<stdio.h>
#include<stdlib.h>
#include<math.h>

int main()
{
 int i;
 double dx=0.2,x=0.0,rs[750]={0.0},rs_b;
 double ra[750]={0.0},ra_b;
 double rb[750]={0.0},rb_b;
 double xi=0.0,xb=0.0,Integ_s=0.0,Integ_a=0.0,Integ_b=0.0,rf_len=2.99;
 double PI=3.14159265359;

 FILE *f1;

 f1 = fopen("rho0.dat","r");
 for(i=0;i<750;i++)
 {
  fscanf(f1,"%lf %lf\n",&x,&rs[i]);
 }
 fclose(f1);


 dx = dx/rf_len;
 xb = 749*dx;
 rs_b = 4.0*PI*rs[749]*pow(xb,3.0)/3.0;

 Integ_s=0.5*rs[749]*xb*xb;
 for(i=1;i<749;i++)
 {
  xi=(double)i*dx;
  Integ_s = Integ_s + rs[i]*xi*xi;
 }
 
 Integ_s = Integ_s*4.0*PI*dx;
 printf("Diff:%.15e\tSystem:%.15e\tBlk:%.15e\n",Integ_s-rs_b,Integ_s,rs_b);

 f1 = fopen("rho1.dat","r");
 for(i=0;i<750;i++)
 {
  fscanf(f1,"%lf %lf\n",&x,&ra[i]);
 }
 fclose(f1);


 xb = 749*dx;
 ra_b = 4.0*PI*ra[749]*pow(xb,3.0)/3.0;
 printf("%.15e\n",ra[749]);
 Integ_a=0.5*ra[749]*xb*xb;
 for(i=1;i<749;i++)
 {
  xi=(double)i*dx;
  Integ_a = Integ_a + ra[i]*xi*xi;
 }

 Integ_a = Integ_a*4.0*PI*dx;
 printf("Diff:%.15e\tSystem:%.15e\tBlk:%.15e\n",Integ_a-ra_b,Integ_a,ra_b);

 f1 = fopen("rho2.dat","r");
 for(i=0;i<750;i++)
 {
  fscanf(f1,"%lf %lf\n",&x,&rb[i]);
 }
 fclose(f1);


 xb = 749*dx;
 rb_b = 4.0*PI*rb[749]*pow(xb,3.0)/3.0;
 printf("%.15e\n",rb[749]);
 Integ_b=0.5*rb[749]*xb*xb;
 for(i=1;i<749;i++)
 {
  xi=(double)i*dx;
  Integ_b = Integ_b + rb[i]*xi*xi;
 }

 Integ_b = Integ_b*4.0*PI*dx;
 printf("Diff:%.15e\tSystem:%.15e\tBlk:%.15e\n",Integ_b-rb_b,Integ_b,rb_b);
}
