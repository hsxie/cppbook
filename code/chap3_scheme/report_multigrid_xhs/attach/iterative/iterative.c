/***************************************************************
    Solved Poisson Equation, Using Gauss-Seidel Mehtod
                      uxx + uyy = -f
      Hua-sheng XIE, huashengxie@gmail.com, 2010-09-12
**************************************************************/

#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<math.h>

#define nx 128
#define ny 64
#define a -2.0
#define b 2.0
#define h ((b-a)/nx)
#define c -1.0
#define d (c+h*ny)
#define eps 1.0e-12
#define pi 3.14159265
#define w 1.0

double f[nx+1][ny+1],u[nx+1][ny+1];
double mxcount=10000,rn;
double qq1=1.0,xx1=-1.0,yy1=0.0,qq2=-1.0,xx2=1.0,yy2=0.0; //charged point(s)

double ffun(double x,double y)  // Real f
{
	double f;
//	f=-( x*x + y*y )*exp( x*y ); // for testing
//	f=2.0*pi*pi*sin(pi*x)*sin(pi*y);

    /* for charged point(s) */
    f=0.0;
	if(fabs(x-xx1)<h&&fabs(y-yy1)<h)
	{
		f+=qq1*(h-fabs(x-xx1))*(h-fabs(y-yy1))/(h*h*h*h);
	}
	if(fabs(x-xx2)<h&&fabs(y-yy2)<h)
	{
		f+=qq2*(h-fabs(x-xx2))*(h-fabs(y-yy2))/(h*h*h*h);
	}
	
	return f;
}

double ufun(double x,double y)  // u
{
	double u;
//	u=exp( x*y );
//	u=sin(pi*x)*sin(pi*y);
	
	u=qq1/(2.0*pi)*log(1.0/sqrt((x-xx1)*(x-xx1)+(y-yy1)*(y-yy1)))+qq2/(2.0*pi)*log(1.0/sqrt((x-xx2)*(x-xx2)+(y-yy2)*(y-yy2))); //for charged point(s)
	
	return u;
}

int main()
{
	int j,k,count;
	double x,y,error,errmx,totalt;
	double tempu[nx+1][ny+1],maxdu,du;
	clock_t t1,t2;
	int jmx,kmx;
	
	FILE *fp1,*fp2,*fp3,*fp4;
	if((fp1=fopen("usolved.dat","w"))==NULL){
		printf("File open error!\n");
		exit(0);
	}
	if((fp2=fopen("ureal.dat","w"))==NULL){
		printf("File open error!\n");
		exit(0);
	}
	if((fp3=fopen("usolved_ureal.dat","w"))==NULL){
		printf("File open error!\n");
		exit(0);
	}
	if((fp4=fopen("for_contour.dat","w"))==NULL){
		printf("File open error!\n");
		exit(0);
	}
	
	/* Initial */
	for(j=1;j<nx;j++)
	{
		for(k=1;k<ny;k++)
		{
			x=a+h*j;
			y=c+h*k;
			u[j][k]=0.0;
			f[j][k]=ffun(x,y);
		}
	}
	for(j=0;j<=nx;j++)
	{
		u[j][0]=ufun(a+j*h,c);
		u[j][ny]=ufun(a+j*h,d);
	}
	for(k=0;k<=ny;k++)
	{
		u[0][k]=ufun(a,c+k*h);
		u[nx][k]=ufun(b,c+k*h);
	}
	
	
	/* Solve */
	t1=clock();
	
	count=0;
	do{
		rn=0.0;
		for(j=1;j<nx;j++)
		{
			for(k=1;k<ny;k++)
			{
				x=a+h*j;
				y=c+h*k;
				maxdu=0.0;
				tempu[j][k]=u[j][k];
				u[j][k]=(1.0-w)*u[j][k]-(-ffun(x,y)*h*h-u[j+1][k]-u[j-1][k]-u[j][k+1]-u[j][k-1])*w/4.0;				
				du=tempu[j][k]-u[j][k];
				rn+=du*du;
			}
		}
		count++;
		printf("count=%d\n",count);
	}while(rn>eps&&count<mxcount);
	
	t2=clock();
	
	/* Output */
	errmx=0.0;
	
	for(k=0;k<=ny;k++)
	{
		for(j=0;j<=nx;j++)
		{
			x=a+h*j;
			y=c+h*k;
			fprintf(fp1,"%f\t%f\t%f\n",x,y,u[j][k]);
			fprintf(fp2,"%f\t%f\t%f\n",x,y,ufun(x,y));
			error=u[j][k]-ufun(x,y);
			fprintf(fp3,"%f\t%f\t%f\n",x,y,error);
			fprintf(fp4,"%f\t",ufun(x,y));
			if(fabs(errmx)<fabs(error))
			{
				jmx=j;
				kmx=k;
				errmx=error;
			}
		}
		fprintf(fp4,"\n");
	}
	
	totalt=((double)(t2-t1))/CLK_TCK;
	printf("Total time used: %f\n",totalt);
	printf("Max error: %8.4e at j=%d k=%d\n",errmx,jmx,kmx);
	
	if(fclose(fp1))
	{
		printf("Can not close the file!\n");
		exit(0);
	}
	if(fclose(fp2))
	{
		printf("Can not close the file!\n");
		exit(0);
	}	
	if(fclose(fp3))
	{
		printf("Can not close the file!\n");
		exit(0);
	}
	if(fclose(fp4))
	{
		printf("Can not close the file!\n");
		exit(0);
	}
}
