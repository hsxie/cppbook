/**************************************************************
      Solved Poisson Equation, Using Multigrid Mehtod
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
#define lq (3*(nx+1)*(ny+1))
#define mm 10

double f[nx+1][ny+1],u[nx+1][ny+1],r[nx+1][ny+1],q[lq];
int nx1=2,ny1=2,nu1=2,nu2=2,ncyc=100;
int m,mx[mm],my[mm],ipu[mm],ipf[mm],level;
double hh[mm],rn;
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

void setgrd()  // Set Grid
{
	int lastq,j,k,msize;
	double dd;
	j=nx;
	k=ny;
	dd=h;
	lastq=0;
	m=0;
	do{
		m=m+1;
		mx[m]=j;
		my[m]=k;
		hh[m]=dd;
		msize=(j+1)*(k+1);
		ipu[m]=lastq+1;
		lastq+=msize;
		ipf[m]=lastq+1;
		lastq+=msize;
		j=j/2;
		k=k/2;
		dd=2*dd;
		printf("m\tmx[m]\tmy[m]\thh[m]\tipu[m]\tipf[m]\n%d\t%d\t%d\t%f\t%d\t%d\n\n",m,mx[m],my[m],hh[m],ipu[m],ipf[m]);
	}while((j%2==0)&&(k%2==0)&&(m<mm));
	level=m;
}

void uftoq() // u[], f[] to q[]
{
	int j,k,i1,i2;
	for(j=0;j<=mx[1];j++)
	{
		for(k=0;k<=my[1];k++)
		{
			i1=ipu[1]+k*(mx[1]+1)+j;
			i2=ipf[1]+k*(mx[1]+1)+j;
			q[i1]=u[j][k];
			q[i2]=f[j][k];
		}
	}
}

void qtouf() // q[] to u[], f[]
{
	int j,k,i1,i2;
	for(j=0;j<=mx[1];j++)
	{
		for(k=0;k<=my[1];k++)
		{
			i1=ipu[1]+k*(mx[1]+1)+j;
			i2=ipf[1]+k*(mx[1]+1)+j;
			u[j][k]=q[i1];
			f[j][k]=q[i2];
		}
	}
}

void mgrlax(int l) // Relax
{
	int j,k,i1,i2;
	double rr,unew;
	rn=0.0;
	for(j=1;j<mx[l];j++)
	{
		for(k=1;k<my[l];k++)
		{
			i1=ipu[l]+k*(mx[l]+1)+j;
			i2=ipf[l]+k*(mx[l]+1)+j;
			unew=0.25*(q[i2]*hh[l]*hh[l]+q[i1+mx[l]+1]+q[i1-mx[l]-1]+q[i1+1]+q[i1-1]);
			rr=unew-q[i1];
			rn+=rr*rr;
			q[i1]=unew;
		}
	}
	rn=sqrt(rn)/hh[l];
}

void mgftoc(int l) // Fine to Coarse
{
	int j,k,jc,kc,i1,i2,i3,i4;
	for(j=1;j<mx[l];j++)
	{
		for(k=1;k<my[l];k++)
		{
			i1=ipu[l]+k*(mx[l]+1)+j;
			i2=ipf[l]+k*(mx[l]+1)+j;
			r[j][k]=q[i2]+(q[i1+mx[l]+1]+q[i1-mx[l]-1]+q[i1+1]+q[i1-1]-4.0*q[i1])/hh[l]/hh[l];
		}
	}
	for(jc=1;jc<mx[l+1];jc++)
	{
		for(kc=1;kc<my[l+1];kc++)
		{
			j=2*jc;
			k=2*kc;
			i3=ipf[l+1]+kc*(mx[l+1]+1)+jc;
			q[i3]=0.25*(r[j][k]+0.5*(r[j][k+1]+r[j][k-1]+r[j+1][k]+r[j-1][k])+0.25*(r[j+1][k+1]+r[j+1][k-1]+r[j-1][k+1]+r[j-1][k-1]));
		}
	}
	for(jc=0;jc<=mx[l+1];jc++)
	{
		for(kc=0;kc<=my[l+1];kc++)
		{
			i4=ipu[l+1]+kc*(mx[l+1]+1)+jc;
			q[i4]=0.0;
		}
	}
}

void mgctof(int l) // Coarse to Fine
{
	int j,k,jc,kc,jf,kf,i1,i2;
	q[ipu[l-1]]+=q[ipu[l]];
	for(kc=1;kc<=my[l];kc++)
	{
		kf=2*kc;
		q[ipu[l-1]+(kf-1)*(mx[l-1]+1)]+=0.5*(q[ipu[l]+(kc-1)*(mx[l]+1)]+q[ipu[l]+kc*(mx[l]+1)]);
		q[ipu[l-1]+kf*(mx[l-1]+1)]+=q[ipu[l]+kc*(mx[l]+1)];
	}
	for(jc=1;jc<=mx[l];jc++)
	{
		jf=2*jc;
		q[ipu[l-1]+jf-1]+=0.5*(q[ipu[l]+jc-1]+q[ipu[l]+jc]);
		q[ipu[l-1]+jf]+=q[ipu[l]+jc];
		for(kc=1;kc<=my[l];kc++)
		{
			kf=2*kc;
			i1=ipu[l-1]+kf*(mx[l-1]+1)+jf;
			i2=ipu[l]+kc*(mx[l]+1)+jc;
			q[i1-mx[l-1]-1-1]+=0.25*(q[i2-mx[l]-1-1]+q[i2-mx[l]-1]+q[i2-1]+q[i2]);
			q[i1-mx[l-1]-1]+=0.5*(q[i2-mx[l]-1]+q[i2]);
			q[i1-1]+=0.5*(q[i2-1]+q[i2]);
			q[i1]+=q[i2];
		}
	}
}

void mgcv(int lf,int lc) // MultiGrid V-Cycle Routine
{
	int l,icycle;
	for(icycle=1;icycle<=ncyc;icycle++)
	{
		for(l=lf;l<lc;l++) //Downward
		{
			mgrlax(l);
			mgftoc(l);
		}
		printf("Cycle %d.\n",icycle);
		l=lc;
		mgrlax(l);
		for(l=lc;l>lf;l--)  //Upward
		{
			mgrlax(l);
			mgctof(l);
		}
		printf("rn=%8.4e\n",rn);
		if(rn<eps) break;
	}
} 

int main() // Main Function
{
	int j,k;
	double x,y,error,errmx,totalt;
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
	
	setgrd();
	uftoq();
	printf("level=%d\n",level);
	
	/* Solve */
	t1=clock();
	mgcv(1,level);
	t2=clock();
	
	qtouf();
	
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
