#include<complex>
#include<iostream>
using namespace std;

/**********Constant Variables************/

const double PI= 3.1415926535897;
const double eps = 1e-12;

/****************************************/



/**********Math Function************/

static const int ncof=28;

const double cof[28] = {-1.3026537197817094, 6.4196979235649026e-1,
1.9476473204185836e-2,-9.561514786808631e-3,-9.46595344482036e-4,
3.66839497852761e-4,4.2523324806907e-5,-2.0278578112534e-5,
-1.624290004647e-6,1.303655835580e-6,1.5626441722e-8,-8.5238095915e-8,
6.529054439e-9,5.059343495e-9,-9.91364156e-10,-2.27365122e-10,
9.6467911e-11, 2.394038e-12,-6.886027e-12,8.94487e-13, 3.13092e-13,
-1.12708e-13,3.81e-16,7.106e-15,-1.523e-15,-9.4e-17,1.21e-16,-2.8e-17};

double erfccheb(double z){
	int j;
	double t,ty,tmp,d=0.,dd=0.;
	if (z < 0.) 
		throw("erfccheb requires nonnegative argument");
	t = 2./(2.+z);
	ty = 4.*t - 2.;
	for (j=ncof-1;j>0;j--) {
		tmp = d;
		d = ty*d - dd + cof[j];
		dd = tmp;
	}
	return t*exp(-z*z + 0.5*(cof[0] + ty*d) - dd);
}

inline double erf(double x) {
	if (x >=0.) 
		return 1.0 - erfccheb(x);
	else 
		return erfccheb(-x) - 1.0;
}

void fft(int isign,const complex<double> *datain,complex<double> *dataout,const int N){
	int r,l,p,k,temp,c;
	r = (int)(log((double)N)/log(2.0));
	complex<double> omega,m,n;
	complex<double> *datatemp = new complex<double>[N];
	for(int i=0;i<N;i++)
		datatemp[i] = datain[i];
	l=0;
//	N=(int) pow(2.0,r);
//	omega.real(cos(-isign*2*PI/N));
//	omega.imag(sin(-isign*2*PI/N));
	omega = complex<double>(cos(-isign*2*PI/N),sin(-isign*2*PI/N));
	if(isign == 1)
		c = 1;
	else if(-isign == 1)
		c = N;
	else
		cout << "error!!!!!" << endl;
	for(int i=0;i<r;i++)
	{
		l++;
		for(int t=0;t<(int)pow(2.0,l);t=t+2)
			for(int s=0;s<N/(int)pow(2.0,l);s++)
			{   
				k=t*N/(int)pow(2.0,l)+s;
				m=datatemp[k];
				n=datatemp[k+N/(int)pow(2.0,l)];
				temp=k>>(r-l);
				p=0;
				for(int j=0;j<r;j++)
				{
					p=p*2+temp%2;
					temp/=2;
				}
				datatemp[k]=m+pow(omega,p)*n;
				datatemp[k+N/(int)pow(2.0,l)]=m-pow(omega,p)*n;
			}
			
	}
	for(int i=0;i<N;i++)
	{
		p=0;
		temp=i>>(r-l);
		for(int j=0;j<r;j++)
		{
			p=p*2+temp%2;
			temp/=2;
		}
		dataout[p]=datatemp[i]/(double)c;
	}
	delete [] datatemp;
}



string num2str(int num)
{
	string res;
	int nn[3];
	for(int i=0;i<3;i++)
	{
		nn[i] = (num/((int)pow(10.0,double(3-i-1))))+48;
		num = num%((int)pow(10.0,double(3-i-1)));
		res += nn[i];
	}
	return res;
}

