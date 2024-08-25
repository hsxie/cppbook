#include<iostream>
#include<complex>
#include<cmath>
#include<fstream>
#include<iomanip>
#include<ctime>   
#include<cstdlib>
#include"parameter.h"
#include"PIC.h"
using namespace std;

//#define PI 3.1415926535897;

/***************Global Variables***********************/


int species_number=2;
Species species[]={Species("electron",N),Species("ion",N)};
double t;
complex<double>E_k[mx];
int it=1;
Particle particle[N][2];
Grid grid(mx,3);
double Epart[N][3],Bpart[N][3];
double Jold[mx+1][3];
const double B_external[dimension]={B_0,0.0,0.0};
double E_t[mx+1];
double v_tilde[N][3][2],v_old[N][3][2];
double x_old[N][2];
double phi[mx+1]={0};
double psi[mx+1]={0};


/***************Global Variables***********************/


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


void Gauss_Seidel1D(double* unew,double* f,const double err,const double dx,const int mx)
{
	double error;
	int k=0;
	double* uold;
	uold = new double[mx+1];
	do
	{
		error =0;
		for(int	i=0;i<=mx;i++)
			uold[i] = unew[i];
		unew[0] = 0.5*(uold[mx-1] + uold[1] - (dx*dx)*f[0]);
		unew[mx] = unew[0];
		for(int	i=1;i<mx;i++)
			unew[i] = 0.5*(unew[i-1] + uold[i+1] - (dx*dx)*f[i]);
		for(int i=0;i<mx;i++)
			if(abs((unew[i]-uold[i])/uold[i])>error)
				error = abs((unew[i]-uold[i])/uold[i]);
		k++;
	}while(error>err);
	delete uold;
	//	cout << "Iteration number :" << k << endl;
//	system("pause");system("pause");
}

void Gauss_Seidel1D_Darwin(double* unew,const double* f,const double err,const double dx,const int mx,const double omega_0,const double* omega)
{
	double error;
	int k=0;
	double* uold;
	uold = new double[mx+1];
	do
	{
		error =0;
		for(int	i=0;i<=mx;i++)
			uold[i] = unew[i];
		unew[0] = (uold[mx-1] + uold[1] - (dx*dx)*f[0]-(dx*dx)*(omega[0]-omega_0)*uold[0])/(2+omega_0*(dx*dx));
		unew[mx] = unew[0];
		for(int	i=1;i<mx;i++)
			unew[i] = (unew[i-1] + uold[i+1] - (dx*dx)*f[i]-(dx*dx)*(omega[i]-omega_0)*uold[i])/(2+omega_0*(dx*dx));
		for(int i=0;i<mx;i++)
			if(abs((unew[i]-uold[i])/uold[i])>error)
				error = abs((unew[i]-uold[i])/uold[i]);
		k++;
	}while(error>err);
	delete uold;
//	cout << "Iteration number :" << k << endl;
//	system("pause");
}
















void clean_file()
{
//MS-DOS	
	system("del *.txt");
    system("del *.dat");
//Linux
/*
	system("rm *.txt");
*/
}
void calc_maxw(double vth,double* v,const int N)
{
//Maxwell distribution
	const int nv=100001;
	double F[nv];
	double vv,vmax;
	double probability;
//vth << v
	vmax = 4*vth;
//  f=1/(vth*sqrt(2*pi))*(exp(-(v/vth)**2/2))
    F[0]=0.0;
	F[nv-1]=1.0;
	for(int i=0;i<nv/2;i++)
	{
		vv=-vmax+2*vmax*double(i-1)/(nv-1);
		F[i]=0.5*(1-erf(-vv/vth));
	}
	for(int i=nv/2;i<nv-1;i++)
	{
		vv=-vmax+2*vmax*double(i-1)/(nv-1);
		F[i]=0.5+0.5*erf(vv/vth);
	}
//  f=1/(2*sqrt(pi))*(exp((v-v0)**2)+exp((v+v0)**2))
	srand((unsigned)time(NULL));
	for(int j=0;j<N;j++)
	{
		probability = (double)rand()/RAND_MAX;
		for(int i=0;i<nv;i++)
			if(probability <= F[i])
			{
				v[j] = -vmax+2*vmax*double(i)/(nv-1);
				break;
			}
	}
}

void distrubution(const double *v,const int N,string name)
{
	const int v_grid = 501;
	int v_index;
	int f[v_grid];
	string file_name;
	double vmax=0.0,vmin=0.0,dv;
	for(int i=0;i<v_grid;i++)
		f[i]=0;
	for(int i=0;i<N;i++)
	{
		vmax = v[i]>vmax?v[i]:vmax;
		vmin = v[i]<vmin?v[i]:vmin;
	}
	dv=(vmax - vmin)/(double)(v_grid-1);
	for(int i=0;i<N;i++)
	{
		v_index = int((v[i]-vmin)/dv);
		f[v_index]++;
	}
	if(name == "electron")
		file_name = "f_e(v).txt";
	else if(name == "ion")
		file_name = "f_i(v).txt";
	ofstream file_f(file_name.c_str(),ios::app);
	for(int i=0;i<v_grid;i++)
		file_f << f[i] << "  " << vmin + i*dv << endl;
}

void distrubution_full(const int N,int it)
{
	const int v_grid = 501;
	int v_index;
	int f[v_grid];
	string file_name;
	double vmax,vmin,dv;
	for(int s=0;s<2;s++)
	{
		vmax=0.0;
		vmin=0.0;
		for(int i=0;i<v_grid;i++)
			f[i]=0;
		for(int i=0;i<N;i++)
		{
			vmax = particle[i][s].v[0]>vmax?particle[i][s].v[0]:vmax;
			vmin = particle[i][s].v[0]<vmin?particle[i][s].v[0]:vmin;
		}
		dv=(vmax - vmin)/(double)(v_grid-1);
		for(int i=0;i<N;i++)
		{
			v_index = int((particle[i][s].v[0]-vmin)/dv);
			f[v_index]++;
		}
		if(s == 0)
			file_name = "f_e(v)_"+num2str(it)+".txt";
		else
			file_name = "f_i(v)_"+num2str(it)+".txt";
		ofstream file_f(file_name.c_str(),ios::app);
		for(int i=0;i<v_grid;i++)
			file_f << f[i] << "  " << vmin + i*dv << endl;
	}
}


void initial(const int condition)
{
	const double x1=Lx/N/10,ks=20*PI;
	double v[N];
	switch(condition)
	{
		case 1:
			for(int s=0;s<species_number;s++)
				for(int i=0;i<N;i++)
				{
					particle[i][s].x = i*Lx/(N-1);
					for(int d=0;d<dimension;d++)
						particle[i][s].v[d] = 0;
				}
			for(int i=0;i<N;i++)
				particle[i][0].x += x1*cos(ks*particle[i][0].x);
			break;
		case 2://ion can mobie
			for(int d=0;d<dimension;d++)
			{
				calc_maxw(1.0,v,N);
				srand((unsigned)time(NULL));
				for(int i=0;i<N;i++)
				{
					particle[i][0].x = (double)rand()/RAND_MAX*Lx;
					particle[i][1].x = i*Lx/(N-1);
					particle[i][0].v[d] = v[i];
					particle[i][1].v[d] = 0;
				}
			}
			distrubution(v,N,"electron");
			break;
		case 3://ion can mobie!
			for(int d=0;d<dimension;d++)
			{
				calc_maxw(sqrt(Te/species[0].mass),v,N);
				distrubution(v,N,"electron");
				srand((unsigned)time(NULL));
				for(int i=0;i<N;i++)
				{
					particle[i][0].x = (double)rand()/RAND_MAX*Lx;
					particle[i][0].v[d] = v[i];
				}
				calc_maxw(sqrt(Ti/species[1].mass),v,N);
				distrubution(v,N,"ion");
				for(int i=0;i<N;i++)
				{
					particle[i][1].x = (double)rand()/RAND_MAX*Lx;
					particle[i][1].v[d] = v[i];
				}
			}
		case 4:
			
			for(int s=0;s<species_number;s++)
				for(int i=0;i<N;i++)
					particle[i][s].x = i*Lx/(N-1);
			srand((unsigned)time(NULL));
			for(int i=0;i<N;i++)
				for(int d=0;d<dimension;d++)
				{
					particle[i][0].v[d] = Te/species[0].mass*(double)rand()/RAND_MAX;
					particle[i][1].v[d] = Ti/species[1].mass*(double)rand()/RAND_MAX;
				}
			for(int i=0;i<N;i++)
				v[i] = particle[i][0].v[0];
			distrubution(v,N,"electron");
			for(int i=0;i<N;i++)
				v[i] = particle[i][1].v[0];
			distrubution(v,N,"ion");
			break;
		case 5://ion can mobie!
			for(int d=0;d<dimension;d++)
			{
				srand((unsigned)time(NULL));
				for(int i=0;i<N;i++)
				{
					particle[i][0].x = (double)rand()/RAND_MAX*Lx;
					particle[i][0].v[d] = 0.00001*cos(2*PI/mx*10*i);
				}
				for(int i=0;i<N;i++)
				{
					particle[i][1].x = (double)rand()/RAND_MAX*Lx;
					particle[i][1].v[d] = 0;
				}
			}
			for(int i=0;i<N;i++)
				v[i] = particle[i][0].v[0];
			distrubution(v,N,"electron");
			for(int i=0;i<N;i++)
				v[i] = particle[i][1].v[0];
			distrubution(v,N,"ion");
			break;
		case 6://ion can mobie!
			for(int s=0;s<species_number;s++)
				for(int i=0;i<N;i++)
					particle[i][s].x = i*Lx/(N-1);
			srand((unsigned)time(NULL));
			for(int d=0;d<dimension;d++)
			{
				
				calc_maxw(sqrt(Te/species[0].mass),v,N);
				for(int i=0;i<N;i++)
					particle[i][0].v[d] = v[i];
				calc_maxw(sqrt(Ti/species[1].mass),v,N);
				for(int i=0;i<N;i++)
					particle[i][1].v[d] = v[i];
			}
			for(int i=0;i<N;i++)
				v[i] = particle[i][0].v[0];
			distrubution(v,N,"electron");
			for(int i=0;i<N;i++)
				v[i] = particle[i][1].v[0];
			distrubution(v,N,"ion");
			break;
		case 7:
			for(int s=0;s<species_number;s++)
				for(int i=0;i<N;i++)
					particle[i][s].x = i*Lx/(N-1);
			for(int i=0;i<N;i++)
			{
				particle[i][0].v[0] = 1.0;
				particle[i][0].v[1] = 5.0;
				particle[i][0].v[2] = 5.0;
			}			
			for(int d=0;d<dimension;d++)
				for(int i=0;i<N;i++)
					particle[i][1].v[d] = 0;
			for(int i=0;i<N;i++)
				v[i] = particle[i][0].v[0];
			distrubution(v,N,"electron");
			for(int i=0;i<N;i++)
				v[i] = particle[i][1].v[0];
			distrubution(v,N,"ion");
			break;


	}
//recycle particles 
	for(int s=0;s<species_number;s++)
		for(int i=0;i<N;i++)
		{
			particle[i][s].x += particle[i][s].v[0]*dt;
			while(particle[i][s].x > Lx)
			{
				particle[i][s].x -= Lx;
			}
			while(particle[i][s].x < 0)
				particle[i][s].x += Lx;
		}
	t = 0.0;
//initial trans_filed
	for(int i=0;i<=mx;i++)
	{
		grid.rho[i]  = 0.0;
		grid.B[i][0] = 0.0;
		grid.B[i][1] = 0.0;
		grid.B[i][2] = 0.0;
		grid.E[i][1] = 0.0;
		grid.E[i][2] = 0.0;
		E_t[i] = 0.0;
	}
}



void particle2grid()
{
	for(int i=0;i<=mx;i++)
		grid.rho[i] = 0.0;

	for(int s=0;s<species_number;s++)
		for(int i=0;i<N;i++)
		{
			int index_left = int(particle[i][s].x/dx);
			int index_right = int(particle[i][s].x/dx)+1;
			double weight = (index_right*dx-particle[i][s].x)/dx;
			if(index_left == mx)
				grid.rho[index_left] += species[s].charge/dx;
			else
			{
				grid.rho[index_left] += weight*species[s].charge/dx;
				grid.rho[index_right] += (1-weight)*species[s].charge/dx;

			}
		}
	grid.rho[0] += grid.rho[mx];
	grid.rho[mx] = grid.rho[0];

}

void current()
{
	for(int s=0;s<species_number;s++)
		for(int d=0;d<dimension;d++)
			for(int i=0;i<N;i++)
				v_tilde[i][d][s] = 1.5*particle[i][s].v[d]-0.5*v_old[i][d][s];

	for(int d=0;d<dimension;d++)
		for(int i=0;i<=mx;i++)
			grid.J[i][d] = 0.0;
	

	for(int s=0;s<species_number;s++)
		for(int d=0;d<dimension;d++)
			for(int i=0;i<N;i++)
			{
				int index_left = int(particle[i][s].x/dx);
				int index_right = int(particle[i][s].x/dx)+1;
				double weight = (index_right*dx-particle[i][s].x)/dx;
				if(index_left == mx)
					grid.J[mx][d] += particle[i][s].v[d];
				else
				{
					grid.J[index_left][d] += weight*v_tilde[i][d][s]*species[s].charge/dx;
					grid.J[index_right][d] += (1-weight)*v_tilde[i][d][s]*species[s].charge/dx;
				}
			}
	for(int d=0;d<dimension;d++)
	{
		grid.J[0][d] += grid.J[mx][d];
		grid.J[mx][d] = grid.J[0][d];
	}
}

void grid2particle()
{

	for(int i=0;i<=mx;i++)
		grid.E[i][0] += E_t[i];
	
	for(int s=0;s<species_number;s++)
		for(int d=0;d<dimension;d++)
			for(int i=0;i<N;i++)
			{
				int index_left = int(particle[i][s].x/dx);
				int index_right = int(particle[i][s].x/dx)+1;
				double weight = (index_right*dx-particle[i][s].x)/dx;
				if(index_left == mx)
				{
					particle[i][s].E[d] = grid.E[index_left][d];
					particle[i][s].B[d] = grid.B[index_left][d] + B_external[d];
				}
				else
				{
                    /*
					if(weight > 0.5)
					{
						particle[i][s].E[d] = grid.E[index_left][d];
						particle[i][s].B[d] = grid.B[index_left][d] + B_external[d] ;
					}
					else
					{
						particle[i][s].E[d] = grid.E[index_right][d];
						particle[i][s].B[d] = grid.B[index_right][d] + B_external[d];
					}
                    */
                    particle[i][s].E[d] =  grid.E[index_left][d]*weight + grid.E[index_right][d]*(1.0-weight);
                    particle[i][s].B[d] = (grid.B[index_left][d]+B_external[d])*weight + (grid.B[index_right][d]+B_external[d])*(1.0-weight);
				}

			}
}

void poisson_eqn(const int method)
{
	complex<double> rho_r[mx],rho_k[mx],E_r[mx],i;
	double kappa;
	i = complex<double>(0,1);
    for(int k=0;k<mx;k++)
		rho_r[k] = complex<double>(grid.rho[k],0);
	fft(1,rho_r,rho_k,mx);
	E_k[0] = 0;

	switch(method)
	{
	case 1:
		for(int k=1;k<mx;k++)
		{
		   kappa = sin(k*2*PI/mx);
		   if(fabs(kappa) < eps)
			   E_k[k] = 0;
		   else
			   E_k[k] = -i*rho_k[k]*dx/kappa/(N/Lx);
		}break;
	case 2:
		for(int k=1;k<=mx/2;k++)
			E_k[k] = rho_k[k]/i/(k*2*PI/Lx)/(N/Lx);
		for(int k=1;k<mx/2;k++)
			E_k[mx/2+k] = rho_k[mx/2+k]/i/((k-mx/2)*2*PI/Lx)/(N/Lx);
		break;
	}
	fft(-1,E_k,E_r,mx);
	for(int j=0;j<mx;j++)
		grid.E[j][0]=E_r[j].real();
	grid.E[mx][0] = grid.E[0][0];
}



void poisson_eqn_GS()
{
	double u[mx+1],f[mx+1];
//	double err=1e-4;
	for(int i=0;i<=mx;i++)
	{
		u[i] = phi[i];
		f[i] = grid.rho[i]/(N/Lx);
	}
	
	Gauss_Seidel1D(u,f,err,dx,mx);
    for(int i=0;i<=mx;i++)
        phi[i] = u[i];
	for(int i=0;i<=mx;i++)
		grid.E[i][0] = (u[i+1]-u[i-1])/2/dx;
	grid.E[0][0] = (u[1]-u[mx-1])/2/dx;
	grid.E[mx][0] = grid.E[0][0];



}


void pusher_x()
{
	for(int s=0;s<species_number;s++)
		for(int i=0;i<N;i++)
			x_old[i][s] = particle[i][s].x;


	for(int s=0;s<species_number;s++)
		for(int i=0;i<N;i++)
		{
			particle[i][s].x += particle[i][s].v[0]*dt;
			while(particle[i][s].x > Lx)
			{
				particle[i][s].x -= Lx;
			}
			while(particle[i][s].x < 0)
				particle[i][s].x += Lx;
		}

}

void pusher_v()
{
//Leap-frog scheme!
// ions imobie
	double *v_prime,*v_plus,*v_sub,*tt,*ss,ttt;
	v_prime = new double[dimension];
	v_plus = new double[dimension];
	v_sub = new double[dimension];
	tt = new double[dimension];
	ss = new double[dimension];
	for(int s=0;s<species_number;s++)
		for(int d=0;d<dimension;d++)
			for(int i=0;i<N;i++)
				v_old[i][d][s] = particle[i][s].v[d];
		
	// specie 0 : electron
	// specie 1 : ion
	for(int s=0;s<species_number;s++)
	{
		double qmt=species[s].charge/species[s].mass*dt/2;
		for(int i=0;i<N;i++)
		{
			ttt = 0;
			for(int d=0;d<dimension;d++)
			{
				tt[d] = qmt*particle[i][s].B[d];
				ttt += pow(tt[d],2);
				v_sub[d] = particle[i][s].v[d] + qmt*particle[i][s].E[d];
			}
			for(int d=0;d<dimension;d++)
				ss[d] = 2*tt[d]/(1+ttt);

			v_prime[0] = v_sub[0] + v_sub[1]*tt[2] - v_sub[2]*tt[1];
			v_prime[1] = v_sub[1] + v_sub[2]*tt[0] - v_sub[0]*tt[2];
			v_prime[2] = v_sub[2] + v_sub[0]*tt[1] - v_sub[1]*tt[0];
			v_plus[0]  = v_sub[0] + v_prime[1]*ss[2] - v_prime[2]*ss[1];
			v_plus[1]  = v_sub[1] + v_prime[2]*ss[0] - v_prime[0]*ss[2];
			v_plus[2]  = v_sub[2] + v_prime[0]*ss[1] - v_prime[1]*ss[0];

			for(int d=0;d<dimension;d++)
				particle[i][s].v[d] = v_plus[d] +qmt*particle[i][s].E[d];
		}
	}
	delete [] v_prime,v_plus,v_sub,tt,ss;
	t=t+dt;
}




void pusher_v_0()
{
	for(int s=0;s<species_number;s++)
	{
		double qmt = species[s].charge/species[s].mass*dt/2;
		for(int i=0;i<N;i++)
			for(int d=0;d<dimension;d++)
				particle[i][s].v[d] -=  qmt*particle[i][s].E[d]/ species[s].mass*dt/2;
	}
}



void diagnosis()
{
	//Ek:kinetic       energy
    //Ee:electrostatic energy
    //Em:magnetic      energy
    double Ek,Ee,Em;
	Ek = 0;;
	Ee = 0;
	Em = 0;
	ofstream file_e("energy.txt",ios::app);
	for(int s=0;s<species_number;s++)
		for(int d=0;d<dimension;d++)
			for(int i=0;i<N;i++)
				Ek += 0.5*(species[s].mass*particle[i][s].v[d]*particle[i][s].v[d]);

	for(int d=1;d<dimension;d++)
		for(int i=0;i<mx;i++)
			Ee += 0.5*grid.E[i][d]*grid.E[i][d]*dx*(N/Lx);

	for(int i=0;i<mx;i++)
		Ee += 0.5*(grid.E[i][0])*(grid.E[i][0])*dx*(N/Lx);

	
	//only magnetic  
	for(int d=0;d<dimension;d++)
	  for(int i=0;i<mx;i++)
		  Em += 0.5*grid.B[i][d]*grid.B[i][d]*dx*(N/Lx)*alpha;
	//	  Ee += 0.5*grid.rho[i]*phi(i)*dx
	//  output energy 
	file_e << setprecision(8) << Ek << "  "<< Ee << "  " << Em << "  " << t << endl;
	ofstream file_rv("rv.txt",ios::app);
	//  output (N/2)th paricle's position and velocity
	file_rv << setprecision(8) << particle[N/2][0].x << "  " << particle[N/2][0].v[0] << "  " <<particle[N/2][0].v[1] << "  " << particle[N/2][0].v[2] << "  " << particle[N/2][1].x << "  " << particle[N/2][1].v[0] << "  " << particle[N/2][1].v[1] << "  " <<particle[N/2][1].v[2] << "  " << t <<endl;
}





void dispersion()
{
//dispersion relation
	static int jkt =1; 
	if(it%perstep == 0 && jkt <= nx)
	{
        ofstream file_Ekt("E_k_t.txt",ios::app);
		ofstream file_Eyxt("Ey_x_t.txt",ios::app);
		ofstream file_Exxt("Ex_x_t.txt",ios::app);
		ofstream file_Ezxt("Ez_x_t.txt",ios::app);
		ofstream file_Byxt("By_x_t.txt",ios::app);
		ofstream file_Bzxt("Bz_x_t.txt",ios::app);
		for(int i=1;i<=mx/2;i++)
			file_Ekt << setprecision(8) << E_k[i].real() << "  " ;
		for(int i=0;i<mx;i++)
		{
			file_Exxt << setprecision(8) << grid.E[i][0] << "  " ;
			file_Eyxt << setprecision(8) << grid.E[i][1] << "  " ;
			file_Ezxt << setprecision(8) << grid.E[i][2] << "  " ;
			file_Byxt << setprecision(8) << grid.B[i][1] << "  " ;
			file_Bzxt << setprecision(8) << grid.B[i][2] << "  " ;
		}
		file_Ekt << endl;
		file_Eyxt << endl;
		file_Exxt << endl;
		file_Ezxt << endl;
		file_Byxt << endl;
		file_Bzxt << endl;
        if(jkt == 1)
        {
            cout << "The diagnosis is started!! the beginning time is : " << t << endl;
            cout << "per how many steps to recrd                      : " << perstep << endl;
        }
        if(jkt == nx)
            cout << "The diagnosis is completed!! the final time is   : " << t << endl;
     jkt++;
	}
}
void record_data(int it)
{
	string str = "particle_info_"+num2str(it)+".txt";
	ofstream file_p(str.c_str(),ios::app);
	str = "grid_info_"+num2str(it)+".txt";
	ofstream file_g(str.c_str(),ios::app);
	file_p << "x    " << "vx    " << "vy    " << "vz    " << "time    " <<endl;
	for(int i=0;i<N;i++)
	{
		file_p <<  particle[i][0].x << "  ";
		for(int d=0;d<dimension;d++)
			file_p << particle[i][0].v[d] << "  ";
		file_p << t << endl;
	}
	file_g << "    Ex    " << "    Ey    " << "    Ez    " << "    Bx    " << "    By    " << "    Bz    " << "    Jx    " << "    Jy    " << "    Jz    " << "    rho    " << "    time" << endl;
	for(int i=0;i<=mx;i++)
	{
		for(int d=0;d<dimension;d++)
			file_g << setprecision(8)<< grid.E[i][d] << "  ";
		for(int d=0;d<dimension;d++)
			file_g << setprecision(8)<< grid.B[i][d] << "  ";
		for(int d=0;d<dimension;d++)
			file_g << setprecision(8)<< grid.J[i][d] << "  ";
		file_g << setprecision(8)<< grid.rho[i] << "  " << t << endl;;
	}
	distrubution_full(N,it);

}



void Darwin_full()
{
	
	double na=(N/Lx)*alpha;
	double* W=new double[mx+1];
	double* omega=new double[mx+1];
	double** omega_u;
	omega_u = new double*[mx+1];
	for(int i=0;i<=mx;i++)
		omega_u[i] = new double[dimension];


	double* f=new double[mx+1];
	double* u=new double[mx+1];

	for(int i=1;i<mx;i++)
		f[i] = (grid.J[i+1][2] - grid.J[i-1][2])/(2*dx)/na; 
	f[0] =  (grid.J[1][2] - grid.J[mx-1][2])/(2*dx)/na;
	f[mx] = f[0];
	for(int i=0;i<=mx;i++)
		u[i] = grid.B[i][1];
	Gauss_Seidel1D(u,f,err,dx,mx);
	for(int i=0;i<=mx;i++)
		grid.B[i][1] = u[i];
	for(int i=1;i<mx;i++)
		f[i] = -(grid.J[i+1][1] - grid.J[i-1][1])/(2*dx)/na; 
	f[0] =  -(grid.J[1][1] - grid.J[mx-1][1])/(2*dx)/na;
	f[mx] = f[0];
	for(int i=0;i<=mx;i++)
		u[i] = grid.B[i][2];
	Gauss_Seidel1D(u,f,err,dx,mx);
	for(int i=0;i<=mx;i++)
		grid.B[i][2] = u[i];
	
	
	for(int i=0;i<=mx;i++)
	{
		W[i] = omega[i] =0.0;
		for(int d=0;d<dimension;d++)
			 omega_u[i][d] = 0.0;
	}
	for(int s=0;s<species_number;s++)
	{
		double q2m = species[s].charge*species[s].charge/species[s].mass/dx;
		for(int i=0;i<N;i++)
		{
			int index_left=int(particle[i][s].x/dx);
			int index_right=int(particle[i][s].x/dx)+1;
			double weight = (index_right*dx-particle[i][s].x)/dx;
			if(weight >0.5)
				weight = 0.0;
			else
				weight =1.0;

			if(index_left == mx)
			{
				W[mx] +=  v_tilde[i][0][s]* v_tilde[i][1][s]*species[s].charge/dx;
				omega[mx] += q2m;
				for(int d=0;d<dimension;d++)
					omega_u[mx][d] += q2m*particle[i][s].v[d];
			}
			else
			{
				W[index_left] += weight*v_tilde[i][0][s]* v_tilde[i][1][s]*species[s].charge/dx;
				W[index_right] += (1-weight)* v_tilde[i][0][s]* v_tilde[i][1][s]*species[s].charge/dx;
				omega[index_left] += weight*q2m;
				omega[index_right] += (1-weight)*q2m;
				for(int d=0;d<dimension;d++)
				{
					omega_u[index_left][d] += weight*q2m*v_tilde[i][d][s];
					omega_u[index_right][d] += (1-weight)*q2m*v_tilde[i][d][s];
				}
			}
		}
	}


	W[0] += W[mx];
	W[mx] = W[0];
	omega[0] += omega[mx];
	omega[mx] = omega[0];
	for(int d=0;d<dimension;d++)
	{
		omega_u[0][d] += omega_u[mx][d];
		omega_u[mx][d] = omega_u[0][d];
	}


	double omega_0,omega_max,omega_min;
	omega_max = 0,0;
	
	for(int i=0;i<=mx;i++)
		omega[i] = omega[i]/na;
	
	
	for(int i=0;i<=mx;i++)
		omega_max = (omega[i]>omega_max)?omega[i]:omega_max;
	omega_min = omega_max;
	for(int i=0;i<=mx;i++)
		omega_min = (omega[i]<omega_min)?omega[i]:omega_min;
	omega_0=0.5*(omega_max+omega_min);

	f[0] = (-(W[1]-W[mx-1])/(2*dx) + omega_u[0][2]*(grid.B[0][0]+B_external[0])-omega_u[0][0]*(grid.B[0][2]+B_external[2]))/na;
	f[mx] = f[0];
	for(int i=1;i<mx;i++)
		f[i] = (-(W[i+1]-W[i-1])/(2*dx) + omega_u[i][2]*(grid.B[i][0]+B_external[0])-omega_u[i][0]*(grid.B[i][2]+B_external[2]))/na;
	for(int i=0;i<=mx;i++)
		u[i] = grid.E[i][1];
	Gauss_Seidel1D_Darwin(u,f,err,dx,mx,omega_0,omega);
	for(int i=0;i<=mx;i++)
		grid.E[i][1] = u[i];



	for(int i=0;i<=mx;i++)
		W[i] = 0.0;
	
	for(int s=0;s<species_number;s++)	
		for(int i=0;i<N;i++)
        {
				int index_left=int(particle[i][s].x/dx);
				int index_right=int(particle[i][s].x/dx)+1;
				double weight = (index_right*dx-particle[i][s].x)/dx;
				if(weight >0.5)
					weight = 0.0;
				else
					weight =1.0;

				if(index_left == mx)
					W[mx] +=  v_tilde[i][0][s]*v_tilde[i][2][s]*species[s].charge/dx;
				else
				{
					W[index_left]  += weight*v_tilde[i][0][s]*v_tilde[i][2][s]*species[s].charge/dx;
					W[index_right] += (1-weight)*v_tilde[i][0][s]*v_tilde[i][2][s]*species[s].charge/dx;
				}
		}
	W[0] += W[mx];
	W[mx] = W[0];


	
	f[0] = (-(W[1]-W[mx-1])/(2*dx) + omega_u[0][0]*(grid.B[0][1]+B_external[1])-omega_u[0][1]*(grid.B[0][0]+B_external[0]))/na;
	f[mx] = f[0];
	for(int i=1;i<mx;i++)
		f[i] = (-(W[i+1]-W[i-1])/(2*dx) + omega_u[i][0]*(grid.B[i][1]+B_external[1])-omega_u[i][1]*(grid.B[i][0]+B_external[0]))/na;
	for(int i=0;i<=mx;i++)
		u[i] = grid.E[i][2];
	Gauss_Seidel1D_Darwin(u,f,err,dx,mx,omega_0,omega);
	for(int i=0;i<=mx;i++)
		grid.E[i][2] = u[i];


	


	for(int i=0;i<=mx;i++)
		W[i] = 0.0;
	
	for(int s=0;s<species_number;s++)	
		for(int i=0;i<N;i++)
		{
			int index_left=int(particle[i][s].x/dx);
			int index_right=int(particle[i][s].x/dx)+1;
			double weight = (index_right*dx-particle[i][s].x)/dx;
			if(weight >0.5)
				weight = 0.0;
			else
				weight =1.0;

			if(index_left == mx)
				W[mx] +=  v_tilde[i][0][s]*v_tilde[i][0][s]*species[s].charge/dx;
			else
			{
				W[index_left] += weight*v_tilde[i][0][s]* v_tilde[i][0][s]*species[s].charge/dx;
				W[index_right] += (1-weight)*v_tilde[i][0][s]* v_tilde[i][0][s]*species[s].charge/dx;
			}
		}
	W[0] += W[mx];
	W[mx] = W[0];


	f[0] = (-(W[1]-W[mx-1])/(2*dx) + omega_u[0][1]*(grid.B[0][2]+B_external[2])-omega_u[0][2]*(grid.B[0][1]+B_external[1])+omega[0]*((N/Lx)*alpha)*grid.E[0][0])/na;
	f[mx] = f[0];
	for(int i=1;i<mx;i++)
		f[i] = ( -(W[i+1]-W[i-1])/(2*dx)
                 + omega[i]*((N/Lx)*alpha)*grid.E[i][0]
                 + omega_u[i][1]*(grid.B[i][2]+B_external[2])-omega_u[i][2]*(grid.B[i][1]+B_external[1])
               )/na;
	
	for(int i=0;i<=mx;i++)
		u[i] = E_t[i];
    Gauss_Seidel1D_Darwin(u,f,err,dx,mx,omega_0,omega);
//	Gauss_Seidel1D_Darwin(u,f,err,dx,mx,omega_0,omega,N,Lx,alpha);
	for(int i=0;i<=mx;i++)
		E_t[i] = u[i];


	for(int i=1;i<mx;i++)
		f[i] = (E_t[i+1] - E_t[i-1])/2/dx;
	f[0] = (E_t[1] - E_t[mx-1])/2/dx;
	f[mx]=f[0];
    for(int i=0;i<=mx;i++)
		u[i] = psi[i];
	Gauss_Seidel1D(u,f,err,dx,mx);
    for(int i=0;i<=mx;i++)
		psi[i] = u[i];
	for(int i=0;i<=mx;i++)
		E_t[i] -= (u[i+1] - u[i-1])/2/dx;
	E_t[0] -= (u[1] - u[mx-1])/2/dx;
	E_t[mx]=E_t[0];




	/*
	double E_t_max,E_t_min;
	E_t_max = u[0];
	E_t_min = u[0];
	for(int i=0;i<=mx;i++)
	{
		E_t_max = (u[i]>E_t_max)?u[i]:E_t_max;
		E_t_min = (u[i]<E_t_min)?u[i]:E_t_min;
	}
	for(int i=0;i<=mx;i++)
		E_t[i] = 0.5*(E_t_max+E_t_min);
*/

//	cout << E_t[20] << "   " << E_t[22] << endl;
//	system("pause");



}

void wave_heating()
{
	double B_0,v_0,omega_A,k_A,v_A;
	B_0 = 0.2*B_external[0];
	v_0 = 0.001;
	k_A = 10;
	v_A = sqrt(B_external[0]/sqrt(species[1].mass)*alpha);
	v_A = 1.0;
	omega_A = v_A*k_A;
//	cout << grid.B[mx][1] << "  " << B_0*cos(k_A*Lx+omega_A*t) << endl;
	grid.B[mx][1] += B_0*cos(k_A*Lx+omega_A*t);
//	for(int i=0;i<N;i++)
//	{
//		particle[i][0].v[1] += v_0*cos(k_A*particle[i][0].x+omega_A*t);
//		particle[i][1].v[1] += v_0*0.01*cos(k_A*particle[i][0].x+omega_A*t);
//	}
	grid.B[mx][1] += B_0*cos(k_A*Lx+omega_A*t);
//	cout << v_A << endl;
}

void parameter_output()
{
    ofstream parameter_output("parameter.dat");
    parameter_output << N << " " << mx << " " << Lx << " " << B_0 << " " << mi_div_me << " " << dt << " " << (nx-1)*perstep*dt << " " << nx << " " << alpha << endl;
    parameter_output.close();
    ofstream parameter_list("parameter_list.dat");
    parameter_list << "------------------------parameter list------------------------" << endl;
    parameter_list << "Number of particles                 : " << N   << endl;
    parameter_list << "Number of grid on x direction       : " << mx  << endl;
    parameter_list << "Length of simulation region         : " << Lx  << endl;
    parameter_list << "Magnitude of guilding field         : " << B_0 << endl;
    parameter_list << "Mass ratio(ion mass/electron mass)  : " << mi_div_me << endl;
    parameter_list << "Time step                           : " << dt << endl;
    parameter_list << "Speed of light(c/v_th)              : " << sqrt(alpha) << endl;
    parameter_list.close();


}
