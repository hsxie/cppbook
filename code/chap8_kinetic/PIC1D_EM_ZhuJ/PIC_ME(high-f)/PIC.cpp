#include<iostream>
#include<complex>
#include<cmath>
#include<fstream>
#include<iomanip>
#include<ctime>   
#include<cstdlib>
#include"PIC.h"
using namespace std;

//#define PI 3.1415926535897;




class Particle{
public:int d;
	   double x;
	   double *v,*E,*B,*B_old;
	   Particle(int dimension):x(0),d(dimension)
	   {
		   v = new double[d];
		   E = new double[d];
		   B = new double[d];
		   B_old = new double[d];
		   for(int i=0;i<d;i++)
			   B_old[i] = B[i] = E[i] = v[i] = 0.0;
	   }
	   Particle():x(0),d(3)
	   {
		   v = new double[d];
		   E = new double[d];
		   B = new double[d];
		   B_old = new double[d];
		   for(int i=0;i<d;i++)
			    B_old[i] = B[i] = E[i] = v[i] = 0.0;
	   }
	   ~Particle()
	   {
		   delete [] v;
		   delete [] B;
		   delete [] E;
		   delete [] B_old;
	   }
};

class Electron:public Particle{
public:const int charge;
	   const double mass;
	   Electron():charge(-1),mass(1.0){}
};

class Ion:public Particle{
public:const int charge;
	   const double mass;
	   Ion():charge(1),mass(1836.0){}
};

class Grid{
private:int m,d;
public:double **E,**B,**J,***v_fluid,*rho;
	   Grid(int mx,int dimension):m(mx),d(dimension)
	   {
		   E = new double*[m+1];
		   B = new double*[m+1];
		   J = new double*[m+1];
		   v_fluid = new double**[m+1];
		   rho = new double[m+1];
		   for(int i=0;i<=m;i++)
		   {
			   E[i] = new double[d];
			   B[i] = new double[d];
			   J[i] = new double[d];
			   v_fluid[i] = new double*[d];
			   for(int j=0;j<d;j++)
				   v_fluid[i][j] = new double[2];

		   }

		   for(int i=0;i<=m;i++)
			   for(int j=0;j<d;j++)
			   {
				   J[i][j]=B[i][j]=E[i][j]=0.0;
				   for(int s=0;s<2;s++)
					   v_fluid[i][j][s] = 0.0;
			   }
	   }
	   ~Grid()
	   {
		   for(int i=0;i<m;i++)
		   {
			   delete [] E[i];
			   delete [] B[i];
			   delete [] J[i];
			   for(int j=0;j<d;j++)
				   delete [] v_fluid[i][j];
		   }

		delete [] E;
		delete [] B;
		delete [] J;
		delete [] v_fluid;
	   }
};









/***************Global Variables***********************/

const int species=1;
const int N=10000;
const int Nspecies[2]={10000,10000};
const int mx=128,my=1,mz=1,dimension=3;
const double Lx=128,dx=Lx/mx,dt=0.1;
double t;
double Ek,Ee,Em;
double Te=1.0,Ti=1.0;
complex<double>E_k[mx];
int it=1;
Electron electron[N];
Ion ion[N];
Grid grid(mx,3);
double Epart[N][3],Bpart[N][3];
double Bpart_old[N][3];
//alpha = (c/vth)^2
double alpha=100;

/***************Global Variables***********************/







void clean_file()
{
//MS-DOS	
	system("del *.txt");
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
	ofstream file_f("f.txt");
	for(int i=0;i<nv;i++)
		file_f << F[i] << endl;
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

void distrubution(const double *v,const int N)
{
	const int v_grid = 501;
	int v_index;
	int f[v_grid];
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
	ofstream file_f("f(v).txt");
	for(int i=0;i<v_grid;i++)
		file_f << f[i] << "  " << vmin + i*dv << endl;
}


void initial(const int condition)
{
	const double x1=Lx/N/10,ks=20*PI;
	double v[N];
	switch(condition)
	{
		case 1:
			for(int i=0;i<N;i++)
			{
				ion[i].x = electron[i].x = i*Lx/(N-1);
				ion[i].v[0] = electron[i].v[0] = 0;
			}
			for(int i=0;i<N;i++)
				electron[i].x += x1*cos(ks*electron[i].x);
			break;
		case 2://ion can mobie
			for(int d=0;d<dimension;d++)
			{
				calc_maxw(1.0,v,N);
				srand((unsigned)time(NULL));
				for(int i=0;i<N;i++)
				{
					ion[i].x = i*Lx/(N-1);
					electron[i].x = (double)rand()/RAND_MAX*Lx;
					ion[i].v[d] = 0;
					electron[i].v[d] = v[i];
				}
			}
			distrubution(v,N);
			break;
		case 3://ion can mobie!
			for(int d=0;d<dimension;d++)
			{
				calc_maxw(sqrt(Te/electron[0].mass),v,N);
				srand((unsigned)time(NULL));
				for(int i=0;i<N;i++)
				{
					electron[i].x = (double)rand()/RAND_MAX*Lx;
					electron[i].v[d] = v[i];
				}
				calc_maxw(sqrt(Ti/ion[0].mass),v,N);
				srand((unsigned)time(NULL));
				for(int i=0;i<N;i++)
				{
					ion[i].x = i*Lx/(N-1);
					ion[i].v[d] = v[i];
				}
			}
			distrubution(v,N);
			break;
	}
//recycle particles 
	for(int i=0;i<N;i++)
	{
	  while(electron[i].x > Lx)//if particle exactly at x=Lx,put it in x=0 
		  electron[i].x -= Lx;
	  while(electron[i].x < 0)
		  electron[i].x += Lx;
	  while(ion[i].x > Lx)//if particle exactly at x=Lx,put it in x=0 
		  ion[i].x -= Lx;
	  while(ion[i].x < 0)
		  ion[i].x += Lx;
	}
	t = 0.0;
//initial trans_filed
	for(int i=0;i<=mx;i++)
	{
		grid.rho[i]  = 0.0;
		grid.B[i][0] = 1.2;//guiding center
		grid.B[i][1] = 0.0;
		grid.B[i][2] = 0.0;
		grid.E[i][1] = 0.0;
		grid.E[i][2] = 0.0;
	}
}



void particle2grid()
{
	for(int i=0;i<=mx+1;i++)
		grid.rho[i] = 0.0;
	for(int i=0;i<=mx;i++)
		for(int d=0;d<dimension;d++)
			for(int s=0;s<2;s++)
				grid.v_fluid[i][d][s] = 0.0;
	for(int i=0;i<N;i++)
	{
		int Xe=int(electron[i].x/dx);
		int Xi=int(ion[i].x/dx);
		if(Xe == mx)
			grid.rho[Xe] += electron[i].charge/dx;
		else if(Xi == mx)
			grid.rho[Xi] += ion[i].charge/dx;
		else
		{
			grid.rho[Xe] += ((Xe+1)*dx-electron[i].x)/dx*electron[i].charge/dx;
			grid.rho[Xe+1] += (electron[i].x-Xe*dx)/dx*electron[i].charge/dx;
			grid.rho[Xi] += ((Xi+1)*dx-ion[i].x)/dx*ion[i].charge/dx;
			grid.rho[Xi+1] += (ion[i].x-Xi*dx)/dx*ion[i].charge/dx;
		}
	}
	for(int d=0;d<dimension;d++)
		for(int i=0;i<N;i++)
		{
			int Xe=int(electron[i].x/dx);
			int Xi=int(ion[i].x/dx);
			if(Xe == mx)
				grid.v_fluid[mx][d][0] += electron[i].v[d];
			else if(Xi == mx)
				grid.v_fluid[mx][d][1] += ion[i].v[d];
			else
			{
				grid.v_fluid[Xe][d][0] += ((Xe+1)*dx-electron[i].x)/dx*electron[i].v[d];
				grid.v_fluid[Xe+1][d][0] += (electron[i].x-Xe*dx)/dx*electron[i].v[d];
				grid.v_fluid[Xi][d][1] += ((Xi+1)*dx-ion[i].x)/dx*ion[i].v[d];
				grid.v_fluid[Xi+1][d][1] += (ion[i].x-Xi*dx)/dx*ion[i].v[d];
			}
		}
	grid.rho[0] += grid.rho[mx];
	grid.rho[mx] = grid.rho[0];
	for(int s=0;s<2;s++)
		for(int d=0;d<dimension;d++)
		{
			grid.v_fluid[0][d][s] += grid.v_fluid[mx][d][s];
			grid.v_fluid[mx][d][s] = grid.v_fluid[0][d][s];
		}
}

void current()
{
	for(int d=0;d<dimension;d++)
		for(int i=0;i<=mx;i++)
			grid.J[i][d] =(grid.v_fluid[i][d][0]*electron[i].charge+grid.v_fluid[i][d][1]*ion[i].charge)/dx;
}

void trans_Efield()
{
	for(int i=1;i<mx;i++)
		{
			grid.E[i][1] += -alpha*dt*(grid.B[i+1][2]-grid.B[i-1][2])/(2*dx)-dt*grid.J[i][1]/(N/Lx);
			grid.E[i][2] +=  alpha*dt*(grid.B[i+1][1]-grid.B[i-1][1])/(2*dx)-dt*grid.J[i][2]/(N/Lx);
		}
	grid.E[0][1] += -alpha*dt*(grid.B[1][2]-grid.B[mx-1][2])/(2*dx)-dt*grid.J[0][1]/(N/Lx);
	grid.E[mx][1] = grid.E[0][1];
	grid.E[0][2] +=  alpha*dt*(grid.B[1][1]-grid.B[mx-1][1])/(2*dx)-dt*grid.J[0][2]/(N/Lx);
	grid.E[mx][2] = grid.E[0][2];
}
void trans_Bfield()
{
	
	for(int i=1;i<mx;i++)
		{
			grid.B[i][1] += dt*(grid.E[i+1][2]-grid.E[i-1][2])/(2*dx);
			grid.B[i][2] += -dt*(grid.E[i+1][1]-grid.E[i-1][1])/(2*dx);
		}
	grid.B[0][1] += dt*(grid.E[1][2]-grid.E[mx-1][2])/(2*dx);
	grid.B[mx][1]=grid.B[0][1];
	grid.B[0][2] += -dt*(grid.E[1][1]-grid.E[mx-1][1])/(2*dx);
	grid.B[mx][2]=grid.B[0][2];
	

}
void trans_Bfield_0()
{
	
	for(int i=1;i<mx;i++)
		{
			grid.B[i][1] += (dt*(grid.E[i+1][2]-grid.E[i-1][2])/(2*dx))/2;
			grid.B[i][2] += -(dt*(grid.E[i+1][1]-grid.E[i-1][1])/(2*dx))/2;
		}
	grid.B[0][1] += (dt*(grid.E[1][2]-grid.E[mx-1][2])/(2*dx))/2;
	grid.B[mx][1]=grid.B[0][1];
	grid.B[0][2] += -(dt*(grid.E[1][1]-grid.E[mx-1][1])/(2*dx))/2;
	grid.B[mx][2]=grid.B[0][2];
	

}
void grid2particle()
{
	for(int d=0;d<dimension;d++)
		for(int i=0;i<N;i++)
		{
			electron[i].B_old[d] = electron[i].B[d];
			ion[i].B_old[d] = ion[i].B[d];
		}
	for(int d=0;d<dimension;d++)
		for(int i=0;i<N;i++)
		{
			int Xe=int(electron[i].x/dx);
			int Xi=int(ion[i].x/dx);
			if(Xe == mx)
			{
				electron[i].E[d] = grid.E[Xe][d];
				electron[i].B[d] = grid.B[Xe][d];
			}
			else if(Xi == mx)
			{
				ion[i].E[d] = grid.E[Xi][d];
				ion[i].B[d] = grid.B[Xi][d];
			}
			else
			{
				electron[i].E[d] = ((Xe+1)*dx-electron[i].x)/dx*grid.E[Xe][d] + (electron[i].x-Xe*dx)/dx*grid.E[Xe+1][d];
				electron[i].B[d]  = ((Xe+1)*dx-electron[i].x)/dx*grid.B[Xe][d] + (electron[i].x-Xe*dx)/dx*grid.B[Xe+1][d];
				ion[i].E[d] = ((Xi+1)*dx-ion[i].x)/dx*grid.E[Xi][d] + (ion[i].x-Xi*dx)/dx*grid.E[Xi+1][d];
				ion[i].B[d]  = ((Xi+1)*dx-ion[i].x)/dx*grid.B[Xi][d] + (ion[i].x-Xi*dx)/dx*grid.B[Xi+1][d];
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


void pusher_0()
{
//set velocity at time = -1/2
	double *v_prime,*v_plus,*v_sub,*tt,*s,ttt;
	v_prime = new double[dimension];
	v_plus = new double[dimension];
	v_sub = new double[dimension];
	tt = new double[dimension];
	s = new double[dimension];
	ttt = 0;
	double q,m,Epart[dimension],Bpart[dimension],Bpart_old[dimension],v[dimension];
	for(int ss=0;ss<species;ss++)
		for(int i=0;i<Nspecies[ss];i++)
		{
			switch(ss)
			{
			case 0:
				q = electron[i].charge;
				m = electron[i].mass;
				for(int j=0;j<dimension;j++)
				{
					v[j] = electron[i].v[j];
					Epart[j] = electron[i].E[j];
					Bpart[j] = electron[i].B[j];
					Bpart_old[j] = electron[i].B_old[j];
				}
				break;
			case 1:
				q = ion[i].charge;
				m = ion[i].mass;
				for(int j=0;j<dimension;j++)
				{
					v[j] = ion[i].v[j];
					Epart[j] = ion[i].E[j];
					Bpart[j] = ion[i].B[j];
					Bpart_old[j] = ion[i].B_old[j];
				}
				break;
			}
			for(int d=0;d<dimension;d++)
			{
				// specie 0 : electron
				// specie 1 : ion
				v_sub[d] = v[d] + (q*Epart[d]/m*dt/2)/2;
				tt[d] = (q*(Bpart[d]+Bpart_old[d])/2/m*dt/2)/2;
				ttt += pow(tt[d],2);
	//			electron[i].v[d] += electron[i].charge*Epart[i][d]/electron[i].mass*dt;
			}
			for(int d=0;d<dimension;d++)
				s[d] = 2*tt[d]/(1+ttt);
			v_prime[0] = v_sub[0] + v_sub[1]*tt[2] - v_sub[2]*tt[1];
			v_prime[1] = v_sub[1] + v_sub[2]*tt[0] - v_sub[0]*tt[2];
			v_prime[2] = v_sub[2] + v_sub[0]*tt[1] - v_sub[1]*tt[0];
			v_plus[0]  = v_sub[0] + v_prime[1]*s[2] - v_prime[2]*s[1];
			v_plus[1]  = v_sub[1] + v_prime[2]*s[0] - v_prime[0]*s[2];
			v_plus[2]  = v_sub[2] + v_prime[0]*s[1] - v_prime[1]*s[0];
			switch(ss)
			{
			case 0:
				for(int d=0;d<dimension;d++)
				{
					electron[i].v[d] = v_plus[d] + (q*Epart[d]/m*dt/2)/2;
				}
				break;
			case 1:
				for(int d=0;d<dimension;d++)
					ion[i].v[d] = v_plus[d] + (q*Epart[d]/m*dt/2)/2;
				break;
			}
			
		}

	delete [] v_prime,v_plus,v_sub,tt,s;
}


void pusher()
{
//Leap-frog scheme!
// ions imobie
	double *v_prime,*v_plus,*v_sub,*tt,*s,ttt;
	v_prime = new double[dimension];
	v_plus = new double[dimension];
	v_sub = new double[dimension];
	tt = new double[dimension];
	s = new double[dimension];
	

	double q,m,Epart[dimension],Bpart[dimension],Bpart_old[dimension],v[dimension];
	for(int ss=0;ss<species;ss++)
		for(int i=0;i<Nspecies[ss];i++)
		{
			switch(ss)
			{
			case 0:
				q = electron[i].charge;
				m = electron[i].mass;
				for(int j=0;j<dimension;j++)
				{
					v[j] = electron[i].v[j];
					Epart[j] = electron[i].E[j];
					Bpart[j] = electron[i].B[j];
					Bpart_old[j] = electron[i].B_old[j];
				}
				break;
			case 1:
				q = ion[i].charge;
				m = ion[i].mass;
				for(int j=0;j<dimension;j++)
				{
					v[j] = ion[i].v[j];
					Epart[j] = ion[i].E[j];
					Bpart[j] = ion[i].B[j];
					Bpart_old[j] = ion[i].B_old[j];
				}
				break;
			}
			ttt = 0;
			for(int d=0;d<dimension;d++)
			{
				// specie 0 : electron
				// specie 1 : ion
				v_sub[d] = v[d] + q*Epart[d]/m*dt/2;
				tt[d] = q*(Bpart[d]+Bpart_old[d])/2/m*dt/2;
				ttt += pow(tt[d],2);
			}
			for(int d=0;d<dimension;d++)
				s[d] = 2*tt[d]/(1+ttt);
			v_prime[0] = v_sub[0] + v_sub[1]*tt[2] - v_sub[2]*tt[1];
			v_prime[1] = v_sub[1] + v_sub[2]*tt[0] - v_sub[0]*tt[2];
			v_prime[2] = v_sub[2] + v_sub[0]*tt[1] - v_sub[1]*tt[0];
			v_plus[0]  = v_sub[0] + v_prime[1]*s[2] - v_prime[2]*s[1];
			v_plus[1]  = v_sub[1] + v_prime[2]*s[0] - v_prime[0]*s[2];
			v_plus[2]  = v_sub[2] + v_prime[0]*s[1] - v_prime[1]*s[0];
			switch(ss)
			{
			case 0:
				for(int d=0;d<dimension;d++)
					electron[i].v[d] = v_plus[d] + q*Epart[d]/m*dt/2;
				break;
			case 1:
				for(int d=0;d<dimension;d++)
					ion[i].v[d] = v_plus[d] + q*Epart[d]/m*dt/2;
				break;
			}
			
		}
	for(int i=0;i<N;i++)
	{
		electron[i].x += electron[i].v[0]*dt;
		while(electron[i].x > Lx)
		{
			electron[i].x -= Lx;
		}
		while(electron[i].x < 0)
			electron[i].x += Lx;
	}

	for(int i=0;i<N;i++)
	{
		ion[i].x += ion[i].v[0]*dt;
		while(ion[i].x > Lx)
		{
			ion[i].x -= Lx;
		}
		while(ion[i].x < 0)
			ion[i].x +=Lx;
	}
	delete [] v_prime,v_plus,v_sub,tt,s;
	t=t+dt;
}


void diagnosis()
{
	//Ek:kinetic ernergy   Ee:electrostatic energy Em:magnetic energy
	Ek = 0;;
	Ee = 0;
	Em = 0;
	ofstream file_e("energy.txt",ios::app);
	for(int d=0;d<dimension;d++)
		for(int i=0;i<N;i++)
			Ek += 0.5*(electron[i].mass*electron[i].v[d]*electron[i].v[d]+ion[i].mass*ion[i].v[d]*ion[i].v[d]);
	for(int d=0;d<dimension;d++)
		for(int i=0;i<mx;i++)
			Ee += 0.5*grid.E[i][d]*grid.E[i][d]*dx*(double (N)/Lx);
	//only magnetic  
	for(int d=1;d<dimension;d++)
	  for(int i=0;i<mx;i++)
		  Em += 0.5*grid.B[i][d]*grid.B[i][d]*dx*(double (N)/Lx)*alpha;
	//	  Ee += 0.5*grid.rho[i]*phi(i)*dx
	//  output energy 
	file_e << setprecision(8) << Ek << "  "<< Ee << "  " << Em << "  " << t << endl;
	ofstream file_rv("rv.txt",ios::app);
	//  output (N/2)th paricle's position and velocity
	file_rv << setprecision(8) << electron[N/2].x << "  " << electron[N/2].v[0] << "  " << electron[N/2].v[1] << "  " << electron[N/2].v[2] << "  " << ion[N/2].x << "  " << ion[N/2].v[0] << "  " << ion[N/2].v[1] << "  " << ion[N/2].v[2] << "  " << t <<endl;
}





void dispersion()
{
//dispersion relation
	static int jkt =1; 
	if(it%1 == 0 && jkt <=1024)
	{
		ofstream file_Ekt("E_k_t.txt",ios::app);
		ofstream file_Eyxt("Ey_x_t.txt",ios::app);
		ofstream file_Exxt("Ex_x_t.txt",ios::app);
		ofstream file_Ezxt("Ez_x_t.txt",ios::app);
		for(int i=1;i<=mx/2;i++)
			file_Ekt << setprecision(8) << E_k[i].real() << "  " ;
		for(int i=0;i<mx;i++)
		{
			file_Eyxt << setprecision(8) << grid.E[i][1] << "  " ;
			file_Exxt << setprecision(8) << grid.E[i][0] << "  " ;
			file_Ezxt << setprecision(8) << grid.E[i][2] << "  " ;
		}
		file_Ekt << endl;
		file_Eyxt << endl;
		file_Exxt << endl;
		file_Ezxt << endl;
     jkt++;
	}
}
void record_data(int it)
{
	string str = "particle.txt_"+num2str(it)+".txt";
	ofstream file_p(str.c_str(),ios::app);
	str = "grid_info_"+num2str(it)+".txt";
	ofstream file_g(str.c_str(),ios::app);
	file_p << "x    " << "vx    " << "vy    " << "vz    " << "time    " <<endl;
	for(int i=0;i<N;i++)
	{
		file_p << electron[i].x << "  ";
		for(int d=0;d<dimension;d++)
			file_p << electron[i].v[d] << "  ";
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

}
	
int main()
{
	clean_file();
	initial(2);
	particle2grid();
	poisson_eqn(1);
	grid2particle();
	pusher_0();
	trans_Bfield_0();
	while(t < 110)
	{
		grid2particle();
		trans_Bfield();
		pusher();
		current();
		particle2grid();
		trans_Efield();
		poisson_eqn(1);
		if( it % 1000 == 0 )
		{
			cout << t << endl;
			record_data(int(it/1000));
		}
		diagnosis();
		dispersion();
		it++;
	}
	system("pause");
	return 0;
}


