#include<complex>
#include<iostream>
using namespace std;


class Particle{
public:int d;
	   double x;
	   double *v,*E,*B;
	   Particle(int dimension):x(0),d(dimension)
	   {
		   v = new double[d];
		   E = new double[d];
		   B = new double[d];		   
		   for(int i=0;i<d;i++)
			   B[i] = E[i] = v[i] = 0.0;
	   }
	   Particle():x(0),d(3)
	   {
		   v = new double[d];
		   E = new double[d];
		   B = new double[d];
		   for(int i=0;i<d;i++)
			    B[i] = E[i] = v[i] = 0.0;
	   }
	   ~Particle()
	   {
		   delete [] v;
		   delete [] B;
		   delete [] E;
	   }
};

class Species{
public:double mass,charge;
	   string label;
	   int N;
	   Species(string name,int number):label(name),N(number)
	   {
		   if(label == "electron")
		   {
			   mass = 1.0;
			   charge = -1.0;
		   }
		   else if(label == "ion")
		   {
			   mass = mi_div_me;
			   charge = 1.0;
			}
	   }
	   Species(string name,double m,double q,int number):label(name),N(number)
	   {
		   if(label == "electron")
		   {
			   mass = 1.0;
			   charge = -1.0;
		   }
		   else if(label == "ion")
		   {
			   mass = mi_div_me;
			   charge = 1.0;
		   }
		   else
		   {
			   mass = m;
			   charge = 1;
		   }
	   }
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
		delete [] rho;
		delete [] E;
		delete [] B;
		delete [] J;
		delete [] v_fluid;
	   }
};










extern const int perstep;
extern const int nx;
extern const double dt;
extern double t;
extern int it;


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

double erfccheb(double z);
inline double erf(double x);
void fft(int isign,const complex<double> *datain,complex<double> *dataout,const int N);
string num2str(int num);
void Gauss_Seidel1D(double* unew,double* f,const double err,const double dx,const int mx);
void Gauss_Seidel1D_Darwin(double* unew,const double* f,const double err,const double dx,const int mx,const double omega_0,const double* omega);





void clean_file();
void calc_maxw(double vth,double* v,const int N);
void distrubution(const double *v,const int N,string name);
void distrubution_full(const int N,int it);
void initial(const int condition);
void particle2grid();
void current();
void grid2particle();
void poisson_eqn(const int method);
void poisson_eqn_GS();
void pusher_x();
void pusher_v();
void pusher_v_0();
void diagnosis();
void dispersion();
void record_data(int it);
void Darwin_full();
void wave_heating();
void parameter_output();