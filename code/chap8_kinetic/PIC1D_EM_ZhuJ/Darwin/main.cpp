#include"parameter.h"
#include"PIC.h"

using namespace std;


int main()
{
	clean_file();
	initial(6);
    parameter_output();
/*
	for(int i=0;i<N;i++)
	{
		double kz=0.003;
		particle[i][0].v[1] += 0.0025*cos(kz*particle[i][0].x);
		particle[i][1].v[1] += -0.00054*0.01*cos(kz*particle[i][0].x);
		particle[i][0].v[2] += -0.027*sin(kz*particle[i][0].x);
		particle[i][1].v[2] += -0.028*sin(kz*particle[i][0].x);
		particle[i][0].v[0] += 0*sin(kz*particle[i][0].x);
		particle[i][1].v[0] += 0.012*sin(kz*particle[i][0].x);
		particle[i][0].x += 0.069*sin(kz*particle[i][0].x);
		particle[i][1].x += -0.071*sin(kz*particle[i][0].x);
	}
*/
	record_data(0);
	particle2grid();
//	poisson_eqn(2);
	poisson_eqn_GS();
	grid2particle();
	pusher_v_0();
	current();
	Darwin_full();
	while(t < totaltime)
	{
//		int timestep = int (totaltime/dt);
//        perstep  = int(timestep/nx);
        

		pusher_v();
		pusher_x();
		particle2grid();
//		poisson_eqn(2);
		poisson_eqn_GS();
		current();//wave_heating();
		Darwin_full();
		grid2particle();
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
