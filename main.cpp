#include <nvector/nvector_serial.h>
#include <sundials/sundials_types.h>
#include <cvode/cvode.h>
#include <sundials/sundials_math.h>
#include "force.h"

void printPoints(realtype[], int);
int main()
{
	
	void *cvode_mem = CVodeCreate(CV_ADAMS, CV_FUNCTIONAL);

	int flag;

	//define two hexagons, with 7 points each (including center point) - two points overlap
	//also define vector containing x, y, and z coordinates
		
	//initialize xyz coordinates for all 12 points
	realtype init[] = { 0.5,0,0.5,1.5,2,1.5,1,3,3.5,3,2,2.5,
				  2,1,0,0,1,2,1,1,2,3,3,2,
				  1,0,1,0,1,0,1,0,1,0,1,0 };
	std::cout << "\nInitial list:";
	printPoints(init, 12);

	N_Vector u = N_VMake_Serial(36, init);
	
	//set tolerances
	realtype reltol = RCONST(1e-9);
	realtype abstol = RCONST(1e-9);

	flag = CVodeInit(cvode_mem, force, RCONST(0), u);
	flag = CVodeSStolerances(cvode_mem, reltol, abstol);

	hexconfig testhex(u);

	flag = CVodeSetUserData(cvode_mem, &testhex);

	realtype tout = 1000;
	realtype t = 1;
	flag = CVode(cvode_mem, tout, u, &t, CV_ONE_STEP);

	realtype final[36];
	for (int i = 0; i < 36; ++i)
	{
		final[i] = NV_Ith_S(u, i);
	}

	
	std::cout << "\nFinal list:";
	printPoints(final,12);
	N_VDestroy(u);
	CVodeFree(&cvode_mem);
}

void printPoints(realtype list[], int numPoints)
{
	std::cout << "\n";
	for (int i = 0; i < (numPoints*3); ++i)
	{
		std::cout << (i < numPoints ? "x" : (i < (2*numPoints) ? "y" : "z")) << 
		(i < numPoints ? i : (i < (2*numPoints) ? i - numPoints : i - (2*numPoints))) << ":\t " << list[i] << "\n";
	}
}