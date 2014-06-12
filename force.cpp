#include <nvector/nvector_serial.h>
#include <sundials/sundials_types.h>
#include <cvode/cvode.h>
#include <sundials/sundials_math.h>
#include <cstdio>
#include <iostream>

#include "force.h"

//define function protoypes
//note: for force calculations, calculate forces between x,y,z coordinates seperately
realtype calcRegularSpring(realtype[],realtype,realtype);
realtype calcTorsionalSpring(realtype,realtype[]);

extern "C" int force(realtype t, N_Vector u, N_Vector udot, void *user_data) 
{
	hexconfig& test = *((hexconfig*)(user_data));

	

	//Number of points for our configuration
	int numPoints = 12, size = numPoints*3;


	//NOTE: assigment of spring points should be in struct; this is a temporary assignment for testing purposes
	/*int espringlist[] = { 0,1,0,1,2,0,2,3,0,3,4,
			       	       0,4,5,0,5,0,0,0,6,1,1,6,1,2,6,1,3,6,1,4,6,1,5,
				       6,1,5,10,0,10,9,0,8,7,0,7,4,0,5,11,1,10,11,1,
				       9,11,1,8,11,1,7,11,1,4,11,1};

	int tspringlist[] = { 1,6,0,5,1,0,1,6,2,1,1,2,6,3,1,2,3,6,4,
						1,3,4,6,5,0,4,5,6,0,1,6,4,5,11,0,4,11,
						5,10,1,5,11,10,9,1,10,11,9,8,1,9,11,
						8,7,1,8,11,7,4,1 };*/
	 
	
	//INITIALIZE UDOT
	for (int i = 0; i < size; ++i)
	{
		NV_Ith_S(udot,i) = 0;
	}

	int count = 0;
	//calculate extensible springs
	for (int i = 0; i < test.espringsize; ++i)
	{
		if (i % test.forceIndexE == 2) {
			int x1index = test.espringlist[i-2];
			int x2index = test.espringlist[i-1];
			int y1index = x1index+numPoints;
			int y2index = x2index+numPoints;
			int z1index = x1index+2*(numPoints);
			int z2index = x2index+2*(numPoints);

			std::cout << x1index << " " << y1index;

			realtype x1 = NV_Ith_S(u,x1index);
			realtype x2 = NV_Ith_S(u,x2index);

			realtype y1 = NV_Ith_S(u,y1index);
			realtype y2 = NV_Ith_S(u,y2index);

			realtype z1 = NV_Ith_S(u,z1index);
			realtype z2 = NV_Ith_S(u,z2index);

			//compute gradients
			realtype dist = sqrt( pow(x1 - x2, 2) + pow(y1 - y2, 2) + pow(z1 - z2, 2));
			realtype xgrad = test.ks*(dist-test.l)*(x1-x2)/(dist);
			realtype ygrad = test.ks*(dist-test.l)*(y1-y2)/(dist);
			realtype zgrad = test.ks*(dist-test.l)*(z1-z2)/(dist);

			NV_Ith_S(udot, x1index) += xgrad;
			NV_Ith_S(udot, x2index) -= xgrad;
			NV_Ith_S(udot, y1index) += ygrad;
			NV_Ith_S(udot, y2index) -= ygrad;
			NV_Ith_S(udot, z1index) += zgrad;
			NV_Ith_S(udot, z2index) -= zgrad;

		} 
	}

	//torsional spring computation
	/*for (int i = 0; i < test.tspringsize; ++i)
	{
		if (i % test.forceIndexT == 4) {
			//get x,y,z coords of points i, j, k, m
			int xiindex = espringlist[i-4];
			int xjindex = espringlist[i-3];
			int xkindex = espringlist[i-2];
			int xmindex = espringlist[i-1];

			int yiindex = xiindex+numPoints;
			int yjindex = xjindex+numPoints;
			int ykindex = xkindex+numPoints;
			int ymindex = xmindex+numPoints;

			int ziindex = xiindex+numPoints*2;
			int zjindex = xjindex+numPoints*2;
			int zkindex = xkindex+numPoints*2;
			int zmindex = xmindex+numPoints*2;

			realtype x1 = NV_Ith_S(u,x1index);
			realtype x2 = NV_Ith_S(u,x2index);

			realtype y1 = NV_Ith_S(u,y1index);
			realtype y2 = NV_Ith_S(u,y2index);

			realtype z1 = NV_Ith_S(u,z1index);
			realtype z2 = NV_Ith_S(u,z2index);
		}
	}*/

	//to be done: van der Waals


}


//returns gradient for x coordinate; to get y or z gradient, just input in different order
realtype calcTorsionalSpring(realtype kb, realtype pointslist[])
{
	realtype x1 = pointslist[0];
	realtype x2 = pointslist[1];
	realtype x3 = pointslist[2];
	realtype x4 = pointslist[3];
	realtype y1 = pointslist[4];
	realtype y2 = pointslist[5];
	realtype y3 = pointslist[6];
	realtype y4 = pointslist[7];
	realtype z1 = pointslist[8];
	realtype z2 = pointslist[9];
	realtype z3 = pointslist[10];
	realtype z4 = pointslist[11];
	realtype normgrad = (x3 + x2 - 2*x4)*(y3*y2+z3*z2)
		+ (x3 - x2)*((y1*y2-y3*y2)+(z1*z2-z3*z1))
		+ (x1 - x3)*(y2*y2 + z2*z2)
		+ (x1 - x2)*(y3*y3 + z3*z3);
	realtype gradient = -kb * normgrad;

	return gradient;
}

