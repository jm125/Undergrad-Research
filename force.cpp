#include <nvector/nvector_serial.h>
#include <sundials/sundials_types.h>
#include <cvode/cvode.h>
#include <sundials/sundials_math.h>

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
	 
	
	//INITIALIZE UDOT
	for (int i = 0; i < size; ++i)
	{
		NV_Ith_S(udot,i) = 0;
	}

	realtype gradient[12];
	for (int i = 0; i < test.espringsize; ++i)
	{
		if (i % test.forceIndexE == 1) {
			int x1index = test.espringlist[i-2], x2index = test.espringlist[i-1];
			int y1index = test.espringlist[i-2]+numPoints, y2index = test.espringlist[i-1]+numPoints;
			int z1index = test.espringlist[i-2]+2*(numPoints), z2index = test.espringlist[i-1]+2*(numPoints);

			realtype x1 = NV_Ith_S(u,x1index);
			realtype x2 = NV_Ith_S(u,x2index);

			realtype y1 = NV_Ith_S(u,y1index);
			realtype y2 = NV_Ith_S(u,y2index);

			realtype z1 = NV_Ith_S(u,z1index);
			realtype z2 = NV_Ith_S(u,z2index);

			//pass in array with point coordinates
			realtype pointslist[] = {x1,x2,y1,y2,z1,z2};
			//re-initialize gradient array to 0
			for (int i = 0; i < 6; ++i)
			{
				gradient[i] = 0;
			}
			gradient[0] = calcRegularSpring(pointslist,test.ks,test.l);
			gradient[1] = calcRegularSpring(pointslist,test.ks,test.l)*-1;



			NV_Ith_S(udot, x1index) += gradient[0];
			NV_Ith_S(udot, x2index) += gradient[1];
			NV_Ith_S(udot, y1index) += gradient[2];
			NV_Ith_S(udot, y2index) += gradient[3];
			NV_Ith_S(udot, z1index) = 0;
			NV_Ith_S(udot, z2index) = 0;
		} 
	}

	//similar for torsional

	//to be done: van der Waals


}


realtype calcRegularSpring(realtype pointslist[], realtype ks, realtype l)
{
	realtype x1 = pointslist[0];
	realtype x2 = pointslist[1];
	realtype y1 = pointslist[2];
	realtype y2 = pointslist[3];
	realtype z1 = pointslist[4];
	realtype z2 = pointslist[5];
	realtype dist = sqrt( pow(x1 - x2, 2) + pow(y1 - y2, 2) + pow(z1 - z2, 2));
	realtype xgrad = ks*(dist-l)*(x1-x2)/dist;
	realtype ygrad = ks*(dist-l)*(y1-y2)/dist;
	realtype zgrad = ks*(dist-l)*(z1-z2)/dist;

	return xgrad;

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

