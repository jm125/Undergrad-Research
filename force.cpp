#include <nvector/nvector_serial.h>
#include <sundials/sundials_types.h>
#include <cvode/cvode.h>
#include <sundials/sundials_math.h>
#include <cstdio>
#include <iostream>

#include "force.h"

void generate_nhbd( N_Vector& u,
                    std::vector<int>& nhbd,
                    std::vector<int>& nhbd_partner,
                    realtype r,
                    int numPoints);

extern "C" int force(realtype t, N_Vector u, N_Vector udot, void *user_data) 
{
	hexconfig& test = *((hexconfig*)(user_data));

	

	//Number of points for our configuration
	int numPoints = 12, size = numPoints*3;


	//NOTE: assigment of spring points should be in struct; this is a temporary assignment for testing purposes
	int espringlist[] = { 0,1,0,1,2,0,2,3,0,3,4,
			       	       0,4,5,0,5,0,0,0,6,1,1,6,1,2,6,1,3,6,1,4,6,1,5,
				       6,1,5,10,0,10,9,0,8,7,0,7,4,0,5,11,1,10,11,1,
				       9,11,1,8,11,1,7,11,1,4,11,1};

	int tspringlist[] = { 1,6,0,5,1,0,1,6,2,1,1,2,6,3,1,2,3,6,4,
						1,3,4,6,5,0,4,5,6,0,1,6,4,5,11,0,4,11,
						5,10,1,5,11,10,9,1,10,11,9,8,1,9,11,
						8,7,1,8,11,7,4,1 };
	 
	
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
			int xiindex = tspringlist[i-4];
			int xjindex = tspringlist[i-3];
			int xkindex = tspringlist[i-2];
			int xmindex = tspringlist[i-1];

			int yiindex = xiindex+numPoints;
			int yjindex = xjindex+numPoints;
			int ykindex = xkindex+numPoints;
			int ymindex = xmindex+numPoints;

			int ziindex = xiindex+numPoints*2;
			int zjindex = xjindex+numPoints*2;
			int zkindex = xkindex+numPoints*2;
			int zmindex = xmindex+numPoints*2;

			realtype xi = NV_Ith_S(u,xiindex);
			realtype xj = NV_Ith_S(u,xjindex);
			realtype xk = NV_Ith_S(u,xkindex);
			realtype xm = NV_Ith_S(u,xmindex);

			realtype yi = NV_Ith_S(u,yiindex);
			realtype yj = NV_Ith_S(u,yjindex);
			realtype yk = NV_Ith_S(u,ykindex);
			realtype ym = NV_Ith_S(u,ymindex);

			realtype zi = NV_Ith_S(u,ziindex);
			realtype zj = NV_Ith_S(u,zjindex);
			realtype zk = NV_Ith_S(u,zkindex);
			realtype zm = NV_Ith_S(u,zmindex);

			realtype normgrad = (xi + xj - 2*xh)*(yi*yi+zi*zj)
		+ (xi - xj)*((yj*yh-yi*yj)+(zj*zh-zi*zh))
		+ (x1 - x3)*(y2*y2 + z2*z2)
		+ (x1 - x2)*(y3*y3 + z3*z3);
		}
	}*/

	//van der Waals
	//generate neighbor list every 10th step
	realtype radius = test.r;
	int curtime;
	curtime = static_cast<int> (t);
	if (curtime % 10 == 0){
		generate_nhbd(u, test.nhbd, test.nhbd_partner, radius, numPoints);
	}

	//compute bending
	int isize = test.nhbd.size();
	realtype eps = test.epsilon;
	realtype del = test.delta;

	for (int i = 0; i < isize; ++i) {
		int x1index = test.nhbd.at(i);
		int x2index = test.nhbd_partner.at(i);
		int y1index = x1index+numPoints;
		int y2index = x2index+numPoints;
		int z1index = x1index+2*(numPoints);
		int z2index = x2index+2*(numPoints);

		realtype x = NV_Ith_S(u,x1index) - NV_Ith_S(u,x2index);
		realtype y = NV_Ith_S(u,y1index) - NV_Ith_S(u,y2index);
		realtype z = NV_Ith_S(u,z1index) - NV_Ith_S(u,z2index);

		realtype LJdist = sqrt( x*x + y*y + z*z );
		realtype LJconst = RCONST(-12)*eps/del;
		realtype LJprime = pow(del/LJdist, 13) - pow(del/LJdist, 7);
		realtype LJx = LJconst*LJprime*x/LJdist;
		realtype LJy = LJconst*LJprime*y/LJdist;
		realtype LJz = LJconst*LJprime*z/LJdist;

		NV_Ith_S(udot, x1index) += LJx;
		NV_Ith_S(udot, x2index) -= LJx;
		NV_Ith_S(udot, y1index) += LJy;
		NV_Ith_S(udot, y2index) -= LJy;
		NV_Ith_S(udot, z1index) += LJz;
		NV_Ith_S(udot, z2index) -= LJz;

	}

}

void generate_nhbd( N_Vector& u,
                    std::vector<int>& nhbd,
                    std::vector<int>& nhbd_partner,
                    realtype r,
                    int numPoints)
{
	nhbd.clear();
	nhbd_partner.clear();

	realtype r2 = r*r;

	for (int i = 0; i < numPoints - 1; ++i)
	{
		realtype x1 = NV_Ith_S(u, i);
		realtype y1 = NV_Ith_S(u, i+numPoints);
		realtype z1 = NV_Ith_S(u, i+2*numPoints);


		for (int j = i+1; j < numPoints; ++j)
		{
			realtype x2 = NV_Ith_S(u, j);
			realtype y2 = NV_Ith_S(u, j+numPoints);
			realtype z2 = NV_Ith_S(u, j+2*numPoints);

			//calculate distance between points
			realtype dist2 = (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) + (z1-z2)*(z1-z2);
			if (dist2 <= r2) {
				nhbd.push_back(i);
				nhbd_partner.push_back(j);
			}
		}
	}

}
