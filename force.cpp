#include <nvector/nvector_serial.h>
#include <sundials/sundials_types.h>
#include <cvode/cvode.h>
#include <sundials/sundials_math.h>
#include <cstdio>
#include <string>
#include <iostream>
#include <fstream>
#include "force.h"

void generate_nhbd( N_Vector& u,
                    std::vector<int>& nhbd,
                    std::vector<int>& nhbd_partner,
                    realtype r,
                    int numPoints);
void writeToXYZFile(std::ofstream outFile,
	int t, N_Vector points, int size);

extern "C" int force(realtype t, N_Vector u, N_Vector udot, void *user_data) 
{
	hexconfig& test = *((hexconfig*)(user_data));
	int size = test.numPoints*3;
	//INITIALIZE UDOT
	for (int i = 0; i < size; i++)
	{
		NV_Ith_S(udot,i) = 0;
	}

	
	int count = 0;
	//calculate extensible springs
	for (int i = 0; i < test.espringsize; ++i)
	{
		if (i % test.forceIndexE == 2) {
			int x1index = test.espringlist.at(i-2);
			int x2index = test.espringlist.at(i-1);
			int y1index = x1index+test.numPoints;
			int y2index = x2index+test.numPoints;
			int z1index = x1index+2*(test.numPoints);
			int z2index = x2index+2*(test.numPoints);

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

			NV_Ith_S(udot, x1index) -= xgrad;
			NV_Ith_S(udot, x2index) += xgrad;
			NV_Ith_S(udot, y1index) -= ygrad;
			NV_Ith_S(udot, y2index) += ygrad;
			NV_Ith_S(udot, z1index) -= zgrad;
			NV_Ith_S(udot, z2index) += zgrad;

		} 
	}

	//torsional spring computation
	//May want to double check that these formulas are correct, until then leave kb = 0
	for (int i = 0; i < test.tspringsize; ++i)
	{
		if (i % test.forceIndexT == 4) {
			//get x,y,z coords of points i, j, k, m
			int xiindex = test.tspringlist[i-4];
			int xjindex = test.tspringlist[i-3];
			int xkindex = test.tspringlist[i-2];
			int xmindex = test.tspringlist[i-1];

			int yiindex = xiindex+test.numPoints;
			int yjindex = xjindex+test.numPoints;
			int ykindex = xkindex+test.numPoints;
			int ymindex = xmindex+test.numPoints;

			int ziindex = xiindex+test.numPoints*2;
			int zjindex = xjindex+test.numPoints*2;
			int zkindex = xkindex+test.numPoints*2;
			int zmindex = xmindex+test.numPoints*2;

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

			realtype xidot = -1*test.kb*(0.75*test.a - (xj-xk)*(xj-xi)*(xm-xj)*(xk-xj)-(xj-xk)*(xk-xj)*(xj-xk)*(xm-xj));
			realtype xjdot = xidot;
			realtype xkdot = -xidot;
			realtype xmdot = -xidot;

			realtype yidot = -1*test.kb*(0.75*test.a - (yj-yk)*(yj-yi)*(ym-yj)*(yk-yj)-(yj-yk)*(yk-yj)*(yj-yk)*(ym-yj));
			realtype yjdot = yidot;
			realtype ykdot = -yidot;
			realtype ymdot = -yidot;

			realtype zidot = -1*test.kb*(0.75*test.a - (zj-zk)*(zj-zi)*(zm-zj)*(zk-zj)-(zj-zk)*(zk-zj)*(zj-zk)*(zm-zj));
			realtype zjdot = zidot;
			realtype zkdot = -zidot;
			realtype zmdot = -zidot;

			NV_Ith_S(udot, xiindex) += xidot;
			NV_Ith_S(udot, xjindex) += xjdot;
			NV_Ith_S(udot, xkindex) += xkdot;
			NV_Ith_S(udot, xmindex) += xmdot;
			NV_Ith_S(udot, yiindex) += yidot;
			NV_Ith_S(udot, yjindex) += yjdot;
			NV_Ith_S(udot, ykindex) += ykdot;
			NV_Ith_S(udot, ymindex) += ymdot;
			NV_Ith_S(udot, ziindex) += zidot;
			NV_Ith_S(udot, zjindex) += zjdot;
			NV_Ith_S(udot, zkindex) += zkdot;
			NV_Ith_S(udot, zmindex) += zmdot;
		}
	}

	//compute vdW
	int isize = test.nhbd.size();
	realtype eps = test.epsilon;
	realtype del = test.delta;

	for (int i = 0; i < isize; ++i) {
		int x1index = test.nhbd.at(i);
		int x2index = test.nhbd_partner.at(i);
		int y1index = x1index+test.numPoints;
		int y2index = x2index+test.numPoints;
		int z1index = x1index+2*(test.numPoints);
		int z2index = x2index+2*(test.numPoints);

		realtype x = NV_Ith_S(u,x1index) - NV_Ith_S(u,x2index);
		realtype y = NV_Ith_S(u,y1index) - NV_Ith_S(u,y2index);
		realtype z = NV_Ith_S(u,z1index) - NV_Ith_S(u,z2index);

		realtype LJdist = sqrt( x*x + y*y + z*z );
		realtype LJconst = RCONST(-12)*eps/del;
		realtype LJprime = pow(del/LJdist, 13) - pow(del/LJdist, 7);
		realtype LJx = LJconst*LJprime*x/LJdist;
		realtype LJy = LJconst*LJprime*y/LJdist;
		realtype LJz = LJconst*LJprime*z/LJdist;

		NV_Ith_S(udot, x1index) -= LJx;
		NV_Ith_S(udot, x2index) += LJx;
		NV_Ith_S(udot, y1index) -= LJy;
		NV_Ith_S(udot, y2index) += LJy;
		NV_Ith_S(udot, z1index) -= LJz;
		NV_Ith_S(udot, z2index) += LJz;

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

	for (int i = 0; i < numPoints; ++i)
	{
		realtype x1 = NV_Ith_S(u, i);
		realtype y1 = NV_Ith_S(u, i+numPoints);
		realtype z1 = NV_Ith_S(u, i+2*numPoints);


		for (int j = i+1; j < numPoints; ++j)
		{
			realtype x2 = NV_Ith_S(u, j);
			realtype y2 = NV_Ith_S(u, j+numPoints);
			realtype z2 = NV_Ith_S(u, j+2*numPoints);

			realtype distSquared = (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) + (z1-z2)*(z1-z2);
			if (distSquared <= r2) {
				nhbd.push_back(i);
				nhbd_partner.push_back(j);
			}
		}
	}

}
void writeToXYZFile(std::ofstream outFile,
	int t, N_Vector u, int size)
{
	outFile.open("output.xyz");


	outFile << size << "\n"
			<< "Frame: t =  " << t << "\n";
	for (int i = 0; i<(size-3); i+=3) {
		outFile << "C  " << NV_Ith_S(u, i) << "  "
			    << NV_Ith_S(u, i+1) << "  " << NV_Ith_S(u, i+2)
			    << "\n";
	}
	outFile.close();
}