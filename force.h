#ifndef FORCE_H
#define FORCE_H

#include <nvector/nvector_serial.h>
#include <sundials/sundials_types.h>

#include <iostream>
#include <cmath>
#include <cstdio>
#include <cassert>
#include <vector>



extern "C" int force(realtype t, N_Vector u, N_Vector udot, void *user_data);
void generate_nhbd( N_Vector& u,
                    std::vector<int>& nhbd,
                    std::vector<int>& nhbd_partner,
                    realtype r,
                    int numPoints);

//Structure for configuration of hexagons. 
struct hexconfig
{
	//Spring constants, length between bonds, and indices to help with looping through lists

	realtype kb, ks, l;
	int forceIndexE, forceIndexT;
	std::vector<int> espringlist;
	int * tspringlist;
	int espringsize, tspringsize; 

	//radius for generating neighbor list
	realtype r;

	//Lennard-Jones constants
	realtype delta, epsilon;

	int numPoints;

	std::vector<int> nhbd;
	std::vector<int> nhbd_partner;

	hexconfig(N_Vector& init) 
	:espringlist({ 0,1,0,1,2,0,2,3,0,3,4,
			       0,4,5,0,5,0,0,0,6,1,1,6,1,2,6,1,3,6,1,4,6,1,5,
				       6,1,5,10,0,10,9,0,8,7,0,7,4,0,5,11,1,10,11,1,
				       9,11,1,8,11,1,7,11,1,4,11,1})
	 {
	 	r = 3;
	 	numPoints = 12;
		kb = 0;
		ks = 10;
		l = 1;
		delta = 1;
		epsilon = 1;

		forceIndexE = 3;
		forceIndexT = 5;
		//regular springs: list two points and bond type
		//bonds: 0 for regular, 1 for "fake" center point
		espringsize = 66;

		//torsional springs: list four points and bond type
		//bonds: 0 for regular, 1 for overlapping
		/*tspringlist = (int[60]){ 1,6,0,5,1,0,1,6,2,1,1,2,6,3,1,2,3,6,4,
						1,3,4,6,5,0,4,5,6,0,1,6,4,5,11,0,4,11,
						5,10,1,5,11,10,9,1,10,11,9,8,1,9,11,
						8,7,1,8,11,7,4,1 };*/
		//stspringsize = sizeof(tspringlist) / sizeof(int);
	}
};

#endif