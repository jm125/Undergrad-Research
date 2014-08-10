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
	//Spring constants, length between bonds
	realtype kb, ks, l;
	int espringsize, tspringsize; 

	//radius for generating neighbor list
	realtype r;

	//Lennard-Jones constants
	realtype delta, epsilon;

	std::vector<int> nhbd;
	std::vector<int> nhbd_partner;
	int forceIndexE, forceIndexT;

	int numPoints;
	std::vector<int> tspringlist;
	std::vector<int> espringlist;
	realtype a;

	hexconfig(N_Vector& init, int npoints, std::vector<int> elist, std::vector<int> tlist)
	 {
	 	//Test data - modify these parameters for different simulations.
	 	r = 3;
		kb = 0;
		ks = 10;
		l = 1;
		delta = 1;
		epsilon = 1;
		a = 0;

		numPoints = npoints;
		for (int i = 0; i < elist.size(); ++i)
		{
			espringlist.push_back(elist.at(i));
		}
		for (int i = 0; i < tlist.size(); ++i)
		{
			tspringlist.push_back(tlist.at(i));
		}
		forceIndexE = 3;
		forceIndexT = 5;
		espringsize = espringlist.size();
		tspringsize = tspringlist.size();
	}
};

#endif