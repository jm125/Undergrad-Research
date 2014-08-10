#include <nvector/nvector_serial.h>
#include <sundials/sundials_types.h>
#include <cvode/cvode.h>
#include <sundials/sundials_math.h>
#include "force.h"
#include <fstream>
#include <vector>

void printPoints(realtype[], int);
void PopulateLists(std::vector<double>&, std::vector<int>&, std::vector<int>&,int&,bool&);
int main()
{
	
	void *cvode_mem = CVodeCreate(CV_ADAMS, CV_FUNCTIONAL);

	int flag;

	//initialize
	int numPoints;
	std::vector<int> extensibleSpringList;
	std::vector<int> torsionalSpringList;
	std::vector<double> init;
	bool returnFlag;
	PopulateLists(init, extensibleSpringList, torsionalSpringList, numPoints, returnFlag);
	if (returnFlag) {
		return(-1);
	}
	N_Vector u = N_VNew_Serial(numPoints*3);
	for (int i = 0; i < init.size(); ++i)
	{
		NV_Ith_S(u,i) = init.at(i);
	}
	
	realtype reltol = RCONST(1e-9);
	realtype abstol = RCONST(1e-9);

	flag = CVodeInit(cvode_mem, force, RCONST(0), u);
	flag = CVodeSStolerances(cvode_mem, reltol, abstol);

	hexconfig testhex(u, numPoints, extensibleSpringList, torsionalSpringList);

	flag = CVodeSetUserData(cvode_mem, &testhex);

	realtype tout = 100;
	realtype t;


	realtype radius = testhex.r;
	generate_nhbd(u, testhex.nhbd, testhex.nhbd_partner, radius, numPoints);

	std::ofstream outFile;
	outFile.open("output.xyz");
	while (t < tout) {
		flag = CVode(cvode_mem, tout, u, &t, CV_ONE_STEP);
		int curtime;
		curtime = static_cast<int> (t);

		//TBD: change to verlet list
		if (curtime % 10 == 0) {
			outFile << numPoints << "\n"
			<< "Frame: t =  " << t << "\n";
			for (int i = 0; i<numPoints; i++) {
				outFile << "C\t" << NV_Ith_S(u,i) << "\t"
				<< NV_Ith_S(u,i+numPoints) << "\t" << NV_Ith_S(u,i+2*numPoints)
				<< "\t\n";
			}
			generate_nhbd(u, testhex.nhbd, testhex.nhbd_partner, radius, numPoints);
			realtype check[numPoints*3];
			for (int i = 0; i < numPoints*3; ++i)
			{
				check[i] = NV_Ith_S(u, i);
				
			}
			std::cout << "\nChecking at point " << t << ":";
			printPoints(check,numPoints);
		}
	}

	outFile.close();

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
void PopulateLists(std::vector<double>& init, 
	std::vector<int>& extensibleSprings, 
	std::vector<int>& torsionalSprings,
	int& numPoints,
	bool& returnFlag)
{
	std::ifstream inFile;
	inFile.open("init.txt");

	double line;

	numPoints = 0;

	returnFlag = false;
	const double ESPRING_INDEX = 1011;
	const double TSPRING_INDEX = 1012;


	if(!inFile) 
	{
		std::cout << "Failed to open input file\n";
		returnFlag = true;
	}
	else
	{
		int listID = 0;
		int initIndex = 0;
		int curPointInt;
		double curPoint;
		while(inFile >> line && !returnFlag)
		{
			if (line == ESPRING_INDEX)
			{
				listID = 1;
			}
			else if (line == TSPRING_INDEX)
			{
				listID = 2;
			}
			else {
				switch(listID) {
					case 0:
						curPoint = line;
						init.push_back(curPoint);
						initIndex++;
						numPoints++;
						break;
					case 1:
						curPointInt = int(line);
						extensibleSprings.push_back(curPointInt);
						break;
					case 2:
						curPointInt = int(line);
						torsionalSprings.push_back(curPointInt);
						break;
					default:
						std::cout << "Error in input text";
						returnFlag = true;
						break;
				}
			}
			
		}
	}
	inFile.close();	
	numPoints = numPoints / 3;
}