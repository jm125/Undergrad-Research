#include <nvector/nvector_serial.h>
#include <sundials/sundials_types.h>
#include <cvode/cvode.h>
#include <sundials/sundials_math.h>
#include "force.h"
#include <fstream>
#include <vector>

void printPoints(realtype[], int);
void PopulateLists(std::vector<double>&, std::vector<int>&, std::vector<int>&,int&);
int main()
{
	
	void *cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);

	int flag;
		
	//initialize
	int numPoints;
	std::vector<int> extensibleSpringList;
	std::vector<int> torsionalSpringList;
	std::vector<double> init;
	PopulateLists(init, extensibleSpringList, torsionalSpringList, numPoints);
	N_Vector u = N_VNew_Serial(numPoints);
	
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
			outFile.close();
			generate_nhbd(u, testhex.nhbd, testhex.nhbd_partner, radius, numPoints);
			realtype check[36];
			for (int i = 0; i < 36; ++i)
			{
				check[i] = NV_Ith_S(u, i);
				
			}
			std::cout << "\nChecking at point " << t << ":";
			printPoints(check,12);
		}
	}

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
	int& numPoints)
{
	std::ifstream inFile;
	inFile.open("init.txt");

	char line[10];
	std::vector<char> fileLines;

	numPoints = 0;

	if(!inFile) 
	{
		std::cout << "Failed to open input file";
	}
	else
	{
		while(!inFile.eof())
		{
			inFile >> line;
			fileLines.push_back(line);
		}
	}
	inFile.close();

	int listID = 0;
	int initIndex = 0;
	int curPointInt;
	double curPoint;
	for (int i = 0; i < fileLines.size(); ++i)
	{
		if (fileLines.at(i) == "e")
		{
			listID = 1;
		}
		else if (fileLines.at(i) == "t")
		{
			listID = 2;
		}
		else {
			switch(listID) {
				case 0:
					curPoint = atof(fileLines.at(i));
					//NV_Ith_S(init, initIndex) = curPoint;
					init.push_back(curPoint);
					initIndex++;
					numPoints++;
					break;
				case 1:
					curPointInt = atoi(fileLines.at(i));
					extensibleSprings.push_back(curPointInt);
					break;
				case 2:
					curPointInt = atoi(fileLines.at(i));
					torsionalSprings.push_back(curPointInt);
					break;
				default:
					std::cout << "Error";
					break;
			}
		}
	}
}