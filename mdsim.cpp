#include<fstream>
#include <cassert>
#include<iostream>
#include <stdlib.h>
#include <string>

double ro;

#include "READ.h"
#include "MD.h"
#include "Timer.h"


using namespace std;




void run_test(ParameterReader readPar){
	string s="-";
	double d;
	readPar.GetParameter("vis_space",s);
	d = stod(s);
	readPar.GetParameter("name",s);
	cerr<<d<<"\n";
	cerr<<s<<"\n";

}

int main(int argc, char *argv[]){
	(void) argc; //to suppress Warnings about unused argc
	assert(argc>0);
	string parData = argv[1];
	string dataFile = argv[2];

	ParameterReader readPar(parData);
	//run_test(readPar);
	//return 0;

	Simulator simulate(readPar);
	simulate.read(dataFile);
	simulate.cellSimulation();

}
