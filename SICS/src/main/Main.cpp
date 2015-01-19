//============================================================================
// Name        : SICS.cpp
// Author      : Juan Liberato
// Version     :
// Copyright   : Undefined
// Description : Hello World in C++, Ansi-style
//============================================================================

#include "Main.h"

void profilerOut(Trace* profile, int type){
	//Types of profiling :
	//	1: Only total time
	//  2: General
	//  3: Verbose
	//  4: Verbose with percents
	if(type==1){
		(*profile)("Profiling file : ",'n');
		(*profile)(profile->filename);
		(*profile)("Total time :",'n');
		(*profile)(profile->dr("total"));
	}
	else if(type==2){
		profilerOut(profile,1);
		(*profile)("Input time :",'n');
		(*profile)(profile->dr("input"));
		(*profile)("Initial time :",'n');
		(*profile)(profile->dr("initial"));
	}
	else if(type==3){
		profilerOut(profile,2);
		(*profile)("Estimation time :",'n');
		(*profile)(profile->dr("estim"));
		//Estep time
		(*profile)("E step time:",'n');
		(*profile)(profile->dr("Et"));
		//Mstep time
		(*profile)("M step time:",'n');
		(*profile)(profile->dr("Mt"));
	}
	else if(type == 4){
		(*profile)("Profiling file : ",'n');
		(*profile)(profile->filename);
		(*profile)("Total time :",'n');
		float total = (float)profile->dr("total");
		float percent = total/100;
		(*profile)((float)profile->dr("total")/percent);
		(*profile)("Input time :",'n');
		(*profile)((float)profile->dr("input")/percent);
		(*profile)("Initial time :",'n');
		(*profile)((float)profile->dr("initial")/percent);
		(*profile)("Estimation time :",'n');
		(*profile)((float)profile->dr("estim")/percent);
		//Estep time
		(*profile)("E step time:",'n');
		(*profile)((float)profile->dr("Et")/percent);
		//Mstep time
		(*profile)("M step time:",'n');
		(*profile)((float)profile->dr("Mt")/percent);
	}
	else{
		(*profile)("No profiling selected please select a profiling mode");
	}
}

void oneRun(char * args) {
	Input input;
	Matrix<double> cuad(41, 2);
	//Create the profiler to profile the program
	Trace* profiler = new Trace("Profile.log");
	profiler->resetTimer("total");
	profiler->startTimer("total");
	profiler->resetTimer("input");
	profiler->startTimer("input");
	input.importCSV((char *) "Cuads.csv", cuad, 1, 0);
	// **** **** Run model complete and ordered process **** ****
	// Create general pars
	Model *model = new Model();
	// Create general model
	ModelFactory *modelFactory = new SICSGeneralModel();
	PatternMatrix *dataSet = new PatternMatrix(0);
	// Load matrix
	input.importCSV(args, *dataSet, 1, 0);
	// set dataset
	dataSet->print();
	//RASCH_A1, RASCH_A_CONSTANT, TWO_PL, THREE_PL
	model->setModel(modelFactory, Constant::TWO_PL);
	//This is where it is decided what model is the test to make
	model->getItemModel()->setDataset(dataSet);		//Sets the dataset.
	// set Theta and weight for the EM Estimation
	Matrix<double> *theta = new Matrix<double>(1, 41);
	Matrix<double> *weight = new Matrix<double>(1, 41);

	for (int k = 0; k < cuad.nR(); k++) {
		(*theta)(0, k) = cuad(k, 0);
		(*weight)(0, k) = cuad(k, 1);
	}


	// build parameter set
	model->getParameterModel()->buildParameterSet(model->getItemModel(),
			model->getDimensionModel());

	// Create estimation
	profiler->stopTimer("input");
	profiler->resetTimer("initial");
	profiler->startTimer("initial");
	EMEstimation *em = new EMEstimation();
	//Here is where quadratures must be set.
	//create the quad nodes
	QuadratureNodes nodes(theta, weight);
	em->setQuadratureNodes(&nodes);
	em->setModel(model);
	cout<<"Setting initial values with option "<<Constant::ANDRADE<<endl;
	cout<<em<<endl;
	em->setInitialValues(Constant::ANDRADE);
	profiler->stopTimer("initial");
	//Pass the profiler to the estimation object so it can be used to profile each step
	em->setProfiler(profiler);
	//Run the estimation
	em->estimate();
	delete modelFactory;
	delete dataSet;
	delete em;
	delete model;
	profiler->stopTimer("total");
	//Out the profiler here
	profilerOut(profiler,4);
	delete profiler;
}
int main(int argc, char *argv[]) {
	double ** ps = new double*[3];
	ps[0] = new double [3];
	ps[1] = new double [2];
	ps[2] = new double [80];
	for (int i = 0 ; i < 3 ; i++){
		ps[0][i]=i;
	}

	Timer tm;
	tm.reset();
	tm.start();
	if (argc < 2) {
		cout << "Please specify an input file" << endl;
		return (0);
	}
	oneRun(argv[1]);
	tm.stop();
	cout<<"time: "<<endl<<tm.totalTime<<endl;
	return (0);
}

