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
		(*profile)("_______________________________");
		(*profile)("Profiling file : ",'n');
		(*profile)(profile->filename);
		(*profile)("Total time :",'n');
		float total = (float)profile->dr("total");
		(*profile)((float)profile->dr("total"));
		(*profile)("Input time :",'n');
		(*profile)((float)profile->dr("input"));
		(*profile)("Initial time :",'n');
		(*profile)((float)profile->dr("initial"));
		(*profile)("Estimation time :",'n');
		(*profile)((float)profile->dr("estim"));
		//Estep time
		(*profile)("E step time:",'n');
		(*profile)((float)profile->dr("Et"));
		(*profile)("First for time:",'n');
		(*profile)((float)profile->dr("for1"));
		(*profile)("Second for time:",'n');
		(*profile)((float)profile->dr("for2"));
		//Mstep time
		(*profile)("M step time:",'n');
		(*profile)((float)profile->dr("Mt"));

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
	profiler->resetTimer("for1");
	profiler->resetTimer("for2");
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
	//RASCH_A1, RASCH_A_CONSTANT, TWO_PL, THREE_PL
	model->setModel(modelFactory, Constant::THREE_PL);
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

void runArgs(char * filename,char * initialValues){
	Input input;
	Matrix<double> cuad(41, 2);
	Trace* profiler = new Trace("Profile.log");
	profiler->resetTimer("total");
	profiler->startTimer("total");
	profiler->resetTimer("input");
	profiler->startTimer("input");
	profiler->resetTimer("for1");
	profiler->resetTimer("for2");
	input.importCSV((char *) "Cuads.csv", cuad, 1, 0);
	// **** **** Run model complete and ordered process **** ****
	// Create general pars

	Model *model = new Model();
	// Create general model
	ModelFactory *modelFactory = new SICSGeneralModel();
	PatternMatrix *dataSet = new PatternMatrix(0);
	// Load matrix
	input.importCSV(filename, *dataSet, 1, 0);
	// set dataset
	//RASCH_A1, RASCH_A_CONSTANT, TWO_PL, THREE_PL
	int model_const =  Constant::THREE_PL;
	model->setModel(modelFactory, model_const);
	//This is where it is decided what model is the test to make
	model->getItemModel()->setDataset(dataSet);//Sets the dataset.
	// set Theta and weight for the EM Estimation
	Matrix<double> *theta = new Matrix<double>(1, 41);
	Matrix<double> *weight = new Matrix<double>(1, 41);

	for (int k = 0; k < cuad.nR(); k++) {
		(*theta)(0, k) = cuad(k, 0);
		(*weight)(0, k) = cuad(k, 1);
	}

	// build parameter set
	model->getParameterModel()->buildParameterSet(model->getItemModel(),model->getDimensionModel());
	profiler->stopTimer("input");
	profiler->resetTimer("initial");
	profiler->startTimer("initial");
	// Create estimation
	EMEstimation *em = new EMEstimation();
	em->setProfiler(profiler);
	//Here is where quadratures must be set.
	//create the quad nodes
	QuadratureNodes nodes(theta,weight);
	em->setQuadratureNodes(&nodes);
	em->setModel(model);

	Matrix<double> mat_initialValues(model->getItemModel()->countItems(), 3);
	input.importCSV(initialValues, mat_initialValues, 1, 0);
	Matrix<double> *a_init = new Matrix<double>(1,(model->getItemModel()->countItems()));
	Matrix<double> *b_init = new Matrix<double>(1,(model->getItemModel()->countItems()));
	Matrix<double> *c_init = new Matrix<double>(1,(model->getItemModel()->countItems()));

	for (int k = 0; k < model->getItemModel()->countItems() ; k++) {
		(*a_init)(0,k) = mat_initialValues(k,0);
		(*b_init)(0,k) = mat_initialValues(k,1);
		(*c_init)(0,k) = mat_initialValues(k,2);
	}
	double *** matrix_initial;

	switch (model_const) {
	case Constant::RASCH_A1:
		//matrix_initial[b] = b_init;
		break;
	case Constant::RASCH_A_CONSTANT:
		//matrix_initial[a] = a_init;
		//matrix_initial[d] = b_init;
		break;
	case Constant::TWO_PL:
		matrix_initial = new double** [2];
		//Pass a
		matrix_initial[0] = new double* [1];
		matrix_initial[0][0] = new double [a_init->nC()];
		for (int var = 0; var < a_init->nC(); ++var) {
			matrix_initial[0][0][var] = (*a_init)(0,var);
		}
		//Pass b
		matrix_initial[1] = new double* [1];
		matrix_initial[1][0] = new double [b_init->nC()];
		for (int var = 0; var < b_init->nC(); ++var) {
			matrix_initial[1][0][var] = (*b_init)(0,var);
		}
		break;
	case Constant::THREE_PL:
		//matrix_initial[a] = a_init;
		//matrix_initial[d] = b_init;
		//matrix_initial[c] = c_init;
		int items = model->getItemModel()->countItems();
		cout<<" My items : "<<items<<endl;

		matrix_initial = new double** [3];
		matrix_initial[0] = new double *[1];
		matrix_initial[1] = new double *[1];
		matrix_initial[2] = new double *[1];

		matrix_initial[0][0] = new double [items];
		matrix_initial[1][0] = new double [items];
		matrix_initial[2][0] = new double [items];
		for (int var = 0; var < a_init->nC(); ++var) {
			matrix_initial[0][0][var] = (*a_init)(0,var);
			cout<<matrix_initial[0][0][var]<<endl;
		}
		for (int var = 0; var < b_init->nC(); ++var) {
			matrix_initial[1][0][var] = (*b_init)(0,var);
			cout<<matrix_initial[1][0][var]<<endl;
		}
		for (int var = 0; var < c_init->nC(); ++var) {
			matrix_initial[2][0][var] = (*c_init)(0,var);
			cout<<matrix_initial[2][0][var]<<endl;
		}
		cout<<"Passing vals"<<endl;
		//Aqui esta el error
		break;
	}



	//Delete the matrix initial aFTER THIS REMEMBER TODO
	em->setInitialValues(matrix_initial);
	// run estimation
	cout<<"Passed vals"<<endl;

	profiler->stopTimer("initial");
	em->setProfiler(profiler);
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
	Timer tm;
	tm.reset();
	tm.start();
	if (argc < 2) {
		cout << "Please specify an input file" << endl;
		return (0);
	}
	//oneRun(argv[1]);
	runArgs(argv[1],argv[2]);
	tm.stop();
	cout<<"time: "<<endl<<tm.totalTime<<endl;
	return (0);
}

