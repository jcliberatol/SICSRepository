#include "interface.h"
#define ESTIMATION_MODEL Constant::RASCH_A_CONSTANT

void estimatingParameters(int ** dataI, int nRowsDataI, int nColumnsDataI, char * modelI, int dimI, char * initValI, double epsilonConvI, int maxIterI, bool verboseI, double *parametersO) {
	Input input;
	Constant::EPSILON = epsilonConvI;
	Matrix<double> cuad(41, 2);
	//Create the profiler to profile the program
	Trace* profiler = new Trace("Profile.log");
	profiler->resetTimer("total");
	profiler->startTimer("total");
	profiler->resetTimer("input");
	profiler->startTimer("input");
	profiler->resetTimer("for1");
	profiler->resetTimer("for2");
	profiler->resetTimer("fyr");
	profiler->resetTimer("optim");
	input.importCSV((char *) "Cuads.csv", cuad, 1, 0);
	// **** **** Run model complete and ordered process **** ****
	// Create general pars
	Model *model = new Model();
	// Create general model
	ModelFactory *modelFactory = new SICSGeneralModel();
	PatternMatrix *dataSet = new PatternMatrix(0);
	// Load matrix
	//input.importCSV(args, *dataSet, 1, 0);

	//copy matrix
	for( int _i_ = 0; _i_ < nRowsDataI; _i_++ )
	{
		vector<char> dset(nColumnsDataI);
		for ( int _j_ = 0; _j_ < nColumnsDataI; _j_++ )
		{
			if ( dataI[_i_][_j_] == 1) dset[nColumnsDataI- _j_ - 1] = true;  // taken of Input.cpp
		}
		dataSet->size = nColumnsDataI;
		dataSet->push(dset);
	}


	//

	// set dataset
	//RASCH_A1, RASCH_A_CONSTANT, TWO_PL, THREE_PL
	model->setModel(modelFactory, ESTIMATION_MODEL);
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
	em->setProfiler(profiler);
	QuadratureNodes nodes(theta, weight);
	em->setQuadratureNodes(&nodes);
	em->setModel(model);
	em->setInitialValues(Constant::ANDRADE); // initvalI here!
	profiler->stopTimer("initial");
	//Pass the profiler to the estimation object so it can be used to profile each step
	em->setProfiler(profiler);
	//Run the estimation
	em->estimate();
	model->parameterModel->getParameters(parametersO);
	delete modelFactory;
	delete dataSet;
	delete em;
	delete model;
	profiler->stopTimer("total");
	//Out the profiler here
	profilerOut(profiler, 4);
	delete profiler;
}

void profilerOut(Trace* profile, int type) {
	//Types of profiling :
	//	1: Only total time
	//  2: General
	//  3: Verbose
	//  4: Verbose with percents
	if (type == 1) {
		(*profile)("Profiling file : ", 'n');
		(*profile)(profile->filename);
		(*profile)("Total time :", 'n');
		(*profile)(profile->dr("total") / 1000000000);

	} else if (type == 2) {
		profilerOut(profile, 1);
		(*profile)("Input time :", 'n');
		(*profile)(profile->dr("input") / 1000000000);
		(*profile)("Initial time :", 'n');
		(*profile)(profile->dr("initial") / 1000000000);
	} else if (type == 3) {
		profilerOut(profile, 2);
		(*profile)("Estimation time :", 'n');
		(*profile)(profile->dr("estim") / 1000000000);
		//Estep time
		(*profile)("E step time:", 'n');
		(*profile)(profile->dr("Et") / 1000000000);
		//Mstep time
		(*profile)("M step time:", 'n');
		(*profile)(profile->dr("Mt") / 1000000000);
	} else if (type == 4) {
		(*profile)("_______________________________");
		(*profile)("Profiling file : ", 'n');
		(*profile)(profile->filename);
		(*profile)("Total time :", 'n');
		float total = (float) profile->dr("total") / 1000000000;
		(*profile)((float) profile->dr("total") / 1000000000);
		(*profile)("Input time :", 'n');
		(*profile)((float) profile->dr("input") / 1000000000);
		(*profile)("Initial time :", 'n');
		(*profile)((float) profile->dr("initial") / 1000000000);
		(*profile)("Estimation time :", 'n');
		(*profile)((float) profile->dr("estim") / 1000000000);
		//Estep time
		(*profile)("E step time:", 'n');
		(*profile)((float) profile->dr("Et") / 1000000000);
		(*profile)("First for time:", 'n');
		(*profile)((float) profile->dr("for1") / 1000000000);
		(*profile)("Second for time:", 'n');
		(*profile)((float) profile->dr("for2") / 1000000000);
		//Mstep time
		(*profile)("M step time:", 'n');
		(*profile)((float) profile->dr("Mt") / 1000000000);
		(*profile)("F and R transference time:", 'n');
		(*profile)((float) profile->dr("fyr") / 1000000000);
		(*profile)("optim time:", 'n');
		(*profile)((float) profile->dr("optim") / 1000000000);

	} else if (type == 5) {
		//Starts pacho profile
		(*profile)("File:", 'n');
		//Filename
		(*profile)(profile->readMessage("filename"));
		//Item parameters
		(*profile)("Item parameters");
		//Pattern Parameters

		//Convergence time
		//Iterations
	} else {
		(*profile)("No profiling selected please select a profiling mode");
	}
}

