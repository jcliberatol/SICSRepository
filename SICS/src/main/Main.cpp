//============================================================================
// Name        : SICS.cpp
// Author      : Juan Liberato
// Version     :
// Copyright   : Undefined
// Description : Hello World in C++, Ansi-style
//============================================================================

#include "Main.h"
#include <stdlib.h>

#define ESTIMATION_MODEL Constant::RASCH_A_CONSTANT

//void profilerOut(Trace* profile, int type){
//	//Types of profiling :
//	//	1: Only total time
//	//  2: General
//	//  3: Verbose
//	//  4: Verbose with percents
//	if(type==1){
//		(*profile)("Profiling file : ",'n');
//		(*profile)(profile->filename);
//		(*profile)("Total time :",'n');
//		(*profile)(profile->dr("total")/1000000000);
//
//	}
//	else if(type==2){
//		profilerOut(profile,1);
//		(*profile)("Input time :",'n');
//		(*profile)(profile->dr("input")/1000000000);
//		(*profile)("Initial time :",'n');
//		(*profile)(profile->dr("initial")/1000000000);
//	}
//	else if(type==3){
//		profilerOut(profile,2);
//		(*profile)("Estimation time :",'n');
//		(*profile)(profile->dr("estim")/1000000000);
//		//Estep time
//		(*profile)("E step time:",'n');
//		(*profile)(profile->dr("Et")/1000000000);
//		//Mstep time
//		(*profile)("M step time:",'n');
//		(*profile)(profile->dr("Mt")/1000000000);
//	}
//	else if(type == 4){
//		(*profile)("_______________________________");
//		(*profile)("Profiling file : ",'n');
//		(*profile)(profile->filename);
//		(*profile)("Total time :",'n');
//		float total = (float)profile->dr("total")/1000000000;
//		(*profile)((float)profile->dr("total")/1000000000);
//		(*profile)("Input time :",'n');
//		(*profile)((float)profile->dr("input")/1000000000);
//		(*profile)("Initial time :",'n');
//		(*profile)((float)profile->dr("initial")/1000000000);
//		(*profile)("Estimation time :",'n');
//		(*profile)((float)profile->dr("estim")/1000000000);
//		//Estep time
//		(*profile)("E step time:",'n');
//		(*profile)((float)profile->dr("Et")/1000000000);
//		(*profile)("First for time:",'n');
//		(*profile)((float)profile->dr("for1")/1000000000);
//		(*profile)("Second for time:",'n');
//		(*profile)((float)profile->dr("for2")/1000000000);
//		//Mstep time
//		(*profile)("M step time:",'n');
//		(*profile)((float)profile->dr("Mt")/1000000000);
//		(*profile)("F and R transference time:",'n');
//		(*profile)((float)profile->dr("fyr")/1000000000);
//		(*profile)("optim time:",'n');
//		(*profile)((float)profile->dr("optim")/1000000000);
//
//	}
//	else if (type == 5){
//		//Starts pacho profile
//		(*profile)("File:",'n');
//		//Filename
//		(*profile)(profile->readMessage("filename"));
//		//Item parameters
//		(*profile)("Item parameters");
//		//Pattern Parameters
//
//		//Convergence time
//		//Iterations
//	}
//	else{
//		(*profile)("No profiling selected please select a profiling mode");
//	}
//}
//
//void oneRun(char * args) {
//	Input input;
//	Matrix<double> cuad(41, 2);
//	//Create the profiler to profile the program
//	Trace* profiler = new Trace("Profile.log");
//	profiler->resetTimer("total");
//	profiler->startTimer("total");
//	profiler->resetTimer("input");
//	profiler->startTimer("input");
//	profiler->resetTimer("for1");
//	profiler->resetTimer("for2");
//	profiler->resetTimer("fyr");
//	profiler->resetTimer("optim");
//	input.importCSV((char *) "Cuads.csv", cuad, 1, 0);
//	// **** **** Run model complete and ordered process **** ****
//	// Create general pars
//	Model *model = new Model();
//	// Create general model
//	ModelFactory *modelFactory = new SICSGeneralModel();
//	PatternMatrix *dataSet = new PatternMatrix(0);
//	// Load matrix
//	input.importCSV(args, *dataSet, 1, 0);
//	// set dataset
//	//RASCH_A1, RASCH_A_CONSTANT, TWO_PL, THREE_PL
//	model->setModel(modelFactory, ESTIMATION_MODEL);
//	//This is where it is decided what model is the test to make
//	model->getItemModel()->setDataset(dataSet);		//Sets the dataset.
//	// set Theta and weight for the EM Estimation
//	Matrix<double> *theta = new Matrix<double>(1, 41);
//	Matrix<double> *weight = new Matrix<double>(1, 41);
//	for (int k = 0; k < cuad.nR(); k++) {
//		(*theta)(0, k) = cuad(k, 0);
//		(*weight)(0, k) = cuad(k, 1);
//	}
//
//
//	// build parameter set
//	model->getParameterModel()->buildParameterSet(model->getItemModel(),
//			model->getDimensionModel());
//
//	// Create estimation
//	profiler->stopTimer("input");
//	profiler->resetTimer("initial");
//	profiler->startTimer("initial");
//	EMEstimation *em = new EMEstimation();
//	//Here is where quadratures must be set.
//	//create the quad nodes
//	em->setProfiler(profiler);
//	QuadratureNodes nodes(theta, weight);
//	em->setQuadratureNodes(&nodes);
//	em->setModel(model);
//	em->setInitialValues(Constant::ANDRADE);
//	profiler->stopTimer("initial");
//	//Pass the profiler to the estimation object so it can be used to profile each step
//	em->setProfiler(profiler);
//	//Run the estimation
//	em->estimate();
//	delete modelFactory;
//	delete dataSet;
//	delete em;
//	delete model;
//	profiler->stopTimer("total");
//	//Out the profiler here
//	profilerOut(profiler,4);
//	delete profiler;
//}
//
//void runArgs(char * filename,char * initialValues){
//	Input input;
//	Matrix<double> cuad(41, 2);
//	Trace* profiler = new Trace("Profile.log");
//	profiler->storeMessage	std::cout<<"hola1"<<endl;("filename",filename);
//	profiler->startCounter("iterations");
//	profiler->resetTimer("total");
//	profiler->startTimer("total");
//	profiler->resetTimer("input");
//	profiler->startTimer("input");
//	profiler->resetTimer("for1");
//	profiler->resetTimer("for2");
//	input.importCSV((char *) "Cuads.csv", cuad, 1, 0);
//	// **** **** Run model complete and ordered process **** ****
//	// Create general pars
//
//	Model *model = new Model();
//	// Create general model
//	ModelFactory *modelFactory = new SICSGeneralModel();
//	PatternMatrix *dataSet = new PatternMatrix(0);
//	// Load matrix
//	input.importCSV(filename, *dataSet, 1, 0);
//	// set dataset
//	//RASCH_A1, RASCH_A_CONSTANT, TWO_PL, THREE_PL
//	int model_const =  ESTIMATION_MODEL;
//	model->setModel(modelFactory, model_const);
//	//This is where it is decided what model is the test to make
//	model->getItemModel()->setDataset(dataSet);//Sets the dataset.
//	// set Theta and weight for the EM Estimation
//	Matrix<double> *theta = new Matrix<double>(1, 41);
//	Matrix<double> *weight = new Matrix<double>(1, 41);
//
//	for (int k = 0; k < cuad.nR(); k++) {
//		(*theta)(0, k) = cuad(k, 0);
//		(*weight)(0, k) = cuad(k, 1);
//	}
//
//	// build parameter set
//	model->getParameterModel()->buildParameterSet(model->getItemModel(),model->getDimensionModel());
//	profiler->stopTimer("input");
//	profiler->resetTimer("initial");
//	profiler->startTimer("initial");
//	// Create estimation
//	EMEstimation *em = new EMEstimation();
//	em->setProfiler(profiler);
//	//Here is where quadratures must be set.
//	//create the quad nodes
//	QuadratureNodes nodes(theta,weight);
//	em->setQuadratureNodes(&nodes);
//	em->setModel(model);
//
//	Matrix<double> mat_initialValues(model->getItemModel()->countItems(), 3);
//	input.importCSV(initialValues, mat_initialValues, 1, 0);
//	Matrix<double> *a_init = new Matrix<double>(1,(model->getItemModel()->countItems()));
//	Matrix<double> *b_init = new Matrix<double>(1,(model->getItemModel()->countItems()));
//	Matrix<double> *c_init = new Matrix<double>(1,(model->getItemModel()->countItems()));
//
//	for (int k = 0; k < model->getItemModel()->countItems() ; k++) {
//		(*a_init)(0,k) = mat_initialValues(k,0);
//		(*b_init)(0,k) = mat_initialValues(k,1);
//		(*c_init)(0,k) = mat_initialValues(k,2);
//	}
//	double *** matrix_initial;
//	int items = model->getItemModel()->countItems();
//	switch (model_const) {
//	case Constant::RASCH_A1:
//		matrix_initial = new double** [1];
//		matrix_initial[0] = new double* [1];
//		matrix_initial[0][0] = new double [items];
//		for (int var = 0; var < items; ++var) {
//			matrix_initial[0][0][var] = (*b_init)(0,var);
//		}
//		break;
//	case Constant::RASCH_A_CONSTANT:
//		matrix_initial = new double** [2];
//		matrix_initial[0] = new double* [1];
//		matrix_initial[1] = new double* [1];
//		matrix_initial[0][0] = new double [1];
//		matrix_initial[1][0] = new double [items];
//		matrix_initial[0][0][0] = (*a_init)(0,0);
//		for (int var = 0; var < items; ++var) {
//			matrix_initial[1][0][var] = (*b_init)(0,var);
//		}
//		break;
//	case Constant::TWO_PL:
//		matrix_initial = new double** [2];
//		matrix_initial[0] = new double* [1];
//		matrix_initial[1] = new double* [1];
//
//		matrix_initial[0][0] = new double [items];
//		matrix_initial[1][0] = new double [items];
//		for (int var = 0; var < items; ++var) {
//			matrix_initial[0][0][var] = (*a_init)(0,var);
//		}
//		for (int var = 0; var < items; ++var) {
//			matrix_initial[1][0][var] = (*b_init)(0,var);
//		}
//		break;
//	case Constant::THREE_PL:
//
//		matrix_initial = new double** [3];
//		matrix_initial[0] = new double *[1];
//		matrix_initial[1] = new double *[1];
//		matrix_initial[2] = new double *[1];
//
//		matrix_initial[0][0] = new double [items];
//		matrix_initial[1][0] = new double [items];
//		matrix_initial[2][0] = new double [items];
//		for (int var = 0; var < items; ++var) {
//			matrix_initial[0][0][var] = (*a_init)(0,var);
//		}
//		for (int var = 0; var < items; ++var) {
//			matrix_initial[1][0][var] = (*b_init)(0,var);
//		}
//		for (int var = 0; var < items; ++var) {
//			matrix_initial[2][0][var] = (*c_init)(0,var);
//		}
//		break;
//	}
//	//delete a_init;
//	//delete b_init;
//	//delete c_init;
//	//Delete the matrix initial aFTER THIS REMEMBER TODO
//	em->setInitialValues(matrix_initial);
//	// run estimation
//	profiler->stopTimer("initial");
//	em->setProfiler(profiler);
//	em->estimate();
//	delete modelFactory;
//	delete dataSet;
//	delete em;
//	delete model;
//	profiler->stopTimer("total");
//	//Out the profiler here
//	profilerOut(profiler,4);
//	delete profiler;
//	switch (model_const) {
//	case Constant::THREE_PL:
//		delete[] matrix_initial[0][0] ;
//		delete[] matrix_initial[1][0] ;
//		delete[] matrix_initial[2][0] ;
//
//		delete[] matrix_initial[0] ;
//		delete[] matrix_initial[1] ;
//		delete[] matrix_initial[2] ;
//
//		delete[] matrix_initial ;
//		break;
//	case Constant::TWO_PL :
//		//delete[] matrix_initial[0][0] ;
//		//delete[] matrix_initial[1][0] ;
//		//delete[] matrix_initial[0] ;
//		//delete[] matrix_initial[1] ;
//		//delete[] matrix_initial ;
//		break;
//	case Constant::RASCH_A1 :
//		//delete[] matrix_initial[0][0] ;
//		//delete[] matrix_initial[0] ;
//		//delete[] matrix_initial ;
//		break;
//	}
//	//delete weight;
//	//delete theta;
//}

int main(int argc, char *argv[]) {
//	Timer tm;
//	tm.reset();
//	tm.start();
//	if (argc < 2) {
//		cout << "Please specify an input file" << endl;
//		return (0);
//	}
//	//oneRun(argv[1]);
//	runArgs(argv[1],argv[2]);
//	//oneRun(argv[1]);
//	tm.stop();
//	//cout<<"time: "<<endl<<tm.totalTime<<endl;

	double *TestHessianaInput;
	double **TestHessianaOutput;
	int sizeHessiana = 5;
	int auxSize = (sizeHessiana*(sizeHessiana+1))/2;
	TestHessianaInput = new double[auxSize];
	TestHessianaOutput= new double *[sizeHessiana];
	for ( int i = 0; i < auxSize; i++ ) TestHessianaInput[i] = i+1;
	for ( int i = 0; i<  sizeHessiana; i++ ) TestHessianaOutput[i] = new double[sizeHessiana];
	UTIL_H_::transformHessiana(TestHessianaInput, TestHessianaOutput, sizeHessiana);
	std::cout<<"Test of transform of Hessiana"<<endl;
	for ( int i = 0; i<  sizeHessiana; i++ )
	{
		for ( int j = 0;j < sizeHessiana; j++ )
		{
			std::cout<<TestHessianaOutput[i][j]<<" ";
		}
		std::cout<<""<<endl;
	}
	std::cout<<"end of Test of transform of Hessiana\n\n"<<endl;
	const int items = 5, peoples = 5;
	int **DataI;
	DataI = new int *[peoples];
	for ( int i = 0; i < peoples; i++ ) DataI[i] = new int[items];
	char  *initValues;
	for ( int i = 0; i < peoples; i++ ) for ( int j = 0; j<  items; j++ ) DataI[i][j] = rand()%2;
	initValues = "ANDRADE";
	double epsilon = 0.001;
	int maxNIteration = 200;
	bool verbose = true;
	double *parameters;
	parameters = new double[items+1];
	int numberOfCycles = -1;
	double logLik = -1;
	double convEp = -1;
	int model = Constant::RASCH_A_CONSTANT;

	estimatingParameters(DataI, peoples, items, model, 1, initValues, epsilon, maxNIteration, verbose, parameters, numberOfCycles, logLik, convEp);
	std::cout<<"print parameters"<<endl;
	for ( int i = 0;i  <= items; i++ )
	{
		std::cout<<parameters[i] <<" ";
	}
	std::cout<<endl;
	std::cout<<"pŕint number of cycles"<<endl;
	std::cout<<numberOfCycles<<endl;
	std::cout<<"print convEp"<<endl;
	std::cout<<convEp<<endl;
	std::cout<<"print loglik"<<endl;
	std::cout<<logLik<<endl;
	return (0);
}
