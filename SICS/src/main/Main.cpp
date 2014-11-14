//============================================================================
// Name        : SICS.cpp
// Author      : Juan Liberato
// Version     :
// Copyright   : Undefined
// Description : Hello World in C++, Ansi-style
//============================================================================

#include "Main.h"


using namespace std;

void oneRun(){
	Input input;
		Matrix<double> cuad(41, 2);
		input.importCSV((char *) "Cuads.csv", cuad, 1, 0);
		// **** **** Run model complete and ordered process **** ****
		// Create general pars
		int It; // Number of items

		Model *model = new Model();
		// Create general model
		ModelFactory *modelFactory = new SICSGeneralModel();
		PatternMatrix *dataSet = new PatternMatrix ();
		// Load matrix
		input.importCSV((char *) "Test_10_1_1000.csv", *dataSet, 1, 0);
		// set dataset
		model->setModel(modelFactory, Constant::THREE_PL);
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

		// Create estimation
		EMEstimation *em = new EMEstimation();
		//Here is where quadratures must be set.
		//create the quad nodes
		QuadratureNodes nodes(theta,weight);
		em->setQuadratureNodes(&nodes);
		em->setModel(model);
		cout << "Model setted" << endl;
		em->setInitialValues("ANDRADE");
		// run estimation
		em->estimate();
		delete modelFactory;
		delete dataSet;
		delete em;
		delete model;
}
int main(int argc, char *argv[]) {

	oneRun();
}


