//============================================================================
// Name        : SICS.cpp
// Author      : Juan Liberato
// Version     :
// Copyright   : Undefined
// Description : Hello World in C++, Ansi-style
//============================================================================

#include "Main.h"

void oneRun(char * filename,char * initialValues){
	Input input;
		Matrix<double> cuad(41, 2);
		input.importCSV((char *) "Cuads.csv", cuad, 1, 0);
		// **** **** Run model complete and ordered process **** ****
		// Create general pars

		Model *model = new Model();
		// Create general model
		ModelFactory *modelFactory = new SICSGeneralModel();
		PatternMatrix *dataSet = new PatternMatrix ();
		// Load matrix
		input.importCSV(filename, *dataSet, 1, 0);
		//input.importCSV((char *) "Test_10_1_1000.csv", *dataSet, 1, 0);
		// set dataset
		//RASCH_A1, RASCH_A_CONSTANT, TWO_PL, THREE_PL
		model->setModel(modelFactory, Constant::RASCH_A1);
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
		//em->setInitialValues(Constant::ANDRADE);

		Matrix<double> mat_initialValues(model->getItemModel()->countItems(), 3);
		input.importCSV(initialValues, mat_initialValues, 1, 0);
		//Matrix<double> *a_init = new Matrix<double>(1,(model->getItemModel()->countItems()));
		Matrix<double> *b_init = new Matrix<double>(1,(model->getItemModel()->countItems()));
		//Matrix<double> *c_init = new Matrix<double>(1,(model->getItemModel()->countItems()));

		for (int k = 0; k < model->getItemModel()->countItems() ; k++) {
			//(*a_init)(0,k) = mat_initialValues(k,0);
			(*b_init)(0,k) = mat_initialValues(k,1);
			//(*c_init)(0,k) = mat_initialValues(k,2);
		}

		map<Parameter, Matrix<double> *> matrix_initial;
		//matrix_initial[a] = a_init;
		matrix_initial[b] = b_init;
		//matrix_initial[c] = c_init;

		em->setInitialValues(matrix_initial);

		// run estimation
		em->estimate();
		delete modelFactory;
		delete dataSet;
		delete em;
		delete model;
}
int main(int argc, char *argv[]) {
	//cout << argv[1] << " " << argv[2] << '\n';
	oneRun(argv[1], argv[2]);
	return (0);
}


