//============================================================================
// Name        : SICS.cpp
// Author      : Juan Liberato
// Version     :
// Copyright   : Undefined
// Description : Hello World in C++, Ansi-style
//============================================================================

#include "Main.h"

using namespace std;

int main() {

	Input input;
	Matrix<double> cuad(41, 2);
	input.importCSV((char *) "Cuads.csv", cuad, 1, 0);
	//cout << cuad;
	Matrix<double> pacho(3,3);
	input.importCSV((char *) "pacho.csv",pacho,0,0);
	cout<<pacho<<endl;
	NCM(pacho);
	cout<<pacho<<endl;

	// **** **** Run model complete and ordered process **** ****

	// Create general pars
	int It; // Number of items

	// Create general model
	Model *model = new Model();
	ModelFactory *modelFactory = new SICSGeneralModel();
	PatternMatrix *dataSet= new PatternMatrix();

	model->setModel(modelFactory);

	// Load matrix
	input.importCSV((char *) "input.csv", *dataSet, 0, 0);
	// set dataset
	model->getItemModel()->setDataset(dataSet);
	// build parameter set
	model->getParameterModel()->buildParameterSet(model->getItemModel(),
			model->getDimensionModel());

	It = model->getItemModel()->countItems();

	// Initial Parameters
	for (int i = 0; i < It; i++) {
		(*model->getParameterModel()->getParameterSet()[a])(0, i) = 0.851;
		(*model->getParameterModel()->getParameterSet()[d])(0, i) = 0.272;
		(*model->getParameterModel()->getParameterSet()[c])(0, i) = 0.2;
	}

	// set Theta and weight
	Matrix<double> *theta = new Matrix<double> (1, 41);
	Matrix<double> *weight = new Matrix<double>(1, 41);

	for (int k = 0; k < cuad.nR(); k++) {
		(*theta) (0,k) = cuad(k,0);
		(*weight) (0,k) = cuad(k,1);
	}
	model->getDimensionModel()->getLatentTraitSet()->setTheta(theta);
	model->getDimensionModel()->getLatentTraitSet()->setWeight(weight);

	// Create estimation
	EMEstimation *em = new EMEstimation();
	em->setModel(model);

	// run estimation
	em->stepE();

	delete modelFactory;
	delete dataSet;
	delete model;
	delete em;

	return (0);
}
