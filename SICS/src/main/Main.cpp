//============================================================================
// Name        : SICS.cpp
// Author      : Juan Liberato
// Version     :
// Copyright   : Undefined
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>
#include <type/Matrix.h>
#include <boost/dynamic_bitset.hpp>
#include <type/PatternMatrix.h>
#include <model/Model.h>
#include <model/ModelFactory.h>
#include <model/SICSGeneralModel.h>
#include <estimation/classical/EMEstimation.h>
#include <input/Input.h>
#include <time.h>
#include <estimation/bayesian/LatentTraitEstimation.h>
#include <type/LatentTraits.h>
#include <util/fitness/ItemFit.h>
#include <util/fitness/PersonFit.h>

//#define ESTIMATION_MODEL Constant::THREE_PL
//RASCH_A1, RASCH_A_CONSTANT, TWO_PL, THREE_PL
//#define ESTIMATION_MODEL Constant::THREE_PL
#define ESTIMATION_MODEL Constant::TWO_PL
//#define ESTIMATION_MODEL Constant::RASCH_A1
//#define ESTIMATION_MODEL Constant::RASCH_A_CONSTANT

void oneRun(char * args) {
	Input input;
	Matrix<double> cuad(41, 2);
	//Create the profiler to profile the program
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
	EMEstimation em;
	//Here is where quadratures must be set.
	//create the quad nodes
	QuadratureNodes nodes(theta, weight);
	em.setQuadratureNodes(&nodes);
	em.setModel(model);
	em.setInitialValues(Constant::ANDRADE);
	//Pass the profiler to the estimation object so it can be used to profile each step
	//Run the estimation
	em.estimate();

	/*
	 * Now we will run the estimation of individual parameter
	 */
	//first create the latenTrait objects
	//LatentTraits * latentTraits;
	//latentTraits = new LatentTraits(dataSet);
	//Now create the estimation
	//LatentTraitEstimation * lte = new LatentTraitEstimation();
	//Pass the model
	//lte->setModel(model);
	//Pass the latent traits
	//lte->setLatentTraits(latentTraits);
	//Pass the quadrature nodes
	//lte->setQuadratureNodes(&nodes);
	//Ready to estimate
	//lte->estimateLatentTraitsEAP();
	//lte->estimateLatentTraitsMAP();
	//finished
	//now read the latent traits but we will do this later
	//lte->getLatentTraits()->print();

	//Matrix<double> data(dataSet->countIndividuals(), dataSet->countItems());
	//input.importCSV(args, data, 1, 0);
	
	//double* itemsf = new double[ data.nC()];
	//itemFit(latentTraits->pm, *(latentTraits->traits), data, model->getParameterModel()->getParameterSet(), model -> type,itemsf);
	//personFit(latentTraits->pm, *(latentTraits->traits), data, model->getParameterModel()->getParameterSet(), model -> type);

	delete modelFactory;
	delete dataSet;
	delete model;
}

int main(int argc, char *argv[])
{
	if (argc < 2) {
		cout << "Please specify an input file" << endl;
		return (0);
	}
	return (0);
}

