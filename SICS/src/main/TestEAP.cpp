/*
 * TestEAP.cpp
 *
 *  Created on: Mar 19, 2015
 *      Author: cristian
 */

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
#include <trace/Trace.h>
#include <estimation/bayesian/LatentTraitEstimation.h>
#include <type/LatentTraits.h>

//#define ESTIMATION_MODEL Constant::THREE_PL
//RASCH_A1, RASCH_A_CONSTANT, TWO_PL, THREE_PL
//#define ESTIMATION_MODEL Constant::THREE_PL
#define ESTIMATION_MODEL Constant::TWO_PL
//#define ESTIMATION_MODEL Constant::RASCH_A1
//#define ESTIMATION_MODEL Constant::RASCH_A_CONSTANT

void oneRun(char * args) {
	Input input;
	Matrix<double> cuad(41, 2);

	input.importCSV((char *) "Cuads.csv", cuad, 1, 0);
	Model *model = new Model();
	ModelFactory *modelFactory = new SICSGeneralModel();
	PatternMatrix *dataSet = new PatternMatrix(0);
	input.importCSV(args, *dataSet, 1, 0);
	model->setModel(modelFactory, ESTIMATION_MODEL);
	model->getItemModel()->setDataset(dataSet);
	Matrix<double> *theta = new Matrix<double>(1, 41);
	Matrix<double> *weight = new Matrix<double>(1, 41);

	for (int k = 0; k < cuad.nR(); k++)
	{
		(*theta)(0, k) = cuad(k, 0);
		(*weight)(0, k) = cuad(k, 1);
	}

	model->getParameterModel()->buildParameterSet(model->getItemModel(),
			model->getDimensionModel());

	EMEstimation *em = new EMEstimation();

	QuadratureNodes nodes(theta, weight);
	em->setQuadratureNodes(&nodes);
	em->setModel(model);
	em->setInitialValues(Constant::ANDRADE);

	em->estimate();

	LatentTraits * latentTraits;
	latentTraits = new LatentTraits(dataSet);
	LatentTraitEstimation * lte = new LatentTraitEstimation();
	lte->setModel(model);
	lte->setLatentTraits(latentTraits);
	lte->setQuadratureNodes(&nodes);
	lte->estimateLatentTraitsEAP();
	//lte->printVectors();
	//lte->getLatentTraits()->print();
	delete modelFactory;
	delete dataSet;
	delete em;
	delete model;
	//Out the profiler here
}

int main(int argc, char *argv[]) {
	if (argc < 2) {
		cout << "Please specify an input file" << endl;
		return (0);
	}
	oneRun(argv[1]);
	return (0);
}

