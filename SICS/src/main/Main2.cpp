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
#include <util/util.h>

#define ESTIMATION_MODEL Constant::THREE_PL

void oneRun(char * args)
{
  Input input;
  Matrix<double> cuad(10, 2);
  cuad.reset();
  Model *model = new Model();
  ModelFactory *modelFactory;
  PatternMatrix *dataSet;
  Matrix<double> *theta;
  Matrix<double> *weight;

  modelFactory = new SICSGeneralModel();
  dataSet = new PatternMatrix(0);
  theta = new Matrix<double>(1, 10);
  weight = new Matrix<double>(1, 10);

  input.importCSV((char *) "quads10.csv", cuad, 1, 0);
 std::cout<<cuad<<std::endl;
  input.importCSV(args, *dataSet, 1, 0);
  //Start by telling the model that it is a multidimensional model.
  int dimstype = 2;
  model->setModel(modelFactory, ESTIMATION_MODEL, dimstype);
  //Here we must set the number of dimensions to estimate in the given model
  model->getDimensionModel()->setDims(2);
  int dims = model->getDimensionModel()->getNumDimensions();
  std::cout<<"Dims used : "<<dims<<std::endl;
  delete modelFactory;

  model->getItemModel()->setDataset(dataSet);

  for (int k = 0; k < cuad.nR(); k++)
  {
    (*theta)(0, k) = cuad(k, 0);
    (*weight)(0, k) = cuad(k, 1);
  }

  // build parameter set
    std::cout << "Building parameterSet" << std::endl;
model->getParameterModel()->buildParameterSet(model->getItemModel(), model->getDimensionModel());

  // Create estimation
std::cout << "Creating EMEstimation" << std::endl;
  EMEstimation em;
  //Here is where quadratures must be set.
  //create the quad nodes
std::cout << "Declaring QuadratureNodes" << std::endl;
  QuadratureNodes nodes(theta, weight);
std::cout << "Setting them" << std::endl;
  em.setQuadratureNodes(&nodes);
std::cout << "Setting Model" << std::endl;
  em.setModel(model);
    std::cout << "Ready to ANDRADE" << std::endl;
  em.setInitialValues(Constant::ANDRADE);

  //Run the estimation
std::cout << "em.estimate" << std::endl;
  em.estimate();

  double fulloglik = em.getLoglik();
      std::cout<<"Ll : "<<fulloglik<<std::endl;
      double* returnpars;
      double* pars;


      returnpars = new double[(dims+2)*dataSet->size];
      pars = new double[(dims+2)*dataSet->size];
      model->parameterModel->getParameters(returnpars);


          // Return in list
          for (int i = 0; i < (dims+2)*dataSet->size; i++){
              pars[i] = returnpars[i];
              std::cout<<pars[i]<<" . "<<std::endl;
      }

  /*
  * Now we will run the estimation of individual parameter
  */
  //Now create the estimation
 // LatentTraitEstimation lte(dataSet);
  //Pass the model
  //lte.setModel(model);
  //Pass the quadrature nodes
  //lte.setQuadratureNodes(&nodes);
  //Ready to estimate
  //lte.estimateLatentTraitsEAP();
  //lte.estimateLatentTraitsMAP();
  //finished
  //now read the latent traits but we will do this later
  //lte.getLatentTraits()->print();

  //Matrix<double> data(dataSet->countIndividuals(), dataSet->countItems());
  //input.importCSV(args, data, 1, 0);

  //double* itemsf = new double[ data.nC()];
  //itemFit(latentTraits->pm, *(latentTraits->traits), data, model->getParameterModel()->getParameterSet(), model -> type,itemsf);
  //personFit(latentTraits->pm, *(latentTraits->traits), data, model->getParameterModel()->getParameterSet(), model -> type);

  delete dataSet;
  delete model;
  delete theta;
  delete weight;
}

int main(int argc, char *argv[]) {
  if (argc < 2) {
    cout << "Please specify an input file" << endl;
    return (0);
  }
  oneRun(argv[1]);
  return (0);
}
