//============================================================================
// Name        : SICS.cpp
// Author      : Juan Liberato
// Version     :
// Copyright   : Undefined
// Description : Hello World in C++, Ansi-style
//============================================================================

#include "Main.h"
#include "util/blasInterface.h"

using namespace std;

double banana(double* args, double* pars, int nargs,
		int npars){
	double x = args[0];
	double y = args[1];
	return ((1-x)*(1-x) + 100*(y-x*x)*(y-x*x));
}

void bananaGradient (double* args, double* pars, int nargs, int npars, double* gradient){
	double h = 0.00000000001;
		//(f(x+h)-f(x))/h
		for(int i = 0 ; i < nargs; i++){
			args[i]=args[i]+h;
			gradient[i]=banana(args,pars,nargs,npars);
			args[i]=args[i]-h;
			gradient[i]-=banana(args,pars,nargs,npars);
			gradient[i]=gradient[i]/h;
		}
}

void rosenbrockTest(){
	//
	//(1-x)2 + 100(y-x2)2
	int nargs = 2;
	int npars = 0;
	double values[2];
	double gradient[2];
	values[0]=3;
	values[1]=4;
	gradient[0]=0;
	gradient[1]=0;
	cout<< "Value of the banana "<<banana(values,gradient,nargs,npars)<<endl;
	bananaGradient(values,gradient,nargs,npars,gradient);
	cout<< "Value of the gradient "<<gradient[0]<<" "<<gradient[1]<<endl;
	double (*fptr)(double*,double*,int,int);
	void (*gptr)(double*,double*,int,int,double*);
	fptr = &banana;
	gptr = &bananaGradient;
	Optimizer* optim;
	optim = new Optimizer();
	optim->searchOptimal(fptr,gptr,gptr,values,gradient,2,2);
	cout<< "Value of the banana "<<banana(values,gradient,nargs,npars)<<endl;
	cout<< "At points  "<<values[0]<<" "<<values[1]<<endl;
	bananaGradient(values,gradient,nargs,npars,gradient);
	cout<< "Value of the gradient "<<gradient[0]<<" "<<gradient[1]<<endl;
	cout<<"Changed values to optimal"<<endl;
	values[0]=1;
	values[1]=1;

	cout<< "Value of the banana "<<banana(values,gradient,nargs,npars)<<endl;
	cout<< "At points  "<<values[0]<<" "<<values[1]<<endl;
	bananaGradient(values,gradient,nargs,npars,gradient);
	cout<< "Value of the gradient "<<gradient[0]<<" "<<gradient[1]<<endl;
}


int main() {
	//rosenbrockTest();
	Input input;
	Matrix<double> pachotest(3,3);
	input.importCSV((char *) "pacho.csv", pachotest, 0, 0);
	ApproximateMatrixInverse(pachotest);
	return (0);
	cout<<1-ThreePLModel::successProbability(-5.1225,1.55162,1.52326,-1.44746)<<" The pii"<<endl;
	Matrix<double> cuad(41, 2);
	input.importCSV((char *) "Cuads.csv", cuad, 1, 0);
	// **** **** Run model complete and ordered process **** ****
	cout << "Imported cuadratures" << cuad<< endl;
	// Create general pars
	int It; // Number of items

	// Create general model
	Model *model = new Model();
	ModelFactory *modelFactory = new SICSGeneralModel();
	PatternMatrix *dataSet = new PatternMatrix();
	cout << "Created model" << endl;
	model->setModel(modelFactory);

	// Load matrix
	input.importCSV((char *) "Test_num_1_1000_individuos.csv", *dataSet, 1, 0);
	// set dataset
	cout << "Dataset size : " << (*dataSet).countItems() << " x "
			<< (*dataSet).countIndividuals() << endl;
	model->getItemModel()->setDataset(dataSet);
	//cout<<*dataSet<<endl;
	// set Theta and weight
	Matrix<double> *theta = new Matrix<double>(1, 41);
	Matrix<double> *weight = new Matrix<double>(1, 41);

	for (int k = 0; k < cuad.nR(); k++) {
		(*theta)(0, k) = cuad(k, 0);
		(*weight)(0, k) = cuad(k, 1);
	}
	model->getDimensionModel()->getLatentTraitSet()->setTheta(theta);
	model->getDimensionModel()->getLatentTraitSet()->setWeight(weight);
	cout << "Quadratures transferred" << endl;

	// build parameter set
	model->getParameterModel()->buildParameterSet(model->getItemModel(),
			model->getDimensionModel());
	cout << "Loaded input matrix" << endl;
	It = model->getItemModel()->countItems();

	// Initial Parameters
	for (int i = 0; i < It; i++) {
		(*model->getParameterModel()->getParameterSet()[a])(0, i) = 0.851;
		(*model->getParameterModel()->getParameterSet()[d])(0, i) = 0.272;
		(*model->getParameterModel()->getParameterSet()[c])(0, i) = 0.15;
	}
	/*
	(*model->getParameterModel()->getParameterSet()[a])(0, 0) = 0.6609297;
	(*model->getParameterModel()->getParameterSet()[a])(0, 1) = 0.8741322;
	(*model->getParameterModel()->getParameterSet()[a])(0, 2) = 0.7340038;
	(*model->getParameterModel()->getParameterSet()[a])(0, 3) = 0.5422908;
	(*model->getParameterModel()->getParameterSet()[a])(0, 4) = 0.9116418;
	(*model->getParameterModel()->getParameterSet()[a])(0, 5) = 0.7105918;
	*/
	(*model->getParameterModel()->getParameterSet()[d])(0, 0) = 1.612291;
	(*model->getParameterModel()->getParameterSet()[d])(0, 1) = 0.4928179;
	(*model->getParameterModel()->getParameterSet()[d])(0, 2) = 0.6055251;
	(*model->getParameterModel()->getParameterSet()[d])(0, 3) = -1.294814;
	(*model->getParameterModel()->getParameterSet()[d])(0, 4) = 0.1676967;
	(*model->getParameterModel()->getParameterSet()[d])(0, 5) = 1.445592;

	cout << "Initial Pars setted" << endl;

	// Create estimation
	EMEstimation *em = new EMEstimation();
	em->setModel(model);
	cout << "Model setted" << endl;
	// run estimation
	em->estimate();

	delete modelFactory;
	delete dataSet;
	delete model;
	delete em;

	return (0);
}


