//============================================================================
// Name        : SICS.cpp
// Author      : Juan Liberato
// Version     :
// Copyright   : Undefined
// Description : Hello World in C++, Ansi-style
//============================================================================

#include "Main.h"


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


void pachotest(){
	Input input;
	Matrix<double> pachotest(3,3);
	input.importCSV((char *) "pacho.csv", pachotest, 0, 0);
	ApproximateMatrixInverse(pachotest);
}

void initTests(char *filename){
	string testConfFile;
	testConfFile.assign(filename);
	//string testConfFile = "/home/mirt/Documentos/tests/config.conf";
	EMTest *emTest = new EMTest (testConfFile );
	emTest->runTest();
}

void oneRun(){
	Input input;
		Matrix<double> cuad(41, 2);
		input.importCSV((char *) "Cuads.csv", cuad, 1, 0);
		// **** **** Run model complete and ordered process **** ****
		cout << "Imported cuadratures"<< cuad << endl;
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

		cout << "Initial Pars setted" << endl;

		// Create estimation
		EMEstimation *em = new EMEstimation();
		em->setModel(model);
		cout << "Model setted" << endl;
		em->setInitialValues("Andradinho");
		// run estimation
		em->estimate();

		delete modelFactory;
		delete dataSet;
		delete model;
		delete em;
}
int main(int argc, char *argv[]) {

	cout<<"The show begins";
	oneRun();
	//initTests(argv[1]);
	//oneRun();
}


