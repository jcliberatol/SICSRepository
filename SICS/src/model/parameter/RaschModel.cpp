/*
 * RaschModel.cpp
 *
 *  Created on: 17 Jun 2014
 *      Author: jlgpisa
 */

#include <model/parameter/RaschModel.h>

RaschModel::RaschModel() {

	parameterSet[b] = NULL;
	probabilityMatrix=NULL;
	nodes = 0;

}

void RaschModel::buildParameterSet(ItemModel* itemModel, DimensionModel* dimensionModel) {

	if (typeid(*itemModel) == typeid(DichotomousModel)) {

		if (typeid(*dimensionModel) == typeid(UnidimensionalModel)) {

			int items = itemModel->countItems();
			parameterSet[b] = new Matrix<double>(1, items);

		}

		else if (typeid(*dimensionModel) == typeid(MultidimensionalModel)) {
			// TODO: Dichotomous Multidimensional
		}

		else if (typeid(*dimensionModel) == typeid(MultiUniDimModel)) {
			// TODO: Dichotomous MultiUniDimensional
		}

	}

	else if (typeid(*dimensionModel) == typeid(PolytomousModel)) {
		// TODO: Polytomous Model for Unidimensional, Multidimensional and MultiUni
	}
}


void RaschModel::setEstimationNodes(QuadratureNodes* n) {
	this->nodes = n;
}

void RaschModel::successProbability(DimensionModel *dimensionModel, QuadratureNodes * quadNodes) {

	int q = 0;
	double b_d, theta_d; // d stands from "double"

	if ( dimensionModel != NULL ) {
		q = quadNodes->size();
	}

	if(typeid(*dimensionModel)==typeid(UnidimensionalModel)) {
		int It = parameterSet[b]->nC();
		if(probabilityMatrix == NULL){
			//Creates the matrix if it is not already created
			probabilityMatrix = new Matrix<double>(q,It);
		}
		for (int k = 0; k < q; k++) {
			for ( int i = 0; i < It; i++ ){
				// Rasch Success Probability Function
				theta_d = (*quadNodes->getTheta())(0,k);
				b_d = (*parameterSet[b])(0,i);
				//c_d = (*parameterSet[c])(0,i);
				double p_d = successProbability ( theta_d, b_d );
				//cout<<a_d<<" "<<d_d<<" "<<c_d<<" "<<theta_d<<" "<<p_d<<" the prox"<<endl;
				(*probabilityMatrix)(k,i) = p_d;

			}
		}
	}
}


double RaschModel::successProbability(double theta, double b) {

	long double exponential = (Constant::NORM_CONST)*(-theta + b);

	if ( exponential > Constant::MAX_EXP ) {
		exponential = Constant::MAX_EXP;
	}

	else if ( exponential < -(Constant::MAX_EXP*1.0) ) {
		exponential = -Constant::MAX_EXP;
	}

	exponential = exp(-exponential) ;

	return (1 / (1 + exponential)) ;

	//return ((ec/(1+ec))+((ec)/((1+ec)*(1+exponential))));
	//return (c + (1.0 - c)/(1.0 + exponential));
}


void RaschModel::setParameterSet(map<Parameter, Matrix<double> *> pair) {
	this->parameterSet = parameterSet;
}


map<Parameter, Matrix<double> *> RaschModel::getParameterSet()  {
	return (this->parameterSet);
}

double RaschModel::getProbability(int node, int item) {
	return (0);
}


double RaschModel::logLikelihood (double* args, double* pars, int nargs,
		int npars) {

	//args
	/*
	 * a[i]
	 * b[i]
	 * c[i]
	 */

	//pars
	/*
	 * q
	 * I
	 * theta[q]
	 * f[q]
	 * r[q*I]
	 */

	int nA = 0;
	int nP = 0;

	int q, It;
	double *theta, *r, *f, *b;


	// Obtain q
	q = pars[nP ++]; // q is obtained and npars is augmented

	// Obtain I
	It = pars[nP ++];

	theta = new double[q];
	r = new double[q*It];
	f = new double[q];
	b = new double[It];



	// Obtain theta
	for (int k=0; k<q; k++) {
		theta[k] = pars[nP ++];
	}

	// Obtain f
	for (int k=0; k<q; k++) {
		f[k] = pars[nP ++];
	}

	// Obtain r
	for (int k=0; k<q; k++) {
		for (int i=0; i<It; i++) {
			r[k*It+i] = pars[nP ++];
		}
	}

	// Obtain b
	for (int i=0; i<It; i++) {
		b[i] = args [nA ++];
	}

	long double tp , tq;
	long double sum = 0;

	for (int k = 0; k < q; ++k) {
		for (unsigned int i = 0; i < It; ++i) {
			tp = (RaschModel::successProbability ( theta[k],  b[i]));
			if (tp==0)tp=1e-08;
			tq = 1-tp;
			if (tq==0)tq=1e-08;
			//suma = suma + (rki*logg(pki)+(fki-rki)*logg(qki))
			sum+=(r[k * It + i]*log(tp))+(f[k]-r[k * It + i])*log(tq);
		}
	}

	//antiLogit(c, I);
	delete[] theta;
	delete[] f;
	delete[] r;
	delete[] b;

	return (-sum);
}



RaschModel::~RaschModel() {
	if (parameterSet[b] != NULL) {
		delete parameterSet[b];
	}
}

