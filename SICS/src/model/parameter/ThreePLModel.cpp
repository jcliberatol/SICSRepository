/*
 * ThreePLModel.cpp
 *
 *  Created on: 17 Jun 2014
 *      Author: jlgpisa
 */

#include <model/parameter/ThreePLModel.h>

ThreePLModel::ThreePLModel() {

	parameterSet[a] = NULL;
	parameterSet[b] = NULL;
	parameterSet[c] = NULL;
	parameterSet[d] = NULL;

}

void ThreePLModel::buildParameterSet(ItemModel* itemModel,
		DimensionModel* dimensionModel) {

	if (typeid(*itemModel) == typeid(DichotomousModel)) {

		if (typeid(*dimensionModel) == typeid(UnidimensionalModel)) {

			int items = itemModel->countItems();

			parameterSet[a] = new Matrix<double>(1, items);
			parameterSet[d] = new Matrix<double>(1, items);
			parameterSet[c] = new Matrix<double>(1, items);

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

void ThreePLModel::successProbability(DimensionModel *dimensionModel) {

	int q = 0;
	double a_d, d_d, c_d, theta_d; // d stands from "double"

	if ( dimensionModel != NULL ) {
		q = dimensionModel->getLatentTraitSet()->getTheta()->nC();
	}

	if(typeid(*dimensionModel)==typeid(UnidimensionalModel)) {
		int I = parameterSet[a]->nC();

		for (int k = 0; k < q; k++) {
			for ( int i = 0; i < I; i++ ){

				// 3PL Success Probability Function
				theta_d = (*dimensionModel->getLatentTraitSet()->getTheta())(0,k);
				a_d = (*parameterSet[a])(0,i);
				d_d = (*parameterSet[d])(0,i);
				c_d = (*parameterSet[c])(0,i);

				double p_d = successProbability ( theta_d, a_d, d_d, c_d );
				(*probabilityMatrix)(k,i) = p_d;

			}
		}
	}

}

map<Parameter, Matrix<double> *> ThreePLModel::getParameterSet()  {
	return (this->parameterSet);
}

void ThreePLModel::setParameterSet(
		map<Parameter, Matrix<double> *> parameterSet) {
	this->parameterSet = parameterSet;
}

double ThreePLModel::successProbability(double theta, double a, double d,
		double c) {
	long double exponential = (Constant::NORM_CONST)*(a*theta+d);

	if ( exponential > Constant::MAX_EXP ) {
		exponential = Constant::MAX_EXP;
	}

	else if ( exponential < -(Constant::MAX_EXP*1.0) ) {
		exponential = -Constant::MAX_EXP;
	}

	exponential = exp(-exponential) ;

	return (c + (1.0 - c)/(1.0 + exponential));
}

double ThreePLModel::getProbability(int node, int item) {
	return ((*probabilityMatrix)(node, item));
}

double ThreePLModel::LogLikelihood (double* args, double* pars, int nargs,
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

	int q, I;
	double *theta, *r, *f, *a, *b, *c;

	// Obtain q
	q = pars[nP ++]; // q is obtained and npars is augmented

	// Obtain I
	I = pars[nP ++];

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
		for (int i=0; i<I; i++) {
			r[k*I+i] = pars[nP ++];
		}
	}

	// Obtain a
	for (int i=0; i<I; i++) {
		a[i] = args [nA ++];
	}

	// Obtain b
	for (int i=0; i<I; i++) {
		b[i] = args [nA ++];
	}

	// Obtain c
	for (int i=0; i<I; i++) {
		c[i] = args [nA ++];
	}

	long double tp , tq;
	long double sum = 0;

	for (int k = 0; k < q; ++k) {
		for (unsigned int i = 0; i < I; ++i) {
			tp = successProbability_cPrime ( theta[k], a[i], b[i], c[i] );
			if (tp==0)tp=1e-08;
			tq = 1-tp;
			if (tq==0)tq=1e-08;
			//suma = suma + (rki*logg(pki)+(fki-rki)*logg(qki))
			sum+=(r[k * I + i]*log(tp))+(f[k]-r[k * I + i])*log(tq);
		}
	}
	//antiLogit(c, I);
	return (-sum);


}

double ThreePLModel::successProbability_cPrime(double theta, double a, double b,
		double c) {
	long double cPrime = exp(c)/(1+exp(c));
	return (successProbability ( theta, a, b, cPrime ));
}

ThreePLModel::~ThreePLModel() {

	if (parameterSet[a] != NULL) {
		delete parameterSet[a];
	}
	if (parameterSet[b] != NULL) {
		delete parameterSet[b];
	}
	if (parameterSet[c] != NULL) {
		delete parameterSet[c];
	}
	if (parameterSet[d] != NULL) {
		delete parameterSet[d];
	}

}

