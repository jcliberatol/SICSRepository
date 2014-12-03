/*
 * OnePLModel.cpp
 *
 *  Created on: Nov 16, 2014
 *      Author: anmrodriguezre
 */

#include <model/parameter/OnePLModel.h>


OnePLModel::OnePLModel() {

	parameterSet[b] = NULL;
	probabilityMatrix = NULL;
	nodes = 0;

}

void OnePLModel::buildParameterSet(ItemModel* itemModel,
		DimensionModel* dimensionModel) {

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

void OnePLModel::setEstimationNodes(QuadratureNodes* n) {
	this->nodes = n;
}

void OnePLModel::successProbability(DimensionModel *dimensionModel,
		QuadratureNodes * quadNodes) {

	int q = 0;
	double b_d, theta_d; // d stands from "double"

	if (dimensionModel != NULL) {
		q = quadNodes->size();
	}

	if (typeid(*dimensionModel) == typeid(UnidimensionalModel)) {
		int It = parameterSet[b]->nC();
		if (probabilityMatrix == NULL) {
			//Creates the matrix if it is not already created
			probabilityMatrix = new Matrix<double>(q, It);
		}
		for (int k = 0; k < q; k++) {
			for (int i = 0; i < It; i++) {
				// Rasch Success Probability Function
				theta_d = (*quadNodes->getTheta())(0, k);
				b_d = (*parameterSet[b])(0, i);
				double p_d = successProbability(theta_d, b_d);
				//cout<<a_d<<" "<<d_d<<" "<<c_d<<" "<<theta_d<<" "<<p_d<<" the prox"<<endl;
				(*probabilityMatrix)(k, i) = p_d;

			}
		}
	}
}

double OnePLModel::successProbability(double theta, double b) {

	long double exponential = (Constant::NORM_CONST) * (-theta + b);

	if (exponential > Constant::MAX_EXP) {
		exponential = Constant::MAX_EXP;
	}

	else if (exponential < -(Constant::MAX_EXP * 1.0)) {
		exponential = -Constant::MAX_EXP;
	}

	exponential = exp(-exponential);

	return (1 / (1 + exponential));
}

void OnePLModel::setParameterSet(map<Parameter, Matrix<double> *> pair) {
	this->parameterSet = parameterSet;
}

map<Parameter, Matrix<double> *> OnePLModel::getParameterSet() {
	return (this->parameterSet);
}

double OnePLModel::getProbability(int node, int item) {
	return ((*probabilityMatrix)(node, item));
}

void OnePLModel::printParameterSet(ostream& out){
	/*out<<"Estimated parameters : "<<endl;*/
	out <<"\"a\" \"b\" \"c\""<<"\n";
	for (int k = 0; k <(*parameterSet[b]).nC(); k++) {
	 out<<1<<" "<<(*parameterSet[b])(0,k)<<" "<<0 <<"\n";
	}
	//out<<*parameterSet[b]<<endl;
}


double OnePLModel::logLikelihood(double* args, double* pars, int nargs,
		int npars) {

	//args
	/*
	 * a = 1
	 * b[i]
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
	q = pars[nP++]; // q is obtained and npars is augmented

	// Obtain I
	It = pars[nP++];

	theta = new double[q];
	r = new double[q * It];
	f = new double[q];
	b = new double[It];

	// Obtain theta
	for (int k = 0; k < q; k++) {
		theta[k] = pars[nP++];
	}

	// Obtain f
	for (int k = 0; k < q; k++) {
		f[k] = pars[nP++];
	}

	// Obtain r
	for (int k = 0; k < q; k++) {
		for (int i = 0; i < It; i++) {
			r[k * It + i] = pars[nP++];
		}
	}

	// Obtain b
	for (int i = 0; i < It; i++) {
		b[i] = args[nA++];
	}

	long double tp, tq;
	long double sum = 0;

	for (int k = 0; k < q; ++k) {
		for (unsigned int i = 0; i < It; ++i) {
			tp = (OnePLModel::successProbability(theta[k], b[i]));
			if (tp == 0)
				tp = 1e-08;
			tq = 1 - tp;
			if (tq == 0)
				tq = 1e-08;
			sum += (r[k * It + i] * log(tp)) + (f[k] - r[k * It + i]) * log(tq);
		}
	}

	delete[] theta;
	delete[] f;
	delete[] r;
	delete[] b;

	return (-sum);
}
/*

void OnePLModel::Hessian(double* args, double* pars, int nargs, int npars,
		double* hessian) {

}
*/


void OnePLModel::gradient(double* args, double* pars, int nargs, int npars,
		double* gradient) {
	/*
	 *
	 * What we need
	 * items
	 * q
	 * theta array
	 * b
	 * f and r
	 */
	int nA = 0;
	int nP = 0;
	int q, items;
	double *theta, *r, *f, *b;

	// Obtain q
	q = pars[nP++]; // q is obtained and npars is augmented
	// Obtain I
	items = pars[nP++];

	theta = new double[q];
	r = new double[q * items];
	f = new double[q];
	b = new double[items];

	// Obtain theta
	for (int k = 0; k < q; k++) {
		theta[k] = pars[nP++];
	}

	// Obtain f
	for (int k = 0; k < q; k++) {
		f[k] = pars[nP++];
	}

	// Obtain r
	for (int k = 0; k < q; k++) {
		for (int i = 0; i < items; i++) {
			r[k * items + i] = pars[nP++];
		}
	}

	// Obtain b
	for (int i = 0; i < items; i++) {
		b[i] = args[nA++];
	}

	long double *h; // Block vector of size I (i.e. I blocks). Each block-element has size of 1
	long double *P;  // Matrix of size q*I
	long double *factor;	  // Matrix of product (r-fP)

	h = new long double[items];
	P = new long double[q * items];
	factor = new long double[q * items];

	for (int k = 0; k < q; k++) {
		for (unsigned int i = 0; i < items; i++) {

			P[k * items + i] = successProbability(theta[k], b[i]);

			factor[k * items + i] =
					(r[k * items + i] - f[k] * P[k * items + i]);

		}
	}
	memset(h, 0, sizeof(long double) * items);
	memset(gradient, 0, sizeof(double) * items);
	for (unsigned int i = 0; i < items; i++) {
		for (int k = 0; k < q; k++) {
			h[i] += factor[k * items + i];
		}
	}

	delete[] P;
	delete[] factor;

	delete[] theta;
	delete[] r;
	delete[] f;
	delete[] b;

	//return h as the gradient
	int hc = 0;
	for (int i = 0; i < items; ++i) {
		gradient[hc++] = -static_cast<double>(h[i]);
	}

	delete[] h;
}

OnePLModel::~OnePLModel() {
	if (parameterSet[b] != NULL) {
		delete parameterSet[b];
	}
}
