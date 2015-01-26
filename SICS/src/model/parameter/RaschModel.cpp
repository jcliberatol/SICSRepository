/*
 * RaschModel.cpp
 *
 *  Created on: 17 Jun 2014
 *      Author: jlgpisa
 */

#include <model/parameter/RaschModel.h>

RaschModel::RaschModel() {

	parameterSet = NULL;
	probabilityMatrix = NULL;
	nodes = 0;

}

string RaschModel::getStringParameters(){
	return ("stringPars");
}

void RaschModel::buildParameterSet(ItemModel* itemModel,
		DimensionModel* dimensionModel) {

	if (typeid(*itemModel) == typeid(DichotomousModel)) {

		if (typeid(*dimensionModel) == typeid(UnidimensionalModel)) {

			int items = itemModel->countItems();
			parameterSet = new double**[1];
			parameterSet[0] = new double*[1];
			parameterSet[0][0] = new double[items];
			//parameterSet[b] = new Matrix<double>(1, items);

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

void RaschModel::successProbability(DimensionModel *dimensionModel,
		QuadratureNodes * quadNodes) {

	int q = 0;
	double b_d, theta_d; // d stands from "double"

	if (dimensionModel != NULL) {
		q = quadNodes->size();
	}

	if (typeid(*dimensionModel) == typeid(UnidimensionalModel)) {
		if (probabilityMatrix == NULL) {
			//Creates the matrix if it is not already created
			probabilityMatrix = new Matrix<double>(q, items);
		}
		for (int k = 0; k < q; k++) {
			for (int i = 0; i < items; i++) {
				// Rasch Success Probability Function
				theta_d = (*quadNodes->getTheta())(0, k);
				//b_d = (*parameterSet[b])(0, i);
				b_d = parameterSet[0][0][i];
				double p_d = successProbability(theta_d, b_d);
				//cout<<a_d<<" "<<d_d<<" "<<c_d<<" "<<theta_d<<" "<<p_d<<" the prox"<<endl;
				(*probabilityMatrix)(k, i) = p_d;

			}
		}
	}
}

double RaschModel::successProbability(double theta, double b) {

	long double exponential = (Constant::NORM_CONST) * (-theta + b);

	if (exponential > Constant::MAX_EXP) {
		exponential = Constant::MAX_EXP;
	}

	else if (exponential < -(Constant::MAX_EXP * 1.0)) {
		exponential = -Constant::MAX_EXP;
	}

	exponential = exp(-exponential);

	return (1 / (1 + exponential));

	//return ((ec/(1+ec))+((ec)/((1+ec)*(1+exponential))));
	//return (c + (1.0 - c)/(1.0 + exponential));
}

void RaschModel::setParameterSet(double*** pair) {
	this->parameterSet = pair;
}

double*** RaschModel::getParameterSet() {
	return (this->parameterSet);
}

double RaschModel::getProbability(int node, int item) {
	return ((*probabilityMatrix)(node, item));
}

double RaschModel::logLikelihood(double* args, double* pars, int nargs,
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
			tp = (RaschModel::successProbability(theta[k], b[i]));
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

void RaschModel::Hessian(double* args, double* pars, int nargs, int npars,
		double* hessian) {

}
*/
void RaschModel::Ngradient(double* args, double* pars, int nargs, int npars,
		double* gradient) {
	//	For each of the gradient thingies increase the args and apply richardsons
	double hh = 0.000001;
	//(f(x+h)-f(x))/h
	for (int i = 0; i < nargs; i++) {
		args[i] = args[i] + hh;
		gradient[i] = logLikelihood(args, pars, nargs, npars);
		args[i] = args[i] - hh;
		gradient[i] -= logLikelihood(args, pars, nargs, npars);
		gradient[i] = gradient[i] / hh;
	}
}

void RaschModel::NHessian(double*args, double* pars, int nargs, int npars,
		double* hessian) {
	//the gradient is composed of items a's items b's items c's take it as such
	double hh = 0.001;
	double gradiente[nargs];
	int items = nargs / 3;
	gradient(args, pars, nargs, npars, gradiente);
	//Calculate the hessian for each item order is a's b's c's in args
	for (int i = 0; i < items; ++i) {
		//Extract the gradient at the points (Put at the hessian)
		for (int j = 0; j < 1; ++j) {
			for (int k = 0; k < 1; ++k) {
				//Now for each of these point change the k  parameter for deriving
				hessian[i * items + j * 3 + k] = gradiente[j * items + i];
			}
		}
	}

	for (int i = 0; i < items; ++i) {
		//Derivate two times for this item // matrix point
		//Extract the gradient at the points (Put at the hessian)
		for (int j = 0; j < 1; ++j) {
			for (int k = 0; k < 1; ++k) {
				memset(gradiente, 0, sizeof(double) * items);
				//Change the argument k in the item i
				args[k * items + i] += hh;
				//Calculate le gradient
				gradient(args, pars, nargs, npars, gradiente);
				//Reset the parameter
				args[k * items + i] -= hh;
				//This is the gradient at the point
				hessian[i * items + j + k] -= gradiente[j * items + i];
				hessian[i * items + j + k] = hessian[i * items + j + k] / hh;
			}
		}
	}

}
void RaschModel::gradient(double* args, double* pars, int nargs, int npars,
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

RaschModel::~RaschModel() {

	if (parameterSet != NULL) {
		delete[] parameterSet[0][0];
		delete[] parameterSet[0];
		delete[] parameterSet;


	}
}

