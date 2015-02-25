/*
 * TwoPLModel.cpp
 *
 *  Created on: 18 Jun 2014
 *      Author: jlgpisa
 */

#include <model/parameter/TwoPLModel.h>

//
TwoPLModel::TwoPLModel() {

	parameterSet = NULL;
	probabilityMatrix = NULL;

}

string TwoPLModel::getStringParameters() {
	return ("stringPars");
}

void TwoPLModel::buildParameterSet(ItemModel* itemModel,
		DimensionModel* dimensionModel) {

	if (typeid(*itemModel) == typeid(DichotomousModel)) {

		if (typeid(*dimensionModel) == typeid(UnidimensionalModel)) {
			parameterSet = new double **[2]; // two parameters

			parameterSet[0] = new double*[1]; // parameter a
			parameterSet[1] = new double*[1]; // parameter b
			items = itemModel->countItems();
			parameterSet[0][0] = new double[items];
			parameterSet[1][0] = new double[items];

			//para<meterSet[a] = new Matrix<double>(1, items);
			//parameterSet[d] = new Matrix<double>(1, items);

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

inline void TwoPLModel::successProbability(DimensionModel *dimensionModel,
		QuadratureNodes *quadNodes) {
	int q = 0;
	double a_d, d_d, theta_d; // d stands from "double"

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
				double* theta = &(*quadNodes->getTheta())(0, k);
				(*probabilityMatrix)(k, i) = successProbability(theta,
						&parameterSet[0][0][i], &parameterSet[1][0][i]);
			}
		}
	}
}

inline double TwoPLModel::successProbability(double *theta, double *a,
		double *d) {

	long double exponential = (Constant::D_CONST * ((*a * *theta) + *d));
	if (exponential > Constant::MAX_EXP) {
		exponential = Constant::MAX_EXP;
	}

	else if (exponential < -(Constant::MAX_EXP)) {
		exponential = -Constant::MAX_EXP;
	}
	return (1 / (1 + exp(-exponential)));
}

double *** TwoPLModel::getParameterSet() {
	return (this->parameterSet);
}

void TwoPLModel::setParameterSet(double *** pair) {
	this->parameterSet = pair;
}

double TwoPLModel::getProbability(int node, int item) {
	return ((*probabilityMatrix)(node, item));
}

void TwoPLModel::gradient(double* args, double* pars, int nargs, int npars,
		double* gradient) {
	/*
	 s	 *
	 * What we need
	 * items
	 * q
	 * D
	 * a, b, c
	 */
	int nA = 0;
	int nP = 0;
	int q, items;
	double *a, *d, *c;

	// Obtain q
	q = pars[nP++]; // q is obtained and npars is augmented
	// Obtain I
	items = pars[nP++];
	a = new double[items];
	d = new double[items];

	// Obtain a
	for (int i = 0; i < items; i++) {
		a[i] = args[nA++];
	}
	// Obtain d
	for (int i = 0; i < items; i++) {
		d[i] = args[nA++];
	}

	long double *h_0;// Block Matrix of size q*I. Each block-element has size of 1*2
	long double *h;	// Block vector of size I (i.e. I blocks). Each block-element has size of 1*2
	long double *P;			// Matrix of size q*I
	long double *factor;	// Matrix of product (r-fP)

	h = new long double[2 * items];
	h_0 = new long double[q * 2 * items];
	factor = new long double[q * items];
	P = new long double[q * items];

	for (int k = 0; k < q; k++) {
		for (unsigned int i = 0; i < items; i++) {
			P[k * items + i] = successProbability(&pars[k + 2], &a[i], &d[i]);
			factor[k * items + i] = (pars[(k * items + i) + 2 + (2 * q)]
					- pars[k + 2 + q] * P[k * items + i]);

			h_0[2 * items * k + 2 * i + 0] = Constant::D_CONST * pars[k + 2];
			h_0[2 * items * k + 2 * i + 1] = Constant::D_CONST;
		}
	}

	memset(h, 0, sizeof(long double) * 2 * items);
	memset(gradient, 0, sizeof(double) * 2 * items);

	for (unsigned int i = 0; i < items; i++) {
		for (int k = 0; k < q; k++) {
			h[2 * i + 0] += factor[k * items + i]
					* h_0[2 * items * k + 2 * i + 0];
			h[2 * i + 1] += factor[k * items + i]
					* h_0[2 * items * k + 2 * i + 1];
		}
	}

	delete[] h_0;
	delete[] P;
	delete[] factor;
	delete[] a;
	delete[] d;
	int hc = 0;
	for (int n = 0; n < 2; ++n) {
		for (int i = 0; i < items; ++i) {
			gradient[hc++] = -static_cast<double>(h[i * 2 + n]);
		}
	}

	delete[] h;
}

void TwoPLModel::Ngradient(double* args, double* pars, int nargs, int npars,
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

void TwoPLModel::NHessian(double*args, double* pars, int nargs, int npars,
		double* hessian) {

	//the gradient is composed of items a's items d's take it as such
	double hh = 0.001;
	double gradiente[nargs];
	int items = nargs / 2;
	gradient(args, pars, nargs, npars, gradiente);
	//Calculate the hessian for each item order is a's b's c's in args
	for (int i = 0; i < items; ++i) {
		//Extract the gradient at the points (Put at the hessian)
		for (int j = 0; j < 2; ++j) {
			for (int k = 0; k < 2; ++k) {
				//Now for each of these point change the k  parameter for deriving
				hessian[i * items + j * 3 + k] = gradiente[j * items + i];
			}
		}
	}

	for (int i = 0; i < items; ++i) {
		//Derivate two times for this item // matrix point
		//Extract the gradient at the points (Put at the hessian)
		for (int j = 0; j < 2; ++j) {
			for (int k = 0; k < 2; ++k) {
				memset(gradiente, 0, sizeof(double) * items * 2);
				//Change the argument k in the item i
				args[k * items + i] += hh;
				//Calculate le gradient
				gradient(args, pars, nargs, npars, gradiente);
				//Reset the parameter
				args[k * items + i] -= hh;
				//This is the gradient at the point
				hessian[i * items + j * 2 + k] -= gradiente[j * items + i];
				hessian[i * items + j * 2 + k] = hessian[i * items + j * 2 + k]
						/ hh;
			}
		}
	}
}

void TwoPLModel::Hessian(double* args, double* pars, int nargs, int npars,
		double* hessian) {

}

double TwoPLModel::logLikelihood(double* args, double* pars, int nargs,
		int npars) {

//args
	/*
	 * a[i]
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
	double *a, *b;
	//double *theta, *r, *f, *a, *b;

// Obtain q
	q = pars[nP++]; // q is obtained and npars is augmented

// Obtain I
	It = pars[nP++];

//	theta = new double[q];
//	r = new double[q * It];
//	f = new double[q];
	a = new double[It];
	b = new double[It];

// Obtain theta
//	for (int k = 0; k < q; k++) {
//		theta[k] = pars[nP++];
//	}

// Obtain f
//	for (int k = 0; k < q; k++) {
//		f[k] = pars[nP++];
//	}

// Obtain r
//	for (int k = 0; k < q; k++) {
//		for (int i = 0; i < It; i++) {
//			r[k * It + i] = pars[nP++];
//		}
//	}

// Obtain a
	for (int i = 0; i < It; i++) {
		a[i] = args[nA++];
	}

// Obtain b
	for (int i = 0; i < It; i++) {
		b[i] = args[nA++];
	}

	long double tp, tq;
	long double sum = 0;

	for (int k = 0; k < q; ++k) {
		for (unsigned int i = 0; i < It; ++i) {
			//tp = (TwoPLModel::successProbability(&theta[k], &a[i], &b[i]));
			tp = (TwoPLModel::successProbability(&pars[k + 2], &a[i], &b[i]));
			if (tp == 0)
				tp = 1e-08;
			tq = 1 - tp;
			if (tq == 0)
				tq = 1e-08;

			//sum += (r[k * It + i] * log(tp)) + (f[k] - r[k * It + i]) * log(tq);
			sum += (pars[(k * It + i) + 2 + (2 * q)] * log(tp))
					+ (pars[k + 2 + q] - pars[(k * It + i) + 2 + (2 * q)])
							* log(tq);
		}
	}

//antiLogit(c, I);
//	delete[] theta;
//	delete[] f;
//	delete[] r;
	delete[] a;
	delete[] b;

	return (-sum);

}

TwoPLModel::~TwoPLModel() {
	if (parameterSet != NULL) {
		delete[] parameterSet[0][0];
		delete[] parameterSet[1][0];

		delete[] parameterSet[0];
		delete[] parameterSet[1];

		delete[] parameterSet;
	}
}

void TwoPLModel::printParameterSet(ostream& out) {
	out << "\"a\" \"b\" \"c\"" << endl;

	for (int i = 0; i < items; i++) {
		out << parameterSet[0][0][i] << " " << parameterSet[1][0][i] << " " << 0
				<< endl;
	}
//	cout<<"2PL Model Parameters :"<<endl;
//	cout<<"Discrimination parameter"<<endl;
//	for (int i = 0; i < items; ++i) {
//		cout<<parameterSet[0][0][i]<<" ";
//	}cout<<endl;
//	cout<<"Dificulty parameter"<<endl;
//	for (int i = 0; i < items; ++i) {
//			cout<<parameterSet[1][0][i]<<" ";
//	}cout<<endl;
}
void TwoPLModel::getParameters(double * parameters)
{
	int i = 0;
	for (int j = 0; i < items; i++) {
		parameters[i++] = parameterSet[0][0][j];
	}
	for (int j = 0; i < items; i++) {
		parameters[i++] = parameterSet[1][0][j];
	}
	for (int j = 0; i < items; i++) {
		parameters[i++] = parameterSet[2][0][j];
	}
}
