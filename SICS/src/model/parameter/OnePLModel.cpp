/*
 * OnePLModel.cpp
 *
 *  Created on: Nov 16, 2014
 *      Author: anmrodriguezre
 */

#include <model/parameter/OnePLModel.h>

OnePLModel::OnePLModel() {

	parameterSet = NULL;
	probabilityMatrix = NULL;
}

string OnePLModel::getStringParameters() {
	return ("stringPars");
}
void OnePLModel::getParameters(double * parameters)
{
	for ( int i = 0; i < items; i++ )
	{
		parameters[i] = parameterSet[1][0][i];
	}
}
inline void OnePLModel::successProbability(DimensionModel *dimensionModel,
		QuadratureNodes * quadNodes) {

	int q = 0;
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
				(*probabilityMatrix)(k, i) = successProbability(
						(*quadNodes->getTheta())(0, k),
						(parameterSet[0][0][i]));
			}
		}
	}
}

inline double OnePLModel::successProbability(double theta, double b) {

	long double exponential = ((theta) - b);

	if (exponential > Constant::MAX_EXP) {
		exponential = Constant::MAX_EXP;
	} else if (exponential < -(Constant::MAX_EXP)) {
		exponential = -Constant::MAX_EXP;
	}

	return (1 / (1.0 + exp(-exponential)));
}

inline double OnePLModel::successProbability(double theta, double * zita) {

	long double exponential = ((theta) - zita[0]);

	if (exponential > Constant::MAX_EXP) {
		exponential = Constant::MAX_EXP;
	} else if (exponential < -(Constant::MAX_EXP)) {
		exponential = -Constant::MAX_EXP;
	}

	return (1 / (1.0 + exp(-exponential)));
}

void OnePLModel::setParameterSet(double*** par) {
	this->parameterSet = par;
}

double*** OnePLModel::getParameterSet() {
	return (this->parameterSet);
}

double OnePLModel::getProbability(int node, int item) {
	return ((*probabilityMatrix)(node, item));
}

void OnePLModel::printParameterSet(ostream& out) {
	out << "\"a\" \"b\" \"c\"" << "\n";
	for (int k = 0; k < items; k++) {
		out << 1 << " " << (parameterSet[0][0][k]) << " " << 0 << "\n";
	}
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

	int k, i, q, It;
	double *b;

	// Obtain q
	q = pars[nP++]; // q is obtained and npars is augmented

	// Obtain I
	It = pars[nP++];

	b = new double[It];

	// Obtain b
	for (i = 0; i < It; i++) {
		b[i] = args[nA++];
	}

	long double tp, tq;
	long double sum = 0;

	for (k = 0; k < q; ++k) {
		for (i = 0; i < It; ++i) {
			tp = (OnePLModel::successProbability(pars[k + 2], b[i]));
			if (tp == 0)
				tp = 1e-08;
			tq = 1 - tp;
			if (tq == 0)
				tq = 1e-08;
			sum += (pars[k * It + i + 2 + (2 * q)] * log(tp))
					+ (pars[k + 2 + q] - pars[k * It + i + 2 + 2 * q])
							* log(tq);
		}
	}

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
	int i, k, q, items;
	double *b;

	// Obtain q
	q = pars[nP++]; // q is obtained and npars is augmented
	// Obtain I
	items = pars[nP++];
	b = new double[items];

	// Obtain b
	for (i = 0; i < items; i++) {
		b[i] = args[nA++];
	}

	long double *h; // Block vector of size I (i.e. I blocks). Each block-element has size of 1
	long double *P;  // Matrix of size q*I
	long double *factor;	  // Matrix of product (r-fP)

	h = new long double[items];
	P = new long double[q * items];
	factor = new long double[q * items];

	for (k = 0; k < q; k++) {
		for (i = 0; i < items; i++) {
			P[k * items + i] = successProbability(pars[k + 2], b[i]);
			factor[k * items + i] = (pars[k * items + i + 2 + 2 * q]
					- pars[k + 2 + q] * P[k * items + i]);
		}
	}
	memset(h, 0, sizeof(long double) * items);
	memset(gradient, 0, sizeof(double) * items);
	for (i = 0; i < items; i++) {
		for (k = 0; k < q; k++) {
			h[i] += factor[k * items + i];
		}
	}

	delete[] P;
	delete[] factor;
	delete[] b;
	for (int i = 0; i < items; ++i) {
		gradient[i] = static_cast<double>(h[i]);
	}
	delete[] h;
}

OnePLModel::~OnePLModel() {

	if (parameterSet != NULL) {
		delete[] parameterSet[0][0];
		delete[] parameterSet[0];
		delete[] parameterSet;
	}
}
