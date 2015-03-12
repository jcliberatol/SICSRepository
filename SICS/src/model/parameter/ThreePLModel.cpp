/*
 * ThreePLModel.cpp
 *
 *  Created on: 17 Jun 2014
 *      Author: jlgpisa
 */

#include <model/parameter/ThreePLModel.h>
#include <string>

ThreePLModel::ThreePLModel() {
	parameterSet = NULL;
	probabilityMatrix=NULL;
	nodes = NULL;
}

string ThreePLModel::getStringParameters(){
	//If the parameter set is null
	if(parameterSet == NULL){
		return ("Error In Parameter set return please move the getStringParameters Method to another point of execution");
	}
	//If the parameter set exists
	string pars = "";
	for (int i = 0; i < items; ++i) {
		pars = pars + "Item : " + std::to_string(i);
		pars = pars + " " + std::to_string(parameterSet[0][0][i]) + "," +
				std::to_string(parameterSet[1][0][i]) + ","+ std::to_string(parameterSet[2][0][i]) + " ";
	}
	parameterSet[0][0] = new double [items];
	parameterSet[1][0] = new double [items];
	parameterSet[2][0] = new double [items];
	return ("stringPars");
}

void ThreePLModel::setEstimationNodes(QuadratureNodes* n) {
	this->nodes = n;
}

void ThreePLModel::successProbability(DimensionModel *dimensionModel, QuadratureNodes * quadNodes) {

	int q = 0;
	double a_d, d_d, c_d, theta_d; // d stands from "double"

	if ( dimensionModel != NULL ) {
		q = quadNodes->size();
	}

	if(typeid(*dimensionModel)==typeid(UnidimensionalModel)) {
		if(probabilityMatrix == NULL){
				//Creates the matrix if it is not already created
				probabilityMatrix = new Matrix<double>(q,items);
			}
		for (int k = 0; k < q; k++) {
			for ( int i = 0; i < items; i++ ){
				// 3PL Success Probability Function
				theta_d = (*quadNodes->getTheta())(0,k);
				a_d = parameterSet[0][0][i];
				d_d = parameterSet[1][0][i];
				c_d = parameterSet[2][0][i];
				double p_d = successProbability ( theta_d, a_d, d_d, c_d );
				(*probabilityMatrix)(k,i) = p_d;

			}
		}
	}

}

double *** ThreePLModel::getParameterSet()  {
	return (this->parameterSet);
}

void ThreePLModel::setParameterSet(double ***) {
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
	double ec = exp(c);
	return ( (ec/(1+ec)) + (1 - (ec/(1+ec))) * (1/(1+exponential)) );
}
void ThreePLModel::getParameters(double * parameters)
{
	int i = 0;
	for (int j = 0; j < items; j++) {
		parameters[i++] = parameterSet[0][0][j];
	}
	for (int j = 0; j < items; j++) {
		parameters[i++] = parameterSet[1][0][j];
	}
	for (int j = 0; j < items; j++) {
		parameters[i++] = parameterSet[2][0][j];
	}
}
double ThreePLModel::successProbability(double theta, double * zita) {

	long double exponential = (Constant::NORM_CONST)*(zita[0]*theta+zita[1]);

	if ( exponential > Constant::MAX_EXP ) {
		exponential = Constant::MAX_EXP;
	}

	else if ( exponential < -(Constant::MAX_EXP*1.0) ) {
		exponential = -Constant::MAX_EXP;
	}

	exponential = exp(-exponential) ;
	double ec = exp(zita[2]);
	return ( (ec/(1+ec)) + (1 - (ec/(1+ec))) * (1/(1+exponential)) );
}

double ThreePLModel::getProbability(int node, int item) {
	return ((*probabilityMatrix)(node, item));
}


void ThreePLModel::gradient (double* args, double* pars, int nargs, int npars, double* gradient){
	/*
	 *
	 * What we need
	 * items
	 * q
	 * theta array
	 * D
	 * a, b, c
	 * f and r
	 */
	int nA = 0;
	int nP = 0;
	int q, items;
	double *theta, *r, *f, *a, *b, *c;

	// Obtain q
	q = pars[nP ++]; // q is obtained and npars is augmented
	// Obtain I
	items = pars[nP ++];

	theta = new double[q];
	r = new double[q*items];
	f = new double[q];
	a = new double[items];
	b = new double[items];
	c = new double[items];
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
		for (int i=0; i<items; i++) {
			r[k*items+i] = pars[nP ++];
		}
	}

	// Obtain a
	for (int i=0; i<items; i++) {
		a[i] = args [nA ++];
	}
	// Obtain b
	for (int i=0; i<items; i++) {
		b[i] = args [nA ++];
	}
	// Obtain c
	for (int i=0; i<items; i++) {
		c[i] = args [nA ++];
	}

	double D = Constant::NORM_CONST;
	long double *h_0; // Block Matrix of size q*I. Each block-element has size of 1*3
	long double *h; // Block vector of size I (i.e. I blocks). Each block-element has size of 1*3
	long double *P_Star, *P;  // Matrix of size q*I
	long double *W;           // Matrix of size q*I
	long double *factor;	  // Matrix of product (r-fP)W
	long double *ec;            // e^c_i
	long double *ecPlus1Inv;	// 1 / (e^c_i + 1)

	h = new long double [3*items];
	h_0 = new long double [q*3*items];
	P = new long double [q*items];
	P_Star = new long double [q*items];
	factor = new long double [q*items];
	W = new long double [q*items];
	ec = new long double [items];
	ecPlus1Inv = new long double [items];

	for( unsigned  int i = 0; i < items; i++ ) {
		ecPlus1Inv[i]=1/(1+exp(c[i]));
		ec[i]=exp(c[i]);
	}
	for ( int k = 0; k < q; k++ ) {
		for ( unsigned  int i = 0; i < items; i++ ) {

			P[k * items + i] = successProbability ( theta[k], a[i], b[i], c[i] );
			//P_Star[k * items + i] = successProbability ( theta[k], a[i], b[i], 0.0 );
			P_Star[k * items + i] = 1/(1+exp(-D*(a[i]*theta[k]+b[i])));

			W[k * items + i] = P_Star[k * items + i] * ( 1 - P_Star[k * items + i] ); // Numerator
			W[k * items + i] /= P[k * items + i] * ( 1 - P[k * items + i] );// Denominator

			factor[k * items + i] = ( r[k * items + i] - f[k]*P[k * items + i] ) * W[k * items + i];

			// h_0 / (P_star*Q_star)
			h_0[3 * items * k + 3 * i + 0] = D * theta[k] * ecPlus1Inv[i];
			h_0[3 * items * k + 3 * i + 1] = D * ecPlus1Inv[i];
			h_0[3 * items * k + 3 * i + 2] = ec[i] * (ecPlus1Inv[i]*ecPlus1Inv[i]) / P_Star[k * items + i];

		}
	}
	memset(h,0,sizeof(long double)*3*items);
	memset(gradient,0,sizeof(double)*3*items);
	for ( unsigned int i = 0; i < items; i++ ) {
		for ( int k = 0; k < q; k++ ) {
			h[3 * i + 0] += factor[k * items + i] * h_0[3 * items * k + 3 * i + 0];
			h[3 * i + 1] += factor[k * items + i] * h_0[3 * items * k + 3 * i + 1];
			h[3 * i + 2] += factor[k * items + i] * h_0[3 * items * k + 3 * i + 2];
		}
	}

	delete [] h_0;
	delete [] P_Star;
	delete [] P;
	delete [] W;
	delete [] factor;
	delete [] ec;
	delete [] ecPlus1Inv;

	delete [] theta;
	delete [] r;
	delete [] f;
	delete [] a;
	delete [] b;
	delete [] c;

//return h as the gradient
	int hc=0;
	for (int n = 0; n < 3; ++n) {
		for(int i = 0 ; i < items ; ++i){
			gradient[hc++]= -static_cast<double>(h[i*3+n]);
		}
	}
	delete [] h;
}
double ThreePLModel::logLikelihood (double* args, double* pars, int nargs,
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
	double *theta, *r, *f, *a, *b, *c;


	// Obtain q
	q = pars[nP ++]; // q is obtained and npars is augmented

	// Obtain I
	It = pars[nP ++];

	theta = new double[q];
	r = new double[q*It];
	f = new double[q];
	a = new double[It];
	b = new double[It];
	c = new double[It];


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

	// Obtain a
	for (int i=0; i<It; i++) {
		a[i] = args [nA ++];
	}

	// Obtain b
	for (int i=0; i<It; i++) {
		b[i] = args [nA ++];
	}

	// Obtain c
	for (int i=0; i<It; i++) {
		c[i] = args [nA ++];
	}

	long double tp , tq;
	long double sum = 0;

	for (int k = 0; k < q; ++k) {
		for (unsigned int i = 0; i < It; ++i) {
			tp = (ThreePLModel::successProbability ( theta[k], a[i], b[i], c[i]));
			if (tp==0)tp=1e-08;
			tq = 1-tp;
			if (tq==0)tq=1e-08;
			sum+=(r[k * It + i]*log(tp))+(f[k]-r[k * It + i])*log(tq);
		}
	}
	delete[] theta;
	delete[] f;
	delete[] r;
	delete[] a;
	delete[] b;
	delete[] c;

	return (-sum);


}

double ThreePLModel::successProbability_cPrime(double theta, double a, double b,
		double c) {
	long double cPrime = exp(c)/(1+exp(c));
	return (successProbability ( theta, a, b, cPrime ));
}

ThreePLModel::~ThreePLModel() {
	if (parameterSet != NULL) {

		delete[] parameterSet[2][0];
		delete[] parameterSet[1][0];
		delete[] parameterSet[0][0];

		delete[] parameterSet[0];
		delete[] parameterSet[1];
		delete[] parameterSet[2];

		delete[] parameterSet;
	}


}

void ThreePLModel::printParameterSet(ostream& out) {
	out << "\"a\" \"b\" \"c\"" << endl;

	for (int i = 0; i < items; i++) {
		out << parameterSet[0][0][i] << " " << parameterSet[1][0][i] << " " << parameterSet[2][0][i]
				<< endl;
	}

}
