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
	probabilityMatrix=NULL;

}

void ThreePLModel::buildParameterSet(ItemModel* itemModel,
		DimensionModel* dimensionModel) {

	if (typeid(*itemModel) == typeid(DichotomousModel)) {

		if (typeid(*dimensionModel) == typeid(UnidimensionalModel)) {

			int items = itemModel->countItems();
			int q = dimensionModel->getLatentTraitSet()->getTheta()->nC();
			probabilityMatrix = new Matrix<double>(q,items);
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
				//cout<<a_d<<" "<<d_d<<" "<<c_d<<" "<<theta_d<<" "<<p_d<<" the prox"<<endl;
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
	double ec = exp(c);
	return ( (ec/(1+ec)) + (1 - (ec/(1+ec))) * (1/(1+exponential)) );

	//return ((ec/(1+ec))+((ec)/((1+ec)*(1+exponential))));
	//return (c + (1.0 - c)/(1.0 + exponential));
}

double ThreePLModel::getProbability(int node, int item) {
	return ((*probabilityMatrix)(node, item));
}

void ThreePLModel::itemHessian(double* args, double* hess, int nargs, int i, double* thess){
	int It = thess[0];
	//fill the hessian for the item
	for (int h = 0 ; h < 9 ; h++){
		thess[h] = hess[i*It+h];
	}
}
void ThreePLModel::itemgradient(double* args, double* grad, int nargs, int i, double* tgrad){
    /*
     * In the item gradient and hessians these procedures take an already calculated hessian
     * and gradient and just access the item that is needed in the optimization step
     * args will then be the ith args for the item
     * pars is going to be the full array of gradients and npars is going to be used to take the ith element
     * the gradient will then be copied to g
     */
	//The first element of tgrad is the number of items in total
	int It = tgrad[0];
	//fill the gradient for the item
	tgrad[0] = grad[i];
	tgrad[1] = grad[It+i];
	tgrad[2] = grad[2*It+i];
}

void ThreePLModel::Hessian(double* args, double* pars, int nargs, int npars, double* hessian){
	//Hessian is composed of 3 * 3 matrix for every item so item * 9 memories

	/*
		 * TODO
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
		//memset(gradient,0,sizeof(double)*3*items);
		for ( unsigned int i = 0; i < items; i++ ) {
			for ( int k = 0; k < q; k++ ) {
				h[3 * i + 0] += factor[k * items + i] * h_0[3 * items * k + 3 * i + 0];
				h[3 * i + 1] += factor[k * items + i] * h_0[3 * items * k + 3 * i + 1];
				h[3 * i + 2] += factor[k * items + i] * h_0[3 * items * k + 3 * i + 2];
			}
		}

		//Now we must start calculating the hessian matrix
		long double *H_0;
		H_0 = new long double [q*9*items];
		memset(hessian,0,sizeof(double)*3*3*items);
		memset(H_0,0,sizeof(long double)*q*3*3*items);
		//#define H_0(k,i,m,n) H_0[9 * I * k + 9 * i + 3 * m + n]
		//#define H(i,m,n) H[9 * i + 3 * m + n]
		//Step 3: Calculate blocks of H_0
			for ( int k = 0; k < q; k++ ) {
				for ( unsigned int i = 0; i < items; i++ ) {
					H_0[9 * items * k + 9 * i + 3 * 0 + 0] = (D*D) * (theta[k]*theta[k]) * ecPlus1Inv[i] * (1-2*P_Star[k * items + i] );
					H_0[9 * items * k + 9 * i + 3 * 0 + 1] = H_0[9 * items * k + 9 * i + 3 * 1 + 0] = (D*D) * theta[k] * ecPlus1Inv[i] * (1-2*P_Star[k * items + i] );
					H_0[9 * items * k + 9 * i + 3 * 1 + 1] = (D*D) * ecPlus1Inv[i] * (1-2*P_Star[k * items + i] );
					H_0[9 * items * k + 9 * i + 3 * 0 + 2] = H_0[9 * items * k + 9 * i + 3 * 2 + 0] = - ec[i] * D * theta[k] * (ecPlus1Inv[i]*ecPlus1Inv[i]);
					H_0[9 * items * k + 9 * i + 3 * 1 + 2] = H_0[9 * items * k + 9 * i + 3 * 2 + 1] = - ec[i] * D * (ecPlus1Inv[i]*ecPlus1Inv[i]);
					H_0[9 * items * k + 9 * i + 3 * 2 + 2] = - ( ec[i] * (-1+ec[i]) * (ecPlus1Inv[i]*ecPlus1Inv[i]*ecPlus1Inv[i]) ) / P_Star[k * items + i] ;
				}
			}

			for ( unsigned int i = 0; i < items; i++ ) {

				// Fill block i of zeros
				for ( int m = 0; m < 3; m++ ) {
					for ( int n = 0; n < 3; n++ ) {
						hessian[9 * i + 3 * m + n]= 0.0;
					} // n
				} // m

				for ( int k = 0; k < q; k++ ) {
					for ( int m = 0; m < 3; m++ ) {
						// Add k-block values to block i
						for ( int n = 0; n < 3; n++ ) {
							//Begins the hessian multiplication
							hessian[9 * i + 3 * m + n] +=
									(
									factor[k * items + i]*H_0[9 * items * k + 9 * i + 3 * m + n]
									- ( r[k * items + i] - 2*P[k * items + i]*r[k * items + i]
									+ f[k]*P[k * items + i]*P[k * items + i] )
									* (W[k * items + i]*W[k * items + i] )*
									(h_0[3 * items * k + 3 * i + m]	* h_0[3 * items * k + 3 * i + n] )
									);
						} // n
					} // m
				} // k

				// Debugging printing
				/*
								for ( int m = 0; m < 3; m++ ) {
									for ( int n = 0; n < 3; n++ ) {
										cout<<9 * i + 3 * m + n<<":"<<hessian[9 * i + 3 * m + n]<<" ";
									} // n
									cout<<endl;
								} // m
								cout<<endl;
								*/

			} // i

		delete [] H_0;

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
}

void ThreePLModel::Ngradient(double* args, double* pars, int nargs, int npars, double* gradient){
	//	For each of the gradient thingies increase the args and apply richardsons
	double hh = 0.000001;
	//(f(x+h)-f(x))/h
	for(int i = 0 ; i < nargs; i++){
		args[i]=args[i]+hh;
		gradient[i]=logLikelihood(args,pars,nargs,npars);
		args[i]=args[i]-hh;
		gradient[i]-=logLikelihood(args,pars,nargs,npars);
		gradient[i]=gradient[i]/hh;
	}
}

void ThreePLModel::NHessian(double*args, double* pars, int nargs, int npars, double* hessian){
	//the gradient is composed of items a's items b's items c's take it as such
	double hh = 0.001;
	double gradiente[nargs];
	int items = nargs/3;
	gradient(args,pars,nargs,npars,gradiente);
	//Calculate the hessian for each item order is a's b's c's in args
	for (int i = 0; i < items; ++i) {
				//Extract the gradient at the points (Put at the hessian)
		for (int j = 0; j < 3; ++j) {
			for (int k = 0; k < 3; ++k) {
				//Now for each of these point change the k  parameter for deriving
				hessian[i*items+j*3+k] = gradiente[j*items+i];
			}
		}
	}

	for (int i = 0; i < items; ++i) {
					//Derivate two times for this item // matrix point
					//Extract the gradient at the points (Put at the hessian)
			for (int j = 0; j < 3; ++j) {
				for (int k = 0; k < 3; ++k) {
					memset(gradiente,0,sizeof(double)*items*3);
					//Change the argument k in the item i
					args[k*items+i]+=hh;
					//Calculate le gradient
					gradient(args,pars,nargs,npars,gradiente);
					//Reset the parameter
					args[k*items+i]-=hh;
					//This is the gradient at the point
					hessian[i*items+j*3+k]-= gradiente[j*items+i];
					hessian[i*items+j*3+k]=hessian[i*items+j*3+k]/hh;
					cout<<i<<" "<<j<<k<<"   ";
				}cout<<endl;
			}cout<<endl;
		}

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

	int q, I;
	double *theta, *r, *f, *a, *b, *c;


	// Obtain q
	q = pars[nP ++]; // q is obtained and npars is augmented

	// Obtain I
	I = pars[nP ++];

	theta = new double[q];
	r = new double[q*I];
	f = new double[q];
	a = new double[I];
	b = new double[I];
	c = new double[I];


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
		//cout<<" "<<c[i];
	}//cout<<endl;

	long double tp , tq;
	long double sum = 0;

	for (int k = 0; k < q; ++k) {
		for (unsigned int i = 0; i < I; ++i) {
			tp = (ThreePLModel::successProbability ( theta[k], a[i], b[i], c[i]));
			if (tp==0)tp=1e-08;
			tq = 1-tp;
			if (tq==0)tq=1e-08;
			//suma = suma + (rki*logg(pki)+(fki-rki)*logg(qki))
			sum+=(r[k * I + i]*log(tp))+(f[k]-r[k * I + i])*log(tq);
		}
	}


	//antiLogit(c, I);
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

