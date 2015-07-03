/*
 * ThreePLModel.cpp
 *
 *  Created on: 17 Jun 2014
 *      Author: jlgpisa
 */

#include <model/parameter/ThreePLModel.h>
#include <string>

ThreePLModel::ThreePLModel()
{
	parameterSet = NULL;
	probabilityMatrix = NULL;
	nodes = NULL;
} 

void ThreePLModel::setEstimationNodes(QuadratureNodes* n) { this->nodes = n; }

void ThreePLModel::successProbability(DimensionModel *dimensionModel, QuadratureNodes * quadNodes)
{
	int q = 0;
	double a_d, d_d, c_d, theta_d; // d stands from "double"

	if ( dimensionModel != NULL )
		q = quadNodes->size();

	if(typeid(*dimensionModel)==typeid(UnidimensionalModel))
	{
		if(probabilityMatrix == NULL)
			//Creates the matrix if it is not already created
			probabilityMatrix = new Matrix<double>(q,items);

		for (int k = 0; k < q; k++)
		{
			for ( int i = 0; i < items; i++ )
			{
				// 3PL Success Probability Function
				theta_d = (*quadNodes->getTheta())(0,k);
				a_d = parameterSet[0][0][i];
				d_d = parameterSet[1][0][i];
				c_d = parameterSet[2][0][i];

				(*probabilityMatrix)(k,i) = successProbability ( theta_d, a_d, d_d, c_d );
			}
		}
	}
}

double *** ThreePLModel::getParameterSet() { return (this->parameterSet); }

void ThreePLModel::setParameterSet(double ***) { this->parameterSet = parameterSet; }

double ThreePLModel::successProbability(double theta, double a, double d, double c)
{
	long double exponential = (Constant::NORM_CONST)*(a*theta+d);

	if ( exponential > Constant::MAX_EXP )
		exponential = Constant::MAX_EXP;

	else if ( exponential < -(Constant::MAX_EXP*1.0) )
		exponential = -Constant::MAX_EXP;

	exponential = exp(-exponential) ;
	double ec = exp(c);

	return ( (ec/(1+ec)) + (1 - (ec/(1+ec))) * (1/(1+exponential)) );
}

void ThreePLModel::getParameters(double * parameters)
{
	int i = 0;

	for (int j = 0; j < items; j++)
		parameters[i++] = parameterSet[0][0][j];
	for (int j = 0; j < items; j++)
		parameters[i++] = parameterSet[1][0][j];
	for (int j = 0; j < items; j++)
		parameters[i++] = parameterSet[2][0][j];
}

void ThreePLModel::setParameters(double * parameters)
{
	int i = 0;

	for (int j = 0; j < items; j++)
		this->parameterSet[0][0][j] = parameters[i++];
	for (int j = 0; j < items; j++)
		this->parameterSet[1][0][j] = parameters[i++];
	for (int j = 0; j < items; j++)
		this->parameterSet[2][0][j] = parameters[i++];
}

double ThreePLModel::successProbability(double theta, double * zita)
{
	long double exponential = (Constant::NORM_CONST)*(zita[0]*theta+zita[1]);

	if ( exponential > Constant::MAX_EXP )
		exponential = Constant::MAX_EXP;

	else if ( exponential < -(Constant::MAX_EXP*1.0) )
		exponential = -Constant::MAX_EXP;

	exponential = exp(-exponential) ;
	double ec = exp(zita[2]);

	return ( (ec/(1+ec)) + (1 - (ec/(1+ec))) * (1/(1+exponential)) );
}

double ThreePLModel::getProbability(int node, int item) { return ((*probabilityMatrix)(node, item)); }

void ThreePLModel::NitemGradient(double* args, double* pars, int nargs, int npars, double* gradient)
{
	itemGradient(args,pars,nargs,npars,gradient);

 	double hh = 1e-08;
 	int items = pars[1];
 	int* idx = new int[3];
 	int index = pars[npars-1];
 	int nA = 0;
 	
 	nA = 0;
 	nA += index;
	idx[0] = nA;
	nA += (items-index);
	nA += index;
	idx[1] = nA;
	nA += (items-index);
	nA += index;
	idx[2] = nA;
 	
 	for(int i = 0 ; i < nargs; i++)
 	{
 		//Add hh to argument
 		args[idx[i]] = args[idx[i]] + hh;
 		//Evaluate
 		gradient[i] = itemLogLik(args,pars,nargs,npars);
 		//Remove hh to argument and evaluate
 		args[idx[i]] = args[idx[i]] - hh;
 		gradient[i] -= itemLogLik(args,pars,nargs,npars);
 		//Substract
 		gradient[i] = gradient[i]/hh;
 	}

 	delete[] idx;
}

void ThreePLModel::itemGradient (double* args, double* pars, int nargs, int npars, double* gradient)
{
	int nA = 0;
	int nP = 0;
	int q, items;
	double *theta, *r, *f;
	double a, b, c;
	double D = Constant::NORM_CONST;
	long double *h_0; // Block Matrix of size q*I. Each block-element has size of 1*3
	long double *h; // Block vector of size I (i.e. I blocks). Each block-element has size of 1*3
	long double *P_Star, *P;  // Matrix of size q*I
	long double *W;           // Matrix of size q*I
	long double *factor;	  // Matrix of product (r-fP)W
	long double ec;            // e^c_i
	long double ecp1i;	// 1 / (e^c_i + 1)
	int index = 0;

	index = pars[npars-1];
	
	// Obtain q
	q = pars[nP ++]; // q is obtained and npars is augmented
	
	// Obtain I
	items = pars[nP ++];
	theta = new double[q];
	r = new double[q];
	f = new double[q];
	
	// Obtain theta
	for (int k=0; k<q; k++)
		theta[k] = pars[nP ++];
	
	// Obtain f
	for (int k=0; k<q; k++)
		f[k] = pars[nP ++];
	
	// Obtain r that becomes a vector
	for (int k=0; k<q; k++)
	{
		nP += index;
		r[k] = pars[nP];
		nP += (items-index); 
	}
	
	// Obtain a, b and c
	nA += index;
	a = args[nA];
	nA += (items-index);
	nA += index;
	b = args[nA];
	nA += (items-index);
	nA += index;
	c = args[nA];

	h = new long double [3];
	h_0 = new long double [q*3];
	P = new long double [q];
	P_Star = new long double [q];
	factor = new long double [q];
	W = new long double [q];

	ecp1i=1/(1+exp(c));
	ec=exp(c);
	
	for ( int k = 0; k < q; k++ )
	{
		P[k] = successProbability ( theta[k], a,b,c);
		P_Star[k] = 1/(1+exp(-D*(a*theta[k]+b)));

		W[k] = P_Star[k] * ( 1 - P_Star[k] ); // Numerator
		W[k] /= P[k] * ( 1 - P[k] );// Denominator

		factor[k] = ( r[k] - f[k]*P[k] ) * W[k];

		// h_0 / (P_star*Q_star)
		h_0[3 * k ] = D * theta[k] * ecp1i;
		h_0[3 * k + 1] = D * ecp1i;
		h_0[3 * k + 2] = ec * (ecp1i*ecp1i) / P_Star[k];
	}

	memset(h,0,sizeof(long double)*3);
	memset(gradient,0,sizeof(double)*3);

	for ( int k = 0; k < q; k++ )
	{
			h[0] += factor[k] * h_0[3 * k + 0];
			h[1] += factor[k] * h_0[3 * k + 1];
			h[2] += factor[k] * h_0[3 * k + 2];
	}

	delete [] h_0;
	delete [] P_Star;
	delete [] P;
	delete [] W;
	delete [] factor;

	delete [] theta;
	delete [] r;
	delete [] f;

	//return h as the gradient
	int hc=0;
	for (int n = 0; n < 3; ++n)
		gradient[hc++]= -static_cast<double>(h[n]);

	delete [] h;
}

void ThreePLModel::itemGradient2 (double* args, double* pars, int nargs, int npars, double* gradient)
{
	int nA = 0;
	int nP = 0;
	int q, items;
	double a, b, c;
	int index = 0;
	double D = Constant::NORM_CONST;
	double *theta, *r, *f;
	long double *h_0; // Block Matrix of size q*I. Each block-element has size of 1*3
	long double *h; // Block vector of size I (i.e. I blocks). Each block-element has size of 1*3
	long double *P_Star, *P;  // Matrix of size q*I
	long double *W;           // Matrix of size q*I
	long double *factor;	  // Matrix of product (r-fP)W
	long double ec;            // e^c_i
	long double ecp1i;	// 1 / (e^c_i + 1)

	index = pars[npars-1];
	
	// Obtain q
	q = pars[nP ++]; // q is obtained and npars is augmented
	
	// Obtain I
	items = pars[nP ++];
	theta = new double[q];
	r = new double[q];
	f = new double[q];
	
	// Obtain theta
	for (int k=0; k<q; k++)
		theta[k] = pars[nP ++];

	// Obtain f
	for (int k=0; k<q; k++)
		f[k] = pars[nP ++];

	// Obtain r that becomes a vector
	for (int k=0; k<q; k++)
	{
		nP += index;
		r[k] = pars[nP];
		nP += (items-index); 
	}

	// Obtain a
	nA += index;
	a = args[nA];
	nA += (items-index);
	nA += index;
	// Obtain b
	b = args[nA];
	nA += (items-index);
	nA += index;
	// Obtain c
	c = args[nA];

	a = args[0];
	b = args[1];
	c = args[2];

	h = new long double [3];
	h_0 = new long double [q*3];
	P = new long double [q];
	P_Star = new long double [q];
	factor = new long double [q];
	W = new long double [q];

	ecp1i=1/(1+exp(c));
	ec=exp(c);
	
	for ( int k = 0; k < q; k++ )
	{
		P[k] = successProbability ( theta[k], a,b,c);
		P_Star[k] = 1/(1+exp(-D*(a*theta[k]+b)));

		W[k] = P_Star[k] * ( 1 - P_Star[k] ); // Numerator
		W[k] /= P[k] * ( 1 - P[k] );// Denominator

		factor[k] = ( r[k] - f[k]*P[k] ) * W[k];

		// h_0 / (P_star*Q_star)
		h_0[3 * k ] = D * theta[k] * ecp1i;
		h_0[3 * k + 1] = D * ecp1i;
		h_0[3 * k + 2] = ec * (ecp1i*ecp1i) / P_Star[k];
	}

	memset(h,0,sizeof(long double)*3);
	memset(gradient,0,sizeof(double)*3);

	for ( int k = 0; k < q; k++ )
	{
			h[0] += factor[k] * h_0[3 * k + 0];
			h[1] += factor[k] * h_0[3 * k + 1];
			h[2] += factor[k] * h_0[3 * k + 2];
	}

	delete [] h_0;
	delete [] P_Star;
	delete [] P;
	delete [] W;
	delete [] factor;

	delete [] theta;
	delete [] r;
	delete [] f;

	//return h as the gradient
	int hc=0;

	for (int n = 0; n < 3; ++n)
			gradient[hc++]= -static_cast<double>(h[n]);

	delete [] h;
}

void ThreePLModel::gradient (double* args, double* pars, int nargs, int npars, double* gradient)
{
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
	unsigned int q, items;
	double D = Constant::NORM_CONST;
	double *theta, *r, *f, *a, *b, *c;
	long double *h_0; // Block Matrix of size q*I. Each block-element has size of 1*3
	long double *h; // Block vector of size I (i.e. I blocks). Each block-element has size of 1*3
	long double *P_Star, *P;  // Matrix of size q*I
	long double *W;           // Matrix of size q*I
	long double *factor;	  // Matrix of product (r-fP)W
	long double *ec;            // e^c_i
	long double *ecPlus1Inv;	// 1 / (e^c_i + 1)

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
	for (unsigned int k = 0; k < q; k++)
		theta[k] = pars[nP ++];

	// Obtain f
	for (unsigned int k = 0; k < q; k++)
		f[k] = pars[nP ++];

	// Obtain r
	for (unsigned int k = 0; k < q; k++)
		for (unsigned int i = 0; i < items; i++)
			r[k*items+i] = pars[nP ++];

	// Obtain a
	for (unsigned int i = 0; i < items; i++)
		a[i] = args [nA ++];

	// Obtain b
	for (unsigned int i = 0; i < items; i++)
		b[i] = args [nA ++];

	// Obtain c
	for (unsigned int i = 0; i < items; i++)
		c[i] = args [nA ++];

	h = new long double [3*items];
	h_0 = new long double [q*3*items];
	P = new long double [q*items];
	P_Star = new long double [q*items];
	factor = new long double [q*items];
	W = new long double [q*items];
	ec = new long double [items];
	ecPlus1Inv = new long double [items];

	for(unsigned  int i = 0; i < items; i++)
	{
		ecPlus1Inv[i]=1/(1+exp(c[i]));
		ec[i]=exp(c[i]);
	}

	for (unsigned int k = 0; k < q; k++)
	{
		for (unsigned  int i = 0; i < items; i++)
		{
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
	
	for (unsigned int i = 0; i < items; i++)
	{
		for (unsigned int k = 0; k < q; k++)
		{
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

	int hc=0;

	for (int n = 0; n < 3; ++n)
		for(unsigned int i = 0 ; i < items ; ++i)
			gradient[hc++]= -static_cast<double>(h[i*3+n]);

	delete [] h;
}

double ThreePLModel::itemLogLik (double* args, double* pars, int nargs, int npars)
{
	int nA = 0;
	int nP = 0;
	int q, items;
	int index = 0;
	double *theta, *r, *f;
	double a, b, c;
	double sum=0;
	long double tp , tq;

	index = pars[npars - 1];
	// Obtain q
	q = pars[nP ++]; // q is obtained and npars is augmented
	// Obtain I
	items = pars[nP ++];
	theta = new double[q];
	r = new double[q];
	f = new double[q];
	
	// Obtain theta
	for (int k=0; k<q; k++)
		theta[k] = pars[nP ++];

	// Obtain f
	for (int k=0; k<q; k++)
		f[k] = pars[nP ++];

	// Obtain r that becomes a vector
	for (int k=0; k<q; k++)
	{
		nP += index;
		r[k] = pars[nP];
		nP += (items-index); 
	}

	// Obtain a
	nA += index;
	a = args[nA];
	nA += (items-index);
	nA += index;
	// Obtain b
	b = args[nA];
	nA += (items-index);
	nA += index;
	// Obtain c
	c = args[nA];

	for (int k = 0; k < q; ++k)
	{
		tp = (ThreePLModel::successProbability (theta[k], a,b,c));

		if (tp<1e-08) tp=1e-08;

		tq = 1-tp;
		
		if (tq<1e-08) tq=1e-08;

		sum+=(r[k]*log(tp))+(f[k]-r[k])*log(tq);
	}

	delete[] theta;
	delete[] f;
	delete[] r;

	return (-sum);
}

double ThreePLModel::itemLogLik2 (double* args, double* pars, int nargs, int npars)
{
	int nP = 0;
	int q, items;
	int index = 0;
	double *theta, *r, *f;
	double a, b, c;
	double sum=0;
	long double tp , tq;
	
	index = pars[npars-1];
	
	// Obtain q
	q = pars[nP ++]; // q is obtained and npars is augmented
	
	// Obtain I
	items = pars[nP ++];
	theta = new double[q];
	r = new double[q];
	f = new double[q];
	
	// Obtain theta
	for (int k=0; k<q; k++)
		theta[k] = pars[nP ++];
	
	// Obtain f
	for (int k=0; k<q; k++)
		f[k] = pars[nP ++];

	// Obtain r that becomes a vector
	for (int k=0; k<q; k++)
	{
		nP += index;
		r[k] = pars[nP];
		nP += (items-index); 
	}
	
	// Obtain a
	a = args[0];
	// Obtain b
	b = args[1];
	// Obtain c
	c = args[2];

	if(abs(a)>5)
		a = 0.851;
	double dd = 0;
	dd = -b/a;
	if(abs(dd)>5)
		b = 0;
	if(abs(c)>5)
		c = 0.3;

	for (int k = 0; k < q; ++k)
	{
		tp = (ThreePLModel::successProbability ( theta[k], a,b,c));
		
		if (tp<1e-08) tp=1e-08;
		
		tq = 1-tp;
		
		if (tq<1e-08) tq=1e-08;
		
		sum+=(r[k]*log(tp))+(f[k]-r[k])*log(tq);
	}

	args[0] = a;
	args[1] = b;
	args[2] = c;

	delete[] theta;
	delete[] f;
	delete[] r;

	return (-sum);
}

double ThreePLModel::logLikelihood (double* args, double* pars, int nargs, int npars)
{
	int nA = 0;
	int nP = 0;
	unsigned int q, It;
	double *theta, *r, *f, *a, *b, *c;
	long double tp , tq;
	long double sum = 0;


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
	for (unsigned int k=0; k<q; k++)
		theta[k] = pars[nP ++];

	// Obtain f
	for (unsigned int k=0; k<q; k++)
		f[k] = pars[nP ++];

	// Obtain r
	for (unsigned int k=0; k<q; k++)
		for (unsigned int i=0; i<It; i++)
			r[k*It+i] = pars[nP ++];

	// Obtain a
	for (unsigned int i=0; i<It; i++)
		a[i] = args [nA ++];

	// Obtain b
	for (unsigned int i=0; i<It; i++)
		b[i] = args [nA ++];

	// Obtain c
	for (unsigned int i=0; i<It; i++) 
		c[i] = args [nA ++];

	for (unsigned int k = 0; k < q; ++k)
	{
		for (unsigned int i = 0; i < It; ++i)
		{
			tp = (ThreePLModel::successProbability ( theta[k], a[i], b[i], c[i]));
			
			if (tp<1e-08) tp=1e-08;

			tq = 1-tp;

			if (tq<1e-08) tq=1e-08;

			sum += (r[k * It + i]*log(tp))+(f[k]-r[k * It + i])*log(tq);
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

double ThreePLModel::successProbability_cPrime(double theta, double a, double b, double c)
{
	long double cPrime = exp(c)/(1+exp(c));
	return (successProbability ( theta, a, b, cPrime ));
}

void ThreePLModel::printParameterSet(ostream& out)
{
	out << "\"a\" \"b\" \"c\"" << endl;

	for (int i = 0; i < items; i++)
		out << parameterSet[0][0][i] << " "
		    << parameterSet[1][0][i] << " "
		    << parameterSet[2][0][i] << endl;
}


double ThreePLModel::banana(double* args, double* pars, int nargs, int npars)
{
	double x1 = args[0];
	double x2 = args[1];

	return (100 * (x2 - x1 * x1) * (x2 - x1 * x1) + (1 - x1) * (1 - x1));
}

double ThreePLModel::Loglik(double* xzita,double* xRmat,double* xf, double* xptcuad, double itemn)
{
	/*
	NumericVector xzita(zita);
	NumericMatrix xRmat(Rmat);
	NumericVector xf(f);
	NumericVector xptcuad(ptcuad);
	*/
	//Printing of the arguments
	const int items = itemn;
	const int numCuads = 41;

	double suma = 0;
	for(int k = 0 ; k < numCuads ; k++)
	{
	    for(int i = 0 ; i < items ; i++)
	    {
	        double pki = Pr(1, xzita[i], xzita[items + i], xzita[2 * items + i], xptcuad[k]);
	        double qki = 1.0 - pki;
	        suma = suma + (xRmat[k*numCuads+i] * log(pki) + (xf[k] - xRmat[k*numCuads+i] * log(qki)));
		}
	}

	return (-suma);
}

double* ThreePLModel::grad(double* xzita,double* xRmat,double* xf, double* xptcuad, double nitems)
{  
	int items = (int)nitems;
	const int numCuads = 41;
	double* h = new double[3];
	double* grad = new double[items * 3];
	Matrix<double>* g = new Matrix<double>(3,items);
	bool *llamadoNear = new bool;

	*llamadoNear = 0;

	for(int i = 0 ; i < items ; i++)
	{
		double *suma = new double[3];
		//Inicializa suma a cero
		for(int j = 0 ; j < 3 ; j++)
			suma[j] = 0;
		for(int k = 0 ; k < numCuads ; k++)
		{
			double p = Pr(1,xzita[i], xzita[items + i] ,xzita[2 * items + i], xptcuad[k]);
			double pm = Prm(1,xzita[i] , xzita[items + i], xptcuad[k]);
			double aux = (xRmat[k*numCuads+i] -  xf[k] * p) * ((pm * (1 - pm)) / (p * (1 - p)));
			
			h[0] = xptcuad[k]*(1/(1+exp(xzita[2 * items + i])));
			h[1] = (1/(1+exp(xzita[2 * items + i])));
			h[2] = pow((1/(1+exp(xzita[2 * items + i]))),2) * exp(xzita[2 * items + i]) / pm;
			
			for(int j = 0 ; j < 3 ; j++)
				suma[j] = suma[j] + (aux * h[j]);
		}

		for(int j = 0 ; j < 3 ; j++)
			(*g)(j,i) = suma[j]; 

		delete[] suma;
	}
	
	for(int j = 0 ; j < 3 ; j++)
		for(int i = 0 ; i < items ; i++)
			grad[j*items + i] = -(*g)(j,i); 

	return grad;
}

double ThreePLModel::Prm(int u,double a, double d, double theta)
{
	double p = 0;
	double z = a * theta + d;
	
	if(abs(z) > 35)
		z = abs(z) /  z * 35; 

	p =  (1/(1 + exp(-1*(z))));

	if(0.0 == p)
		p = sqrt(2.2e-16);
	else if(1.0 == p)
		p = 1 - 1e-6;

	if(1 == u)
		return(p);
	else
		return(1-p);
}

double ThreePLModel::Pr(int u,double a, double d, double c, double theta)
{
	double p = 0;
	double z = a * theta + d;

	if(abs(z) > 35)
		z = abs(z) /  z * 35; 

	p = exp(c) / (1+exp(c)) + (1-(exp(c) / (1+exp(c)))) * (1/(1 + exp(-1*(z))));

	if(0.0 == p)
		p = sqrt(2.2e-16);
	else if(1.0 == p)
		p = 1 - 1e-6;

	if(1 == u)
		return(p);
	else
		return(1-p);
}

ThreePLModel::~ThreePLModel()
{
	if(nodes != NULL)
	{
		delete nodes;
		nodes = NULL;
	}
	if(profiler != NULL)
	{
		delete profiler;
		profiler = NULL;
	}
	if(probabilityMatrix != NULL)
	{
		delete probabilityMatrix;
		probabilityMatrix = NULL;
	}
	if (parameterSet != NULL)
	{
		delete[] parameterSet[2][0];
		delete[] parameterSet[1][0];
		delete[] parameterSet[0][0];

		delete[] parameterSet[0];
		delete[] parameterSet[1];
		delete[] parameterSet[2];

		delete[] parameterSet;
		parameterSet = NULL;
	}
}
