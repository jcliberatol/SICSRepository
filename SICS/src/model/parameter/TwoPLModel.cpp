/*
 * TwoPLModel.cpp
 *
 *  Created on: 18 Jun 2014
 *      Author: cesandovalp
 */

#include <model/parameter/TwoPLModel.h>

TwoPLModel::TwoPLModel()
{
	parameterSet = NULL;
	probabilityMatrix = NULL;
}

void TwoPLModel::untransform()
{
    for (  int i = 0; i < itemModel->getDataset()->countItems(); ++i)
        parameterSet[1][0][i] /= -parameterSet[0][0][i];
}

void TwoPLModel::getParameters(double * parameters)
{
      int i = 0;

    for (  int j = 0; j < items; j++)
        parameters[i++] = parameterSet[0][0][j];
    for (  int j = 0; j < items; j++)
        parameters[i++] = parameterSet[1][0][j];
}

void TwoPLModel::setParameters(double * parameters)
{
      int i = 0;

    for (  int j = 0; j < items; j++)
        this->parameterSet[0][0][j] = parameters[i++];
    for (  int j = 0; j < items; j++)
        this->parameterSet[1][0][j] = parameters[i++];
}

inline void TwoPLModel::successProbability(DimensionModel *dimensionModel, QuadratureNodes *quadNodes)
{
	  int q = 0;
	double a_d, d_d, theta_d; // d stands from "double"

	if (dimensionModel != NULL)
		q = quadNodes->size();

	if (typeid(*dimensionModel) == typeid(UnidimensionalModel))
	{
		if (probabilityMatrix == NULL)
			//Creates the matrix if it is not already created
			probabilityMatrix = new Matrix<double>(q, items);
		
		for (  int k = 0; k < q; k++)
		{
			for (  int i = 0; i < items; i++)
			{
				theta_d = (*quadNodes->getTheta())(0, k);
				a_d = parameterSet[0][0][i];
				d_d = parameterSet[1][0][i];
				(*probabilityMatrix)(k, i) = successProbability(theta_d, a_d, d_d);
			}
		}
	}
}

inline double TwoPLModel::successProbability(double theta, double a, double d)
{
	long double exponential = (Constant::D_CONST * ((a * theta) + d));

	if (exponential > Constant::MAX_EXP)
		exponential = Constant::MAX_EXP;
	else if (exponential < -(Constant::MAX_EXP))
		exponential = -Constant::MAX_EXP;

	return (1 / (1 + exp(-exponential)));
}

inline double TwoPLModel::successProbability(double theta, double * zita) { return successProbability(theta, zita[0], zita[1]); }

double *** TwoPLModel::getParameterSet() { return (this->parameterSet); }

void TwoPLModel::setParameterSet(double *** pair) { this->parameterSet = pair; }

double TwoPLModel::getProbability(int node, int item) { return ((*probabilityMatrix)(node, item)); }

void TwoPLModel::itemGradient (double* args, double* pars, int nargs, int npars, double* gradient)
{
	int nP, q, items, index, hc;
	double *theta, *r, *f;
	double a, b, D;
	long double *h_0; // Block Matrix of size q*I. Each block-element has size of 1*2
	long double *h; // Block vector of size I (i.e. I blocks). Each block-element has size of 1*2
	long double *P;  // Matrix of size q*I
	long double *factor;	  // Matrix of product (r-fP)W

	index = pars[npars-1];
	hc = nP = 0;
	D = Constant::NORM_CONST;
	a = args[0];
	b = args[1];
        
	q = pars[nP ++];
	items = pars[nP ++];
	theta = new double[q];
	r = new double[q];
	f = new double[q];
	h = new long double [2];
	h_0 = new long double [q*2];
	P = new long double [q];
	factor = new long double [q];
	
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

	for ( int k = 0; k < q; k++ )
	{
		P[k] = successProbability ( theta[k], a,b);
		factor[k] = ( r[k] - f[k]*P[k] );
		h_0[2 * k ] = D * theta[k];
		h_0[2 * k + 1] = D;
	}

	memset(h,0,sizeof(long double)*2);
	memset(gradient,0,sizeof(double)*2);

	for ( int k = 0; k < q; k++ )
	{
		h[0] += factor[k] * h_0[2 * k + 0];
		h[1] += factor[k] * h_0[2 * k + 1];
	}

	delete [] h_0;
	delete [] P;
	delete [] factor;

	delete [] theta;
	delete [] r;
	delete [] f;

	//return h as the gradient
	for (int n = 0; n < 2; ++n) 
		gradient[hc++]= -static_cast<double>(h[n]);

	delete [] h;
}

double TwoPLModel::itemLogLik (double* args, double* pars, int nargs, int npars)
{
	int nP = 0;
	  int q, items;
	int index = 0;
	double *theta, *r, *f;
	double a, b;
	double sum=0;
	long double tp , tq;
	
	index = pars[npars-1];
	a = args[0];
	b = args[1];

	// Obtain q
	q = pars[nP ++]; // q is obtained and npars is augmented

	// Obtain I
	items = pars[nP ++];
	theta = new double[q];
	r = new double[q];
	f = new double[q];

	// Obtain theta
	for (  int k=0; k<q; k++)
		theta[k] = pars[nP ++];

	// Obtain f
	for (  int k=0; k<q; k++)
		f[k] = pars[nP ++];

	// Obtain r that becomes a vector
	for (  int k=0; k<q; k++)
	{
		nP += index;
		r[k] = pars[nP];
		nP += (items-index);
	}
	
	if(abs(a)>5)
		a = 0.851;

	double dd = 0;
	dd = -b/a;
	
	if(abs(dd)>5)
		b = 0;

	for (  int k = 0; k < q; ++k)
	{
		tp = (TwoPLModel::successProbability ( theta[k], a,b));
		
		if (tp < 1e-08) tp = 1e-08;

		tq = 1-tp;

		if (tq < 1e-08) tq = 1e-08;

		sum += (r[k]*log(tp))+(f[k]-r[k])*log(tq);
	}

	args[0] = a;
	args[1] = b;

	delete[] theta;
	delete[] f;
	delete[] r;

	return (-sum);
}

void TwoPLModel::printParameterSet(ostream& out)
{
	cout << "\"a\" \"b\" \"c\"" << endl;
	
	for (  int i = 0; i < items; i++)
		cout << parameterSet[0][0][i] << " "
		     << parameterSet[1][0][i] << " "
		     << 0 << endl;
}

