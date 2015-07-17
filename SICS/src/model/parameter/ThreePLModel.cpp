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

void ThreePLModel::transform()
{
	for (int i = 0; i < itemModel->countItems(); ++i)
	{
		double qc = parameterSet[2][0][i];
		parameterSet[2][0][i] = log(qc / (1 - qc));
	}
}

void ThreePLModel::untransform()
{
	for (int i = 0; i < itemModel->getDataset()->countItems(); ++i)
	{
		double qa = parameterSet[0][0][i];
		double qb = parameterSet[1][0][i];
		double qc = parameterSet[2][0][i];
		double ec = exp(qc);
		parameterSet[2][0][i] = ec / (1 + ec);
		parameterSet[1][0][i] = -qb / qa; //Transformacion del B   d=-b/a
	}
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

double ThreePLModel::successProbability(double theta, double * zita) { return successProbability(theta, zita[0], zita[1], zita[2]); }

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

double ThreePLModel::getProbability(int node, int item) { return ((*probabilityMatrix)(node, item)); }

void ThreePLModel::itemGradient(double* args, double* pars, int nargs, int npars, double* gradient)
{
	int nP, q, items, index;
	double a, b, c;
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
	nP = 0;
	a = args[0];
	b = args[1];
	c = args[2];

	q = pars[nP ++];
	items = pars[nP ++];
	theta = new double[q];
	r = new double[q];
	f = new double[q];
	h = new long double [3];
	h_0 = new long double [q*3];
	P = new long double [q];
	P_Star = new long double [q];
	factor = new long double [q];
	W = new long double [q];
	
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

double ThreePLModel::itemLogLik(double* args, double* pars, int nargs, int npars)
{
	double *theta, *r, *f;
	unsigned int nP, q, items, index;
	double a, b, c, sum;
	//long double tp , tq;
	double tp , tq;
	sum = nP = index = 0;
	
	a = args[0];
	b = args[1];
	c = args[2];

        q = pars[nP ++]; // q is obtained and npars is augmented
        items = pars[nP ++];
	index = pars[npars-1];
	
	theta = new double[q];
	r = new double[q];
	f = new double[q];
	
	// Obtain theta
	for (unsigned int k=0; k<q; k++)
		theta[k] = pars[nP ++];
	
	// Obtain f
	for (unsigned int k=0; k<q; k++)
		f[k] = pars[nP ++];

	// Obtain r that becomes a vector
	for (unsigned int k=0; k<q; k++)
	{
		nP += index;
		r[k] = pars[nP];
		nP += (items-index); 
	}

	if(abs(a) > 5)
		a = 0.851;
	double dd = -b/a;
	if(abs(dd)>5)
		b = 0;
	if(abs(c)>5)
		c = 0.1;

	for (unsigned int k = 0; k < q; ++k)
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
