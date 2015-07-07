/*
 * OnePLModel.cpp
 *
 *  Created on: Nov 16, 2014
 *      Author: anmrodriguezre
 */

#include <model/parameter/OnePLModel.h>

OnePLModel::OnePLModel()
{
	parameterSet = NULL;
	probabilityMatrix = NULL;
}

inline void OnePLModel::successProbability(DimensionModel *dimensionModel, QuadratureNodes * quadNodes)
{
	int q = 0;

	if (dimensionModel != NULL)
		q = quadNodes->size();

	if (typeid(*dimensionModel) == typeid(UnidimensionalModel))
	{
		if (probabilityMatrix == NULL)
			//Creates the matrix if it is not already created
			probabilityMatrix = new Matrix<double>(q, items);

		for (int k = 0; k < q; k++)
			for (int i = 0; i < items; i++)
				// Rasch Success Probability Function
				(*probabilityMatrix)(k, i) =
					successProbability((*quadNodes->getTheta())(0, k), (parameterSet[0][0][i]));
	}
}

inline double OnePLModel::successProbability(double theta, double b)
{
	long double exponential = ((theta) - b);

	if (exponential > Constant::MAX_EXP)
		exponential = Constant::MAX_EXP;
	else if (exponential < -(Constant::MAX_EXP))
		exponential = -Constant::MAX_EXP;

	return (1 / (1.0 + exp(-exponential)));
}

inline double OnePLModel::successProbability(double theta, double * zita)
{
	long double exponential = ((theta) - zita[0]);

	if (exponential > Constant::MAX_EXP)
		exponential = Constant::MAX_EXP;
	else if (exponential < -(Constant::MAX_EXP))
		exponential = -Constant::MAX_EXP;

	return (1 / (1.0 + exp(-exponential)));
}

void OnePLModel::setParameterSet(double*** par) { this->parameterSet = par; }

double*** OnePLModel::getParameterSet() { return (this->parameterSet); }

void OnePLModel::getParameters(double * parameters)
{
	for (int i = 0; i < items; i++)
		parameters[i] = parameterSet[0][0][i];
}

void OnePLModel::setParameters(double * parameters)
{
	for (int i = 0; i < items; i++)
		this->parameterSet[0][0][i] = parameters[i];
}

double OnePLModel::getProbability(int node, int item) { return ((*probabilityMatrix)(node, item)); }

void OnePLModel::printParameterSet(ostream& out)
{
	out << "\"a\" \"b\" \"c\"" << "\n";

	for (int k = 0; k < items; k++)
		out << 1 << " " << (parameterSet[0][0][k]) << " " << 0 << "\n";
}

double OnePLModel::logLikelihood(double* args, double* pars, int nargs, int npars)
{
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
	long double tp, tq;
	long double sum = 0;

	// Obtain q
	q = pars[nP++]; // q is obtained and npars is augmented

	// Obtain I
	It = pars[nP++];

	b = new double[It];

	// Obtain b
	for (i = 0; i < It; i++)
		b[i] = args[nA++];

	for (k = 0; k < q; ++k)
	{
		for (i = 0; i < It; ++i)
		{
			tp = (OnePLModel::successProbability(pars[k + 2], b[i]));

			if (tp == 0)
				tp = 1e-08;
			tq = 1 - tp;
			if (tq == 0)
				tq = 1e-08;
			
			sum += (pars[k * It + i + 2 + (2 * q)] * log(tp))
					+ (pars[k + 2 + q] - pars[k * It + i + 2 + 2 * q]) * log(tq);
		}
	}

	delete[] b;

	return (-sum);
}

double OnePLModel::itemLogLik (double* args, double* pars, int nargs, int npars)
{
	double *theta, *r, *f;
	unsigned int nP, q, items, index;
	double a, sum;
	double tp , tq;
	sum = nP = index = 0;
	
	a = args[0];

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

	for (unsigned int k = 0; k < q; ++k)
	{
		tp = (OnePLModel::successProbability ( theta[k], a));
		
		if (tp<1e-08) tp=1e-08;
		
		tq = 1-tp;
		
		if (tq<1e-08) tq=1e-08;
		
		sum+=(r[k]*log(tp))+(f[k]-r[k])*log(tq);
	}

	args[0] = a;

	delete[] theta;
	delete[] f;
	delete[] r;

	return (-sum);
}

void OnePLModel::gradient(double* args, double* pars, int nargs, int npars, double* gradient)
{
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

void OnePLModel::itemGradient (double* args, double* pars, int nargs, int npars, double* gradient)
{
    double a;
    double *theta, *r, *f;
	int nP, q, items, index;
	long double h;
	long double *P;  // Matrix of size q*I
	long double *factor;	  // Matrix of product (r-fP)W

	nP = index = h = 0;
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

	a = args[0];

	P = new long double [q];
	factor = new long double [q];
	
	for ( int k = 0; k < q; k++ )
	{
		P[k] = successProbability (theta[k], a);
		factor[k] = ( r[k] - f[k]*P[k] );
	}
	
	for ( int k = 0; k < q; k++ )
		h += factor[k];

	memset(gradient,0,sizeof(double));

	delete [] P;
	delete [] factor;

	delete [] theta;
	delete [] r;
	delete [] f;

	//return h as the gradient
	gradient[0]= static_cast<double>(h);
}

void OnePLModel::NitemGradient (double* args, double* pars, int nargs, int npars, double* gradient)
{
 	double hh = 1e-08;
 	
	for(int i = 0 ; i < nargs; i++)
 	{
 		//Add hh to argument
 		args[i] = args[i] + hh;
 		//Evaluate
 		gradient[i] = itemLogLik(args,pars,nargs,npars);
 		//Remove hh to argument and evaluate
 		args[i] = args[i] - hh;
 		gradient[i] -= itemLogLik(args,pars,nargs,npars);
 		//Substract
 		gradient[i] = gradient[i]/hh;
 	}
 	
 	cout << "h → " << gradient[0] << endl;
}
