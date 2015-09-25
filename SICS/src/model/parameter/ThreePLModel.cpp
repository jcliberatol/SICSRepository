/*
 * ThreePLModel.cpp
 *
 *  Created on: 17 Jun 2014
 *      Author: jlgpisa
 */

#include <model/parameter/ThreePLModel.h>
#include <string>
#include <math.h>
#include <util/util.h>

ThreePLModel::ThreePLModel()
{
	parameterSet = NULL;
	probabilityMatrix = NULL;
	nodes = NULL;
	multiweights = NULL;
}

//Same method in uni and multidimensional
void ThreePLModel::transform()
{
	for (unsigned int i = 0; i < itemModel->countItems(); ++i)
	{
		double qc = parameterSet[2][0][i];
		parameterSet[2][0][i] = log(qc / (1 - qc));
	}
}

void ThreePLModel::untransform()
{
	for (unsigned int i = 0; i < itemModel->getDataset()->countItems(); ++i)
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


//remember that model class wraps this so only quad nodes is needed to call this.
void ThreePLModel::successProbability(DimensionModel *dimensionModel, QuadratureNodes * quadNodes)
{
	unsigned int q = 0;
	double a_d, d_d, c_d, theta_d; // d stands from "double"

	if ( dimensionModel != NULL )
		q = quadNodes->size();

	if(typeid(*dimensionModel)==typeid(UnidimensionalModel))
	{
		if(probabilityMatrix == NULL)
			//Creates the matrix if it is not already created
			probabilityMatrix = new Matrix<double>(q,items);

		for (unsigned int k = 0; k < q; k++)
		{
			for (unsigned int i = 0; i < items; i++ )
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

	if(typeid(*dimensionModel)==typeid(MultidimensionalModel))
	{
		//std::cout << "Multidim prob matrix lolololo" << std::endl;
		int dims = dimensionModel->getNumDimensions();
		int totalNodes = pow(q,dims);
		if(probabilityMatrix == NULL){
			//Creates the matrix if it is not already created
			//This matrix has many more columns because in the multidim case it stores

			probabilityMatrix = new Matrix<double>(totalNodes,items);}
			//Filling the matrix
			//Perform times selections in
			//std::cout << "total nodes in matrix : " <<totalNodes<< std::endl;
			//For each node every place in the probmatrix must be filled.
			//LEAKS
			double * theta = new double[dims];
			int * theta_index = new int[dims];
			multiweights = new double[totalNodes];
			for (int k = 0; k < totalNodes; k++) {
				multiweights[k] = 1;
			}
			for (int k = 0; k < q; k++) {
				theta_index[k] = 0;
			}
			for (int k = 0; k < totalNodes; k++) {
				//Calculate theta index
				fullpermutations(dims,q,k,theta_index);
				//Index the theta array at the theta_index
				for (int j = 0; j < dims; j++) {
					theta[j] = (*quadNodes->getTheta())(0,theta_index[j]);
					multiweights[k] *= (*quadNodes->getWeight())(0,theta_index[j]);
				}

				//Now calculate the probability for each item using the theta array.
				// a alias :   parameterSet[0][0] *
				//std::cout << "Must enter the universe" << std::endl;
				for (int i = 0; i < items; i++ )
				{
					// 3PL Success Probability Function
					d_d = parameterSet[1][0][i];
					c_d = parameterSet[2][0][i];
					//std::cout<<"d : "<< d_d << "th : "<< theta[0]<<"  "<< theta[1]<<"	"<<std::endl;
					(*probabilityMatrix)(k,i) = successProbabilityMD ( theta,
					parameterSet[0][0]+(dims*i)  //This is a array passed directly
					, d_d, c_d , dims );
				}
			}

	}
}

double *** ThreePLModel::getParameterSet() { return (this->parameterSet); }

void ThreePLModel::setParameterSet(double ***) { this->parameterSet = parameterSet; }

double ThreePLModel::successProbabilityMD(double * theta, double * a , double d , double c , int dims ){
	double exponential = 0;
	for (int i = 0; i < dims; i++) {
		exponential += theta[i] * a[i];
	}
	exponential += d;
	if ( exponential > Constant::MAX_EXP )
		exponential = Constant::MAX_EXP;

	else if ( exponential < -(Constant::MAX_EXP*1.0) )
		exponential = -Constant::MAX_EXP;

	exponential = exp(-exponential) ;
	double ec = exp(c);

	return ( (ec/(1+ec)) + (1 - (ec/(1+ec))) * (1/(1+exponential)) );
}

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
	unsigned int i = 0;

	for (unsigned int j = 0; j < items; j++)
		parameters[i++] = parameterSet[0][0][j];
	for (unsigned int j = 0; j < items; j++)
		parameters[i++] = parameterSet[1][0][j];
	for (unsigned int j = 0; j < items; j++)
		parameters[i++] = parameterSet[2][0][j];
}

void ThreePLModel::setParameters(double * parameters)
{
	unsigned int i = 0;

	for (unsigned int j = 0; j < items; j++)
		this->parameterSet[0][0][j] = parameters[i++];
	for (unsigned int j = 0; j < items; j++)
		this->parameterSet[1][0][j] = parameters[i++];
	for (unsigned int j = 0; j < items; j++)
		this->parameterSet[2][0][j] = parameters[i++];
}

double ThreePLModel::getProbability(int node, int item) { return ((*probabilityMatrix)(node, item)); }

void ThreePLModel::destroyWeights(){
	delete [] multiweights;
}

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

double ThreePLModel::itemLogLikMultiDim(double* args, double* pars, int nargs, int npars){
	// args are a array of dims size, b and c , singles ...
	// nargs thus is dims + 2.
	int dims = nargs - 2;
	double tp , tq; //Probabilities
	double * a = new double [dims];
	double * thetafull , *r , *f;
 	int nP = 0;
	double sum = 0;
	int q = pars[nP ++ ]; // totalNodes ( 10 nodes at 3 dimensions  = 1000 nodes for instance)
	int qs  = pars[nP ++]; //Small Nodes (10 nodes for instance)

	// pars are : number of quads  (dims), f and r, theta. and thats  it.
	//Pars contain number of thetas,

	thetafull = new double[qs]; // qs is number of quadrature nodes
	r = new double[q];
	f = new double[q];

	// Obtain theta
	for (unsigned int k=0; k<qs; k++)
		thetafull[k] = pars[nP ++];

	// Obtain f
	for (unsigned int k=0; k<q; k++)
		f[k] = pars[nP ++];

	// Obtain r that becomes a vector
	for (unsigned int k=0; k<q; k++)
	{
		r[k] = pars[nP++];
	}
	//Restart
	nP = 0;
	for (unsigned int k=0; k<dims; k++){
		a[k] = args[nP ++];
		if(a[k] > 5 ) a[k] = 0.851;
	}

	double d = args[nP ++];
	double c = args[nP ++];


	if(abs(d)>5)
		d = 0;
	if(abs(c)>5)
		c = 0.1;

	//Here things change because the two thetas must be send to the optimizing function
	//Ergo the same permutations must occur here
	double * theta = new double[dims];
	int * theta_index = new int[dims];

	for (int k = 0; k < dims; k++) {
		theta_index[k] = 0;
	}
	for (int k = 0; k < q; k++) {
		//Calculate theta index
		fullpermutations(dims,qs,k,theta_index);
		//Index the theta array at the theta_index
		for (int j = 0; j < dims; j++) {
			theta[j]  = thetafull[theta_index[k]];
			tp = ThreePLModel::successProbabilityMD(theta, a, d, c , dims);
			if (tp<1e-08) tp=1e-08;
			tq = 1-tp;
			if (tq<1e-08) tq=1e-08;
			sum+=(r[k]*log(tp))+(f[k]-r[k])*log(tq);
		}
	}

	delete[] theta;
	delete[] f;
	delete[] r;
	delete[] a;
	delete[] thetafull;

	return -sum;
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

	for (unsigned int i = 0; i < items; i++)
		out << parameterSet[0][0][i] << " "
		    << parameterSet[1][0][i] << " "
		    << parameterSet[2][0][i] << endl;
}
