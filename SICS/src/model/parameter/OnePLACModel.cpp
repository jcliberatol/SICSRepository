#include <model/parameter/OnePLACModel.h>

OnePLACModel::OnePLACModel()
{
	items = 0;
	parameterSet = NULL;
	probabilityMatrix = NULL;
	nodes = 0;
}

void OnePLACModel::setEstimationNodes(QuadratureNodes* n) { this->nodes = n; }

void OnePLACModel::successProbability(DimensionModel *dimensionModel, QuadratureNodes * quadNodes)
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

		a_d = parameterSet[0][0][0];

		for (int k = 0; k < q; k++)
		{
			for (int i = 0; i < items; i++)
			{
				theta_d = (*quadNodes->getTheta())(0, k);
				d_d = parameterSet[1][0][i];
				(*probabilityMatrix)(k, i) = successProbability(theta_d, a_d, d_d);
			}
		}
	}
}

double *** OnePLACModel::getParameterSet() { return (this->parameterSet); }

void OnePLACModel::getParameters(double * parameters)
{
	parameters[0] = parameterSet[0][0][0];

	for ( int i = 0; i < items; i++ )
		parameters[i+1] = parameterSet[1][0][i];
}

void OnePLACModel::setParameters(double * parameters)
{
	this->parameterSet[0][0][0] = parameters[0];

	for ( int i = 0; i < items; i++ )
		 this->parameterSet[1][0][i] = parameters[i+1];
}

void OnePLACModel::setParameterSet(double *** parameterSet) { this->parameterSet = parameterSet; }

double OnePLACModel::successProbability(double theta, double a, double d)
{
	double exponential = (Constant::NORM_CONST) * (a * theta + d);
	
	if (exponential > Constant::MAX_EXP)
		exponential = Constant::MAX_EXP;
	else if (exponential < -(Constant::MAX_EXP * 1.0))
		exponential = -Constant::MAX_EXP;

	return (1.0 / (1.0 + exp(-exponential)));
}

double OnePLACModel::successProbability(double theta, double * zita)
{
	double exponential = (Constant::NORM_CONST) * (zita[0] * theta + zita[1]);

	if (exponential > Constant::MAX_EXP)
		exponential = Constant::MAX_EXP;
	else if (exponential < -(Constant::MAX_EXP * 1.0))
		exponential = -Constant::MAX_EXP;

	return (1.0 / (1.0 + exp(-exponential)));
}

double OnePLACModel::getProbability(int node, int item) { return ((*probabilityMatrix)(node, item)); }

void OnePLACModel::gradient(double* args, double* pars, int nargs, int npars, double* gradient)
{
	int nA = 0;
	int nP = 0;
	int q, items;
	double *theta, *r, *f, *a, *d;
	double sumTA = 0.0;
	double sumTBs = 0.0;
	double aux;

	// Obtain q
	q = pars[nP++]; // q is obtained and npars is augmented
	// Obtain I
	items = pars[nP++];

	theta = new double[q];
	r = new double[q * items];
	f = new double[q];
	a = new double[1];
	d = new double[items];

	// Obtain theta
	for (int k = 0; k < q; k++)
		theta[k] = pars[nP++];

	// Obtain f
	for (int k = 0; k < q; k++)
		f[k] = pars[nP++];

	// Obtain r
	for (int k = 0; k < q; k++)
		for (int i = 0; i < items; i++)
			r[k * items + i] = pars[nP++];

	// Obtain a
	a[0] = args[nA++];
	
	// Obtain d
	for (int i = 0; i < items; i++)
		d[i] = args[nA++];

	for (int i = 0; i < items; i++)
	{
		sumTBs = 0.0;
		
		for (int k = 0; k < q; k++)
		{
			aux = (r[k * items + i] - f[k] * successProbability(theta[k], a[0], d[i]));
			sumTBs += aux;
			sumTA += theta[k] * aux;
		}

		gradient[1 + i] = -sumTBs;
	}

	gradient[0] = -sumTA;

	delete[] theta;
	delete[] r;
	delete[] f;
	delete[] a;
	delete[] d;
}

double OnePLACModel::logLikelihood(double* args, double* pars, int nargs, int npars)
{
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
	double *theta, *r, *f, *a, *d;
	double tp, tq;
	double sum = 0;

	// Obtain q
	q = pars[nP++]; // q is obtained and npars is augmented

	// Obtain I
	It = pars[nP++];

	theta = new double[q];
	r = new double[q * It];
	f = new double[q];
	a = new double[It];
	d = new double[It];

	// Obtain theta
	for (int k = 0; k < q; k++)
		theta[k] = pars[nP++];

	// Obtain f
	for (int k = 0; k < q; k++)
		f[k] = pars[nP++];

	// Obtain r
	for (int k = 0; k < q; k++)
		for (int i = 0; i < It; i++)
			r[k * It + i] = pars[nP++];

	// Obtain a
	a[0] = args[nA++];

	// Obtain b
	for (int i = 0; i < It; i++)
		d[i] = args[nA++];

	for (int k = 0; k < q; ++k)
	{
		for (int i = 0; i < It; ++i)
		{
			tp = (OnePLACModel::successProbability(theta[k], a[0], d[i]));

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
	delete[] a;
	delete[] d;

	return (-sum);
}

void OnePLACModel::printParameterSet(ostream& out)
{
	out << "\"a\" \"b\" \"c\"" << endl;

	for (int _i = 0; _i < items; _i++)
		out << parameterSet[0][0][0] <<
		" " << parameterSet[1][0][_i] <<
		" 0.25" << endl;
}
