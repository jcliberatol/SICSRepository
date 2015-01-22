#include <model/parameter/OnePLACModel.h>

OnePLACModel::OnePLACModel() {

	parameterSet[a] = NULL;
	parameterSet[b] = NULL;
	parameterSet[c] = NULL;
	parameterSet[d] = NULL;
	probabilityMatrix = NULL;
	nodes = 0;
}

void OnePLACModel::setEstimationNodes(QuadratureNodes* n) {
	this->nodes = n;
}
//Nodes must be brought as a parameter of the dimension model
void OnePLACModel::buildParameterSet(ItemModel* itemModel,
		DimensionModel* dimensionModel) {

	if (typeid(*itemModel) == typeid(DichotomousModel)) {

		if (typeid(*dimensionModel) == typeid(UnidimensionalModel)) {

			int items = itemModel->countItems();
			parameterSet[a] = new Matrix<double>(1, 1);
			parameterSet[d] = new Matrix<double>(1, items);

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

void OnePLACModel::successProbability(DimensionModel *dimensionModel,
		QuadratureNodes * quadNodes) {

	int q = 0;
	double a_d, d_d, theta_d; // d stands from "double"

	if (dimensionModel != NULL) {
		q = quadNodes->size();
	}

	if (typeid(*dimensionModel) == typeid(UnidimensionalModel)) {
		int It = parameterSet[d]->nC();
		if (probabilityMatrix == NULL) {
			//Creates the matrix if it is not already created
			probabilityMatrix = new Matrix<double>(q, It);
		}
		a_d = (*parameterSet[a])(0, 0);
		for (int k = 0; k < q; k++) {
			for (int i = 0; i < It; i++) {
				// 3PL Success Probability Function
				theta_d = (*quadNodes->getTheta())(0, k);
				d_d = (*parameterSet[d])(0, i);
				double p_d = successProbability(theta_d, a_d, d_d);
				(*probabilityMatrix)(k, i) = p_d;

			}
		}
	}

}

map<Parameter, Matrix<double> *> OnePLACModel::getParameterSet() {
	return (this->parameterSet);
}

void OnePLACModel::setParameterSet(
		map<Parameter, Matrix<double> *> parameterSet) {
	this->parameterSet = parameterSet;
}

double OnePLACModel::successProbability(double theta, double a, double d) {
	long double exponential = (Constant::NORM_CONST) * (a * theta + d);
	if (exponential > Constant::MAX_EXP)
		exponential = Constant::MAX_EXP;
	else if (exponential < -(Constant::MAX_EXP * 1.0))
		exponential = -Constant::MAX_EXP;
	return (1.0 / (1.0 + exp(-exponential)));
}

double OnePLACModel::getProbability(int node, int item) {
	return ((*probabilityMatrix)(node, item));
}

void OnePLACModel::gradient(double* args, double* pars, int nargs, int npars,
		double* gradient) {
//	double numerico;
//	double normal;
//	double hh = 0.000001;
//	double grad = logLikelihood(args, pars, nargs, npars);
//	cout<<"numerico"<<endl;
//	for (int i = 0; i < nargs; i++) {
//		args[i] = args[i] + hh;
//		gradient[i] = logLikelihood(args, pars, nargs, npars);
//		args[i] = args[i] - hh;
//		gradient[i] -=  grad;
//		gradient[i] = gradient[i] / hh;
//		cout<<gradient[i]<<endl;
//	}

	int nA = 0;
	int nP = 0;
	int q, items;
	double *theta, *r, *f, *a, *d, *c;

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
	for (int k = 0; k < q; k++) {
		theta[k] = pars[nP++];
	}

	// Obtain f
	for (int k = 0; k < q; k++) {
		f[k] = pars[nP++];
	}

	// Obtain r
	for (int k = 0; k < q; k++) {
		for (int i = 0; i < items; i++) {
			r[k * items + i] = pars[nP++];
		}
	}

	// Obtain a
	a[0] = args[nA++];
	// Obtain d
	for (int i = 0; i < items; i++) {
		d[i] = args[nA++];
	}
	double sumTA = 0.0;
	double sumTBs = 0.0;
	double aux;
	for ( int i = 0; i < items; i++ )
	{
		sumTBs = 0.0;
		for ( int k = 0; k < q; k++ )
		{
			aux = (r[k*items+i]-f[k]*successProbability(theta[k], a[0], d[i]));
			sumTBs += aux;
			sumTA += theta[k]*aux;
		}
		gradient[1+i] = -sumTBs;
	}

    gradient[0] = -sumTA;
//    cout<<"normal"<<endl;
//    for ( int i = 0; i < nargs; i++)
//    	cout<<gradient[i]<<endl;
	delete[] theta;
	delete[] r;
	delete[] f;
	delete[] a;
	delete[] d;

}
double OnePLACModel::logLikelihood(double* args, double* pars, int nargs,
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
	double *theta, *r, *f, *a, *d;

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
	for (int k = 0; k < q; k++) {
		theta[k] = pars[nP++];
	}

	// Obtain f
	for (int k = 0; k < q; k++) {
		f[k] = pars[nP++];
	}

	// Obtain r
	for (int k = 0; k < q; k++) {
		for (int i = 0; i < It; i++) {
			r[k * It + i] = pars[nP++];
		}
	}

	// Obtain a
	a[0] = args[nA++];

	// Obtain b
	for (int i = 0; i < It; i++) {
		d[i] = args[nA++];
	}

	long double tp, tq;
	long double sum = 0;

	for (int k = 0; k < q; ++k) {
		for (unsigned int i = 0; i < It; ++i) {
			tp = (OnePLACModel::successProbability(theta[k], a[0], d[i]));
			if (tp == 0)
				tp = 1e-08;
			tq = 1 - tp;
			if (tq == 0)
				tq = 1e-08;
			sum += (r[k * It + i] * log(tp)) + (f[k] - r[k * It + i]) * log(tq);
		}
	}
	//antiLogit(c, I);
	delete[] theta;
	delete[] f;
	delete[] r;
	delete[] a;
	delete[] d;
	return (-sum);
}

OnePLACModel::~OnePLACModel() {

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

void OnePLACModel::printParameterSet(ostream& out) {
	out << "\"a\" \"b\" \"c\"" << endl;
	int sizeParameter = (*parameterSet[d]).nC();
	for (int _i = 0; _i < sizeParameter; _i++) {
		out << (*parameterSet[a])(0, 0) << " " << (*parameterSet[d])(0, _i)
				<< " 0.25" << endl;
	}
}
