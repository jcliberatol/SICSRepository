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

        for (  int k = 0; k < q; k++)
            for (  int i = 0; i < items; i++)
                // Rasch Success Probability Function
                (*probabilityMatrix)(k, i) = successProbability((*quadNodes->getTheta())(0, k), (parameterSet[0][0][i]));
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

inline double OnePLModel::successProbability(double theta, double * zita) { return successProbability(theta, zita[0]); }

void OnePLModel::setParameterSet(double*** par) { this->parameterSet = par; }

double*** OnePLModel::getParameterSet() { return (this->parameterSet); }

void OnePLModel::getParameters(double * parameters)
{
    for(  int i = 0; i < items; i++)
        parameters[i] = 1;
    for(  int i = items; i < 2*items; i++)
        parameters[i] = parameterSet[0][0][i - items];
}

void OnePLModel::setParameters(double * parameters)
{
    for (  int i = 0; i < items; i++)
        this->parameterSet[0][0][i] = parameters[i];
}

double OnePLModel::getProbability(int node, int item) { return ((*probabilityMatrix)(node, item)); }

void OnePLModel::printParameterSet(ostream& out)
{
    out << "\"a\" \"b\" \"c\"" << "\n";

    for (  int k = 0; k < items; k++)
        out << 1 << " " << (parameterSet[0][0][k]) << " " << 0 << "\n";
}

double OnePLModel::itemLogLik (double* args, double* pars, int nargs, int npars)
{
    double *theta, *r, *f;
      int nP, q, items, index;
    double a, sum;
    double tp , tq;

    sum = nP = index = 0;
    a = args[0];
    q = pars[nP ++];
    items = pars[nP ++];
    index = pars[npars-1];
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

    if(abs(a) > 5)
        a = 0.851;

    for (  int k = 0; k < q; ++k)
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

void OnePLModel::itemGradient (double* args, double* pars, int nargs, int npars, double* gradient)
{
    double a;
    double *theta, *r, *f;
    int nP, q, items, index;
    long double h;
    long double *P;  // Matrix of size q*I
    long double *factor;

    nP = index = h = 0;
    index = pars[npars-1];
    a = args[0];

    q = pars[nP ++];
    P = new long double [q];
    factor = new long double [q];
    theta = new double[q];
    r = new double[q];
    f = new double[q];

    // Obtain I
    items = pars[nP ++];

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
