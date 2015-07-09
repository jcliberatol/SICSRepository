#ifndef EM1PLAC_H_
#define EM1PLAC_H_
#include <estimation/classical/EMEstimators/EMEstimator.h>
#include <model/parameter/OnePLACModel.h>
class EM1PLAC: public EMEstimator
{

private:

public:

	virtual void setInitialValues(double *** pset, Model* m) { m->getParameterModel()->setParameterSet(pset); }

	virtual void setInitialValues(int method, Model* m)
	{
		int items = m->getParameterModel()->items;
		pset = m->getParameterModel()->getParameterSet();
		if (method == Constant::RANDOM) {
			std::srand(std::time(0));
			// use current time as seed for random generator
			pset[0][0][0] = randomd() * 2;
			for (int i = 0; i < items; i++) {

				//fill b
				pset[1][0][i] = randomd() * 4 - 2;
			}
		}
		if (method == Constant::ANDRADE)
		{
			double * result = Andrade();
			int ifault;
			
			for (int i = 0; i < items; i++)
			{
				if (!i) pset[0][0][0] = std::sqrt((result[1] * result[1]) / (1.0 - result[1] * result[1]));
			    pset[1][0][i] = -(ppnd(result[0], &ifault)) / result[1] ;
			}

			delete [] result;
		}
	}

	EM1PLAC(Model* m, QuadratureNodes* nodes, Matrix<double>* f, Matrix<double>* r)
	{
		this->nodes = nodes;
		this->m = m;
		this->f = f;
		this->r = r;
		sum = 0.0;
		data = m->getItemModel()->getDataset();
		pm = m->getParameterModel();
		q = this->nodes->size();
		faux = new long double[q];
		weights = this->nodes->getWeight();
		items = data->countItems();
		fptr = &OnePLACModel::logLikelihood;
		gptr = &OnePLACModel::gradient;
		hptr = NULL;

		bitset_list = data->getBitsetList();
		frequency_list = data->getFrequencyList();

		size = data->matrix.size();

	}

	virtual void stepM(double *** parameters, int * nargs) {
		int It = m->getItemModel()->getDataset()->countItems();
		int q = nodes->size();
		*nargs = It + 1;
		int npars = 2 + 2 * q + q * It;
		double args[*nargs];
		double pars[npars];

		//filling args
		int nA = 0;
		// Obtain a
		//A Matrix
		double *** pset = m->getParameterModel()->getParameterSet();
		double** A = pset[0];
		double** D = pset[1];
		Matrix<double> DA(A, 1, 1);
		Matrix<double> DD(D, 1, items);

		args[nA++] = A[0][0];

		// Obtain d
		for (int i = 0; i < It; i++) {
			args[nA] = D[0][i];
			nA++;
		}

		//Filling pars
		int nP = 0;
		// Obtain q
		pars[nP++] = q;
		// Obtain I
		pars[nP++] = It;
		// Obtain theta
		//Thetas

		Matrix<double>* thetas = nodes->getTheta();
		for (int k = 0; k < q; k++) {
			pars[nP++] = (*thetas)(0, k);//TODO correct indexing on this and nearby matrices
		}
		// Obtain f
		for (int k = 0; k < q; k++) {
			pars[nP++] = (*f)(0, k);
		}
		// Obtain r
		for (int k = 0; k < q; k++) {
			for (int i = 0; i < It; i++) {
				pars[nP++] = (*r)(k, i);
			}
		}
		*nargs = nA;
		npars = nP;

		//BFGS
		Optimizer* optim;
		optim = new Optimizer();
		optim->searchOptimal(fptr, gptr, hptr, args, pars, *nargs, npars);

		std::copy(&((*parameters)[1][0]), (&((*parameters)[1][0])) + *nargs,
				&((*parameters)[0][0]));
		std::copy(&((*parameters)[2][0]), (&((*parameters)[2][0])) + *nargs,
				&((*parameters)[1][0]));
		std::copy(&args[0], &args[0] + *nargs, &((*parameters)[2][0]));

		nA = 0;

		// Obtain a
		A[0][0] = args[nA++];
		if (fabs(A[0][0]) > abs(10)) {
			//cout << "A reset." <<endl;
			A[0][0] = 0.851;
		}

		// Obtain d
		for (int i = 0; i < It; i++) {
			D[0][i] = args[nA++];
			if (fabs(D[0][i]) > abs(-50)) {
				D[0][i] = 0.5;
				//cout << "D reset." << endl;
			}
		}

		//Perform substracts
		DA(0, 0) = DA(0, 0) - A[0][0];
		double maxDelta = DA(0, 0);
		for (int v1 = 0; v1 < It; ++v1) {

			DD(0, v1) = (DD(0, v1) - D[0][v1]);
			if (fabs(DD(0, v1)) > maxDelta) {
				maxDelta = fabs(DD(0, v1));
			}
		}
		
		Constant::EPSILONC = maxDelta;
		if (maxDelta < Constant::CONVERGENCE_DELTA) {
			m->itemParametersEstimated = true;
		}
		//And set the parameter sets
		double*** parSet;
		parSet = m->getParameterModel()->getParameterSet();
		parSet[0] = A;
		parSet[1] = D;

		// llenar las tres matrices
		m->getParameterModel()->setParameterSet(parSet);

	}

	virtual void stepRamsay(double *** parameters, int * nargs, int t_size,
			bool continue_flag) {

		if (continue_flag) {
			ramsay(parameters, *nargs);
			double *** parSet = m->getParameterModel()->getParameterSet();

			std::copy(&((*parameters)[2][0]),
					&((*parameters)[2][0]) + (t_size / 3), &(parSet[0][0][0]));

			m->getParameterModel()->setParameterSet(parSet);
		}
	}

};
#endif /* EM2PL_H_ */
