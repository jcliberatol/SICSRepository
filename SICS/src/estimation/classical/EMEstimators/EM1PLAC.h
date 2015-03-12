#ifndef EM1PLAC_H_
#define EM1PLAC_H_
#include <estimation/classical/EMEstimators/EMEstimator.h>
#include <model/parameter/OnePLACModel.h>
class EM1PLAC: public EMEstimator {
private:

public:
	virtual ~EM1PLAC() {
	}

	virtual void transform() {

	}

	virtual void untransform() {
		double *** pset = m->getParameterModel()->getParameterSet();
		double qa = pset[0][0][0];
		for (int i = 0; i < m->getItemModel()->countItems(); ++i) {
			pset[1][0][i] = -pset[1][0][i] / qa;
		}
	}

	virtual void setInitialValues(double *** pset, Model* m) {
		m->getParameterModel()->setParameterSet(pset);
	}

	virtual void setInitialValues(int method, Model* m) {
		//TODO MOVE ALGORITHMS TO ANOTHER FILE
		/*TODO
		 * Possible methods
		 * ANDRADE
		 * OSPINA
		 * RANDOM
		 *
		 * The default method is OSPINA
		 */
		int items = m->getParameterModel()->items;
		pset = m->getParameterModel()->getParameterSet();
		if (!method == Constant::RANDOM) {
			std::srand(std::time(0));
			// use current time as seed for random generator
			pset[0][0][0] = randomd() * 2;
			for (int i = 0; i < items; i++) {

				//fill b
				pset[1][0][i] = randomd() * 4 - 2;
			}
		}
		if (method == Constant::ANDRADE) {
			int pSize = 0;
			int ifault;
			PatternMatrix* data =
					dynamic_cast<PatternMatrix *>(m->getItemModel()->getDataset());
			double Ni = data->countIndividuals();
			double PII;
			double frequencyV;
			double mT;
			double mU;
			double mTU;
			double mUU;
			double covar;
			double sdU;
			double sdT;
			double corr;
			double result;

			pSize = data->matrix.size();

			double *T = new double[pSize];
			double *U = new double[pSize];
			double *TU = new double[pSize];
			double *UU = new double[pSize];
			double *Tm = new double[pSize];
			double *Um = new double[pSize];

			for (int i = 0; i < items; i++) {
				PII = 0;
				mT = mU = mTU = mUU = 0.0;
				for (int index = 0; index < size; index++) {
					frequencyV = frequency_list[index];

					T[index] = 0;
					T[index] = data->countBitSet(bitset_list[index], index);
					PII += frequencyV * bitset_list[index][i];
					U[index] = bitset_list[index][i];
					TU[index] = T[index] * U[index];
					UU[index] = U[index] * U[index];
					mT += frequencyV * T[index];
					mU += frequencyV * U[index];
					mTU += frequencyV * TU[index];
					mUU += frequencyV * UU[index];
				}

				PII /= Ni;
				mT /= Ni;
				mU /= Ni;
				mTU /= Ni;
				mUU /= Ni;
				covar = mTU - mU * mT;
				sdT = 0.0;
				sdU = 0.0;

				for (int index = 0; index < size; index++) {
					frequencyV = frequency_list[index];
					Tm[index] = T[index] - mT;
					Um[index] = U[index] - mU;
					sdT += frequencyV * Tm[index] * Tm[index];
					sdU += frequencyV * Um[index] * Um[index];
				}

				sdT = std::sqrt(sdT / (Ni - 1.0));
				sdU = std::sqrt(sdU / (Ni - 1.0));
				corr = covar / (sdT * sdU);
				if (!i)
					pset[0][0][i] = std::sqrt(
							(corr * corr) / (1.0 - corr * corr));
				pset[1][0][i] = -(ppnd(PII, &ifault)) / corr;
			}
		}
	}

	EM1PLAC(Model* m, QuadratureNodes* nodes, Matrix<double>* f,
			Matrix<double>* r) {
		profiler = NULL;
		this->nodes = nodes;
		this->m = m;
		this->f = f;
		this->r = r;
		sum = 0.0;
		data = dynamic_cast<PatternMatrix *>(m->getItemModel()->getDataset());
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

	virtual void stepM(double *** parameters) {
		int It = m->getItemModel()->getDataset()->countItems();
		int q = nodes->size();
		int nargs = It + 1;
		int npars = 2 + 2 * q + q * It;
		double args[nargs];
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
		nargs = nA;
		npars = nP;

		//BFGS
		Optimizer* optim;
		optim = new Optimizer();
		optim->searchOptimal(fptr, gptr, hptr, args, pars, nargs, npars);

		std::copy(&((*parameters)[1][0]), (&((*parameters)[1][0])) + nargs, &((*parameters)[0][0]));
		std::copy(&((*parameters)[2][0]), (&((*parameters)[2][0])) + nargs, &((*parameters)[1][0]));
		std::copy(&args[0], &args[0] + nargs, &((*parameters)[2][0]));

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
		Constant::LOGLIKO = fptr(args, pars, nargs, npars);
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
	;

};
#endif /* EM2PL_H_ */
