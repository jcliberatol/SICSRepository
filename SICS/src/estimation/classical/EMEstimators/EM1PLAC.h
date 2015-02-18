#ifndef EM1PLAC_H_
#define EM1PLAC_H_
#include <estimation/classical/EMEstimators/EMEstimator.h>
#include <model/parameter/OnePLACModel.h>
class EM1PLAC: public EMEstimator {
private:
	PatternMatrix* data;
	Model* m;
	int items;
	ParameterModel* pm;
	QuadratureNodes* nodes;
	int q;
	Matrix<double>* weights;
	long double * faux;
	long double sum;
	Matrix<double>* f;
	Matrix<double>* r;
	double (*fptr)(double*, double*, int, int);
	void (*gptr)(double*, double*, int, int, double*);
	void (*hptr)(double*, double*, int, int, double*);
	bool** bitset_list;
	int size;
	int * frequency_list;

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
		double *** pset = m->getParameterModel()->getParameterSet();
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
			int numeroDePatrones = 0, iter, ifault;
			PatternMatrix* data =
					dynamic_cast<PatternMatrix *>(m->getItemModel()->getDataset());
			double Ni = data->countIndividuals(), PII, frequencyV, mT, mU, mTU,
					mUU, covar, sdU, sdT, corr, result;
			for (data->resetIterator(); !data->checkEnd(); data->iterate())
				numeroDePatrones++; // esto se debe poder hacer de una forma mas optima! en patternMatrix tener el tamaño!
			double *T = new double[numeroDePatrones], *U =
					new double[numeroDePatrones], *TU =
					new double[numeroDePatrones], *UU =
					new double[numeroDePatrones], *Tm =
					new double[numeroDePatrones], *Um =
					new double[numeroDePatrones];
			for (int i = 0; i < items; i++) {

				iter = 0;
				PII = 0;
				mT = mU = mTU = mUU = 0.0;
				for (data->resetIterator(); !data->checkEnd();
						data->iterate()) {
					frequencyV = data->getCurrentFrequency();
					T[iter] = 0; // do this in a function
					for (int i_ = 0; i_ < data->size; i_++) {
						if (data->getCurrentBitSet()[i_])
							T[iter]++;
					}
					PII += frequencyV * data->getCurrentBitSet()[items - i - 1];
					U[iter] = data->getCurrentBitSet()[items - i - 1];
					TU[iter] = T[iter] * U[iter];
					UU[iter] = U[iter] * U[iter];
					mT += frequencyV * T[iter];
					mU += frequencyV * U[iter];
					mTU += frequencyV * TU[iter];
					mUU += frequencyV * UU[iter];
					iter++;
				}
				PII /= Ni;
				mT /= Ni;
				mU /= Ni;
				mTU /= Ni;
				mUU /= Ni;
				covar = mTU - mU * mT;
				iter = 0;
				sdT = 0.0;
				sdU = 0.0;
				for (data->resetIterator(); !data->checkEnd();
						data->iterate()) {
					frequencyV = data->getCurrentFrequency();
					Tm[iter] = T[iter] - mT;
					Um[iter] = U[iter] - mU;
					sdT += frequencyV * Tm[iter] * Tm[iter];
					sdU += frequencyV * Um[iter] * Um[iter];
					iter++;
				}
				sdT = std::sqrt(sdT / (Ni - 1.0));
				sdU = std::sqrt(sdU / (Ni - 1.0));
				corr = covar / (sdT * sdU);
				if (!i)
					pset[0][0][0] = std::sqrt(
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

		map<vector<char>, int>::const_iterator it;
		map<vector<char>, int>::const_iterator begin = data->matrix.begin();
		map<vector<char>, int>::const_iterator end = data->matrix.end();

		bitset_list = new bool*[data->matrix.size()];
		for (int j = 0; j < data->matrix.size(); j++) {
			bitset_list[j] = new bool[data->size];
		}

		size = data->matrix.size();

		frequency_list = new int[size];

		int counter = 0;
		for (it = begin; it != end; ++it, ++counter) {
			copy(it->first.begin(), it->first.end(), bitset_list[counter]);
			frequency_list[counter] = it->second;
		}
	}
	virtual void stepE() {
		sum = 0.0;
		f->reset();
		r->reset();
		//Calculates the success probability for all the nodes.
		m->successProbability(nodes);

		int k, i;
		double prob;
		double prob_matrix[q][(int) items];

		for (k = 0; k < q; k++) {
			for (i = 0; i < items; i++) {
				prob_matrix[k][i] = pm->getProbability(k, i);
			}
		}

		//TODO CAREFULLY PARALLELIZE FOR
		for (int index = 0; index < size; index++) {
			sum = 0.0;
			//Calculate g*(k) for all the k's
			//first calculate the P for each k and store it in the array f aux

			int counter_temp[items];
			for (int p = 0; p < items; ++p) {
				counter_temp[p] = 0;
			}
			//profiler->startTimer("for1");
			for (k = 0; k < q; k++) {
				faux[k] = (*weights)(0, k);
				//Calculate the p (iterate over the items in the productory)
				int counter_set = 0;
				for (i = 0; i < items; i++) {
					if (bitset_list[index][items - i - 1]) {
						counter_temp[counter_set++] = i + 1;
						prob = prob_matrix[k][i];
					} else {
						prob = 1 - prob_matrix[k][i];
					}
					faux[k] = faux[k] * prob;
				}
				//At this point the productory is calculated and faux[k] is equivalent to p(u_j,theta_k)
				//Now multiply by the weight
				sum += faux[k];
			}
			//profiler->stopTimer("for1");
			//profiler->startTimer("for2");
			for (k = 0; k < q; k++) {
				faux[k] *= frequency_list[index] / sum; //This is g*_j_k
				(*f)(0, k) += faux[k];
				//Now selectively add the faux to the r
				for (i = 0; i < items; i++) {
					if (counter_temp[i] == 0)
						break;
					(*r)(k, counter_temp[i] - 1) += faux[k];
				} // for
			} // for
			  //profiler->stopTimer("for2");
		}

	}
	virtual void stepM() {
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
		m->Hessiana = new double[(nargs*(nargs+1))/2];
		optim->searchOptimal(fptr, gptr, hptr, args, pars, nargs, npars, m->Hessiana);
		if (Constant::ITER > 3) {
			if (Constant::ITER % 3 == 1) {
				m->back_2 = new double[nargs];
				for ( int ii = 0; ii < nargs; ii++ ) m->back_2[ii] = args[ii];

			} else if (Constant::ITER % 3 == 2) {
				m->back_1 = new double[nargs];
				for ( int ii = 0; ii < nargs; ii++ ) m->back_1[ii] = args[ii];

			} else {
				ramsay(args, m->back_1, m->back_2, nargs);
			}
		}
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
