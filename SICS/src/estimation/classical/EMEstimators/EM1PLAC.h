/*
 * EM2PL.h
 *
 *  Created on: Nov 14, 2014
 *      Author: jliberato
 */

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
	boost::dynamic_bitset<> current_bitset;
	double (*fptr)(double*, double*, int, int);
	void (*gptr)(double*, double*, int, int, double*);
	void (*hptr)(double*, double*, int, int, double*);
	//boost::dynamic_bitset<> * bitset_list;
	bool** bitset_list;
	int size;
	int * frequency_list;

public:
	virtual ~EM1PLAC() {
	}
	//useless ( used for 3 things )
	virtual void transform() {

	}

	virtual void untransform() {
		double qa = (*m->getParameterModel()->getParameterSet()[a])(0, 0);
		for (int i = 0; i < m->getItemModel()->countItems(); ++i) {
			(*m->getParameterModel()->getParameterSet()[d])(0, i) =
					-(*m->getParameterModel()->getParameterSet()[d])(0, i) / qa;
		}
	}

	virtual void setInitialValues(map<Parameter, Matrix<double>*> parameterSet,
			Model* m) {
		m->getParameterModel()->setParameterSet(parameterSet);
	}

	virtual void setInitialValues(int method, Model* m) {
		cout << "EM1PLAC.hSetInitialVAlues" << endl;
		//TODO MOVE ALGORITHMS TO ANOTHER FILE
		/*TODO
		 * Possible methods
		 * ANDRADE
		 * OSPINA
		 * RANDOM
		 *
		 * The default method is OSPINA
		 */
		if (method == Constant::RANDOM) {

			std::srand(std::time(0)); // use current time as seed for random generator
			int items = m->getParameterModel()->getParameterSet()[b]->nC();
			for (int i = 0; i < items; i++) {
				//fill b
				(*m->getParameterModel()->getParameterSet()[b])(0, i) =
						randomd() * 4 - 2;
			}
		}

		//ANDRADE O( items * numberOfPattern )
		if (method == Constant::ANDRADE) {
			int items = m->getParameterModel()->getParameterSet()[d]->nC(),
					numeroDePatrones = 0, iter, ifault;
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
					(*m->getParameterModel()->getParameterSet()[a])(0, 0) =
							std::sqrt((corr * corr) / (1.0 - corr * corr));

				(*m->getParameterModel()->getParameterSet()[d])(0, i) = -(ppnd(
						PII, &ifault)) / corr;
			}
		}
	}

	EM1PLAC(Model* m, QuadratureNodes* nodes, Matrix<double>* f,
			Matrix<double>* r) {
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

		//map<boost::dynamic_bitset<>, int>::const_iterator it;
		map<bool*, int>::const_iterator it;
		//map<boost::dynamic_bitset<>, int>::const_iterator begin =
		map<bool*, int>::const_iterator begin = data->matrix.begin();
		//map<boost::dynamic_bitset<>, int>::const_iterator end =
		map<bool*, int>::const_iterator end = data->matrix.end();

		//bitset_list = new boost::dynamic_bitset<>[data->matrix.size()];
		bitset_list = new bool*[data->matrix.size()];
		for (int j = 0; j < data->matrix.size(); j++) {
			bitset_list[j] = new bool[data->size];
		}

		size = data->matrix.size();

		frequency_list = new int[size];

		int counter = 0;
		for (it = begin; it != end; ++it, ++counter) {
			bitset_list[counter] = it->first;
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
				for (k = 0; k < q; k++) {
					faux[k] = (*weights)(0, k);
					//Calculate the p (iterate over the items in the productory)
					for (i = 0; i < items; i++) {
						prob = prob_matrix[k][i];
						//if (!current_bitset.test(items - i - 1)) {
						if (!bitset_list[index][items - i - 1]) {
							prob = 1 - prob;
						}
						faux[k] = faux[k] * prob;
					}
					//At this point the productory is calculated and faux[k] is equivalent to p(u_j,theta_k)
					//Now multiply by the weight
					sum += faux[k];
				}

				for (k = 0; k < q; k++) {
					faux[k] = faux[k] / sum; //This is g*_j_k
					//Multiply the f to the frequency of the pattern
					faux[k] = ((long double) frequency_list[index]) * faux[k];
					(*f)(0, k) += faux[k];
					//Now selectively add the faux to the r
					for (i = 0; i < items; i++) {
						if (bitset_list[index][items - i - 1]) {
							(*r)(k, i) += faux[k];
						} // if
					} // for
				} // for
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
		Matrix<double>* A = m->getParameterModel()->getParameterSet()[a];

		//B Matrix
		Matrix<double>* D = m->getParameterModel()->getParameterSet()[d];

		Matrix<double> DA(*A);
		Matrix<double> DD(*D);

		args[nA++] = (*A)(0, 0);

		// Obtain d
		for (int i = 0; i < It; i++) {
			args[nA] = (*D)(0, i);
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
		cout<<"inicio optimizador"<<endl;
		optim->searchOptimal(fptr, gptr, hptr, args, pars, nargs, npars);
		cout<<"fin optimizador"<<endl;
		nA = 0;

		// Obtain a
		(*A)(0, 0) = args[nA++];
		if (fabs((*A)(0, 0)) > abs(10)) {
			//cout << "A reset." <<endl;
			(*A)(0, 0) = 0.851;
		}

		// Obtain d
		for (int i = 0; i < It; i++) {
			(*D)(0, i) = args[nA++];
			if (fabs((*D)(0, i)) > abs(-50)) {
				(*D)(0, i) = 0.5;
				//cout << "D reset." << endl;
			}
		}

		//Perform substracts
		DA(0, 0) = DA(0, 0) - (*A)(0, 0);
		double maxDelta = DA(0, 0);
		for (int v1 = 0; v1 < It; ++v1) {

			DD(0, v1) = (DD(0, v1) - (*D)(0, v1));
			if (fabs(DD(0, v1)) > maxDelta) {
				maxDelta = fabs(DD(0, v1));
			}
		}
		if (maxDelta < 0.001) {
			m->itemParametersEstimated = true;
		}
		//And set the parameter sets
		map<Parameter, Matrix<double> *> parSet;
		parSet[a] = A;
		parSet[d] = D;
		// llenar las tres matrices
		m->getParameterModel()->setParameterSet(parSet);

	}
	;

};
#endif /* EM2PL_H_ */
