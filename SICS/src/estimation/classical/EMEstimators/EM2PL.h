/*
 * EM2PL.h
 *
 *  Created on: Nov 14, 2014
 *      Author: jliberato
 */

#ifndef EM2PL_H_
#define EM2PL_H_
#include <estimation/classical/EMEstimators/EMEstimator.h>
#include <model/parameter/TwoPLModel.h>
class EM2PL: public EMEstimator {
private:
	PatternMatrix* data;
	//Trace* profiler;
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
	virtual ~EM2PL() {
	}
	//useless ( used for 3 things )
	virtual void transform() {

	}

	virtual void untransform() {
		double *** pset = m->getParameterModel()->getParameterSet();
		for (int i = 0; i < m->getItemModel()->countItems(); ++i) {
			pset[1][0][i] /= -pset[0][0][i];
		}

		for (int i = 0; i < items; i++) {
						cout<<pset[0][0][i]<<"				";
						cout<<pset[1][0][i]<<"				"<<endl;
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
			for (int i = 0; i < items; i++) {
				pset[0][0][i] = randomd() * 2;
				//fill b
				pset[1][0][i] = randomd() * 4 - 2;
			}
		}

		if (method == Constant::ANDRADE) {
			int numeroDePatrones = 0;
			int iter, ifault;
			PatternMatrix* data =
					dynamic_cast<PatternMatrix *>(m->getItemModel()->getDataset());
			double Ni = data->countIndividuals();
			double PII, frequencyV, mT, mU, mTU, mUU, covar, sdU, sdT, corr,
					result;
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

					T[iter] = 0;
					for (int i_ = 0; i_ < data->size; i_++) {
						if (data->getCurrentBitSet()[i_])
							T[iter]++;
					}
					//T[iter] = data->getCurrentBitSet().count();
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
				pset[0][0][i] = std::sqrt((corr * corr) / (1.0 - corr * corr));
				pset[1][0][i] = -(ppnd(PII, &ifault)) / corr;
			}
		}
	}

	EM2PL(Model* m, QuadratureNodes* nodes, Matrix<double>* f,
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
		fptr = &TwoPLModel::logLikelihood;
		gptr = &TwoPLModel::gradient;
		hptr = NULL;

		//map<boost::dynamic_bitset<>, int>::const_iterator it;
		map<vector<char>, int>::const_iterator it;
		//map<boost::dynamic_bitset<>, int>::const_iterator begin =
		map<vector<char>, int>::const_iterator begin = data->matrix.begin();
		//map<boost::dynamic_bitset<>, int>::const_iterator end =
		map<vector<char>, int>::const_iterator end = data->matrix.end();

		//bitset_list = new boost::dynamic_bitset<>[data->matrix.size()];
		bitset_list = new bool*[data->matrix.size()];
		for (int j = 0; j < data->matrix.size(); j++) {
			bitset_list[j] = new bool[data->size];
		}

		size = data->matrix.size();

		frequency_list = new int[size];

		int counter = 0;
		for (it = begin; it != end; ++it, ++counter) {
			copy(it->first.begin(), it->first.end(), bitset_list[counter]);
			//bitset_list[counter] = it->first;
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
				counter_temp[p]=0;
			}
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
			for (k = 0; k < q; k++) {
				faux[k] = faux[k] / sum; //This is g*_j_k
				//Multiply the f to the frequency of the pattern
				faux[k] = ((long double) frequency_list[index]) * faux[k];
				(*f)(0, k) += faux[k];
				//Now selectively add the faux to the r
				for (i = 0; i < items; i++) {
					if (counter_temp[i] == 0)
						break;
					(*r)(k, counter_temp[i] - 1) += faux[k];
				} // for
			} // for
		}

	}

	virtual void stepM() {
		/*
		 */
		//Step M implementation using the BFGS Algorithm
		/*
		 * What we need
		 * fptr the pointer to loglik
		 * gprt the pointer to gradient
		 * hessptr the pointer to hessian matrix calculatrix
		 * args the a,b, and c
		 * pars, the other parameters q and stuff
		 * nargs, npars, sizes.
		 */
		//fptr
		int It = m->getItemModel()->getDataset()->countItems();
		int q = nodes->size();
		double args[3 * It]; //TODO not 2 * it here ?
		double pars[2 + 2 * q + q * It];
		int nargs = 3 * It;
		int npars = 2 + 2 * q + q * It;
		//filling args
		int nA = 0;
		// Obtain a
		//A Matrix
		double *** pset = m->getParameterModel()->getParameterSet();
		double** A = pset[0];
		double** B = pset[1];
		//Matrix<double>* A = m->getParameterModel()->getParameterSet()[a];
		//B Matrix
		//Matrix<double>* B = m->getParameterModel()->getParameterSet()[d];
		//C Matrix
		//Matrix<double>* C = model->getParameterModel()->getParameterSet()[c];

		Matrix<double> DA(A, 1, items);
		Matrix<double> DB(B, 1, items);
		//Matrix<double> DC(*C);

		for (int i = 0; i < It; i++) {
			args[nA] = A[0][i];		//(*A)(0, i);
			nA++;
		}

		// Obtain b
		for (int i = 0; i < It; i++) {
			args[nA] = B[0][i];		//(*B)(0, i);
			nA++;
		}

		//Filling pars
		int nP = 0;
		// Obtain q
		pars[nP] = q;
		nP++;
		// Obtain I
		pars[nP] = It;
		nP++;
		// Obtain theta
		//Thetas

		Matrix<double>* thetas = nodes->getTheta();
		for (int k = 0; k < q; k++) {
			pars[nP] = (*thetas)(0, k);	//TODO correct indexing on this and nearby matrices
			nP++;
		}
		// Obtain f
		for (int k = 0; k < q; k++) {
			pars[nP] = (*f)(0, k);
			nP++;
		}
		// Obtain r
		for (int k = 0; k < q; k++) {
			for (int i = 0; i < It; i++) {
				pars[nP] = (*r)(k, i);
				nP++;
			}
		}
		nargs = nA;
		npars = nP;
		/*
		 * Chooses the method
		 * method 1 is NR
		 * method 2 is BFGS
		 */
		Optimizer* optim;
		optim = new Optimizer();
		optim->searchOptimal(fptr, gptr, hptr, args, pars, nargs, npars);
		delete optim;
		// Now pass the optimals to the Arrays.

		nA = 0;
		// Obtain a
		for (int i = 0; i < It; i++) {
			A[0][i] = args[nA++];
			if (fabs(A[0][i]) > abs(10)) {
				A[0][i] = 0.851;
				//			cout << "A reset." << endl;
			}

		}
		// Obtain b
		for (int i = 0; i < It; i++) {
			B[0][i] = args[nA++];
			if (fabs(B[0][i]) > abs(-50)) {
				B[0][i] = 0.5;
				//			cout << "B reset." << endl;
			}
		}

		//Boundary regularize the arguments
		//B = 0.5;
		//A = 0.851

		//Obtain the deltas
		//Perform substracts
		double maxDelta = 0;
		double meanDelta = 0;
		int DeltaC = 0;
		for (int v1 = 0; v1 < It; ++v1) {
			DA(0, v1) = DA(0, v1) - A[0][v1];
			DB(0, v1) = DB(0, v1) - B[0][v1];
			meanDelta = +fabs(DA(0, v1));
			meanDelta = +fabs(DB(0, v1));
			DeltaC += 3;
			if (fabs(DA(0, v1)) > maxDelta) {
				maxDelta = fabs(DA(0, v1));
			}
			if (fabs(DB(0, v1)) > maxDelta) {
				maxDelta = fabs(DB(0, v1));
			}
		}
		meanDelta = meanDelta / DeltaC;
		if (meanDelta < 0.0001 and maxDelta < 0.001) {
			m->itemParametersEstimated = true;
		}
		//And set the parameter sets
		double*** parSet;
		//Must set the parset equal to the original memory in the parameter set
		parSet = m->getParameterModel()->getParameterSet();
		parSet[0] = A;
		parSet[1] = B;
		// llenar las tres matrices
		m->getParameterModel()->setParameterSet(parSet);

	}
	;

};
#endif /* EM2PL_H_ */
