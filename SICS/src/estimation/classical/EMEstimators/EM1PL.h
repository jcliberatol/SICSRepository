/*
 * EM1PL.h
 *
 *  Created on: Nov 16, 2014
 *      Author: jcliberatol
 */

#ifndef EM1PL_H_
#define EM1PL_H_
#include <estimation/classical/EMEstimators/EMEstimator.h>
#include <model/parameter/OnePLModel.h>
class EM1PL: public EMEstimator {
private:

public:
	virtual ~EM1PL() {
	}

	virtual void transform() {
	}

	virtual void untransform() {
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
		if (method == Constant::RANDOM) {
			std::srand(std::time(0));
			// use current time as seed for random generator
			for (int i = 0; i < items; i++) {
				//fill b
				pset[0][0][i] = randomd() * 4 - 2;
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
				pset[0][0][i] = -(ppnd(PII, &ifault)) / corr;
			}
		}
	}

	EM1PL(Model* m, QuadratureNodes* nodes, Matrix<double>* f,
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
		fptr = &OnePLModel::logLikelihood;
		gptr = &OnePLModel::gradient;
		hptr = NULL;

		bitset_list = data->getBitsetList();
		frequency_list = data->getFrequencyList();

		size = data->matrix.size();
	}

	virtual void stepM(double *** parameters) {
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
		double args[It];
		double pars[2 + 2 * q + q * It];
		int nargs = It;
		int npars = 2 + 2 * q + q * It;
		//filling args
		int nA = 0;
		double *** pset = m->getParameterModel()->getParameterSet();

		double** B = pset[0];
		//Matrix<double>* B = m->getParameterModel()->getParameterSet()[d];

		Matrix<double> DB(B, 1, items);

		// Obtain b
		for (int i = 0; i < It; i++) {
			args[nA++] = B[0][i];		//(*B)(0, i);
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

		std::copy(&((*parameters)[1][0]), (&((*parameters)[1][0])) + nargs, &((*parameters)[0][0]));
		std::copy(&((*parameters)[2][0]), (&((*parameters)[2][0])) + nargs, &((*parameters)[1][0]));
		std::copy(&args[0], &args[0] + nargs, &((*parameters)[2][0]));

		// Now pass the optimals to the Arrays.
		nA = 0;
		// Obtain b
		for (int i = 0; i < It; i++) {
			B[0][i] = args[nA++];
			if (fabs(B[0][i]) > abs(-50)) {
				B[0][i] = 0.5;
				cout << "B reset." << endl;
			}
		}

		//Boundary regularize the arguments
		//B = 0.5;

		//Obtain the deltas
		//Perform substracts
		double maxDelta = 0;
		double meanDelta = 0;
		for (int v1 = 0; v1 < It; ++v1) {
			DB(0, v1) = DB(0, v1) - B[0][v1];
			if (fabs(DB(0, v1)) > maxDelta) {
				maxDelta = fabs(DB(0, v1));
			}
		}
		Constant::EPSILONC = maxDelta;
		Constant::LOGLIKO = fptr(args, pars, nargs, npars);
		if (maxDelta < Constant::CONVERGENCE_DELTA) {
			m->itemParametersEstimated = true;
		}

		//And set the parameter sets
		double*** parSet;
		//Must set the parset equal to the original memory in the parameter set
		parSet = m->getParameterModel()->getParameterSet();
		parSet[0] = B;
		// llenar las tres matrices
		m->getParameterModel()->setParameterSet(parSet);

	}
	;

};

#endif /* EM1PL_H_ */
