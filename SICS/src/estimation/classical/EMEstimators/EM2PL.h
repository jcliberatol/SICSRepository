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

public:
	virtual ~EM2PL() {
	}

	virtual void transform() {

	}

	virtual void untransform() {
		double *** pset = m->getParameterModel()->getParameterSet();
		for (int i = 0; i < m->getItemModel()->countItems(); ++i) {
			pset[1][0][i] /= -pset[0][0][i];
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
			for (int i = 0; i < items; i++) {
				pset[0][0][i] = randomd() * 2;
				//fill b
				pset[1][0][i] = randomd() * 4 - 2;
			}
		}

		if (method == Constant::ANDRADE) {
			Andrade();
			for (int i = 0; i < items; i++) {
				pset[2][0][i] = 0;
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
		double args[2 * It];
		double pars[2 + 2 * q + q * It];
		int nargs = 2 * It;
		int npars = 2 + 2 * q + q * It;
		//filling args
		int nA = 0;
		// Obtain a
		//A Matrix
		double *** pset = m->getParameterModel()->getParameterSet();
		double** A = pset[0];
		double** B = pset[1];

		Matrix<double> DA(A, 1, items);
		Matrix<double> DB(B, 1, items);

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

		std::copy(&((*parameters)[1][0]), (&((*parameters)[1][0])) + nargs, &((*parameters)[0][0]));
		std::copy(&((*parameters)[2][0]), (&((*parameters)[2][0])) + nargs, &((*parameters)[1][0]));
		std::copy(&args[0], &args[0] + nargs, &((*parameters)[2][0]));

		delete optim;
		// Now pass the optimals to the Arrays.

		for (int i = 0; i < It; i++) {
			A[0][i] = args[i];
			B[0][i] = args[It + i];
		}

		//Boundary regularize the arguments
		//B = 0.5;
		//A = 0.851

		//Obtain the deltas
		//Perform substracts
		double maxDelta = 0;
		double meanDelta = 0;
		for (int v1 = 0; v1 < It; ++v1) {
			DA(0, v1) = DA(0, v1) - A[0][v1];
			DB(0, v1) = DB(0, v1) - B[0][v1];

			if (fabs(DA(0, v1)) > maxDelta) {
				maxDelta = fabs(DA(0, v1));
			}
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
		parSet[0] = A;
		parSet[1] = B;
		// llenar las tres matrices
		m->getParameterModel()->setParameterSet(parSet);

	}
};
#endif /* EM2PL_H_ */
