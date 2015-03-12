/*
 * EM3PL.h
 *
 *  Created on: Nov 14, 2014
 *      Author: jliberato
 */

#ifndef EM3PL_H_
#define EM3PL_H_
#include <estimation/classical/EMEstimators/EMEstimator.h>
#include <model/parameter/ThreePLModel.h>
class EM3PL: public EMEstimator {
public:
	virtual ~EM3PL() {
	}

	virtual void transform() {
		for (int i = 0; i < m->getItemModel()->countItems(); ++i) {
			double *** pset = m->getParameterModel()->getParameterSet();
			double qa = pset[0][0][i];
			double qb = pset[1][0][i];
			double qc = pset[2][0][i];
			pset[2][0][i] = log(qc / (1 - qc));
		}
	}

	virtual void untransform() {
		for (int i = 0; i < m->getItemModel()->countItems(); ++i) {
			double *** pset = m->getParameterModel()->getParameterSet();
			double qa = pset[0][0][i];
			double qb = pset[1][0][i];
			double qc = pset[2][0][i];
			double ec = exp(qc);
			pset[2][0][i] = ec / (1 + ec);
			pset[1][0][i] = -qb / qa; //Transformacion del B   d=-b/a
		}
	}

	virtual void setInitialValues(double *** npset, Model* m) {
		double *** pset = m->getParameterModel()->getParameterSet();
		items = m->getParameterModel()->items;
		for (int i = 0; i < items; i++) {
			pset[0][0][i] = npset[0][0][i];
			pset[1][0][i] = npset[1][0][i];
			pset[2][0][i] = npset[2][0][i];
		}
	}

	virtual void setInitialValues(int method, Model* m) {
		items = m->getParameterModel()->items;

		pset = m->getParameterModel()->getParameterSet();
		for (int i = 0; i < items; i++) {
			pset[0][0][i] = 0;
			pset[1][0][i] = 0;
			pset[2][0][i] = 0;
		}
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
			for (int i = 0; i < items; i++) {
				pset[0][0][i] = randomd() * 2;
				//fill b
				pset[1][0][i] = randomd() * 4 - 2;
				//fill c
				int cassualOptions = 4;
				pset[2][0][i] = randomd() * (2 / (double) cassualOptions);
			}
		}
		//ANDRADE O( items * numberOfPattern )
		if (method == Constant::ANDRADE) {
			Andrade();
			for (int i = 0; i < items; i++) {
				pset[2][0][i] = 0.2;
			}
		}
	}

	EM3PL(Model* m, QuadratureNodes* nodes, Matrix<double>* f,
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
		fptr = &ThreePLModel::logLikelihood;
		gptr = &ThreePLModel::gradient;
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
		double args[3 * It];
		double pars[2 + 2 * q + q * It];
		int nargs = 3 * It;
		int npars = 2 + 2 * q + q * It;
		//filling args
		int nA = 0;
		double *** pset = m->getParameterModel()->getParameterSet();
		double** A = pset[0];
		double** B = pset[1];
		double** C = pset[2];

		Matrix<double> DA(A, 1, items);
		Matrix<double> DB(B, 1, items);
		Matrix<double> DC(B, 1, items);

		for (int i = 0; i < It; i++) {
			DA(0, i) = A[0][i];
			args[nA] = A[0][i];
			nA++;
		}

		// Obtain b
		for (int i = 0; i < It; i++) {
			DB(0, i) = B[0][i];
			args[nA] = B[0][i];
			nA++;
		}
		// Obtain c
		for (int i = 0; i < It; i++) {
			DC(0, i) = C[0][i];
			args[nA] = C[0][i];
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
		// Now pass the optimals to the Arrays.

		nA = 0;
		// Obtain a
		for (int i = 0; i < It; i++) {
			A[0][i] = args[nA++];
			if (fabs(A[0][i]) > abs(5)) { //5
				A[0][i] = 0.851;
				//cout<<"A";
			}

		}
		// Obtain b
		for (int i = 0; i < It; i++) {
			B[0][i] = args[nA++];
			double a = A[0][i];
			double d = B[0][i];
			double b = -d / a;
			if (fabs(b) > abs(5)) {
				B[0][i] = 0;
				//cout<<"B";
			}
		}

		for (int i = 0; i < It; i++) {
			C[0][i] = args[nA++];
			if (fabs(C[0][i]) > abs(20)) {
				C[0][i] = 0.5;
				//cout<<"C";
			}
		}

		//Boundary regularize the arguments
		//B = 0.5;
		//A = 0.851

		//Obtain the deltas
		//Perform substracts
		double maxDelta = 0;
		int DeltaC = 0;
		for (int v1 = 0; v1 < It; ++v1) {
			DA(0, v1) = DA(0, v1) - A[0][v1];
			DB(0, v1) = DB(0, v1) - B[0][v1];
			DC(0, v1) = DC(0, v1) - C[0][v1];
			DeltaC += 3;
			if (fabs(DA(0, v1)) > maxDelta) {
				maxDelta = fabs(DA(0, v1));

			}
			if (fabs(DB(0, v1)) > maxDelta) {
				maxDelta = fabs(DB(0, v1));
			}
			if (fabs(DC(0, v1)) > maxDelta) {
				maxDelta = fabs(DC(0, v1));
			}

		}
		//TODO change by constant file
		Constant::EPSILONC = maxDelta;
		Constant::LOGLIKO = fptr(args, pars, nargs, npars);
		if (maxDelta < Constant::CONVERGENCE_DELTA) {
			m->itemParametersEstimated = true;
		}
		//cout<<maxDelta<<endl;
		//cout<<maxDelta<<endl;
		//And set the parameter sets
		double *** parSet = m->getParameterModel()->getParameterSet();
		parSet[0] = A;
		parSet[1] = B;
		parSet[2] = C;
		// llenar las tres matrices
		m->getParameterModel()->setParameterSet(parSet);

	}
};

#endif /* EM3PL_H_ */
