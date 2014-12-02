/*
 * EM1PL.h
 *
 *  Created on: Nov 16, 2014
 *      Author: anmrodriguezre
 */

#ifndef EM1PL_H_
#define EM1PL_H_
#include <estimation/classical/EMEstimators/EMEstimator.h>
#include <model/parameter/OnePLModel.h>
class EM1PL : public EMEstimator {
public:
	EM1PL(){}
	virtual ~EM1PL(){}

	virtual void transform(Model* m){
		for (int i = 0; i < m->getItemModel()->countItems(); ++i) {
			double qb = (*m->getParameterModel()->getParameterSet()[b])(0, i);
		}
	}

	virtual void untransform(Model* m){
		for (int i = 0; i < m->getItemModel()->countItems(); ++i) {
			(*m->getParameterModel()->getParameterSet()[b])(0, i) = -(*m->getParameterModel()->getParameterSet()[b])(0, i);
		}
	}

	virtual void setInitialValues(map<Parameter, Matrix<double>*> parameterSet, Model* m) {
		m->getParameterModel()->setParameterSet(parameterSet);
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
		if(!method == Constant::RANDOM){
			std::srand(std::time(0)); // use current time as seed for random generator
			int items = m->getParameterModel()->getParameterSet()[b]->nC();
			for (int i = 0 ; i < items ; i++){
				//fill b
				(*m->getParameterModel()->getParameterSet()[b])(0, i)= randomd()*4-2 ;
			}
		}


		//ANDRADE O( items * numberOfPattern )
		if (method == Constant::ANDRADE) {
			int items = m->getParameterModel()->getParameterSet()[b]->nC(), numeroDePatrones = 0 , iter, ifault;
			PatternMatrix* data =
					dynamic_cast<PatternMatrix *>(m->getItemModel()->getDataset());
			double Ni = data->countIndividuals(), PII, frequencyV, mT, mU, mTU, mUU, covar, sdU, sdT, corr, result;
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
					T[iter] = data->getCurrentBitSet().count();
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
				for (data->resetIterator(); !data->checkEnd(); data->iterate()) {
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
				(*m->getParameterModel()->getParameterSet()[b])(0, i) = -(ppnd(PII, &ifault)) / corr;
			}
		}
	}

	virtual void stepE(Model* model, Matrix<double>* f, Matrix<double>* r,  QuadratureNodes* nodes){
			//Dataset by patterns
			PatternMatrix* data = dynamic_cast<PatternMatrix *>(model->getItemModel()->getDataset());
			//Pattern iterator is data->iterator
			//Item number
			const double items = data->countItems();
			//Success probability matrix is obtained via pm->getProbability(int,int)
			ParameterModel* pm = model->getParameterModel();
			//Thetas
			Matrix<double>* thetas = nodes->getTheta();
			//Amount of nodes
			const int q = nodes->size();
			//Weights
			Matrix<double>* weights =nodes->getWeight();
			//B Matrix
			Matrix<double>* B = model->getParameterModel()->getParameterSet()[b];
			//Auxiliar array for the nodes
			long double faux[q];
			long double sum = 0.0;
			//Restart f and r to zero
			f->reset();
			r->reset();
			//Calculates the success probability for all the nodes.
			model->successProbability(nodes);

			//TODO CAREFULLY PARALLELIZE FOR
			for (data->resetIterator(); !data->checkEnd(); data->iterate()) {
				//Initialize faux in 1 to later calculate the productory
				for (int k = 0; k < q; k++) {
					faux[k] = 1;
				}
				//Calculate g*(k) for all the k's
				//first calculate the P for each k and store it in the array f aux
				for (int k = 0; k < q; k++) {
					//Calculate the p (iterate over the items in the productory)
					for (unsigned int i = 0; i < items; i++) {
						double prob = pm->getProbability(k, i);
						if (!data->getCurrentBitSet()[items - i - 1]) {
							prob = 1 - prob;
						}
						faux[k] = faux[k] * prob;
					}
					//At this point the productory is calculated and faux[k] is equivalent to p(u_j,theta_k)
					//Now multiply by the weight
					faux[k] = faux[k] * (*weights)(0, k);
				}
				//compute the total of the p*a' (denominator of the g*)
				sum = 0.0;
				for (int k = 0; k < q; k++) {
					sum += faux[k];
				}
				for (int k = 0; k < q; k++) {
					faux[k] = faux[k] / sum;	//This is g*_j_k
					//Multiply the f to the frequency of the pattern
					faux[k] = ((long double) data->getCurrentFrequency()) * faux[k];

					(*f)(0, k) += faux[k];
					//Now selectively add the faux to the r
					for (unsigned int i = 0; i < items; i++) {
						if (data->getCurrentBitSet()[items - i - 1]) {
							(*r)(k, i) = (*r)(k, i) + faux[k];
						} // if
					} // for
				} // for

			}
	}


	virtual void stepM(Model* m, Matrix<double>* f, Matrix<double>* r,  QuadratureNodes* nodes){

		//Step M implementation using the BFGS Algorithm
		//Function pointers to represent the loglikelihood, gradient and hessian
		double (*fptr)(double*, double*, int, int);
		void (*gptr)(double*, double*, int, int, double*);
		void (*hptr)(double*, double*, int, int, double*);
		fptr = &OnePLModel::logLikelihood;
		gptr = &OnePLModel::gradient;
		hptr = NULL;
		int It = m->getItemModel()->getDataset()->countItems();
		int q = nodes->size();
		//TODO Array is not deleted at the end of the method, find memory for the array
		double args[It];
		//TODO Find memory
		double pars[2 + 2 * q + q * It];
		int nargs = It;
		int npars = 2 + 2 * q + q * It;
		//filling args
		int nA = 0;
		//B Matrix
		Matrix<double>* B = m->getParameterModel()->getParameterSet()[b];

		//TODO Remove delta creations , find a fast way.
		Matrix<double> DB(*B);

		// Obtain b
		for (int i = 0; i < It; i++) {
			args[nA] = (*B)(0, i);
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
		//Thetas

		Matrix<double>* thetas =nodes->getTheta();
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
		//BFGS
		//TODO Convert optimizer to static function
		Optimizer* optim;
		optim = new Optimizer();
		optim->searchOptimal(fptr, gptr, hptr, args, pars, nargs, npars);

		// Now pass the optimals to the Arrays.

		nA = 0;

		// Obtain b
		for (int i = 0; i < It; i++) {
			(*B)(0, i) = args[nA++];
			if (fabs((*B)(0, i)) > abs(-50)) {
				//(*C)(0, i) = 0.5;
			}
		}
		//Boundary regularize the arguments
		//C = -1.7346
		//B = 0.5;
		//A = 0.851

		//Obtain the deltas
		//Perform substracts
		double maxDelta = 0;
		double meanDelta = 0;
		int DeltaC = 0;
		for (int v1 = 0; v1 < It; ++v1) {
			DB(0, v1) = DB(0, v1) - (*B)(0, v1);
			meanDelta += fabs(DB(0, v1));
			DeltaC ++;
			if (fabs(DB(0, v1)) > maxDelta) {
				maxDelta = fabs(DB(0, v1));
			}
		}
		meanDelta = meanDelta / DeltaC;
		if (meanDelta < 0.0001 and maxDelta < 0.001) {
			m->itemParametersEstimated = true;
		}
		//And set the parameter sets
		map<Parameter, Matrix<double> *> parSet;
		parSet[b] = B;
		// llenar las tres matrices
		m->getParameterModel()->setParameterSet(parSet);
	}

};


#endif /* EM3PL_H_ */
