/*
 * EM3PL.h
 *
 *  Created on: Nov 14, 2014
 *      Author: jliberato
 */

#ifndef EM3PL_H_
#define EM3PL_H_
#include <estimation/classical/EMEstimators/EMEstimator.h>

class EM3PL : public EMEstimator {
public:
	EM3PL(){}
	virtual ~EM3PL(){}

	virtual void transform(Model* m){
		for (int i = 0; i < m->getItemModel()->countItems(); ++i) {
				double qa = (*m->getParameterModel()->getParameterSet()[a])(0, i);
				double qb = (*m->getParameterModel()->getParameterSet()[d])(0, i);
				double qc = (*m->getParameterModel()->getParameterSet()[c])(0, i);
				(*m->getParameterModel()->getParameterSet()[c])(0, i) = log(
						qc / (1 - qc));
			}
	}

	virtual void untransform(Model* m){
		for (int i = 0; i < m->getItemModel()->countItems(); ++i) {
				double qa = (*m->getParameterModel()->getParameterSet()[a])(0, i);
				double qb = (*m->getParameterModel()->getParameterSet()[d])(0, i);
				double qc = (*m->getParameterModel()->getParameterSet()[c])(0, i);
				//(*model->getParameterModel()->getParameterSet()[d])(0,i)= -qb/qa;
				double ec = exp(qc);
				(*m->getParameterModel()->getParameterSet()[c])(0, i) = ec
						/ (1 + ec);
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
			int items = m->getParameterModel()->getParameterSet()[a]->nC();
			for (int i = 0 ; i < items ; i++){
				(*m->getParameterModel()->getParameterSet()[a])(0, i)= randomd()*2;
				//fill b
				(*m->getParameterModel()->getParameterSet()[d])(0, i)= randomd()*4-2 ;
				//fill c
				int cassualOptions = 4;
				(*m->getParameterModel()->getParameterSet()[c])(0, i)= randomd()*(2/(double)cassualOptions);

			}
		}

		if(!method == Constant::ANDRADE){
			//Andrade method
			int items = m->getParameterModel()->getParameterSet()[a]->nC();
			//sums of the patterns
			int totalscores = 0 ;
			int *itemscores = new int [items];
			memset(itemscores,0,sizeof(int)*items);
			double *covariances = new double [items];
			memset(covariances,0,sizeof(double)*items);
			double variance = 0;
			PatternMatrix* data = dynamic_cast<PatternMatrix *>(m->getItemModel()->getDataset());
			double Ni = (double)data->countIndividuals();
			for (data->resetIterator(); !data->checkEnd(); data->iterate()) {
				double df = (double)data->getCurrentFrequency();
				double bs = (double)data->getCurrentBitSet().count();
				for (int i = 0 ; i < items ; i++){
					if(data->getCurrentBitSet()[i]){
						itemscores[i]+=df;
					}
				}
				totalscores += bs*df;
			}
			cout<<"Score total : "<<totalscores<<endl;
			//calculate variances and covariances
			for (data->resetIterator(); !data->checkEnd(); data->iterate()){
				double df = (double)data->getCurrentFrequency();
				double bs = (double)data->getCurrentBitSet().count();
				for (int i = 0 ; i < items ; i++){
					if(data->getCurrentBitSet()[i]){
						covariances[i]+=((1-itemscores[i]/Ni)*(1-bs/items))*df;
					}
				}
				variance+=((bs-((double)totalscores/Ni))*(bs-((double)totalscores/Ni)))*df;
			}
			variance /= Ni;
			for (int i = 0 ; i < items ; i++){
				covariances[i] /= (Ni-1);
				cout<<"cov : "<<i<<" "<<covariances[i]<<endl;
			}
			//Now calculate the standard deviations for the sums and the items
			long double*stddevs = new long double [items];
			memset(stddevs,0,sizeof(long double)*items);
			long double*pearson = new long double [items];
			memset(pearson,0,sizeof(long double)*items);
			long double*pis = new long double [items];
			memset(pis,0,sizeof(long double)*items);
			for (int i = 0 ; i < items ; i++){
				double avg= totalscores/Ni;
				stddevs[i]=stdDev_bin(itemscores[i],Ni,avg);
				pis[i]=itemscores[i]/Ni;
				pearson[i]=(covariances[i]/(stddevs[i]*std::sqrt(variance)));
				//fill a sqrt(pCoef * pCoef / (1.0 - pCoef * pCoef));
				(*m->getParameterModel()->getParameterSet()[a])(0, i)= std::sqrt((pearson[i]*pearson[i])/(1/pearson[i]*pearson[i]));
				//fill b
				(*m->getParameterModel()->getParameterSet()[d])(0, i)=normalInverse(pis[i]);
				//fill c
				int cassualOptions = 4;

				(*m->getParameterModel()->getParameterSet()[c])(0, i)= 1 / (double)cassualOptions;//TODO CHANGE BY CONSTANT FROM CONST.H FILE
			}
		}
	}

	virtual void stepE(Model* m, Matrix<double>* f, Matrix<double>* r,  QuadratureNodes* nodes){

		//Dataset by patterns
		PatternMatrix* data = dynamic_cast<PatternMatrix *>(m->getItemModel()->getDataset());
		//Pattern iterator is data->iterator
		//Item number
		const double items = data->countItems();
		//Success probability matrix is obtained via pm->getProbability(int,int)
		ParameterModel* pm = m->getParameterModel();
		//Thetas
		Matrix<double>* thetas = nodes->getTheta();
		//Amount of nodes
		const int q = nodes->size();
		//Weights
		Matrix<double>* weights =nodes->getWeight();
		//A Matrix
		Matrix<double>* A = m->getParameterModel()->getParameterSet()[a];
		//B Matrix
		Matrix<double>* B = m->getParameterModel()->getParameterSet()[d];
		//C Matrix
		Matrix<double>* C = m->getParameterModel()->getParameterSet()[c];
		//Auxiliar array for the nodes
		long double faux[q];
		long double sum = 0.0;
		//Restart f and r to zero
		f->reset();
		r->reset();
		//Calculates the success probability for all the nodes.
		m->successProbability(nodes);

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
		fptr = &ThreePLModel::logLikelihood;
		gptr = &ThreePLModel::gradient;
		hptr = NULL;
		int It = m->getItemModel()->getDataset()->countItems();
		int q = nodes->size();
		//TODO Array is not deleted at the end of the method, find memory for the array
		double args[3 * It];
		//TODO Find memory
		double pars[2 + 2 * q + q * It];
		int nargs = 3 * It;
		int npars = 2 + 2 * q + q * It;
		//filling args
		int nA = 0;
		// Obtain a
		//A Matrix
		Matrix<double>* A = m->getParameterModel()->getParameterSet()[a];
		//B Matrix
		Matrix<double>* B = m->getParameterModel()->getParameterSet()[d];
		//C Matrix
		Matrix<double>* C = m->getParameterModel()->getParameterSet()[c];
		//TODO Remove delta creations , find a fast way.
		Matrix<double> DA(*A);
		Matrix<double> DB(*B);
		Matrix<double> DC(*C);

		for (int i = 0; i < It; i++) {
			args[nA] = (*A)(0, i);
			nA++;
		}

		// Obtain b
		for (int i = 0; i < It; i++) {
			args[nA] = (*B)(0, i);
			nA++;
		}

		// Obtain c
		for (int i = 0; i < It; i++) {
			args[nA] = (*C)(0, i);
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
		// Obtain a
		for (int i = 0; i < It; i++) {
			(*A)(0, i) = args[nA++];
			if (fabs((*A)(0, i)) > abs(10)) {
				(*C)(0, i) = 0.852;
			}

		}
		// Obtain b
		for (int i = 0; i < It; i++) {
			(*B)(0, i) = args[nA++];
			if (fabs((*B)(0, i)) > abs(-50)) {
				(*C)(0, i) = 0.5;
			}
		}
		// Obtain c
		for (int i = 0; i < It; i++) {
			(*C)(0, i) = args[nA++];
			if (fabs((*C)(0, i)) > abs(600)) {
				(*C)(0, i) = -1.7364;
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
			DA(0, v1) = DA(0, v1) - (*A)(0, v1);
			DB(0, v1) = DB(0, v1) - (*B)(0, v1);
			DC(0, v1) = DC(0, v1) - (*C)(0, v1);
			meanDelta = +fabs(DA(0, v1));
			meanDelta = +fabs(DB(0, v1));
			meanDelta = +fabs(DC(0, v1));
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
		meanDelta = meanDelta / DeltaC;
		if (meanDelta < 0.0001 and maxDelta < 0.001) {
			m->itemParametersEstimated = true;
		}
		//And set the parameter sets
		map<Parameter, Matrix<double> *> parSet;
		parSet[a] = A;
		parSet[d] = B;
		parSet[c] = C;
		// llenar las tres matrices
		m->getParameterModel()->setParameterSet(parSet);
	}

};


#endif /* EM3PL_H_ */
