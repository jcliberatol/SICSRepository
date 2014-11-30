#ifndef EM1PLAC_H_
#define EM1PLAC_H_



#include <estimation/classical/EMEstimators/EMEstimator.h>
#include <model/parameter/OnePLACModel.h>
class EM1PLAC : public EMEstimator {
public:
	EM1PLAC(){}
	virtual ~EM1PLAC(){}

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

		if(!method == Constant::ANDRADE){
			int items = m->getParameterModel()->getParameterSet()[d]->nC();
							for (int i = 0; i < items; i++) {

									PatternMatrix* data = dynamic_cast<PatternMatrix *>(m->getItemModel()->getDataset());
									double Ni = data->countIndividuals();
									int numeroDePatrones = 0, iter = 0;
									for (data->resetIterator(); !data->checkEnd(); data->iterate()) numeroDePatrones++; // esto se debe poder hacer de una forma mas optima! en patternMatrix tener el tamaño!
								    double *T = new double[numeroDePatrones],
								    	   *U = new double[numeroDePatrones],
								    	   *TU = new double[numeroDePatrones],
								           *UU = new double[numeroDePatrones],
								    	   *Tm = new double[numeroDePatrones],
								    	   *Um = new double[numeroDePatrones];
								    double PII = 0.0;
								    double frequencyV;
								    cout<<i<<endl;
									for (data->resetIterator(); !data->checkEnd(); data->iterate())
									{
										    frequencyV = data->getCurrentFrequency();
											T[iter] = data->getCurrentBitSet().count();
											PII += frequencyV*data->getCurrentBitSet()[items - i - 1];
											U[iter] = data->getCurrentBitSet()[items - i - 1];
											TU[iter] = T[iter]*U[iter];
											UU[iter] = U[iter]*U[iter];
											iter++;
									}
									PII /= Ni;

									double mT, mU, mTU, mUU, covar, sdU, sdT, corr, result;
									mT = mU = mTU = mUU = 0.0;
									iter = 0;
									mU = 0.0;
									for ( data->resetIterator(); !data->checkEnd(); data->iterate())
									{
										frequencyV = data->getCurrentFrequency();
										mT += frequencyV*T[iter];
										mU += frequencyV*U[iter];
										mTU += frequencyV*TU[iter];
										mUU += frequencyV*UU[iter];
										iter++;
									}
									mT /= Ni;
									mU /= Ni;
									mTU /= Ni;
									mUU /= Ni;
									covar = mTU - mU*mT;
									iter = 0;
									for ( data->resetIterator(); !data->checkEnd(); data->iterate())
									{
										Tm[iter] = T[iter] - mT;
										Um[iter] = U[iter] - mU;
										iter++;
									}

									iter = 0;
									sdT = 0.0; sdU = 0.0;
									for ( data->resetIterator(); !data->checkEnd(); data->iterate())
									{
										frequencyV = data->getCurrentFrequency();
										sdT += frequencyV*Tm[iter]*Tm[iter];
										sdU += frequencyV*Um[iter]*Um[iter];
										iter++;
									}
							        sdT = std::sqrt(sdT/(Ni-1.0));
							        sdU = std::sqrt(sdU/(Ni-1.0));

							        corr = covar/(sdT*sdU);
							        result = std::sqrt((corr*corr)/(1.0-corr*corr)); // este es A
							        // aqui
							        cout <<result<<endl;
							        int ifault;
							        cout<<(ppnd(PII, &ifault))/corr<<endl; // este es B
									//double T = new double[];
								//(*model->getParameterModel()->getParameterSet()[a])(0, 0) = result;
								//fill b
								(*m->getParameterModel()->getParameterSet()[d])(0, i) = (ppnd(PII, &ifault))/corr;
				}
		}
	}

	virtual void stepE(Model* model, Matrix<double>* f, Matrix<double>* r,  QuadratureNodes* nodes){
		PatternMatrix* data =
					dynamic_cast<PatternMatrix *>(model->getItemModel()->getDataset());
			//Pattern iterator is data->iterator
			//Item number
			const double items = data->countItems();
			//Success probability matrix is obtained via pm->getProbability(int,int)
			ParameterModel* pm = model->getParameterModel();
			//Amount of nodes
			const int q = nodes->size();
			//Weights
			Matrix<double>* weights = nodes->getWeight();
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

		double (*fptr)(double*, double*, int, int);
			void (*gptr)(double*, double*, int, int, double*);
			void (*hptr)(double*, double*, int, int, double*);
			fptr = &OnePLACModel::logLikelihood;
			gptr = &OnePLACModel::gradient;
			hptr = NULL;
			//
			//cout<<"Address : "<<&gptr<<" "<<&hptr<<endl;
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
				pars[nP++] = (*thetas)(0, k);	//TODO correct indexing on this and nearby matrices
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

			nA = 0;

			// Obtain a
			(*A)(0,0) = args[nA++];
			if (fabs((*A)(0, 0)) > abs(10)) {
				//cout << "A reset." <<endl;
				(*A)(0,0) = 0.851;
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
			DA(0,0) = DA(0, 0) - (*A)(0,0);
			double maxDelta = DA(0,0);
			for (int v1 = 0; v1 < It; ++v1) {

				DD(0, v1) = (DD(0, v1) - (*D)(0, v1));
				if (fabs(DD(0, v1)) > maxDelta) {
					maxDelta = fabs(DD(0, v1));
				}
			}
			if ( maxDelta < 0.001) {
				m->itemParametersEstimated = true;
			}
			//And set the parameter sets
			map<Parameter, Matrix<double> *> parSet;
			parSet[a] = A;
			parSet[d] = D;
			// llenar las tres matrices
			m->getParameterModel()->setParameterSet(parSet);
	}

};


#endif
