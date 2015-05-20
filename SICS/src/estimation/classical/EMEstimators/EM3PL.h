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
		fptr = &ThreePLModel::itemLogLik;
		gptr = &ThreePLModel::itemGradient;
		hptr = NULL;

		bitset_list = data->getBitsetList();
		frequency_list = data->getFrequencyList();

		size = data->matrix.size();

	}
	virtual void calculateHessiana()
	{
		int It = m->getItemModel()->getDataset()->countItems();
		int q = nodes->size();
		double args[3 * It];
		double pars[2 + 2 * q + q * It + 1]; //Plus one for item carry.
		int nargs = 3; //Now calculated in item per item basis.
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
		pars[nP++] = 0; //For holding the item.
		npars = nP;
		/*
		 * Chooses the method
		 * method 1 is NR
		 * method 2 is BFGS
		 */
		Optimizer* optim;
		optim = new Optimizer();
		fptr = &ThreePLModel::itemLogLik2;
		double* iargs = new double[3];
		double *xH = new double[3];
		double ** HH;
		HH = new double*[3];
		for( int iterI = 0; iterI < 3; iterI++ ) HH[iterI] = new double[3];
		double ** DeltaH;
		DeltaH = new double*[3];
		for( int iterI = 0; iterI < 3; iterI++ ) DeltaH[iterI] = new double[3];
		double ** xdH;
		xdH = new double*[4];
		for ( int iterI = 0; iterI < 4; iterI++ ) xdH[iterI] = new double[3];
		double *ffH;
		ffH = new double[3];
		for (int i = 0; i < It; i++)
		{

					pars[npars - 1] = i; //Passing the item to fu nctions
					int dims=3;



					int nA = 0;
					nA += i;
					iargs[0] = args[nA];
					int parA = nA;
					nA += (items-i);
					nA += i;
					iargs[1] = args[nA];
					int parB = nA;
					nA += (items-i);
					nA += i;
					iargs[2] = args[nA]; //just nA;
					if(abs(iargs[0])>5){
					iargs[0] = 0.851;

					}

					double dd = 0;
					dd = -iargs[1]/iargs[0];
					if(abs(dd)>5){
					iargs[1] = 0;

					}
					if(abs(iargs[2])>5){
					iargs[2] = 0.3;
					}
					//cout<<"hola2"<<endl;
					fptr(iargs, pars, dims, npars);
					//optim->searchOptimal(fptr, gptr, hptr, iargs, pars, dims, npars);
					//args[parA] = iargs[0];args[parB] = iargs[1];args[nA] = iargs[2];
					//cout<<"hola1";
					int dimH = dims;

					xH[0] = iargs[0];
					xH[1] = iargs[1];
					xH[2] = iargs[2];
					double deltaH = 0.001;


					for ( int iterI = 0; iterI < dimH; iterI++ )
					{
						for ( int iterJ = 0; iterJ < dimH; iterJ++ )
						{
							if ( iterI == iterJ )
							{
								DeltaH[iterI][iterJ] = deltaH;
							}
							else
							{
								DeltaH[iterI][iterJ] = 0.0;
							}
						}
					}

					for ( int iterI = 0; iterI < dimH; iterI++) ffH[iterI] = 0;
					for( int iterI = 0; iterI < dimH; iterI++ )
					{
						for ( int iterJ = 0; iterJ < dimH; iterJ++ )
						{
							for( int iterK = 0; iterK < dimH; iterK++ )
							{
								xdH[0][iterK] = xH[iterK] + DeltaH[iterI][iterK] + DeltaH[iterJ][iterK];
								xdH[1][iterK] = xH[iterK] + DeltaH[iterI][iterK] - DeltaH[iterJ][iterK];
								xdH[2][iterK] = xH[iterK] - DeltaH[iterI][iterK] + DeltaH[iterJ][iterK];
								xdH[3][iterK] = xH[iterK] - DeltaH[iterI][iterK] - DeltaH[iterJ][iterK];

								ffH[0] =  fptr( xdH[0], pars, dims, npars);
								ffH[1] =  fptr( xdH[1], pars, dims, npars);
								ffH[2] =  fptr( xdH[2], pars, dims, npars);
								ffH[3] =  fptr( xdH[3], pars, dims, npars);

								HH[iterI][iterJ] = ( ffH[0] - ffH[1] - ffH[2] + ffH[3]) / (4.0*deltaH*deltaH);
							}
						}
					}
					cout<<"hessiana item:"<<i<<endl;
					for ( int iterI = 0; iterI < dimH; iterI++ )
					{
						for( int iterJ = 0; iterJ < dimH; iterJ++ )
						{
							cout<<HH[iterI][iterJ]<<" ";
						}
						cout<<endl;
					}
					cout<<endl;


		}
		//delete xH;
		//delete[] HH;
		//delete[] DeltaH;
		//delete[] xdH;
		//delete ffH;
		//delete iargs;
	}
	virtual void stepM(double *** parameters, int * nargs) {
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
		double pars[2 + 2 * q + q * It + 1]; //Plus one for item carry.
		*nargs = 3; //Now calculated in item per item basis.
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
		*nargs = nA;
		pars[nP++] = 0; //For holding the item.
		npars = nP;
		/*
		 * Chooses the method
		 * method 1 is NR
		 * method 2 is BFGS
		 */
		Optimizer* optim;
		optim = new Optimizer();
		fptr = &ThreePLModel::itemLogLik2;
		gptr = &ThreePLModel::itemGradient2;
		for (int i = 0; i < It; i++) {
			pars[npars - 1] = i; //Passing the item to fu nctions
			int dims=3;
			double* iargs = new double[3];
			int nA = 0;
			nA += i;
			iargs[0] = args[nA];
			int parA = nA;
			nA += (items-i);
			nA += i;
			iargs[1] = args[nA];
			int parB = nA;
			nA += (items-i);
			nA += i;
			iargs[2] = args[nA]; //just nA;
			if(abs(iargs[0])>5){
			iargs[0] = 0.851;
			}
			double dd = 0;
			dd = -iargs[1]/iargs[0];
			if(abs(dd)>5){
			iargs[1] = 0;

			}
			if(abs(iargs[2])>5){
			iargs[2] = 0.3;
			}
			optim->searchOptimal(fptr, gptr, hptr, iargs, pars, dims, npars);
			args[parA] = iargs[0];args[parB] = iargs[1];args[nA] = iargs[2];
		}
		//Here we call the optimizer for each of the items

		std::copy(&((*parameters)[1][0]), (&((*parameters)[1][0])) + *nargs,
				&((*parameters)[0][0]));
		std::copy(&((*parameters)[2][0]), (&((*parameters)[2][0])) + *nargs,
				&((*parameters)[1][0]));
		std::copy(&args[0], &args[0] + *nargs, &((*parameters)[2][0]));
		// Now pass the optimals to the Arrays.

		nA = 0;
		// Obtain a
		for (int i = 0; i < It; i++) {
			A[0][i] = args[nA++];
			if (fabs(A[0][i]) > abs(5)) { //5
				A[0][i] = 0.851;
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
			}
		}

		for (int i = 0; i < It; i++) {
			C[0][i] = args[nA++];
			if (fabs(C[0][i]) > abs(20)) {
				C[0][i] = 0.5;
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
		if (maxDelta < Constant::CONVERGENCE_DELTA) {
			m->itemParametersEstimated = true;
		}
		//cout<<maxDelta<<endl;
		//And set the parameter sets
		double *** parSet = m->getParameterModel()->getParameterSet();
		parSet[0] = A;
		parSet[1] = B;
		parSet[2] = C;
		// llenar las tres matrices
		m->getParameterModel()->setParameterSet(parSet);

	}

	virtual void stepRamsay(double *** parameters, int * nargs, int t_size,
			bool continue_flag) {
		if (continue_flag) {
			ramsay(parameters, *nargs);
			double *** parSet = m->getParameterModel()->getParameterSet();

			std::copy(&((*parameters)[2][0]),
					&((*parameters)[2][0]) + (t_size / 3), &(parSet[0][0][0]));
			std::copy(&((*parameters)[2][0]) + (t_size / 3),
					&((*parameters)[2][0]) + (2 * (t_size / 3)),
					&(parSet[1][0][0]));
			std::copy(&((*parameters)[2][0]) + (2 * (t_size / 3)),
					&((*parameters)[2][0]) + (3 * (t_size / 3)),
					&(parSet[2][0][0]));

			m->getParameterModel()->setParameterSet(parSet);
		}
	}
};

#endif /* EM3PL_H_ */
