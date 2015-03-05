/*
 * EMEstimator.h
 *
 *  Created on: Nov 14, 2014
 *      Author: jliberato
 */

#ifndef EMESTIMATOR_H_
#define EMESTIMATOR_H_
#include <string>
#include <trace/Trace.h>
#include <type/Matrix.h>
#include <type/PatternMatrix.h>
#include <type/QuadratureNodes.h>
#include <optimizer/Optimizer.h>
#include <util/util.h>
#include <type/Constant.h>

class EMEstimator {
public:

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
	bool ** bitset_list;
	int size;
	int * frequency_list;
	double *** pset;

	Trace* profiler = 0;
	EMEstimator() {
	}
	EMEstimator(Model* m, QuadratureNodes* nodes, Matrix<double>* f,
			Matrix<double>* r) {
	}
	//Transforms the parameters before starting an estimation process
	virtual void transform() = 0;
	//Transforms back the parameters after estimating them
	virtual void untransform() = 0;
	//Step E needs the model , the f and r, and the thetas, besides from the data.
	void stepE() {
		sum = 0.0;
		f->reset();
		r->reset();
		//Calculates the success probability for all the nodes.
		m->successProbability(nodes);

		int k, i;
		double prob;
		double prob_matrix[q][(int) items];

		for (k = 0; k < q; ++k) {
			for (i = 0; i < items; ++i) {
				prob_matrix[k][i] = pm->getProbability(k, i);
			}
		}

		int counter_temp[items];
		int counter_set;

		//TODO CAREFULLY PARALLELIZE FOR
		for (int index = 0; index < size; index++) {
			sum = 0.0;
			//Calculate g*(k) for all the k's
			//first calculate the P for each k and store it in the array f aux
			for (k = 0; k < q; k++) {
				faux[k] = (*weights)(0, k);
				//Calculate the p (iterate over the items in the productory)
				counter_set = 0;
				for (i = 0; i < items; i++) {
					if (bitset_list[index][i]) {
						counter_temp[counter_set++] = i + 1;
						prob = prob_matrix[k][i];
					} else {
						prob = 1 - prob_matrix[k][i];
					}
					faux[k] *= prob;
				}
				//At this point the productory is calculated and faux[k] is equivalent to p(u_j,theta_k)
				//Now multiply by the weight
				sum += faux[k];
			}

			for (k = 0; k < q; k++) {
				faux[k] *= frequency_list[index] / sum; //This is g*_j_k
				(*f)(0, k) += faux[k];
				for (i = 0; i < counter_set; i++)
					(*r)(k, counter_temp[i] - 1) += faux[k];
			} // for
		}
	}
	//Step M also needs the model, quad nodes, f and r
	virtual void stepM(double *** parameters) = 0;
	//in this cases the model are needed to be filled
	virtual void setInitialValues(int, Model*) = 0;
	virtual void setInitialValues(double***, Model*) = 0;

	void Andrade() {
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
			pset[0][0][i] = std::sqrt((corr * corr) / (1.0 - corr * corr));
			pset[1][0][i] = -(ppnd(PII, &ifault)) / corr;
			pset[2][0][i] = 0.2;
		}
	}
	void setProfiler(Trace* t) {
		profiler = t;
	}
	virtual ~EMEstimator() {
	}
};

#endif /* EMESTIMATOR_H_ */
