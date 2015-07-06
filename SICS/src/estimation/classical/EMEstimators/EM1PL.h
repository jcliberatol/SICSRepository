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

class EM1PL: public EMEstimator
{

private:

public:

	virtual void transform() {}

	virtual void untransform() {}

	virtual void setInitialValues(double *** pset, Model* m) { m->getParameterModel()->setParameterSet(pset); }

	virtual void setInitialValues(int method, Model* m)
	{
		int items = m->getParameterModel()->items;
		pset = m->getParameterModel()->getParameterSet();
		
		if (method == Constant::RANDOM)
		{
			std::srand(std::time(0));
			for (int i = 0; i < items; i++)
				pset[0][0][i] = randomd() * 4 - 2;
		}

		if (method == Constant::ANDRADE)
		{
			double * result = Andrade();
			int ifault;
			for (int i = 0; i < items; i++)
			{
				pset[0][0][i] = -(ppnd(result[0], &ifault)) / result[1];
			}
		}
	}

	EM1PL(Model* m, QuadratureNodes* nodes, Matrix<double>* f, Matrix<double>* r)
	{
		this->nodes = nodes;
		this->m = m;
		this->f = f;
		this->r = r;
		this->dims = 1;
		this->sum = 0.0;
		this->data = m->getItemModel()->getDataset();
		this->pm = m->getParameterModel();
		this->q = this->nodes->size();
		this->faux = new long double[q];
		this->weights = this->nodes->getWeight();
		this->items = data->countItems();
		this->fptr = &OnePLModel::logLikelihood;
		this->gptr = &OnePLModel::gradient;
		this->hptr = NULL;

		this->bitset_list = data->getBitsetList();
		this->frequency_list = data->getFrequencyList();

		this->size = data->matrix.size();
	}

	virtual void stepRamsay(double *** parameters, int * nargs, int t_size, bool continue_flag)
	{
		if (continue_flag)
		{
			ramsay(parameters, *nargs);
			double *** parSet = m->getParameterModel()->getParameterSet();

			std::copy(&((*parameters)[2][0]), &((*parameters)[2][0]) + (t_size / 3), &(parSet[0][0][0]));

			m->getParameterModel()->setParameterSet(parSet);
		}
	}

};

#endif /* EM1PL_H_ */
