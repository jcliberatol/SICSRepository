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
				pset[0][0][i] = -(ppnd(result[0], &ifault)) / result[1];

			delete [] result;
		}
	}

	EM1PL(Model* m, QuadratureNodes* nodes, Matrix<double>* f, Matrix<double>* r) : EMEstimator(m, nodes, f, r)
	{
		this->fptr = &OnePLModel::itemLogLik;
        	this->gptr = &OnePLModel::itemGradient;
		this->dims = 1;
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
