/*
 * EM3PL.h
 *
 *  Created on: Nov 14, 2014
 *      Author: jliberato
 *      Updated by: cesandovalp
 */

#ifndef EM3PL_H_
#define EM3PL_H_
#include <estimation/classical/EMEstimators/EMEstimator.h>
#include <model/parameter/ThreePLModel.h>

class EM3PL: public EMEstimator
{

public:

	virtual void setInitialValues(double *** npset, Model* m)
	{
		double *** pset = m->getParameterModel()->getParameterSet();

		items = m->getParameterModel()->items;

		for (int i = 0; i < items; i++)
		{
			pset[0][0][i] = npset[0][0][i];
			pset[1][0][i] = npset[1][0][i];
			pset[2][0][i] = npset[2][0][i];
		}
	}

	virtual void setInitialValues(int method, Model* m)
	{
		items = m->getParameterModel()->items;

		pset = m->getParameterModel()->getParameterSet();

		if (method == Constant::RANDOM)
		{
			std::srand(std::time(0));
			for (int i = 0; i < items; i++)
			{
				pset[0][0][i] = randomd() * 2;
				pset[1][0][i] = randomd() * 4 - 2;
				int cassualOptions = 4;
				pset[2][0][i] = randomd() * (2 / (double) cassualOptions);
			}
		}

		if (method == Constant::ANDRADE)
		{
			double * result = Andrade();
			int ifault;
			
			for (int i = 0; i < items; i++)
			{
				pset[0][0][i] = std::sqrt((result[1] * result[1]) / (1.0 - result[1] * result[1]));
				pset[1][0][i] = -(ppnd(result[0], &ifault)) / result[1];
				pset[2][0][i] = 0.2;
			}

			delete [] result;
		}
	}

	EM3PL(Model* m, QuadratureNodes* nodes, Matrix<double>* f, Matrix<double>* r) : EMEstimator(m, nodes, f, r)
    {
        this->fptr = &ThreePLModel::itemLogLik;
        this->gptr = &ThreePLModel::itemGradient;
        this->dims = 3;
    }

	virtual void stepRamsay(double *** parameters, int * nargs, int t_size, bool continue_flag)
	{
		if (continue_flag)
		{
			ramsay(parameters, *nargs);
			double *** parSet = m->getParameterModel()->getParameterSet();

			std::copy(&((*parameters)[2][0]), &((*parameters)[2][0]) + (t_size / 3), &(parSet[0][0][0]));
			std::copy(&((*parameters)[2][0]) + (t_size / 3), &((*parameters)[2][0]) + (2 * (t_size / 3)),
						&(parSet[1][0][0]));
			std::copy(&((*parameters)[2][0]) + (2 * (t_size / 3)), &((*parameters)[2][0]) + (3 * (t_size / 3)),
						&(parSet[2][0][0]));

			m->getParameterModel()->setParameterSet(parSet);
		}
	}
};

#endif /* EM3PL_H_ */
