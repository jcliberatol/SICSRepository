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
class EM2PL: public EMEstimator
{

public:

    virtual void transform() {}

    virtual void untransform()
    {
        double *** pset = m->getParameterModel()->getParameterSet();

        for (int i = 0; i < m->getItemModel()->getDataset()->countItems(); ++i)
            pset[1][0][i] /= -pset[0][0][i];
    }

    virtual void setInitialValues(double *** pset, Model* m) { m->getParameterModel()->setParameterSet(pset); }

    virtual void setInitialValues(int method, Model* m)
    {
        int items = m->getParameterModel()->items;

        pset = m->getParameterModel()->getParameterSet();

        if (method == Constant::RANDOM)
        {
            std::srand(std::time(0));

            for (int i = 0; i < items; i++)
            {
                pset[0][0][i] = randomd() * 2;
                pset[1][0][i] = randomd() * 4 - 2;
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
            }

            delete [] result;
        }
    }

    EM2PL(Model* m, QuadratureNodes* nodes, Matrix<double>* f, Matrix<double>* r)
    {
        this->nodes = nodes;
        this->m = m;
        this->f = f;
        this->r = r;
        this->dims = 2;
        this->sum = 0.0;
        this->data = m->getItemModel()->getDataset();
        this->pm = m->getParameterModel();
        this->q = this->nodes->size();
        this->faux = new long double[q];
        this->weights = this->nodes->getWeight();
        this->items = data->countItems();
        this->fptr = &TwoPLModel::itemLogLik;
        this->gptr = &TwoPLModel::itemGradient;
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
            std::copy(&((*parameters)[2][0]) + (t_size / 3), &((*parameters)[2][0]) + (2 * (t_size / 3)),
                        &(parSet[1][0][0]));
            std::copy(&((*parameters)[2][0]) + (2 * (t_size / 3)), &((*parameters)[2][0]) + (3 * (t_size / 3)),
                        &(parSet[2][0][0]));

            m->getParameterModel()->setParameterSet(parSet);
        }
    }
};

#endif /* EM2PL_H_ */
