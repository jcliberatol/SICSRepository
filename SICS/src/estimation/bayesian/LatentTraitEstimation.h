/*
 * LatentTraitEstimation.h
 *
 *      Author: cesandovalp
 */

#ifndef LATENTTRAITESTIMATION_H_
#define LATENTTRAITESTIMATION_H_
#include <model/Model.h>
#include <type/LatentTraits.h>
#include <sstream>
#include <type/Constant.h>
#include <optimizer/Brent_fmin.h>
#include <cstdlib>

#define t_zita t_model->parameterModel->parameterSet
#define gg t_model->parameterModel->successProbability
#define probability_matrix (*model->parameterModel->probabilityMatrix)
#define _A 0][0][i
#define _B 1][0][i
#define _C 2][0][i
#define _DELTA 0.0001220703

class LatentTraitEstimation
{

public:

    Model* model;
    QuadratureNodes * quadNodes;
    LatentTraits * lt;

    typedef struct
    {
        bool * pattern;
        int size;
        int node;
        Model *model;
    } parameters_logL;

    LatentTraitEstimation() {}

    LatentTraitEstimation(PatternMatrix * dataSet) { lt = new LatentTraits(dataSet); }

    virtual ~LatentTraitEstimation() { delete lt; }

    inline double probabilities(bool * pattern, int size, int node)
    {
        double p = 1;

        for (int i = 0; i < size; i++)
        {
            if (pattern[i])
                p *= probability_matrix(node, i);
            else
                p *= 1 - probability_matrix(node, i);
        }

        return (p);
    }

    static double probabilities(double theta, bool * pattern, int size, int node, Model * t_model)
    {
        double p = 1;
        double temp_zita[3];

        for (int i = 0; i < size; i++)
        {
            temp_zita[0] = t_zita[_A];
            temp_zita[1] = t_zita[_B];
            temp_zita[2] = t_zita[_C];

            if (pattern[i])
                p *= gg(theta, temp_zita);
            else

                p *= 1 - gg(theta, temp_zita);
        }

        return (p);
    }

    static double probabilities(double theta, bool * pattern, int size, int node, Model * t_model, double *** parSet)
    {
        double p = 1;
        double temp_zita[3];

        for (int i = 0; i < size; i++)
        {
            temp_zita[0] = parSet[_A];
            temp_zita[1] = parSet[_B];
            temp_zita[2] = parSet[_C];

            if (pattern[i])
                p *= gg(theta, temp_zita);
            else
                p *= 1 - gg(theta, temp_zita);
        }

        return (p);
    }

    LatentTraits * getLatentTraits(){ return (lt); }

    void setLatentTraits(LatentTraits * ltt) { lt = ltt; }

    void setModel(Model* m) { model = m; }

    void estimateLatentTraitsEAP(double *** parSet)
    {
        for (int i = 0; i < lt->pm->size; ++i)
        {
            double qc = parSet[2][0][i];
            parSet[2][0][i] = log(qc / (1 - qc));
            parSet[1][0][i] = -parSet[1][0][i]*parSet[0][0][i];
        }

        bool ** list = lt->pm->getBitsetList();
        int size = lt->pm->matrix.size();

        int counter = 0, i;
        double sum_num, sum_den, pp;
        double temp_zita[3];

        for (int pattern = 0; pattern < size; pattern++, ++counter)
        {
            sum_num = 0;
            sum_den = 0;

            for (i = 0; i < quadNodes->size(); ++i)
            {
                pp = probabilities((*quadNodes->getTheta())(0, i), list[pattern], lt->pm->size, i, this->model, parSet);

                sum_num += (*quadNodes->getTheta())(0, i) * ((*quadNodes->getWeight())(0, i)) * pp;
                sum_den += (*quadNodes->getWeight())(0, i) * pp;
            }

            (*lt->traits)(counter, lt->dim - 1) = sum_num / sum_den;
        }
    }

    void estimateLatentTraitsEAP()
    {
        bool ** list = lt->pm->getBitsetList();
        int size = lt->pm->matrix.size();

        int counter = 0, i;
        double sum_num, sum_den, pp;

        for (int pattern = 0; pattern < size; pattern++, ++counter)
        {
            sum_num = 0;
            sum_den = 0;

            for (i = 0; i < quadNodes->size(); ++i)
            {
                pp = probabilities(list[pattern], lt->pm->size, i);

                sum_num += (*quadNodes->getTheta())(0, i) * ((*quadNodes->getWeight())(0, i)) * pp;
                sum_den += (*quadNodes->getWeight())(0, i) * pp;
            }

            (*lt->traits)(counter, lt->dim - 1) = sum_num / sum_den;
        }
    }

    static double logL(double theta, bool * pattern, int size, int node, Model *model)
    { return (-(log(probabilities(theta, pattern, size, node, model)) - ((theta * theta) / 2))); }

    static double logLP(double theta, bool * pattern, int size, int node, Model *model, double *** parSet)
    { return (-(log(probabilities(theta, pattern, size, node, model, parSet)) - ((theta * theta) / 2))); }

    void estimateLatentTraitsMAP(double *** parSet)
    {
        double BOUNDS[] = {-5,5};
        bool ** pattern_list = lt->pm->getBitsetList();
        int size = lt->pm->matrix.size();

        int counter = 0;

        for (int index = 0; index < size; index++, ++counter)
            (*lt->traits)(counter, lt->dim - 1) = Brent_fmin(BOUNDS, _DELTA, &logLP, pattern_list[index],
                lt->pm->size, counter, this->model, parSet, 1);
    }

    void setQuadratureNodes(QuadratureNodes *nodes) { quadNodes = nodes; }
};

#endif /* LATENTTRAITESTIMATION_H_ */
