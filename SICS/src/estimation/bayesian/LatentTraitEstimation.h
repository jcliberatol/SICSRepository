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
            if (pattern[i])
                p *= probability_matrix(node, i);
            else
                p *= 1 - probability_matrix(node, i);

        return (p);
    }

    static double probabilities(double theta, bool * pattern, int size, int node, Model * t_model)
    {
        double p = 1;

        for (int i = 0; i < size; i++)
            if (pattern[i] > 0)
                p *= gg(theta, new double[3] { t_zita[_A], t_zita[_B], t_zita[_C] });
            else
                p *= 1 - gg(theta, new double[3] { t_zita[_A], t_zita[_B], t_zita[_C] });

        return (p);
    }

    static double probabilities(double theta, bool * pattern, int size, int node, Model * t_model, double *** parSet)
    {
        double p = 1;

        for (int i = 0; i < size; i++)
            if (pattern[i] > 0)
                p *= gg(theta, new double[3] { parSet[_A], parSet[_B], parSet[_C] });
            else
                p *= 1 - gg(theta, new double[3] { parSet[_A], parSet[_B], parSet[_C] });

        return (p);
    }

    LatentTraits * getLatentTraits(){ return (lt); }

    void setLatentTraits(LatentTraits * ltt) { lt = ltt; }

    void setModel(Model* m) { model = m; }

    void estimateLatentTraitsEAP()
    {
        bool ** list = lt->pm->getBitsetList();
        int size = lt->pm->matrix.size();

        int counter = 0;

        for (int pattern = 0; pattern < size; pattern++, ++counter)
        {
            double sum_num = 0;
            double sum_den = 0;

            for (int i = 0; i < quadNodes->size(); ++i)
            {
                double pp = probabilities(list[pattern], lt->pm->size, i);

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

    static double logLR(double theta, void * params)
    {
        parameters_logL * temp = (parameters_logL*) params;
        return (-(log(probabilities(theta, temp->pattern, temp->size, temp->node, temp->model)) - ((theta * theta) / 2)));
    }

    void printVectors()
    {
        bool ** pattern_list = lt->pm->getBitsetList();

        int counter = 0;

        for (double i = -3.3; i <= 2.8; i += .01)
            cout << i << "," << logL(i, pattern_list[0], lt->pm->size, counter, this->model) << endl;
    }

    void evaluate_theta(double theta, int pattern)
    {
        bool ** pattern_list = lt->pm->getBitsetList();

        int counter = 0;

        cout << theta << "," << logL(theta, pattern_list[pattern], lt->pm->size, counter, this->model) << endl;
    }

    void estimateLatentTraitsMAP()
    {
        bool ** pattern_list = lt->pm->getBitsetList();
        int size = lt->pm->matrix.size();

        int counter = 0;

        for (int index = 0; index < size; index++, ++counter)
            (*lt->traits)(counter, lt->dim - 1) =
            Brent_fmin(Constant::BOUNDS, _DELTA, &logL, pattern_list[index], lt->pm->size, counter, this->model, 1);
    }

    void estimateLatentTraitsMAP(double *** parSet)
    {
        bool ** pattern_list = lt->pm->getBitsetList();
        int size = lt->pm->matrix.size();

        int counter = 0;

        for (int index = 0; index < size; index++, ++counter)
            (*lt->traits)(counter, lt->dim - 1) = Brent_fmin(Constant::BOUNDS, _DELTA, &logLP, pattern_list[index],
                lt->pm->size, counter, this->model, parSet, 1);
    }

    void estimateLatentTraitsMAP_R()
    {
        bool ** pattern_list = lt->pm->getBitsetList();
        int size = lt->pm->matrix.size();
        parameters_logL temp;

        int counter = 0;

        for (int index = 0; index < size; index++, ++counter)
        {
            temp.model = model;
            temp.node = counter;
            temp.pattern = pattern_list[index];
            temp.size = lt->pm->size;

            (*lt->traits)(counter, lt->dim - 1) = Brent_fmin(-5, 5, &logLR, (void*) &temp, 0.0001220703);
        }
    }

    void setQuadratureNodes(QuadratureNodes *nodes)
    {
        quadNodes = nodes;
        model->successProbability(quadNodes);
    }
};

#endif /* LATENTTRAITESTIMATION_H_ */
