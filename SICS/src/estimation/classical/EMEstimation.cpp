/*
 * EM.cpp
 *
 *  Created on: May 29, 2014
 *      Author: mirt
 */

#include <estimation/classical/EMEstimation.h>
#include <util/util.h>
#include <ctime>

EMEstimation::EMEstimation()
{
    iterations = 0;
    model = NULL;
    f = NULL;
    r = NULL;
    convergenceSignal = false;
    quadNodes = NULL;
    estimator = NULL;
}

EMEstimation::~EMEstimation()
{
    delete estimator;
    delete f;
    delete r;
}

/**
 * Sets the model to be estimated, currently only supports 3PL model
 */
void EMEstimation::setModel(Model * model)
{    
    unsigned int q;
    unsigned int It;
    unsigned int d = 1;

    this->model = model;
    q = quadNodes->size();
    It = this->model->getItemModel()->getDataset()->countItems();

    this->f = new Matrix<double>(d, q);
    this->r = new Matrix<double>(q, It);

    //Discriminate by models
    if (this->model->Modeltype() == Constant::THREE_PL)
    {
        estimator = new EM3PL(this->model, quadNodes, f, r);
        return;
    }

    if (this->model->Modeltype() == Constant::ONE_PL)
    {
        estimator = new EM1PL(this->model, quadNodes, f, r);
        return;
    }

    if (this->model->Modeltype() == Constant::TWO_PL)
    {
        estimator = new EM2PL(this->model, quadNodes, f, r);
        return;
    }

    if (this->model->Modeltype() == Constant::RASCH)
    {
        estimator = new EM1PLAC(this->model, quadNodes, f, r);
        return;
    }
}

/**
 * Sets the initial values for the estimation, use this for inputting a matrix as initial values
 */
void EMEstimation::setInitialValues(double*** parameterSet) { estimator->setInitialValues(parameterSet, model); }

/**
 * Sets the initial values according to a method of calculating the values
 * Possible methods :
 * ANDRADE,
 * OSPINA,
 * RANDOM,
 *
 */
void EMEstimation::setInitialValues(int method) { estimator->setInitialValues(method, model); }

/**
 * Main loop of the EM estimation
 * orchestrates the parameters for the estimation, and holds the estimation
 * for the iterations.
 *
 * TODO : read maxiterations as a input parameter , idea : calculate the max iterations depending on the items
 * TODO : Output last estimation onto the json for recovery in the program.
 */
void ** EMEstimation::estimate()
{
    double ** args_hist;
    int nargs;
    int size;

    // [1] -> iterations
    // [2] -> convergenceSignal
    void ** return_list = new void*[3];

    iterations = 0;
    size = 3 * model->getItemModel()->getDataset()->countItems();

    (args_hist) = new double*[3];
    (args_hist)[0] = new double[size];
    (args_hist)[1] = new double[size];
    (args_hist)[2] = new double[size];

    estimator->pm->transform();

    for (int i = 0; i < 3; i++)
        for (int j = 0; j < size; j++)
            args_hist[i][j] = 0;

    for (;!(iterations++ > Constant::MAX_EM_ITERS || convergenceSignal);)
    {
    //    cout << iterations << endl;
        estimator->stepE();
        estimator->stepM(&args_hist, &nargs);
        estimator->stepRamsay(&args_hist, &nargs, size, iterations > 5 && (iterations) % 3 == 0);

        convergenceSignal = model->itemParametersEstimated;
    }

    return_list[0] = new int(iterations);
    return_list[1] = new bool(convergenceSignal);
    return_list[2] = model->parameterModel->probabilityMatrix;

    estimator->pm->untransform();
    //model->printParameterSet(cout);

    delete [] (args_hist)[0];
    delete [] (args_hist)[1];
    delete [] (args_hist)[2];
    delete [] (args_hist);

    return return_list;
}

/**Returns the iterations that took the estimation to obtain an answer*/
int EMEstimation::getIterations() const { return (iterations); }

QuadratureNodes* EMEstimation::getQuadratureNodes() const { return (quadNodes); }

void EMEstimation::setQuadratureNodes(QuadratureNodes* nodes) { this->quadNodes = nodes; }
