/*
 * EMEstimator.h
 *
 *  Created on: Nov 14, 2014
 *      Author: jliberato
 *      Update by: cesandovalp
 */

#ifndef EMESTIMATOR_H_
#define EMESTIMATOR_H_
#include <string>
#include <type/Matrix.h>
#include <type/PatternMatrix.h>
#include <type/QuadratureNodes.h>
#include <optimizer/Optimizer.h>
#include <util/util.h>
#include <type/Constant.h>
#include <ctime>

class EMEstimator
{

public:

    QuadratureNodes *nodes;
    PatternMatrix *data;
    ParameterModel *pm;
    Model *m;
    int items;
    int q;
    int size;
    int dims;
    int * frequency_list;
    long double sum;
    long double * faux;
    double *** pset;
    Matrix<double>* weights;
    Matrix<double>* f;
    Matrix<double>* r;
    double (*fptr)(double*, double*, int, int);
    void (*gptr)(double*, double*, int, int, double*);
    void (*hptr)(double*, double*, int, int, double*);
    bool ** bitset_list;

    EMEstimator() {}

    EMEstimator(Model* m, QuadratureNodes* nodes, Matrix<double>* f, Matrix<double>* r) {}

    //Reduce iterations
    virtual void stepRamsay(double *** parameters, int * nargs, int t_size, bool continue_flag) = 0;

    //In this cases the model are needed to be filled
    virtual void setInitialValues(int, Model *) = 0;
    virtual void setInitialValues(double ***, Model *) = 0;

    //Step E needs the model , the f and r, and the thetas, besides from the data.
    void stepE()
    {
        double prob;
        double prob_matrix[q][(int) items];
        int k, i;
        int counter_temp[items];
        int counter_set;

        f->reset();
        r->reset();

        //Calculates the success probability for all the nodes.
        m->successProbability(nodes);

        for (k = 0; k < q; ++k)
            for (i = 0; i < items; ++i)
                prob_matrix[k][i] = pm->getProbability(k, i);

        for (int index = 0; index < size; index++)
        {
            sum = 0.0;
            //Calculate g*(k) for all the k's
            //first calculate the P for each k and store it in the array f aux
            for (k = 0; k < q; k++)
            {
                faux[k] = (*weights)(0, k);
                //Calculate the p (iterate over the items in the productory)
                counter_set = 0;

                for (i = 0; i < items; i++) 
                {
                    if (bitset_list[index][i])
                    {
                        counter_temp[counter_set++] = i + 1;
                        prob = prob_matrix[k][i];
                    } else 
                        prob = 1 - prob_matrix[k][i];

                    faux[k] *= prob;
                }
                //At this point the productory is calculated and faux[k] is equivalent to p(u_j,theta_k)
                //Now multiply by the weight
                sum += faux[k];
            }

            for (k = 0; k < q; k++)
            {
                faux[k] *= frequency_list[index] / sum; //This is g*_j_k
                (*f)(0, k) += faux[k];

                for (i = 0; i < counter_set; i++)
                    (*r)(k, counter_temp[i] - 1) += faux[k];
            }
        }
    }

    //Step M also needs the model, quad nodes, f and r
    void stepM(double *** parameters, int * nargs)
    {
        Optimizer optim;
        Matrix<double> ** tri;
        Matrix<double> * thetas;
        double *** pset;
        double * args;
        double * pars;

        int It;
        int q;
        int npars;
        int nA;
        int nP;
        int for_counter;

        *nargs = dims;
        It = m->getItemModel()->getDataset()->countItems();
        q = nodes->size();
        args = new double[dims * It];
        pars = new double[2 + 2 * q + q * It + 1]; 
        npars = 2 + 2 * q + q * It;
        
        //filling args
        nA = 0;
        pset = m->getParameterModel()->getParameterSet();

        tri = new Matrix<double>*[dims];
        
        for(int i = 0; i < dims; i++)
            tri[i] = new Matrix<double>(pset[i], 1, items);

        for(int j = 0; j < dims; j++)
            for (int i = 0; i < It; i++, nA++)
                args[nA] = pset[j][0][i];

        //Filling pars
        nP = 0;
        
        // Obtain q
        pars[nP] = q;
        nP++;

        // Obtain I
        pars[nP] = It;
        nP++;

        // Obtain theta
        thetas = nodes->getTheta();
        for (int k = 0; k < q; k++, nP++)
            pars[nP] = (*thetas)(0, k);
        
        // Obtain f
        for (int k = 0; k < q; k++, nP++)
            pars[nP] = (*f)(0, k);
        
        // Obtain r
        for (int k = 0; k < q; k++)
            for (int i = 0; i < It; i++, nP++)
                pars[nP] = (*r)(k, i);

        *nargs = nA;
        pars[nP++] = 0; //For holding the item.
        npars = nP;

        for_counter = 0;

        for (int i = 0; i < It; i++)
        {
            double* iargs;
            int par_index[dims];
            
            pars[npars - 1] = i;
            iargs = new double[dims];

            for(for_counter = 0; for_counter < dims; for_counter++)
            {
                par_index[for_counter] = i + (i*for_counter) + (for_counter*(items - i));
                iargs[for_counter] = args[par_index[for_counter]];
            }

            if(abs(iargs[0]) > 5)
                iargs[0] = 0.851;

            if(dims > 1)
            {
                if(abs(-iargs[1]/iargs[0]) > 5)
                    iargs[1] = 0;

                if(dims > 2)
                    if(abs(iargs[2]) > 5)
                        iargs[2] = 0.3;
            }

            optim.searchOptimal(fptr, gptr, hptr, iargs, pars, dims, npars);

            for(for_counter = 0; for_counter < dims; for_counter++)
                args[par_index[for_counter]] = iargs[for_counter];

            delete [] iargs;
        }

        std::copy(&((*parameters)[1][0]), (&((*parameters)[1][0])) + *nargs, &((*parameters)[0][0]));
        std::copy(&((*parameters)[2][0]), (&((*parameters)[2][0])) + *nargs, &((*parameters)[1][0]));
        std::copy(&args[0], &args[0] + *nargs, &((*parameters)[2][0]));

        // Now pass the optimals to the Arrays.
        nA = 0;

        // Obtain a
        for (int i = 0; i < It; i++)
        {
            pset[0][0][i] = args[nA++];
            if (fabs(pset[0][0][i]) > abs(5))
                pset[0][0][i] = 0.851;

        }

        // Obtain b
        if(dims > 1)
        for (int i = 0; i < It; i++)
        {
            pset[1][0][i] = args[nA++];
            if (fabs(-pset[1][0][i] / pset[0][0][i]) > abs(5))
                pset[1][0][i] = 0;
        }

        // Obtain c
        if(dims > 2)
        for (int i = 0; i < It; i++)
        {
            pset[2][0][i] = args[nA++];
            if (fabs(pset[2][0][i]) > abs(20))
                pset[2][0][i] = 0.5;
        }

        //Obtain the deltas
        double maxDelta = 0;
        for (int v1 = 0; v1 < It; ++v1)
        {
            for(int j = 0; j < dims; j++)
                tri[j]->setIndex(0, v1, tri[j]->getIndex(0, v1) - pset[j][0][v1]);

            for(int j = 0; j < dims; j++)
                if (fabs(tri[j]->getIndex(0, v1)) > maxDelta)
                    maxDelta = fabs(tri[j]->getIndex(0, v1));
        }
        //TODO change by constant file
        Constant::EPSILONC = maxDelta;
        if (maxDelta < Constant::CONVERGENCE_DELTA)
        {
            m->itemParametersEstimated = true;
        }

        //And set the parameter sets
        double *** parSet = m->getParameterModel()->getParameterSet();
        for(int i = 0; i < dims; i++)
            parSet[i] = pset[i];
        // llenar las tres matrices
        m->getParameterModel()->setParameterSet(parSet);

        for(int i = 0; i < dims; i++)
            delete tri[i];

        delete [] tri;
        delete [] args;
        delete [] pars;
    }

    //Initial values
    double * Andrade()
    {
        PatternMatrix* data;
        int pSize;
        int index;
        double Ni;
        double frequencyV;
        double PII;
        double corr;
        double mT;
        double mU;
        double mTU;
        double mUU;
        double covar;
        double sdU;
        double sdT;
        double *T;
        double *U;
        double *TU;
        double *UU;
        double *Tm;
        double *Um;
        double *result;

        data = m->getItemModel()->getDataset();
        pSize = data->matrix.size();
        Ni = data->countIndividuals();

        T = new double[pSize];
        U = new double[pSize];
        TU = new double[pSize];
        UU = new double[pSize];
        Tm = new double[pSize];
        Um = new double[pSize];

        for (int i = 0; i < items; i++)
        {
            PII = mT = mU = mTU = mUU = 0.0;

            for (index = 0; index < size; index++)
            {
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
            sdT = sdU = 0.0;

            for (index = 0; index < size; index++)
            {
                frequencyV = frequency_list[index];
                Tm[index] = T[index] - mT;
                Um[index] = U[index] - mU;
                sdT += frequencyV * Tm[index] * Tm[index];
                sdU += frequencyV * Um[index] * Um[index];
            }

            sdT = std::sqrt(sdT / (Ni - 1.0));
            sdU = std::sqrt(sdU / (Ni - 1.0));
            corr = covar / (sdT * sdU);
        }

        result = new double[2]{PII, corr};

        delete []T;
        delete []U;
        delete []TU;
        delete []UU;
        delete []Tm;
        delete []Um;

        return (result);
    }

    virtual ~EMEstimator() { delete [] faux; }
};

#endif /* EMESTIMATOR_H_ */
