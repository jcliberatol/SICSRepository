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

    EMEstimator(Model* m, QuadratureNodes* nodes, Matrix<double>* f, Matrix<double>* r)
    {
        this->nodes = nodes;
        this->m = m;
        this->f = f;
        this->r = r;
        this->sum = 0.0;
        this->data = m->getItemModel()->getDataset();
        this->pm = m->getParameterModel();
        this->items = data->countItems();
        //Here q must change =? in the multidim case ?
        if(m->getDimensionModel()->getNumDimensions() == 1){
           this->q = this->nodes->size();
        }
        else{
            this->q = (int) pow(items,m->getDimensionModel()->getNumDimensions());
        }
        this->faux = new long double[q];
        this->weights = this->nodes->getWeight();
        this->hptr = NULL;
        this->bitset_list = data->getBitsetList();
        this->frequency_list = data->getFrequencyList();

        this->size = data->matrix.size();
    }

    //Reduce iterations
    virtual void stepRamsay(double *** parameters, int * nargs, int t_size, bool continue_flag) = 0;

    //In this cases the model are needed to be filled
    virtual void setInitialValues(int, Model *) = 0;
    virtual void setInitialValues(double ***, Model *) = 0;

    //Step E needs the model , the f and r, and the thetas, besides from the data.
    //Unidim step E
    void stepEUnidim(){
        double prob_matrix[q][(int) items];
        int k, i;
        int counter_temp[items];
        int counter_set;

        f->reset();
        r->reset();

        //Calculates the success probability for all the nodes.
        //The model obj calculates the successProbability
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
                        faux[k] *= prob_matrix[k][i];
                    }
                    else
                    faux[k] *= 1 - prob_matrix[k][i];
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

    void stepEMultidim(){
            //Things to do in EStep multidim

                    /* code */
                    int counter_temp[items];
                    int counter_set;

                    //calling to m function successProbability also needs the dimensional model now.

                    m->successProbability(nodes);

                    cout<<*(m->getParameterModel()->probabilityMatrix);
                    cout<<m->getParameterModel()->probabilityMatrix->nC()<<std::endl;
                    //cout<<m->getParameterModel()->probabilityMatrix->nR()<<std::endl;
                    int totalNodes = m->getParameterModel()->probabilityMatrix->nR();
                    //With this two indexes now we can compute the matrix
                    double prob_matrix[totalNodes][(int) items];
                    //And copy that bitch!
                    for (int k = 0; k < totalNodes; ++k)
                    for (int i = 0; i < items; ++i)
                    prob_matrix[k][i] = pm->getProbability(k, i);
                    //The new node array must be forged with the permutations so some way to forge it must pass
                    //Bringing the weights array.
                    double * mweights = m->getParameterModel()->multiweights;
                   /* for (int i = 0; i < totalNodes; i++) {
                            std::cout<<mweights[i]<<", ";
                    }*/

                    //std::cout<<std::endl;

                    for (int index = 0; index < size; index++)
                    {
                        sum = 0.0;
                        //Calculate g*(k) for all the k's
                        //first calculate the P for each k and store it in the array f aux
                        for (int k = 0; k < q; k++)
                        {
                            faux[k] = mweights[k];
                            //Calculate the p (iterate over the items in the productory)
                            counter_set = 0;

                            for (int i = 0; i < items; i++)
                            {
                                if (bitset_list[index][i])
                                {
                                    counter_temp[counter_set++] = i + 1;
                                    faux[k] *= prob_matrix[k][i];
                                }
                                else
                                faux[k] *= 1 - prob_matrix[k][i];
                            }
                            //At this point the productory is calculated and faux[k] is equivalent to p(u_j,theta_k)
                            //Now multiply by the weight
                            sum += faux[k];
                        }

                        for (int k = 0; k < q; k++)
                        {
                            faux[k] *= frequency_list[index] / sum; //This is g*_j_k
                            (*f)(0, k) += faux[k];

                            for (int i = 0; i < counter_set; i++)
                            (*r)(k, counter_temp[i] - 1) += faux[k];
                        }
                    }

            // 1. Create quadnodes in a wavy fashion
            // 2 . Bring probabilityMatrix
            //3. ?
            //Start Now

    }

    void stepE()
    {
        if(m->getDimensionModel()->getNumDimensions() == 1){
            stepEUnidim();
        }
        else{
            stepEMultidim();
        }
    }


    void stepMUnidim(double *** parameters, int * nargs){

        std::cout<<"Starting , this is going to be implauslbe"<<std::endl;
        std::cout<<"Dims  "<<dims<<std::endl;
        int par_index[dims];
        Optimizer optim;
        Matrix<double> ** tri;
        Matrix<double> * thetas;
        double *** pset;
        double * args, * pars, * iargs;

        int It, q, npars, nA, nP, for_counter;

        *nargs = dims;
        It = m->getItemModel()->getDataset()->countItems();
        q = nodes->size();
        iargs = new double[dims];
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
            pars[npars - 1] = i;

            for(for_counter = 0; for_counter < dims; for_counter++)
            {
                par_index[for_counter] = i + (i*for_counter) + (for_counter*(items - i));
                iargs[for_counter] = args[par_index[for_counter]];
            }

            if(iargs[0]<0)
            iargs[0] = 0.851;

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
        {
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
        m->itemParametersEstimated = true;

        //And set the parameter sets
        double *** parSet = m->getParameterModel()->getParameterSet();
        for(int i = 0; i < dims; i++)
        parSet[i] = pset[i];

        // llenar las tres matrices
        m->getParameterModel()->setParameterSet(parSet);

        for(int i = 0; i < dims; i++)
        delete tri[i];

        delete [] tri;
        delete [] iargs;
        delete [] args;
        delete [] pars;

    }
        void stepMMultidim(double *** parameters, int * nargs){
                Matrix<double> * thetas;
                std::cout<<"Starting , this is going to be implauslbe"<<std::endl;
                dims = m->getDimensionModel()->getNumDimensions();
                std::cout<<"Dims  "<<dims<<std::endl;

                //We need all item pars and all f's and r's nothing more, except for the QuadratureNodes

                Optimizer optim;

                //Number of items
                int It = m->getItemModel()->getDataset()->countItems();
                //Number of quadnodes
                q = nodes->size();
                //std::cout<<" the pset is ;: "<<pset<<std::endl;
                std::cout<<" Dims ;: "<<dims<<std::endl;

                pset = m->getParameterModel()->getParameterSet();
                //Pset is accessed with tripple pointers
                cout<<endl;
                cout<<"A ;"<<std::endl;
                for (size_t yy = 0; yy < dims * It; yy++) {
                        std::cout<<"  : "<<pset[0][0][yy]<<" ";
                }
                cout<<endl;
                cout<<"D ;"<<std::endl;
                for (size_t yy = 0; yy < It; yy++) {
                        std::cout<<"  : "<<pset[1][0][yy]<<" ";
                }
                cout<<endl;
                cout<<"C ;"<<std::endl;
                for (size_t yy = 0; yy <  It; yy++) {
                        std::cout<<"  : "<<pset[2][0][yy]<<" ";
                }

                //Now the f's and R's
                 //(*f); 1 x 100
                 //(*r) 100 x 10;   k, i
                std::cout<<" the f  is ;: "<<"R"<<std::endl<<f->nR()<<"   "<<f->nC()<<std::endl<<"R"<<"  "<<std::endl<<r->nR()<<"   "<<r->nC()<<std::endl;


                //Now the motherfreaking nodules.


                thetas = nodes->getTheta();
                std::cout<<*(thetas)<<"     "<<thetas->nR()<<" "<<thetas->nC()<<std::endl;
                //


                //Now the algorithm for every item is take the item parameters according to the function loglik specified,
                //Loglik receives

                double * args = new double [dims+2];
                // 2 parametreras one for total nodes , this is f->nC() , other for small nodes, this is thetas nC()
                double * pars = new double[2 + (f->nC() ) * 2 + thetas->nC()];
                int nP = 2;
                pars[0] = f->nC();
                pars[1] = thetas->nC();

                for (int oo = 0; oo < thetas->nC(); oo++) {
                        pars[nP ++] = (*thetas)(0,oo);
                }

                for (int oo = 0; oo < f->nC(); oo++) {
                        pars[nP ++ ] = (*f)(0,oo);
                }
                 int nPbk = nP;
                //Now the item filling
                double maxDelta = 0;

                for (int i = 0; i < It; i++) {
                        nP = nPbk;
                        //Item optimization
                        //Fill the R-
                        for (int oo = 0; oo < f->nC(); oo++) {
                                pars[nP ++ ] = (*r)(oo,i);
                        }
                        //Fill the a's
                        for (int oo = 0; oo < dims; oo++) {
                                args[oo] = pset[0][0][i*dims+oo];
                        }
                        //Fill the b's and c's
                        nP = dims;
                        args[nP ++ ] = pset[1][0][i];
                        args[nP ++] = pset[2][0][i];

                        //Arrays seem to be complete-
                        int npars = 2 + thetas -> nC() + f->nC() *  2;
                        int numargs = dims + 2 ;
                        std::cout<<"Ready to optimize , nargs and npars  : "<<numargs<<"  "<<npars<<"  In item : "<<i<<std::endl;
                        for (int arr = 0; arr < numargs; arr++) {
                                std::cout<<"  "<<args[arr]<<" ";
                        }std::cout<<std::endl;

                        optim.searchOptimal(fptr, gptr, hptr, args, pars, numargs, npars);
                        for (int arr = 0; arr < numargs; arr++) {
                                std::cout<<"  "<<args[arr]<<" ";
                        }std::cout<<std::endl;
                        double delta = 0;
                        //Copy the optimal to the pset.
                        for (int oo = 0; oo < dims; oo++) {
                                delta =abs( args[oo] - pset[0][0][i*dims+oo] );
                                if(delta > maxDelta){maxDelta = delta;}
                                pset[0][0][i*dims+oo] = args[oo];
                        }

                        nP = dims;
                        delta =abs( args[nP] - pset[1][0][i] );
                        if(delta > maxDelta){maxDelta = delta;}
                        pset[1][0][i] = args[nP ++ ];
                        delta =abs( args[nP] - pset[2][0][i] );
                        if(delta > maxDelta){maxDelta = delta;}
                        pset[2][0][i] = args[nP ++ ];
                }
                cout<<"Max Delta "<<maxDelta<<std::endl;

                delete [] args;
                delete [] pars;
                /*


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
                m->itemParametersEstimated = true;

                //And set the parameter sets
                double *** parSet = m->getParameterModel()->getParameterSet();
                for(int i = 0; i < dims; i++)
                parSet[i] = pset[i];

                // llenar las tres matrices
                m->getParameterModel()->setParameterSet(parSet);

                for(int i = 0; i < dims; i++)
                delete tri[i];

                delete [] tri;
                delete [] iargs;
                delete [] args;
                delete [] pars;
                */
        }
    //Step M also needs the model, quad nodes, f and r
    void stepM(double *** parameters, int * nargs)
    {
            if(m->getDimensionModel()->getNumDimensions() == 1){
                stepMUnidim(parameters, nargs);
            }
            else{
                stepMMultidim(parameters, nargs);
            }

    }

    //Initial values
    double * Andrade()
    {
        PatternMatrix* data;
        int pSize, index;
        double Ni, frequencyV, PII, corr, mT;
        double sdT, sdU, covar, mUU, mTU, mU;
        double *T, *U, *TU, *UU, *Tm, *Um;
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

        PII = corr = 0;

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
