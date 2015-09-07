/*
 * LatentTraits.h
 *
 *  Created on: Feb 2, 2015
 *      Author: cristian
 */

#ifndef TYPE_LATENTTRAITS_H_
#define TYPE_LATENTTRAITS_H_
#include <type/PatternMatrix.h>

class LatentTraits
{

public:

    unsigned int dim;

    LatentTraits(PatternMatrix * p, const unsigned int dims = 1)
    {
        unsigned int rows;
        
        pm = p;
        rows = pm->matrix.size();
        traits = new Matrix<double>(rows, dims);
        dim = dims;
    }

    double ** getListPatternTheta()
    {
        double ** result;
        bool ** pattern_list;
        unsigned int items;

        pattern_list = pm->getBitsetList();
        result = new double*[pm->matrix.size()];
        items = pm->countItems();

        for(unsigned int i = 0; i < pm->matrix.size(); i++)
        {
            result[i] = new double[items + 1];

            for(unsigned int j = 0; j < items; j++)
                result[i][j] = pattern_list[i][j];
            
            result[i][items] = (*traits)(i,0);
        }

        return result;
    }

    void deleteListPatternTheta(double ** p)
    {
        for(unsigned int i = 0; i < pm->matrix.size(); i++)
            delete [] p[i];
        delete [] p;
    }

    virtual ~LatentTraits() { delete traits; };

    void print()
    {
        bool ** pattern_list = pm->getBitsetList();

        for(unsigned int i = 0; i < pm->matrix.size(); i++)
        {
            for(unsigned int j = 0; j < pm->countItems(); j++)
                cout << pattern_list[i][j] << " ";
            cout<<(*traits)(i,0)<<endl;
        }
    }

    PatternMatrix *pm;
    Matrix<double> * traits;
};

#endif /* TYPE_LATENTTRAITS_H_ */
