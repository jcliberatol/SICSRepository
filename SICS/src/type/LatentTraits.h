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

	int dim;

	LatentTraits(PatternMatrix * p, const int dims = 1)
	{
		int rows;
		
		pm = p;
		rows = pm->matrix.size();
		traits = new Matrix<double>(rows, dims);
		dim = dims;
	}

	virtual ~LatentTraits() { delete traits; };

	void print()
	{
		bool ** pattern_list = pm->getBitsetList();

		for(unsigned int i = 0; i < pm->matrix.size(); i++)
		{
			for(int j = 0; j < pm->countItems(); j++)
				cout << pattern_list[i][j] << " ";
			cout<<(*traits)(i,0)<<endl;
		}
	}

	PatternMatrix *pm;
	Matrix<double> * traits;
};

#endif /* TYPE_LATENTTRAITS_H_ */
