/*
 * LatentTraits.h
 *
 *  Created on: Feb 2, 2015
 *      Author: cristian
 */

#ifndef TYPE_LATENTTRAITS_H_
#define TYPE_LATENTTRAITS_H_
#include <type/PatternMatrix.h>

class LatentTraits {
public:
	int dim;

	LatentTraits(PatternMatrix * p, const int dims = 1){
		pm = p;
		int rows = pm->matrix.size();
		traits = new Matrix<double>(rows, dims);
		dim = dims;
	}
	virtual ~LatentTraits(){
		delete traits;
	};

	void print(){
		cout<<(*traits);
	}
	PatternMatrix *pm;
	Matrix<double> * traits;

};

#endif /* TYPE_LATENTTRAITS_H_ */
