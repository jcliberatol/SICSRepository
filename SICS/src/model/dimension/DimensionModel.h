/*
 * DimensionModel.h
 *
 *  Created on: 17 Jun 2014
 *      Author: jlgpisa
 */

#ifndef DIMENSIONMODEL_H_
#define DIMENSIONMODEL_H_

#include <vector>
#include <type/LatentTraitSet.h>

using namespace std;

class DimensionModel {
protected:
	LatentTraitSet *latentTraitSet;
public:
	// Methods
	virtual int getNumDimensions () = 0;
	virtual vector<double> getDimVector() = 0;

	// Getters and Setters
	virtual const LatentTraitSet* getLatentTraitSet() const = 0;
	virtual void setLatentTraitSet(LatentTraitSet* latentTraitSet) = 0;

	// Destructor
	virtual ~DimensionModel();
};

#endif /* DIMENSIONMODEL_H_ */
