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
	LatentTraitSet *latentTraitSet;//Latent traits in use for the estimation , must be moved to estimation procedures
public:
	// Methods
	virtual int getNumDimensions () = 0;
	virtual vector<int> getDimVector() = 0;

	// Getters and Setters
	virtual LatentTraitSet* getLatentTraitSet() const = 0;
	virtual void setLatentTraitSet(LatentTraitSet* latentTraitSet) = 0;//Only settable by estimation methods

	// Destructor
	virtual ~DimensionModel();
};

#endif /* DIMENSIONMODEL_H_ */
