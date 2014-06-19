/*
 * RaschModel.h
 *
 *  Created on: 17 Jun 2014
 *      Author: jlgpisa
 */

#ifndef RASCHMODEL_H_
#define RASCHMODEL_H_

#include <model/parameter/ParameterModel.h>

class RaschModel: public ParameterModel {
public:
	// Constructor
	RaschModel();

	// Methods
	void buildParameterSet(ItemModel *, DimensionModel *);
	void setInitialPars(map<Parameter, Matrix<double> > *);
	void calculateInitialPars();
	void successProbability();

	// Getters and Setters
	const map<Parameter, Matrix<double> >* getParameterSet() const;
	void setParameterSet(map<Parameter, Matrix<double> >*);

	// Destructor
	virtual ~RaschModel();
};

#endif /* RASCHMODEL_H_ */
