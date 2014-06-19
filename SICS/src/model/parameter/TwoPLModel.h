/*
 * TwoPLModel.h
 *
 *  Created on: 18 Jun 2014
 *      Author: jlgpisa
 */

#ifndef TWOPLMODEL_H_
#define TWOPLMODEL_H_

#include <model/parameter/ParameterModel.h>

class TwoPLModel: public ParameterModel {
public:
	// Constructor
	TwoPLModel();

	// Methods
	void buildParameterSet(ItemModel *, DimensionModel *);
	void setInitialPars(map<Parameter, Matrix<double> > *);
	void calculateInitialPars();
	void successProbability();

	// Getters and Setters
	const map<Parameter, Matrix<double> >* getParameterSet() const;
	void setParameterSet(map<Parameter, Matrix<double> >*);

	// Destructor
	virtual ~TwoPLModel();
};

#endif /* TWOPLMODEL_H_ */
