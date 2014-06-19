/*
 * ThreePLModel.h
 *
 *  Created on: 17 Jun 2014
 *      Author: jlgpisa
 */

#ifndef THREEPLMODEL_H_
#define THREEPLMODEL_H_

#include <model/parameter/ParameterModel.h>

class ThreePLModel : public ParameterModel {
public:
	// Constructor
	ThreePLModel();

	// Methods
	void buildParameterSet ( ItemModel *, DimensionModel * );
	void setInitialPars ( map <Parameter, Matrix<double> > * );
	void calculateInitialPars ();
	void successProbability ();

	// Getters and Setters
	const map<Parameter, Matrix<double> >* getParameterSet() const;
	void setParameterSet(map<Parameter, Matrix<double> >*);

	virtual ~ThreePLModel();
};

#endif /* THREEPLMODEL_H_ */
