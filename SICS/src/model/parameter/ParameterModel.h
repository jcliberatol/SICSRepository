/*
 * ParameterModel.h
 *
 *  Created on: 17 Jun 2014
 *      Author: jlgpisa
 */

#ifndef PARAMETERMODEL_H_
#define PARAMETERMODEL_H_

#include <map>
#include <type/Matrix.h>
#include <model/item/ItemModel.h>
#include <model/dimension/DimensionModel.h>

using namespace std;

enum Parameter {a,b,c,d};

class ParameterModel {
protected:
	map <Parameter, Matrix<double> > *parameterSet;
public:
	// Methods
	virtual void buildParameterSet ( ItemModel *, DimensionModel * ) = 0;
	virtual void setInitialPars ( map <Parameter, Matrix<double> > * ) = 0;
	virtual void calculateInitialPars () = 0;
	virtual void successProbability () = 0;

	// Getters and Setters
	virtual const map<Parameter, Matrix<double> >* getParameterSet() const = 0;
	virtual void setParameterSet(map<Parameter, Matrix<double> >*) = 0;

	// Destructor
	virtual ~ParameterModel();

};

#endif /* PARAMETERMODEL_H_ */
