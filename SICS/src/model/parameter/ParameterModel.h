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
#include <type/DataSet.h>

using namespace std;

enum Parameter {a,b,c,d};

class ParameterModel {
protected:
	Matrix<double> * probabilityMatrix;
	map <Parameter, Matrix<double> * > parameterSet;
public:
	// Methods
	virtual void buildParameterSet ( ItemModel *, DimensionModel * ) = 0;
	virtual void successProbability (DimensionModel *) = 0;


	// Getters and Setters
	virtual map<Parameter, Matrix<double> *> getParameterSet()  = 0;
	virtual void setParameterSet(map<Parameter, Matrix<double> *> parameterSet) = 0;
	virtual double getProbability (int, int) = 0;

	// Destructor
	virtual ~ParameterModel();
};

#endif /* PARAMETERMODEL_H_ */
