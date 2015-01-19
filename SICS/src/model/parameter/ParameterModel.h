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
#include <type/QuadratureNodes.h>

using namespace std;


class ParameterModel {
protected:

	//map <Parameter, Matrix<double> * > parameterSet;
	double *** parameterSet;
public:
	int items;
	Matrix<double> * probabilityMatrix;
	// Methods
	virtual void buildParameterSet ( ItemModel *, DimensionModel *) = 0;
	virtual void successProbability (DimensionModel *, QuadratureNodes *) = 0;


	// Getters and Setters
	virtual double *** getParameterSet()  = 0;
	virtual void setParameterSet(double *** parameterSet) = 0;
	virtual double getProbability (int, int) = 0;
	virtual void printParameterSet(ostream&)=0;

	// Destructor
	virtual ~ParameterModel();
};

#endif /* PARAMETERMODEL_H_ */
