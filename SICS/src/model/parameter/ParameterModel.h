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
#include <model/item/DichotomousModel.h>
#include <model/item/PolytomousModel.h>
#include <model/dimension/UnidimensionalModel.h>
#include <model/dimension/MultidimensionalModel.h>
#include <model/dimension/MultiUniDimModel.h>
#include <type/DataSet.h>
#include <type/QuadratureNodes.h>
#include <trace/Trace.h>
#include <string>
#include <typeinfo>

using namespace std;

class ParameterModel {
public:
	//map <Parameter, Matrix<double> * > parameterSet;
	double *** parameterSet;
	int items;
	Trace* profiler = 0;
	Matrix<double> * probabilityMatrix;
	// Methods
	virtual void buildParameterSet(ItemModel * itemModel,
			DimensionModel * dimensionModel) {
		if (typeid(*itemModel) == typeid(DichotomousModel)) {

			if (typeid(*dimensionModel) == typeid(UnidimensionalModel)) {

				items = itemModel->countItems();

				parameterSet = new double**[3];
				parameterSet[0] = new double *[1];
				parameterSet[1] = new double *[1];
				parameterSet[2] = new double *[1];

				parameterSet[0][0] = new double[items];
				parameterSet[1][0] = new double[items];
				parameterSet[2][0] = new double[items];

			}

			else if (typeid(*dimensionModel) == typeid(MultidimensionalModel)) {
				// TODO: Dichotomous Multidimensional
			}

			else if (typeid(*dimensionModel) == typeid(MultiUniDimModel)) {
				// TODO: Dichotomous MultiUniDimensional
			}

		}

		else if (typeid(*dimensionModel) == typeid(PolytomousModel)) {
			// TODO: Polytomous Model for Unidimensional, Multidimensional and MultiUni
		}
	}
	;
	virtual void successProbability(DimensionModel *, QuadratureNodes *) = 0;
	virtual double successProbability(double, double*) = 0;

	// Getters and Setters
	virtual double *** getParameterSet() = 0;
	virtual void setParameterSet(double *** parameterSet) = 0;
	virtual void getParameters(double * parameters) = 0;
	virtual double getProbability(int, int) = 0;
	virtual void printParameterSet(ostream&)=0;
	void setProfiler(Trace* t) {
		profiler = t;
	}
	virtual string getStringParameters() = 0;
	// Destructor
	virtual ~ParameterModel();
};

#endif /* PARAMETERMODEL_H_ */
