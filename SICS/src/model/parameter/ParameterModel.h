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
#include <string>
#include <typeinfo>

using namespace std;

class ParameterModel
{

public:

	double *** parameterSet;
	int items;
	Matrix<double> * probabilityMatrix;
	ItemModel * itemModel;
	DimensionModel * dimensionModel;
	
	ParameterModel() {}

	virtual void buildParameterSet(ItemModel * itemModel, DimensionModel * dimensionModel)
	{
		this->itemModel = itemModel;
		this->dimensionModel = dimensionModel;

		if (typeid(*itemModel) == typeid(DichotomousModel))
		{
			if (typeid(*dimensionModel) == typeid(UnidimensionalModel))
			{
				items = itemModel->countItems();

				parameterSet = new double**[3];
				parameterSet[0] = new double *[1];
				parameterSet[1] = new double *[1];
				parameterSet[2] = new double *[1];

				parameterSet[0][0] = new double[items];
				parameterSet[1][0] = new double[items];
				parameterSet[2][0] = new double[items];
			}
			else if (typeid(*dimensionModel) == typeid(MultidimensionalModel))
			{
				// TODO: Dichotomous Multidimensional
			}
			else if (typeid(*dimensionModel) == typeid(MultiUniDimModel))
			{
				// TODO: Dichotomous MultiUniDimensional
			}

		}
		else if (typeid(*dimensionModel) == typeid(PolytomousModel))
		{
			// TODO: Polytomous Model for Unidimensional, Multidimensional and MultiUni
		}
	};

	virtual void successProbability(DimensionModel *, QuadratureNodes *) = 0;
	virtual double successProbability(double, double*) = 0;

	//Transforms the parameters before starting an estimation process
    virtual void transform() = 0;
    //Transforms back the parameters after estimating them
    virtual void untransform() = 0;

	// Getters and Setters
	virtual double *** getParameterSet() = 0;
	virtual void setParameterSet(double *** parameterSet) = 0;
	virtual void getParameters(double * parameters) = 0;
	virtual void setParameters(double * parameters) = 0;
	virtual double getProbability(int, int) = 0;
	virtual void printParameterSet(ostream&)=0;
	
	// Destructor
	virtual ~ParameterModel()
	{
		//delete probabilityMatrix;
		
		if (typeid(*itemModel) == typeid(DichotomousModel))
		{
			if (typeid(*dimensionModel) == typeid(UnidimensionalModel))
			{
				for(int i = 0; i < 3; i++)
					delete [] parameterSet[i][0];

				for(int i = 0; i < 3; i++)
					delete [] parameterSet[i];

				delete [] parameterSet;
			}
			else if (typeid(*dimensionModel) == typeid(MultidimensionalModel))
			{
				// TODO: Dichotomous Multidimensional
			}
			else if (typeid(*dimensionModel) == typeid(MultiUniDimModel))
			{
				// TODO: Dichotomous MultiUniDimensional
			}

		}
		else if (typeid(*dimensionModel) == typeid(PolytomousModel))
		{
			// TODO: Polytomous Model for Unidimensional, Multidimensional and MultiUni
		}
	}
};

#endif /* PARAMETERMODEL_H_ */
