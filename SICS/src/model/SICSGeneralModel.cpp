/*
 * SICSGeneralModel.cpp
 *
 *  Created on: 17 Jun 2014
 *      Author: jlgpisa
 */

#include <model/SICSGeneralModel.h>

SICSGeneralModel::SICSGeneralModel() {}

ParameterModel* SICSGeneralModel::createParameterModel(int model)
{
	if (model == Constant::THREE_PL)
		return (new ThreePLModel());

	if (model == Constant::TWO_PL)
		return (new TwoPLModel());

	if (model == Constant::ONE_PL)
		return (new OnePLModel());

	if (model == Constant::RASCH)
		return (new OnePLACModel());

	return (new ThreePLModel());
}

ItemModel* SICSGeneralModel::createItemModel()
{
	//TODO FIX FOR DECISION OF NEW MODELS
	return (new DichotomousModel());
}

DimensionModel* SICSGeneralModel::createDimensionModel(int dimstype)
{

	if (dimstype == 1) {
		std::cout<<"Estimating unidimensional model"<<std::endl;
		return (new UnidimensionalModel());
	}
	else if (dimstype == 2){
		std::cout<<"Estimating multidimensional model"<<std::endl;
		return (new MultidimensionalModel());
	}
	else if (dimstype == 3){
		std::cout<<"Estimating Multi-unidimensional model"<<std::endl;
		return ((new MultiUniDimModel()));
	}
	else{
		//Dummy return for wrong values
		std::cout<<"WARNING : WRONG dimstype when creating dimensionalmodel , using UnidimensionalModel"<<std::endl;
		std::cout<<"DIMSTYPE used : "<<dimstype<<std::endl;
		return (new UnidimensionalModel());
	}

}

SICSGeneralModel::~SICSGeneralModel() {}
