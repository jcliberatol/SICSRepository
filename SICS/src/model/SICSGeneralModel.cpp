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
	//TODO FIX FOR DECISION OF NEW MODELS
	if (dims == 1) {
		return (new UnidimensionalModel());
	}
	if (dims == 2){
		return (new MultidimensionalModel());
	}
	if (dimstype == 3){
		return ((new MultiUniDimModel()));
	}
}

SICSGeneralModel::~SICSGeneralModel() {}
