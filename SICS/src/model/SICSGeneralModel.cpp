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

	if (model == Constant::RASCH_A1)
		return (new OnePLModel());

	if (model == Constant::RASCH_A_CONSTANT)
		return (new OnePLACModel());
	
	return (new ThreePLModel());
}

ItemModel* SICSGeneralModel::createItemModel()
{
	//TODO FIX FOR DECISION OF NEW MODELS
	return (new DichotomousModel());
}

DimensionModel* SICSGeneralModel::createDimensionModel()
{
	//TODO FIX FOR DECISION OF NEW MODELS
	return (new UnidimensionalModel());
}

SICSGeneralModel::~SICSGeneralModel() {}
