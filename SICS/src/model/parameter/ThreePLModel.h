/*
 * ThreePLModel.h
 *
 *  Created on: 17 Jun 2014
 *      Author: jlgpisa
 */

#ifndef THREEPLMODEL_H_
#define THREEPLMODEL_H_

#include <typeinfo>
#include <model/parameter/ParameterModel.h>
#include <model/item/DichotomousModel.h>
#include <model/item/PolytomousModel.h>
#include <model/dimension/UnidimensionalModel.h>
#include <model/dimension/MultidimensionalModel.h>
#include <model/dimension/MultiUniDimModel.h>
#include <type/PatternMatrix.h>
#include <type/Constant.h>
#include <cmath>

class ThreePLModel: public ParameterModel {
public:

	static double successProbability(double, double, double, double);
	static double successProbability_cPrime (double, double, double, double);
	// Constructor
	ThreePLModel();

	// Methods
	void buildParameterSet(ItemModel *, DimensionModel *);
	void successProbability(DimensionModel *);
	static double logLikelihood(double*, double*, int, int);
	static void gradient(double*,double*,int,int,double*);
	static void Ngradient(double* args, double* pars, int nargs, int npars, double* gradient);
	static void Hessian(double* args, double* pars, int nargs, int npars, double* Hessian);
	static void NHessian(double* args, double* pars, int nargs, int npars, double* Hessian);
	static void itemHessian(double* args, double* pars, int nargs, int npars, double* Hessian);
	static void itemgradient(double*,double*,int,int,double*);
	// Getters and Setters
	map<Parameter, Matrix<double> *> getParameterSet() ;
	void setParameterSet(map<Parameter, Matrix<double> *>);
	double getProbability(int, int);

	// Destructor
	virtual ~ThreePLModel();
};

#endif /* THREEPLMODEL_H_ */
