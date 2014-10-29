/*
 * TwoPLModel.h
 *
 *  Created on: 18 Jun 2014
 *      Author: jlgpisa
 */

#ifndef TWOPLMODEL_H_
#define TWOPLMODEL_H_

#include <model/parameter/ParameterModel.h>
#include <type/Constant.h>
#include <cmath>
#include <typeinfo>
#include <model/dimension/UnidimensionalModel.h>
#include <model/dimension/MultidimensionalModel.h>
#include <model/dimension/MultiUniDimModel.h>
#include <model/item/PolytomousModel.h>
#include <model/item/DichotomousModel.h>

class TwoPLModel: public ParameterModel {
public:
	// Constructor
	TwoPLModel();

	// Methods
	void buildParameterSet(ItemModel *, DimensionModel *);
	void successProbability(DimensionModel *, QuadratureNodes *);
	static double successProbability(double, double, double);
	static double logLikelihood(double*, double*, int, int);
	static double patternProbability();
	static void gradientAux(long double tp, long double tq,
			long double * grandient);
	static void gradient(double*, double*, int, int, double*);
	static void Ngradient(double* args, double* pars, int nargs, int npars,
			double* gradient);
	static void Hessian(double* args, double* pars, int nargs, int npars,
			double* Hessian);
	static void NHessian(double* args, double* pars, int nargs, int npars,
			double* Hessian);
	static void itemHessian(double* args, double* pars, int nargs, int npars,
			double* Hessian);
	static void itemgradient(double*, double*, int, int, double*);

	// Getters and Setters
	map<Parameter, Matrix<double> *> getParameterSet();
	void setParameterSet(map<Parameter, Matrix<double> *>);
	double getProbability(int, int);

	// Destructor
	virtual ~TwoPLModel();
};

#endif /* TWOPLMODEL_H_ */
