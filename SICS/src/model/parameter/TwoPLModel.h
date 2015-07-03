/*
 * TwoPLModel.h
 *
 *  Created on: 18 Jun 2014
 *      Author: cesandovalp
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

class TwoPLModel: public ParameterModel
{

public:

	TwoPLModel();

	// Methods
	inline void successProbability(DimensionModel *, QuadratureNodes *);
	inline static double successProbability(double , double , double );
	double successProbability(double, double*);
	static double logLikelihood(double*, double*, int, int);
	static double patternProbability();
	static void gradientAux(long double, long double, long double * );
	static void gradient(double*, double*, int, int, double*);
	static void Ngradient(double* , double*, int, int, double*);
	static void Hessian(double*, double*, int, int, double*);
	static void NHessian(double*, double*, int, int, double*);
	static void itemHessian(double*, double*, int, int, double*);
	static void itemgradient(double*, double*, int, int, double*);
	static void itemGradient (double*, double*, int, int, double*);
	static double itemLogLik (double*, double*, int, int);

	double *** getParameterSet();
	double getProbability(int, int);
	void getParameters(double *);
	void setParameterSet(double***);
	void setParameters(double *);
	void printParameterSet(ostream&);

	virtual ~TwoPLModel();
};

#endif /* TWOPLMODEL_H_ */
