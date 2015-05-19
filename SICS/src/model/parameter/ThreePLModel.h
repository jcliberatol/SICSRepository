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
#include <type/QuadratureNodes.h>
/**
 * Model for the 3pl model, uses parameters a d c.
 * unidimensional
 * */
class ThreePLModel: public ParameterModel {

private:
	QuadratureNodes* nodes;
public:
	static double banana(double*, double* , int , int);
	static double Prm(int u,double a, double d, double theta);
	static double Pr(int u,double a, double d, double c, double theta);
	static double* grad(double* xzita,double* xRmat,double* xf, double* ptcuad, double nitems);
	static double Loglik(double* xzita,double* xRmat,double* xf, double* xptcuad, double itemn);

	static double itemLogLik(double*, double* , int, int);
	static void itemGradient(double*, double*, int, int, double*);
	static double itemLogLik2(double*, double* , int, int);
	static void itemGradient2(double*, double*, int, int, double*);
	static void NitemGradient(double*, double*, int, int, double*); 
	static double successProbability(double, double, double, double);
	double successProbability(double, double *);
	static double successProbability_cPrime (double, double, double, double);
	// Constructor
	ThreePLModel();

	// Methods
	void setEstimationNodes(QuadratureNodes * );
	void successProbability(DimensionModel *, QuadratureNodes *);
	static double logLikelihood(double*, double*, int, int);
	static void gradient(double*,double*,int,int,double*);
	// Getters and Setters
	double*** getParameterSet() ;
	void getParameters(double *);
	void setParameters(double *);
	void setParameterSet(double ***);
	void printParameterSet(ostream&);
	string getStringParameters();
	double getProbability(int, int);

	// Destructor
	virtual ~ThreePLModel();
};

#endif /* THREEPLMODEL_H_ */
