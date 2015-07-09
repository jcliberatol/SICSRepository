/*
 * Constant.cpp
 *
 *  Created on: 21 Jul 2014
 *      Author: jlgpisa
 */

#include <type/Constant.h>

double Constant::CONVERGENCE_DELTA = 0.0002;
double Constant::NORM_CONST = 1;
double Constant::MAX_EXP = 35.0;
double Constant::INFINITE = 1e30;
double Constant::EPSILON = 1e-20;
double Constant::stepredn =	0.020;
double Constant::acctol	=	0.0001;
double Constant::reltest	=	10.0;
double Constant::abstol  =    0.00001;
double Constant::reltol  =    1e-8;
double Constant::D_CONST = 1;
double Constant::EPSILONC = 10000;
double Constant::LOGLIKO = 0;
int Constant::ITER = 0;
int Constant::MAX_EM_ITERS = 200;
bool Constant::CAPTURE_HESSIANA = true;
string Constant::INITIAL_VALUE_METHOD = "ANDRADE";
double * Constant::BOUNDS = new double[2] { -5, 5 };