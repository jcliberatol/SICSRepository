/*
 * Constant.h
 *
 *  Created on: 21 Jul 2014
 *      Author: jlgpisa
 */

#ifndef CONSTANT_H_
#define CONSTANT_H_

#include <string>

using namespace std;
/**
 * Defines constants used in the SICS library
 * import this class when using a constant
 * TODO : Config file modification for constants
 */
class Constant
{

public:

	enum model_type{RASCH_A1, RASCH_A_CONSTANT, TWO_PL, THREE_PL};
	enum initial_value_type{ANDRADE,RANDOM,OSPINA,FIXED};

	static double CONVERGENCE_DELTA;
	static int MAX_EM_ITERS;
	static double NORM_CONST;
	static double MAX_EXP;
	static double INFINITE;
	static double EPSILON;
	static string INITIAL_VALUE_METHOD;
	//BFGS METHOD SPECIFIC CONSTANTS
	static double stepredn ;
	static double acctol ;
	static double reltest ;
	static double abstol ;
	static double reltol ;
	static double D_CONST ;
	static int ITER;
	static bool CAPTURE_HESSIANA;
	static double EPSILONC;
	static double LOGLIKO;
	static double * BOUNDS;
};

#endif /* CONSTANT_H_ */
