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

class Constant {
public:
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
};

#endif /* CONSTANT_H_ */
