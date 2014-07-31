/*
 * TestConstant.h
 *
 *  Created on: 25 Jul 2014
 *      Author: jlgpisa
 */

#ifndef TESTCONSTANT_H_
#define TESTCONSTANT_H_

#include <string>

using namespace std;

class TestConstant {
public:
	/**
	 * Configuration Files constants
	 */
	static string kNewConf;
	static string kFileType;
	static string kFullPath;
	static string kHeader;
	static string kSep;
	static string kNodeCount;

	/**
	 * Input File types
	 */
	static string kDataset;
	static string kInitialValue;
	static string kConvergence;
	static string kPob;

	/**
	 * Thetas file constants
	 */
	static int q;
};

#endif /* TESTCONSTANT_H_ */
