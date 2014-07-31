/*
 * EMTest.h
 *
 *  Created on: 25 Jul 2014
 *      Author: jlgpisa
 */

#ifndef EMTEST_H_
#define EMTEST_H_

#include <map>
#include <iostream>
#include <fstream>
#include <string>

#include <../test/TestConfigFile.h>
#include <../test/TestInput.h>
#include <../test/TestReport.h>

#include <model/Model.h>
#include <model/ModelFactory.h>
#include <model/SICSGeneralModel.h>
#include <estimation/Estimation.h>
#include <estimation/classical/EMEstimation.h>

#include "TestConstant.h"

using namespace std;

class EMTest {
	string configFile;
	map<InTestFiletype, TestConfigFile *> inputConfig;

	PatternMatrix *pM;
	Matrix<double> *initialValues;
	Matrix<double> *convergence;
	Matrix<double> *pob;

public:
	// Constructor
	EMTest();
	EMTest(string);

	// Methods
	void loadConfiguration();
	void loadInput(map<InTestFiletype, string>);
	void runTest();

	// Getters and Setters
	const string& getConfigFile() const;
	void setConfigFile(const string& configFile);
	map<InTestFiletype, TestConfigFile *> getInputConfig();
	void setInputConfig(map<InTestFiletype, TestConfigFile *> inputConfig);

	// Destructor
	virtual ~EMTest();

};

#endif /* EMTEST_H_ */
