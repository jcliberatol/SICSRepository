/*
 * TestReport.h
 *
 *  Created on: 30 Jul 2014
 *      Author: jlgpisa
 */

#ifndef TESTREPORT_H_
#define TESTREPORT_H_

#include <string>
#include <cmath>

#include <type/Matrix.h>
#include <model/item/ItemModel.h>
#include <model/parameter/ParameterModel.h>
#include <trace/Trace.h>

using namespace std;

class TestReport {
	Trace *trace;
	Matrix<double> * maxDif;
public:
	//Constructor
	TestReport(string);

	// Methods
	void reportDif (Matrix<double> *,  ParameterModel *);
	void reportMaxDif ();

	// Getters and Setters
	Matrix<double>* getMaxDif();
	void setMaxDif(Matrix<double>* maxDif);
	Trace* getTrace();
	void setTrace(Trace* trace);

	// Destructor
	virtual ~TestReport();
};

#endif /* TESTREPORT_H_ */
