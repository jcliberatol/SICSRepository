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
#include <ctime>

#include <type/Matrix.h>
#include <trace/Trace.h>
#include <model/Model.h>
#include <estimation/classical/EMEstimation.h>

#include "json.h"

using namespace std;

using namespace std;

class TestReport {
	Model *model;
	EMEstimation *emEstimation;
	Trace *trace;
	Matrix<double> * maxDif;
	string reportName;

	json::Object jsonReport;
	json::Array jsonResults;

	clock_t start;
	double timeSpent;
public:
	//Constructor
	TestReport(string);

	// Methods
	void report (Matrix<double> *,  Matrix<double> *);
	void addJsonResult (Matrix<double> *,  Matrix<double> *);
	void reportJSON ();
	Matrix<double> * reportDif (Matrix<double> *,  ParameterModel *);
	Matrix<double> * reportMatrixDif (Matrix<double> *,  Matrix<double> *);
	void reportMaxDif ();
	void startTime();
	void endTime();

	// Getters and Setters
	Matrix<double>* getMaxDif();
	void setMaxDif(Matrix<double>* maxDif);
	Trace* getTrace();
	void setTrace(Trace* trace);
	EMEstimation* getEmEstimation();
	void setEmEstimation(EMEstimation* emEstimation);
	Model* getModel();
	void setModel(Model* model);
	const string& getReportName() const;
	void setReportName(const string& reportName);

	// Destructor
	virtual ~TestReport();

};

#endif /* TESTREPORT_H_ */
