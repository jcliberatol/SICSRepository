/*
 * TestReport.cpp
 *
 *  Created on: 30 Jul 2014
 *      Author: jlgpisa
 */

#include "TestReport.h"

TestReport::TestReport(string reportFile) {
	trace = new Trace((const char *) reportFile.c_str());
	maxDif = new Matrix<double> (1,3);
	emEstimation = NULL;
	model = NULL;
}

void TestReport::report(Matrix<double> * convergence,  Matrix<double> * pob) {

	(*trace)(reportName);

	(*trace)("TIME");
	(*trace)(timeSpent);

	(*trace)("ITERATIONS");
	(*trace)(emEstimation->getIterations() );

	(*trace)("MIRT DIFFERENCES");
	reportDif(convergence, model->getParameterModel());

	(*trace)("POBLATIONAL DIFFERENCES:");
	reportDif(pob, model->getParameterModel());

	(*trace)("MIRT - POB DIFFERENCES:");
	reportMatrixDif(pob, convergence);
	(*trace)("------------------------------------------------");

}


void TestReport::reportDif(Matrix<double>* theoric, ParameterModel* paramModel) {
	int items = paramModel->getParameterSet()[a]->nC();
	Matrix<double> diff(3, items);

	for (int i=0;i<items;i++){

		diff(0,i) = (*paramModel->getParameterSet()[a])(0,i);
		diff(1,i) = (*paramModel->getParameterSet()[d])(0,i);
		diff(2,i) = (*paramModel->getParameterSet()[c])(0,i);

		diff(0,i) -= (*theoric)(0,i);
		diff(1,i) -= (*theoric)(1,i);
		diff(2,i) -= (*theoric)(2,i);

		diff(0,i) = abs(diff(0,i));
		diff(1,i) = abs(diff(1,i));
		diff(2,i) = abs(diff(2,i));

		if ( diff(0,i) > (*maxDif)(0,0) ) {
			(*maxDif)(0,0) = diff(0,i);
		}

		if ( diff(1,i) > (*maxDif)(0,1) ) {
			(*maxDif)(0,1) = diff(1,i);
		}

		if ( diff(2,i) > (*maxDif)(0,2) ) {
			(*maxDif)(0,2) = diff(2,i);
		}

	}

	// Matrix of diff is printed
	diff.del = ',';
	(*trace)(diff);
	(*trace)("");
}

void TestReport::reportMatrixDif(Matrix<double>* theoric1, Matrix<double>* theoric2) {
	int items = theoric1->nC();
	Matrix<double> diff(3, items);

	for ( int j=0; j<3; j++ ) {
		for (int i=0;i<items;i++){

			diff(j,i) = (*theoric1)(j,i);
			diff(j,i) -= (*theoric2)(j,i);
			diff(j,i) = abs(diff(j,i));

		}
	}

	// Matrix of diff is printed
	(*trace)(diff);
	(*trace)("");
}


void TestReport::reportMaxDif() {
	(*trace)("\nMaximum differences reported were:");
	(*trace)(*maxDif);
}

Matrix<double>* TestReport::getMaxDif() {
	return (maxDif);
}

void TestReport::startTime() {
	start = clock();
}

void TestReport::endTime() {
	clock_t end = clock();
	clock_t timediff = end - start;
	timeSpent = ((float)timediff) / CLOCKS_PER_SEC;
}

TestReport::~TestReport() {
	delete trace;
	delete maxDif;
	delete model;
	delete emEstimation;
}

void TestReport::setMaxDif(Matrix<double>* maxDif) {
	this->maxDif = maxDif;
}

Trace* TestReport::getTrace() {
	return (trace);
}

void TestReport::setTrace(Trace* trace) {
	this->trace = trace;
}

EMEstimation* TestReport::getEmEstimation() {
	return (emEstimation);
}

void TestReport::setEmEstimation(EMEstimation* emEstimation) {
	this->emEstimation = emEstimation;
}

Model* TestReport::getModel()  {
	return (model);
}

void TestReport::setModel(Model* model) {
	this->model = model;
}

const string& TestReport::getReportName() const {
	return (reportName);
}

void TestReport::setReportName(const string& reportName) {
	this->reportName = reportName;
}

