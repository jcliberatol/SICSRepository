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

void TestReport::setMaxDif(Matrix<double>* maxDif) {
	this->maxDif = maxDif;
}

Trace* TestReport::getTrace() {
	return (trace);
}

void TestReport::setTrace(Trace* trace) {
	this->trace = trace;
}

TestReport::~TestReport() {
	delete trace;
	delete maxDif;
}

