/*
 * TestInput.cpp
 *
 *  Created on: 25 Jul 2014
 *      Author: jlgpisa
 */

#include "TestInput.h"

TestInput::TestInput() {
	separator = ' ';
	headerLines = 0;
	beginCol = 0;
}

int TestInput::getHeaderLines() const {
	return (headerLines);
}

void TestInput::setHeaderLines(int headerLines) {
	this->headerLines = headerLines;
}

const string& TestInput::getPath() const {
	return (path);
}

void TestInput::setPath(const string& path) {
	this->path = path;
}

char TestInput::getSeparator() const {
	return (separator);
}

void TestInput::setSeparator(char separator) {
	this->separator = separator;
}

int TestInput::getBeginCol() const {
	return (beginCol);
}

void TestInput::setConfig(TestConfigFile *testConfigFile) {
	separator = testConfigFile->getSeparator();
	headerLines = testConfigFile->getHeader();
	beginCol = testConfigFile->getBeginCol();
}

Matrix<double>* TestInput::loadMatrix(int rows, int cols) {
	Matrix<double> *matrix = new Matrix<double>(rows, cols);
	Input *input = new Input();

	input->setDel(separator);

	char * filename = const_cast<char *>(path.c_str());

	input->importCSV(filename, *matrix, headerLines, beginCol);

	delete input;
	// TODO: delete filename
	return (matrix);
}

Matrix<double>* TestInput::loadMatrixTransp(int rows, int cols) {
	Matrix<double> *matrix = new Matrix<double>(rows, cols);
	Matrix<double> *matrixT = new Matrix<double>(cols, rows);

	Input *input = new Input();

	input->setDel(separator);

	char * filename = const_cast<char *>(path.c_str());

	input->importCSV(filename, *matrix, headerLines, beginCol);

	for (int i = 0; i < rows; i++) {
		for(int j=0; j<cols; j++) {
			(*matrixT) (j,i) = (*matrix) (i,j);
		}
	}

	delete input;
	delete matrix;
	// TODO: delete filename
	return (matrixT);
}

PatternMatrix* TestInput::loadPattern() {
PatternMatrix *pM = new PatternMatrix();
Input *input = new Input();

input->setDel(separator);

char * filename = const_cast<char *>(path.c_str());

input->importCSV(filename, *pM, headerLines, beginCol);

delete input;
// TODO: delete filename
return (pM);
}

void TestInput::setBeginCol(int beginCol) {
this->beginCol = beginCol;
}

TestInput::~TestInput() {
// TODO Auto-generated destructor stub
}

