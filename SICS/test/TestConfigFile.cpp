/*
 * TestConfigFile.cpp
 *
 *  Created on: 25 Jul 2014
 *      Author: jlgpisa
 */

#include "TestConfigFile.h"

TestConfigFile::TestConfigFile() {
	header = 0;
	filetype = DATASET;
	separator = ' ';
	beginCol = 0;
}

InTestFiletype TestConfigFile::getFiletype(string fileTypeStr) {
	InTestFiletype ret;

	if (fileTypeStr == TestConstant::kDataset) {
		ret = DATASET;
	}
	else if (fileTypeStr == TestConstant::kInitialValue) {
		ret = INITIAL_VALUE;
	}
	else if (fileTypeStr == TestConstant::kConvergence) {
		ret = CONVERGENCE;
	}
	else if (fileTypeStr == TestConstant::kPob) {
		ret = POB;
	}

	return (ret);
}

InTestFiletype TestConfigFile::getFiletype() const {
	return (filetype);
}

void TestConfigFile::setFiletype(InTestFiletype filetype) {
	this->filetype = filetype;
}

const string& TestConfigFile::getFullPath() const {
	return (fullPath);
}

void TestConfigFile::setFullPath(const string& fullPath) {
	this->fullPath = fullPath;
}

char TestConfigFile::getSeparator() const {
	return (separator);
}

void TestConfigFile::setSeparator(char separator) {
	this->separator = separator;
}

int TestConfigFile::getHeader() const {
	return (header);
}

void TestConfigFile::setHeader(int header) {
	this->header = header;
}

int TestConfigFile::getBeginCol() const {
	return (beginCol);
}

void TestConfigFile::setBeginCol(int beginCol) {
	this->beginCol = beginCol;
}

TestConfigFile::~TestConfigFile() {
	// TODO Auto-generated destructor stub
}


