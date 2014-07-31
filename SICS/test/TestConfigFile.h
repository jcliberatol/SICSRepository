/*
 * TestConfigFile.h
 *
 *  Created on: 25 Jul 2014
 *      Author: jlgpisa
 */

#ifndef TESTCONFIGFILE_H_
#define TESTCONFIGFILE_H_

#include <string>
#include "TestConstant.h"

using namespace std;

enum InTestFiletype { DATASET, INITIAL_VALUE, CONVERGENCE, POB };

class TestConfigFile {
	InTestFiletype filetype;
	string fullPath;
	int header;
	char separator;
	int beginCol;
public:
	// Constructor
	TestConfigFile();

	// Methods
	static InTestFiletype getFiletype(string);

	// Getters and setters
	InTestFiletype getFiletype() const;
	void setFiletype(InTestFiletype filetype);
	const string& getFullPath() const;
	void setFullPath(const string& fullPath);
	char getSeparator() const;
	void setSeparator(char separator);
	int getHeader() const;
	void setHeader(int header);

	// Destructor
	virtual ~TestConfigFile();
	int getBeginCol() const;
	void setBeginCol(int beginCol);
};

#endif /* TESTCONFIGFILE_H_ */
