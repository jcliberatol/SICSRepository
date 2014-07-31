/*
 * TestInput.h
 *
 *  Created on: 25 Jul 2014
 *      Author: jlgpisa
 */

#ifndef TESTINPUT_H_
#define TESTINPUT_H_

#include <string>
#include <vector>

#include <type/Matrix.h>
#include <type/PatternMatrix.h>
#include <input/Input.h>
#include <../test/TestConfigFile.h>

using namespace std;

class TestInput {
	string path;
	char separator;
	int headerLines;
	int beginCol;
public:
	// Constructor
	TestInput();

	// Methods
	void setConfig (TestConfigFile *);
	Matrix<double> * loadMatrix (int, int);
	Matrix<double> * loadMatrixTransp (int, int);
	PatternMatrix * loadPattern ();

	// Getters and Setters
	const vector<string>& getData() const;
	void setData(const vector<string>& data);
	int getHeaderLines() const;
	void setHeaderLines(int headerLines);
	const string& getPath() const;
	void setPath(const string& path);
	char getSeparator() const;
	void setSeparator(char separator);

	// Destructor
	virtual ~TestInput();
	int getBeginCol() const;
	void setBeginCol(int beginCol);
};

#endif /* TESTINPUT_H_ */
