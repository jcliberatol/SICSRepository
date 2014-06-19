/*
 * input.cpp
 *
 *  Created on: May 27, 2014
 *      Author: mirt
 */

#include "Input.h"

Input::Input() {
	// TODO Auto-generated constructor stub

}

Input::~Input() {
	// TODO Auto-generated destructor stub
}


bool Input::importCSV( char* filename, PatternMatrix& M, unsigned int rowIdx, unsigned int colIdx ) {

	ifstream inFile;
	Trace trace("inputErrors.log");
	bool eof=false;

	trace("Input Errors: ");

	inFile.open ( filename, std::ifstream::in );
	if(!inFile.good()) {
		trace( "File does not exists:" );
		trace( filename );
	}

	// currentLine holds the characters of the current line to read
	string currentLine;

	// Header lines are ignored
	for ( unsigned int i=0; i<rowIdx; i++ ) {
		getline ( inFile, currentLine );
	}


	while ( !eof ) {

		getline ( inFile, currentLine );

		eof = inFile.eof();

		const char* processLine = currentLine.c_str();

		//Clean string of the ignored columns
		for (int k = 0; k < colIdx; ++k) {
			processLine = strchr ( processLine, del );
			processLine = &processLine[1]; //Skip one character
		}

		//Count the binary characters
		int i = 0;
		int dlen = 0;

		while( processLine[i]!='\0' ) {
			if( processLine[i]=='0'||processLine[i]=='1' ) dlen ++;
			i++;
		}

		//bitset that holds the bits of a row
		boost::dynamic_bitset<> dset (dlen);

		i = 0;
		int chars = 1;

		while( processLine[i]!= '\0' ) {

			if( processLine[i]=='0') {
				chars++;
			}
			if( processLine[i]=='1') {
				dset.set(dlen-chars);
				chars++;
			}

			i++;
		}
		//Bitset is now filled
		M.push(dset);
	}

	inFile.close();

	return 1;
}

char Input::getDel() const {
	return del;
}

void Input::setDel(char del) {
	this->del = del;
}
