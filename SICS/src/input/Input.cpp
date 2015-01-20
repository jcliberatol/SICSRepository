/*
 * input.cpp
 *
 *  Created on: May 27, 2014
 *      Author: mirt
 */

#include "Input.h"

Input::Input() {
	// TODO Auto-generated constructor stub
	del = ',';
}

Input::~Input() {
	// TODO Auto-generated destructor stub
}
bool Input::importCSV( char* filename, Matrix<double>& M, unsigned int rowIdx, unsigned int colIdx){
	ifstream inFile;
		Trace trace("inputErrors.log");
		bool eof=false;
		int row=0;
		int col=0;

		inFile.open ( filename, std::ifstream::in );
		if(!inFile.good()) {
			trace( "File does not exists:" );
			trace( filename );
		}

		// currentLine holds the characters of the current line to read
		string currentLine;

		// Header lines are ignored
		for ( unsigned int i=0; i<rowIdx; i++ ) {
			trace("Ignored a header line : ");
			getline ( inFile, currentLine );
			trace(currentLine);
		}
		while ( !eof ) {
			getline ( inFile, currentLine );

			eof = inFile.eof();

			const char* processLine = currentLine.c_str();
			if(strlen(processLine)==0){
				eof=true;
				trace("Unproper end of file, read cancelled");
				break;
			}
			//Clean string of the ignored columns
			for (unsigned int k = 0; k < colIdx; ++k) {
				processLine = strchr ( processLine, del );
				processLine = &processLine[1]; //Skip one character
			}
			/*
			 * Now we must parse the line into the doubles.
			 *for this we take the current address and hold it, then we find the next address and hold it , pass the
			 *double into a string very carefully until the memory addresses are the same
			 */
			col = 0;
			while(strlen(processLine)>0){
				char * auxptr;
				M(row,col) = strtod(processLine,&auxptr);
				col++;
				processLine = auxptr;
				if(strlen(processLine)>0){
				processLine = &processLine[1];
				}
			}


			row++;
		}

		inFile.close();
		return (1);
}
/*
 * Imports a CSV file whose elements repeat and generally is composed of only zeroes and ones.
 */
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

	int line = 0;
	int linelen = 0;
	while ( !eof ) {
		getline ( inFile, currentLine );
		if(line==0)
			{
			linelen=currentLine.length();
			}
		else{
			if(linelen!=currentLine.length()){
				trace("Inconsistent line length , stopped importing");
				break;
			}

		}
		line++;
		eof = inFile.eof();

		const char* processLine = currentLine.c_str();
		//Clean string of the ignored columns
		for (unsigned int k = 0; k < colIdx; ++k) {
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

		M.size = dlen;
		//bitset that holds the bits of a row
		//boost::dynamic_bitset<> dset (dlen);S
		//bool * dset = new bool[dlen];
		vector<char> dset(dlen);
		//cout<<&dset[0]<<" "<<dset[0]<<" ";
		i = 0;
		int chars = 1;

		while( processLine[i]!= '\0' ) {

			if( processLine[i]=='0') {
				chars++;
			}
			if( processLine[i]=='1') {
				dset[dlen-chars] = true;
				chars++;
			}

			i++;
		}
		//Bitset is now filled

		M.push(dset);
	}

	inFile.close();

	return (1);
}

char Input::getDel() const {
	return (del);
}

void Input::setDel(char del) {
	this->del = del;
}
