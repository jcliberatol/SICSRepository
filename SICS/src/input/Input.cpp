/*
 * input.cpp
 *
 *  Created on: May 27, 2014
 *      Author: mirt
 */

#include "Input.h"

Input::Input() { del = ','; }

Input::~Input() {}

bool Input::importCSV(char* filename, Matrix<double>& M,   int rowIdx,   int colIdx)
{
	bool eof = false;
	int row = 0;
	int col = 0;

	ifstream inFile;

	inFile.open(filename, std::ifstream::in);
	if (!inFile.good())
	{
		cout << "File does not exists:" << filename << endl;
	}

	// currentLine holds the characters of the current line to read
	string currentLine;

	// Header lines are ignored
	for (  int i = 0; i < rowIdx; i++)
	{
		//cout << "Ignored a header line : " << endl;
		getline(inFile, currentLine);
		//cout << currentLine << endl;
	}
	
	while (!eof)
	{
		getline(inFile, currentLine);

		eof = inFile.eof();

		const char* processLine = currentLine.c_str();
		
		if (strlen(processLine) == 0)
		{
			eof = true;
			//cout << "Unproper end of file, read cancelled" << endl;
			break;
		}
		
		//Clean string of the ignored columns
		for (  int k = 0; k < colIdx; ++k)
		{
			processLine = strchr(processLine, del);
			processLine = &processLine[1]; //Skip one character
		}
		/*
		 * Now we must parse the line into the doubles.
		 *for this we take the current address and hold it, then we find the next address and hold it , pass the
		 *double into a string very carefully until the memory addresses are the same
		 */
		col = 0;
		while (strlen(processLine) > 0)
		{
			char * auxptr;
			M(row, col) = strtod(processLine, &auxptr);
			col++;
			processLine = auxptr;
			if (strlen(processLine) > 0)
				processLine = &processLine[1];
		}

		row++;
	}

	inFile.close();
	return (1);
}
/*
 * Imports a CSV file whose elements repeat and generally is composed of only zeroes and ones.
 */
bool Input::importCSV(char* filename, PatternMatrix& M,   int rowIdx,   int colIdx)
{
	ifstream inFile;
	bool eof;
	string currentLine; // currentLine holds the characters of the current line to read
	  int linelen;
	int line;
	eof = false;

	inFile.open(filename, std::ifstream::in);

	if (!inFile.good())
	{
		cout << "File does not exists:" << filename << endl;
	}

	// Header lines are ignored
	for (  int i = 0; i < rowIdx; i++)
		getline(inFile, currentLine);

	line = 0;
	linelen = 0;
	
	while (!eof)
	{
		int i = 0;
		int dlen = 0;
		int chars = 0;

		getline(inFile, currentLine);
		if (line == 0)
			linelen = currentLine.length();
		else
			if (linelen != currentLine.length())
			{
				//cout << "Inconsistent line length , stopped importing" << endl;
				break;
			}

		line++;
		eof = inFile.eof();

		const char* processLine = currentLine.c_str();
		
		//Clean string of the ignored columns
		for (  int k = 0; k < colIdx; ++k)
		{
			processLine = strchr(processLine, del);
			processLine = &processLine[1]; //Skip one character
		}

		//Count the binary characters
		while (processLine[i] != '\0')
		{
			if (processLine[i] == '0' || processLine[i] == '1')
				dlen++;
			i++;
		}

		vector<char> dset(dlen);
		M.size = dlen;
		i = 0;
		chars = 0;

		while (processLine[i] != '\0')
		{
			if (processLine[i] == '0')
				chars++;
			if (processLine[i] == '1')
			{
				dset[chars] = true;
				chars++;
			}

			i++;
		}
		//Bitset is now filled

		M.push(dset);
	}

	inFile.close();

	return (0);
}

char Input::getDel() const { return (del); }

void Input::setDel(char del) { this->del = del; }
