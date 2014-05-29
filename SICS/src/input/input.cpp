/*
 * input.cpp
 *
 *  Created on: May 27, 2014
 *      Author: mirt
 */

#include "input.h"

input::input() {
	// TODO Auto-generated constructor stub

}

input::~input() {
	// TODO Auto-generated destructor stub
}
/*
bool input::importCSV(char* p, GeMatrix<FullStorage<int, RowMajor> >&mat, int xOff, int yOff){

	bool retval = true;
	if (!exists(p))retval = false;

	FILE *ff = freopen (p, "r", stdin );

		char eol = 10;
		char c;
		unsigned int j = 0; 		//index of row
		unsigned int i = 0; 		//index of col

		while ( !feof(ff) ) {

			//first  character in line is read
			c = getchar();

			// Second while loops over all characters in the current line
			while ( c != eol && !feof(ff) ) {
				if ( isdigit(c) ) {
					//Is a number
					i++;
				}
				else if ( i >= xOff && j >= yOff ) {
					//Character is processed.
					mat(i,j) = c - '0';
				}
				//Next character is read
				c = getchar();
			}
			j++;
			i = 0;

		}
		fclose ( ff );
	return retval;
}
*/
