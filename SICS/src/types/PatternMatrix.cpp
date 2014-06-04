/*
 * PatternMatrix.cpp
 *
 *  Created on: May 30, 2014
 *      Author: mirt
 */

#include "PatternMatrix.h"


PatternMatrix::PatternMatrix() {
	// TODO Auto-generated constructor stub

}

PatternMatrix::~PatternMatrix() {
	// TODO Auto-generated destructor stub
}

void PatternMatrix::push(boost::dynamic_bitset<> n){
	bitset[n]++;
}

void PatternMatrix::push(boost::dynamic_bitset<> n, int k){
	bitset[n]=bitset[n]+k;
}

long int & PatternMatrix::operator()(boost::dynamic_bitset<> n)
{
    return bitset[n];
}

std::ostream& operator<< (std::ostream & out, PatternMatrix & pm){
	for ( pm.iterator = pm.bitset.begin(); pm.iterator != pm.bitset.end(); ++pm.iterator ) {
			//cout << iterator->first << " " << iterator->second << endl;
			for ( unsigned int var = 0; var < pm.iterator->first.size(); ++var ) {
				int k = pm.iterator->first[var];
				out<<k;
			} out<<" "<<pm.iterator->second<<std::endl;
		}
	return out;
}
