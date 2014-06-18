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
	matrix[n]++;
}

void PatternMatrix::push(boost::dynamic_bitset<> n, int k){
	matrix[n]=matrix[n]+k;
}

void PatternMatrix::flush(){
	matrix.clear();
}

long int & PatternMatrix::operator()(boost::dynamic_bitset<> n)
{
    return matrix[n];
}

std::ostream& operator<< (std::ostream & out, PatternMatrix & pm){
	for ( pm.iterator = pm.matrix.begin(); pm.iterator != pm.matrix.end(); ++pm.iterator ) {
			for ( unsigned int var = 0; var < pm.iterator->first.size(); ++var ) {
				int k = pm.iterator->first[var];
				out<<k;
			} out<<" "<<pm.iterator->second<<std::endl;
		}
	return out;
}

int PatternMatrix::getBitsetLength() const {
	//Checks if map is empty if not returns the size of the first bitset
	return ( matrix.empty()) ? 0 : matrix.begin()->first.size();
}
