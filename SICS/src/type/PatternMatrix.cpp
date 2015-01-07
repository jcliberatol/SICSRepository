/*
 * PatternMatrix.cpp
 *
 *  Created on: May 30, 2014
 *      Author: mirt
 */

#include "PatternMatrix.h"

PatternMatrix::PatternMatrix(int size) {
	this->size = size;

}

int PatternMatrix::countItems() const {
	//Checks if map is empty if not returns the size of the first bitset
	return ((matrix.empty()) ? 0 : size);
}

int PatternMatrix::countIndividuals() const {

	//map<boost::dynamic_bitset<>, int>::const_iterator it;
	map<bool*, int>::const_iterator it;
	int counter = 0;

	for (it = matrix.begin(); it != matrix.end(); ++it) {
		counter += it->second;
	}

	return (counter);
}

PatternMatrix::~PatternMatrix() {
	// TODO Auto-generated destructor stub
}

//void PatternMatrix::push(boost::dynamic_bitset<> n) {
void PatternMatrix::push(bool* n) {
	matrix[n]++;
}

//void PatternMatrix::push(boost::dynamic_bitset<> n, int k) {
void PatternMatrix::push(bool* n, int k) {
	matrix[n] = matrix[n] + k;
}

void PatternMatrix::flush() {
	matrix.clear();
}

//int & PatternMatrix::operator()(boost::dynamic_bitset<> n) {
//int & PatternMatrix::operator()(bool* n) {
//	return (matrix[n]);
//}

//std::ostream& operator<<(std::ostream & out, PatternMatrix & pm) {
//	for (pm.iterator = pm.matrix.begin(); pm.iterator != pm.matrix.end();
//			++pm.iterator) {
//		for (unsigned int var = 0; var < size; ++var) {
//			int k = pm.iterator->first[var];
//			out << k;
//		}
//		out << " " << pm.iterator->second << std::endl;
//	}
//	return (out);
//}
