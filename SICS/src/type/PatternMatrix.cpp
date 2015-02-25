/*
 * PatternMatrix.cpp
 *
 *  Created on: May 30, 2014
 *      Author: mirt
 */

#include "PatternMatrix.h"

PatternMatrix::PatternMatrix(int size) {
	this->size = size;
	count_set_bits = NULL;
	bitset_list = NULL;
	frequency_list = NULL;
}

int PatternMatrix::countItems() const {
	return ((matrix.empty()) ? 0 : size);
}
int PatternMatrix::freq(vector<char> bitset) {
	return (matrix[bitset]);
}
int PatternMatrix::countIndividuals() const {

	map<vector<char>, int>::const_iterator it;
	int counter = 0;

	for (it = matrix.begin(); it != matrix.end(); ++it) {
		counter += it->second;
	}

	return (counter);
}

PatternMatrix::~PatternMatrix() {
	// TODO Auto-generated destructor stub
}

void PatternMatrix::push(vector<char> n) {
	matrix[n]++;
}

void PatternMatrix::push(vector<char> n, int k) {
	matrix[n] = matrix[n] + k;
}

void PatternMatrix::flush() {
	matrix.clear();
}

void PatternMatrix::print() {
	for (iterator = matrix.begin(); iterator != matrix.end(); ++iterator) {
		for (unsigned int var = 0; var < size; ++var) {
			int k = iterator->first[var];
			cout << k;
		}
		cout << " " << iterator->second << std::endl;
	}
}

int PatternMatrix::countBitSet(bool * bitset, int index) {
	if (count_set_bits == NULL) {
		count_set_bits = new int[matrix.size()];
		for (int i = 0; i < matrix.size(); i++)
			count_set_bits[i] = -1;
	}

	if (count_set_bits[index] == -1) {
		count_set_bits[index] = 0;
		for (int i = 0; i < size; i++)
			if (bitset[i])
				count_set_bits[index]++;
	}

	return (count_set_bits[index]);
}

bool ** PatternMatrix::getBitsetList() {

	if (bitset_list == NULL) {

		map<vector<char>, int>::const_iterator it;
		map<vector<char>, int>::const_iterator begin = matrix.begin();
		map<vector<char>, int>::const_iterator end = matrix.end();

		bitset_list = new bool*[matrix.size()];

		for (int j = 0; j < matrix.size(); j++) {
			bitset_list[j] = new bool[size];
		}

		frequency_list = new int[matrix.size()];
		//frequency_list = data->getFrequencyList();

		int counter = 0;
		for (it = begin; it != end; ++it, ++counter) {
			copy(it->first.begin(), it->first.end(), bitset_list[counter]);
			frequency_list[counter] = it->second;
		}
	}

	return bitset_list;
}

int * PatternMatrix::getFrequencyList() {

	if (bitset_list == NULL)
		getBitsetList();
	return frequency_list;
}
