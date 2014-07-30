/*
 * Matrix.h
 *
 *  Created on: May 28, 2014
 *      Author: mirt
 */
#ifndef MATRIX_H_
#define MATRIX_H_

#include <iostream>
#include <stdio.h>
#include <cstring>
using namespace std;

template<typename T>
class Matrix;

template<typename T>
ostream& operator<<(ostream &, Matrix<T> &);

template<class T>
class Matrix {
private:
	int nCol;
	int nRow;
	T *memory;
public:
	Matrix(); //Empty object
	Matrix(int, int); //Two dimensional Matrix Constructor allocates memory
	Matrix(Matrix<T>&); //Copy constructor
	void reset();
	int nR(); // Returns number of rows
	int nC(); //Returns number of columns
	T sum(); // Returns the sum of all objects
	T & operator()(const int nCol, const int nRow); //Accessing operator for a element
	friend ostream& operator<<<T>(ostream &, Matrix<T> &); //Output operator
	//Matrix<T>& operator+=(const Matrix<T>& rhs);
	//Matrix<T> operator+(Matrix<T> lhs , const Matrix<T>rhs);
	virtual ~Matrix();
};

template<class T>
T Matrix<T>::sum() {
	T sum = 0;

	for (int i = 0; i < nRow * nCol; i++) {
		sum += memory[i];
	}

	return (sum);
}

template<class T>
int Matrix<T>::nR() {
	return (nRow);
}

template<class T>
int Matrix<T>::nC() {
	return (nCol);
}

template<class T>
Matrix<T>::Matrix() {
	// TODO Auto-generated constructor stub
	nCol = 0;
	nRow = 0;
	memory = NULL;
}

template<class T>
Matrix<T>::Matrix(Matrix<T>& a) {
	// TODO Auto-generated constructor stub
	nCol = a.nCol;
	nRow = a.nRow;
	memory = new T[nCol * nRow];
	memcpy(memory,a.memory,sizeof(T)*nCol*nRow);
}

template<class T>
Matrix<T>::Matrix(int r, int c) {
	nCol = c;
	nRow = r;
	memory = new T[c * r];
}
template<class T>
T & Matrix<T>::operator()(const int r, const int c) {
	return (memory[nCol * r + c]);
}
template<class T>
Matrix<T>::~Matrix() {
	if (memory != NULL) {
		delete[] memory;
	}
}

template<class T>
void Matrix<T>::reset() {
	memset(memory, 0.0, (nCol * nRow) * sizeof(T));
}

template<class T>
ostream& operator<<(ostream &out, Matrix<T> &M) {
	// Since operator<< is a friend of the Point class, we can access
	// Point's members directly.
	for (int i = 0; i < M.nR(); ++i) {
		for (int j = 0; j < M.nC(); j++) {
			out << M(i, j) << " ";
		}
		out << endl;
	}
	return (out);
}

#endif /* MATRIX_H_ */
