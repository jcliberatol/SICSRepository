/*
 * Matrix.h
 *
 *  Created on: May 28, 2014
 *      Author: mirt
 */
#ifndef MATRIX_H_
#define MATRIX_H_

#include <iostream>
using namespace std;

template <typename T>
class Matrix;

template <typename T>
ostream& operator<< (ostream &, Matrix<T> &);

template <class T>
class Matrix {
private:
	int nCol = 0;
	int nRow = 0;
	T *memory;
public:
	Matrix();//Empty object
	Matrix(int,int);//Two dimensional Matrix Constructor allocates memory
	Matrix(Matrix&);//Copy constructor
	int nR(); // Returns number of rows
	int nC(); //Returns number of columns
	double& operator()(const int nCol, const int nRow); //Accessing operator for a element
	double& operator()(const int nCol);//Accessing operator for a row
	friend ostream& operator<< <T>(ostream &, Matrix<T> &);//Output operator


	virtual ~Matrix();
};

template <class T>
int Matrix<T>::nR(){
	return nRow;
}

template <class T>
int Matrix<T>::nC(){
	return nCol;
}

template <class T>
Matrix<T>::Matrix() {
	// TODO Auto-generated constructor stub
	nCol = 0;
	nRow = 0;
	memory = NULL;
}

template <class T>
Matrix<T>::Matrix(int r, int c){
	nCol = c;
	nRow = r;
	memory = new T [c*r];
}
template <class T>
double & Matrix<T>::operator()(const int c, const int r)
{
    return memory[nRow][nCol];
}
template <class T>
Matrix<T>::~Matrix() {
	// TODO Auto-generated destructor stub
	delete[] memory;
}
template <class T>
ostream& operator<< (ostream &out, Matrix<T> &M)
{
    // Since operator<< is a friend of the Point class, we can access
    // Point's members directly.
    for (int i = 0; i < M.nR(); ++i) {
		for (int j = 0 ; j < M.nC(); j++){
			out<<M(i,j)<<" ";
		}
		out<<endl;
	}
    return out;
}

#endif /* MATRIX_H_ */
