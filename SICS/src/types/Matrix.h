/*
 * Matrix.h
 *
 *  Created on: May 28, 2014
 *      Author: mirt
 */

#ifndef MATRIX_H_
#define MATRIX_H_

class Matrix {
public:
	Matrix();//Empty object
	Matrix(int,int);//Dim * Dim
	Matrix(Matrix&);


	virtual ~Matrix();
};

#endif /* MATRIX_H_ */
