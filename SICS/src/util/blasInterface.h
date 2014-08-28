/*
 * blasInterface.h
 *
 *  Created on: Aug 21, 2014
 *      Author: jliberato
 */

#ifndef BLASINTERFACE_H_
#define BLASINTERFACE_H_
#include "type/Matrix.h"
#include <openblas/cblas.h>
#include <openblas/lapacke/lapacke.h>

/*
 * Functions in this file
 * matrixMultiply : doubles
 * matrixInverse : doubles
 */



//Double general matrix multiplication C = A.B
int matrixMultiply(Matrix<double> A , Matrix<double> B , Matrix<double> &C){
	//Operation C = A.B
	//Determine Sizes
	//A is m by k , B is k by n so C will be m by n
	int k,m,n;
	m = A.nR();
	n = B.nC();
	k = A.nC();
	//Check for bad conditioned C matrix
	if(C.nC()!= n or C.nR()!= m){
		cout<<"BAD CONDS";
		return (1); // Badly conditioned Output Matrix
	}
	int alpha = 1;
	//Determine transposition of each matrix A and B
	CBLAS_TRANSPOSE at,bt;
	at = CblasNoTrans;
	bt=at;
	if(A.transposed){
		at = CblasTrans;
	}
	if(B.transposed){
		bt = CblasTrans;
	}
	cblas_dgemm(CblasRowMajor,at,bt,m,n,k,alpha,A.memory,k,B.memory,n,alpha,C.memory,n);
	return(0);
}


int ApproximateMatrixInverse(Matrix<double> M){
	//This procedure relies heavily on blas to determine if any eigenvalues are negative or near zero
	//For the eigenvalues that arent  between the boundaries, the spectral composition is not performed

	//First step, create symetrical array
	/*
	 * Calling arguments explanation
	 */
	char jobs = 'V'; //compute both eigenvalues and eigenvectors
	char range = 'V'; //Only find the eigenvalues and eigenvector in the interval given
	char uplo = 'L'; //Lower representation of the symm matrix
	int n = M.nC();
	//a is M memory
	//lda is n
	//vl and vu lower and uper bounds of the interval
	double vl = 1e-30;
	double vu = 1e100;
	double abstol = 1e-10;
	//il and iu are not used since our search is interval based
	int ilu = 0;
	int m = 0 ; //Place holder for the number of eigenvalues found
	//It seems better to create matrix structures for eigenvalues and eigenvectors
	Matrix<double> evals(1,n);
	Matrix<double> evecs(n,n);
	evals.reset();evecs.reset();
	Matrix<int> isuppz(1,n);
	isuppz.reset();//Indices of the found eigenvectors
	int info=0; //information about convergence
	int lwork = -1;
	int liwork = -1;
	double flw = 0;
	int fliw = 0;
	int dummy=0;//Dummy integer
	/*
	 * int dsyevr_(char *jobz, char *range, char *uplo, integer *n,
	doublereal *a, integer *lda, doublereal *vl, doublereal *vu, integer *
	il, integer *iu, doublereal *abstol, integer *m, doublereal *w,
	doublereal *z__, integer *ldz, integer *isuppz, doublereal *work,
	integer *lwork, integer *iwork, integer *liwork, integer *info)
	 */
	cout<<"Ready to start decompose"<<endl;
	//First queries the function for the size of the array that are used to work
	dsyevr_(&jobs,&range,&uplo,&n,M.memory,&n,&vl,&vu,&ilu,&ilu,&abstol,&m,evals.memory,evecs.memory,&n,isuppz.memory,&flw,&lwork,&fliw,&liwork,&info);
	cout<<endl<<"work size and lwork size : "<<flw<<" "<<fliw<<endl;
	//Now create memory for the work spaces
	lwork = (int)flw;
	liwork = fliw;
	Matrix<double> work(1,lwork);
	Matrix<int> iwork(1,liwork);
	//Call the dsyevr procedure
	dsyevr_(&jobs,&range,&uplo,&n,M.memory,&n,&vl,&vu,&ilu,&ilu,&abstol,&m,evals.memory,evecs.memory,&n,isuppz.memory,work.memory,&lwork,iwork.memory,&liwork,&info);
	cout<<"Matrice after eigendecompose"<<M<<endl;
	cout<<"Eigenvalues"<<endl;
	cout<<"N : "<<n<<"  "<<"W : "<<endl<<evals<<endl<<evecs<<endl;

	//Desde aqui empezo jose
	vector <double> eigenValV;
	vector <vector <double> > eigenVectV;
	for (int i=0; i<evals.nC(); i++ ){
		if (evals(0,i)!=0) {
			eigenValV.push_back(evals(0,i));

			vector<double> eVTemp;
			for (int j=0;j<evecs.nC();j++){
				eVTemp.push_back(evecs(i,j));
			}
			eigenVectV.push_back(eVTemp);
		}
	}

	// If all eigens values are zeroes, function is finished
	if (eigenValV.size() == 0) {
		return (3);
	}

	// final Eigen vectors and eigen values
	Matrix<double> eigenvalues (1, eigenValV.size());
	Matrix<double> eigenvectors ((int)eigenValV.size(), 3);
	Matrix<double> identity ('I', eigenValV.size());


	for (int i=0; i<eigenValV.size(); i++ ){
		identity(i,i) = 1/eigenValV[i];
		for (int j=0; j<eigenvectors.nC();j++){
					eigenvectors(i,j) = 1/eigenVectV[i][j];
				}
	}
	Matrix<double> tail (identity.nR(),eigenvectors.nC());
	// Desde aqui son las dudas

	matrixMultiply(identity, eigenvectors, tail);
	Matrix<double> inverse (eigenvectors.nC(),eigenvectors.nC());
	eigenvectors.transpose();
	matrixMultiply(eigenvectors, tail, inverse);
	cout << "identity\n" << identity;
	cout << "eigenvectors\n" << eigenvectors;
	cout << "eigenInv\n" << inverse;


	return (0);
}





#endif /* BLASINTERFACE_H_ */
