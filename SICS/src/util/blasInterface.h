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
#include <vector>


/*
 * Functions in this file
 * matrixMultiply : doubles
 * matrixInverse : doubles
 */

inline void arrayPrint(double * t, int m){
	for (int var = 0; var < m; ++var) {
		cout<<t[var]<<" ";
	}cout<<endl;
}
inline void arrayPrint(long double * t, int m){
	for (int var = 0; var < m; ++var) {
		cout<<t[var]<<" ";
	}cout<<endl;
}

//Double general matrix multiplication C = A.B
inline int matrixMultiply(Matrix<double> &A , Matrix<double> &B , Matrix<double> &C){
	//Operation C = A.B
	//Determine Sizes
	//A is m by k , B is k by n so C will be m by n
	int k,m,n;
	m = A.nR();
	n = B.nC();
	k = A.nC();
	//cout<<"Multiplying 2 matrices of "<<m<<" by "<<k<<" and "<<B.nR()<<" by "<<n<<endl<<"with result of "<<C.nR()<<" by "<<C.nC()<<endl;
	//Check for bad conditioned C matrix
	if(C.nC()!= n or C.nR()!= m){
		cout<<"BAD CONDS";//TODO CHANGE TO LOGGER
		return (1); // Badly conditioned Output Matrix
	}
	double alpha = 1.0, beta = 0.0;
	//Determine transposition of each matrix A and B
	CBLAS_TRANSPOSE at,bt;
	at = CblasNoTrans;
	bt=at;
	if(A.transposed){
		at = CblasTrans;
		//cout<<"atrans"<<endl;
	}
	if(B.transposed){
		bt = CblasTrans;
		//cout<<"b trans"<<endl;
	}
	//cout<<"Leading dims"<<k<<" "<<n<<" "<<n<<" "<<endl;
	cblas_dgemm(CblasRowMajor,at,bt,m,n,k,alpha,A.memory,A.ld,B.memory,B.ld,beta,C.memory,C.ld);
	//cout<<"Matrice of output : "<<C.nR()<<" "<<C.nC()<<endl<<C<<endl;
	return(0);
}


inline int ApproximateMatrixInverse(Matrix<double> &M){
	//Backup M
	Matrix<double> cM(M);
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
	double vl = 1e-100;
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
	//cout<<"Ready to start decompose"<<endl;
	//First queries the function for the size of the array that are used to work
	dsyevr_(&jobs,&range,&uplo,&n,M.memory,&n,&vl,&vu,&ilu,&ilu,&abstol,&m,evals.memory,evecs.memory,&n,isuppz.memory,&flw,&lwork,&fliw,&liwork,&info);
	//cout<<endl<<"work size and lwork size : "<<flw<<" "<<fliw<<endl;
	//Now create memory for the work spaces
	lwork = (int)flw;
	liwork = fliw;
	Matrix<double> work(1,lwork);
	Matrix<int> iwork(1,liwork);
	//Call the dsyevr procedure
	dsyevr_(&jobs,&range,&uplo,&n,M.memory,&n,&vl,&vu,&ilu,&ilu,&abstol,&m,evals.memory,evecs.memory,&n,isuppz.memory,work.memory,&lwork,iwork.memory,&liwork,&info);
	//cout<<"Matrix after eigendecompose"<<M<<endl;
	//cout<<"Eigenvalues"<<endl;//TODO CHANGE TO LOGGER
	//cout<<"N : "<<n<<"  "<<"W : "<<endl<<evals<<endl<<evecs<<endl;//TODO CHANGE TO LOGGER

	//Desde aqui empezo jose
	vector <double> eigenValV;
	Matrix<double> eigenVectV (evecs.nR(),evecs.nC());
	for (int i=0; i<evals.nC(); i++ ){
		if (evals(i)!=0) {
			eigenValV.push_back(evals(i));
			for (int j=0;j<evecs.nC();j++){
				eigenVectV(i,j) = evecs(i,j);
			}
		}
	}
	// If all eigens values are zeroes, function is finished
	if (eigenValV.size() == 0) {
	 //restore original matrix
		memcpy(M.memory,cM.memory,sizeof(double)*M.nC()*M.nR());
		return (3);
	}

	// final Eigen vectors and eigen values
	Matrix<double> eigenvalues (1, eigenValV.size());
	Matrix<double> eigenvectors ((int)eigenValV.size(), n);
	Matrix<double> identity ('I', eigenValV.size());

	for (int i=0; i<eigenValV.size(); i++ ){
		identity(i,i) = 1/eigenValV[i];
		for (int j=0; j<eigenvectors.nC();j++){
					eigenvectors(i,j) = eigenVectV(i,j);
				}
	}

	// Desde aqui son las dudas
	eigenvectors.transpose();
	//cout << "eigenvectors\n" << eigenvectors;
	//cout << "identity\n" << identity;
	Matrix<double> head (eigenvectors.nR(),identity.nC());
	matrixMultiply(eigenvectors, identity, head);
	//cout << "head\n" << head;
	eigenvectors.transpose();
	Matrix<double> inverse (eigenvectors.nC(),eigenvectors.nC());
	matrixMultiply(head,eigenvectors, inverse);
	//cout << "eigenvectors\n" << eigenvectors;
	//cout << "eigenInv\n" << inverse;
	M.copy(inverse);
	int x;
	return (0);
}





#endif /* BLASINTERFACE_H_ */
