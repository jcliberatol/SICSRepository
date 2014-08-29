/*
 * NewtonOptimizer.h
 *
 *  Created on: May 29, 2014
 *      Author: mirt
 */

#ifndef NEWTONOPTIMIZER_H_
#define NEWTONOPTIMIZER_H_

#include <math.h>
#include <type/Matrix.h>
#include <util/blasInterface.h>


//TODO OPTIMIZE NEWTONS USE OF NEWTON SIMPLE METHOD, PASS GRADIENT FOR INDEPENDANT POINTS
//
/*
 * Returns :
 *  0 Success
 *  1 Approximate result, had to use spectral decomposition for inverse
 *  2 Very probably bad conditioned matrix, no eigenvalue is positive , use other optimization method
 */
int static newton(void (*&gradient)(double*, double*, int, int, double*),
        void (*&hessian)(double*, double*, int, int, double*), double * args,
        double * hpars,double * gpars, int nvars, int npars, int maxiter, double * g , double * h){
    //IN the three PL MODEL this newton procedure is called for each item.
    int n = nvars;
    //double *g,*h;
    //g = new double[n];
    //h = new double[n*n];

    //Step one, calculating the gradient of the function
    (*gradient)(args, gpars, nvars, npars, g);
    //Step two , calculating the hessian of the function
    (*hessian)(args, hpars, nvars, npars, h);

    Matrix<double> G(n,1);
    Matrix<double> H(n,n);
    Matrix<double> delta(n,1);

    //Pass gradient
    for(int i=0;i<n;++i){
        G(i,1) = g[i];

        //Pass hessian
        for (int j=0;j<n;++j){
        H(i,j) = h[i*n+j];
        }
    }

    //Set hessian as symm
    //H.symmetry(true);
    H.setSymmetric(true);
    //Invert hessian
    int inversionStatus = 0;
    inversionStatus = ApproximateMatrixInverse(H);
    if (inversionStatus == 3){
        cout<<"Badly Conditioned hessian matrix, all eigenvalues are near zero or negative"<<endl;
        return (2);
    }
    //make the deltas
    matrixMultiply(H,G,delta);
    //Substract
    for(int i=0;i<n;++i){
            args[i]= args[i]-delta(i,1);
        }
    return (inversionStatus);
}

int static multinewton(void (*&gradient)(double*, double*, int, int, double*),
		void (*&hessian)(double*, double*, int, int, double*), double * args,
		double * pars, int nvars, int npars, int maxiter, int hessianSize, int hessianCount) {

	//Steps of the newton rhapson
	int n  = nvars;
	double *g , *h;
	g = new double[n];
	h = new double[hessianSize];//n is the a's plus the d's

	//Step one, calculating the gradient of the function
	(*gradient)(args, pars, nvars, npars, g);
	//Step two , calculating the hessian of the function
	(*hessian)(args, pars, nvars, npars, h);
	//Step three , Inverting the hessian matrix
		//For each mini hessian matrix
		//Break the hessian matrices pointers for each of the individual hessians
	int hSize = (int) sqrt((double)(hessianSize/hessianCount));//In 3pl model this value must count to three
	int invertResponse = 0;
	Matrix<double> m(hSize,hSize); //TODO delete this matrix
	Matrix<double> gr(hSize,1); //First derivative matrix at the point
	Matrix<double> delta(nvars/hessianCount,1);//This is the number of derivative units (nvars/hessiancount)
	//TODO PARALLEL FOR
	/*
	 * @Dependable Matrix structures m, gr, delta (Clone for each)
	 */
	for (int i = 0; i < hessianCount; ++i) {
		//Pass the matrix
		for(int k = 0 ; k < hSize ; ++k){
			for (int j = 0 ; j < hSize ; ++j){
				m(k,j) = h[i*hSize*hSize+k*hSize+j];
			}
		}
		//Step 3.1 Invert Matrix

		 //Invert the matrix here

		//Return to fast matrix representation
		for(int k = 0 ; k < hSize ; ++k){
					for (int j = 0 ; j < hSize ; ++j){
						h[i*hSize*hSize+k*hSize+j] = m(k,j);
					}
				}
		// Pass the call to BFGS
		if(invertResponse == 2) {
				return (2);
			}
		//Step four, Producing and Substracting deltas from the current values
		//Production step, we must multiply the inverse of the hessian matrix by the vector of the gradient.
		//This production step is done item by item to extract each hessian matrix deltae

		//Take the ith hessian inverted matrix (this is m)
		//Take the ith gradient point

		//Multiply them, this yields a 3x1 delta vector for each hessian counter.


			//Substraction step : the deltas are substracted from each of the element array, this is the returned array.
			//Substraction must be done carefully, taking into account the order of the hessian and the gradient vector
	}

	// Pass the Call  to BFGS

return (0);

}
#endif /* NEWTONOPTIMIZER_H_ */
