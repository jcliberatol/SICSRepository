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
/**
 * Newton optimizer for generic SICS function, uses the gradient and the hessian matrix
 * for multivariable optimization.
 * In item optimization, the newton method can be parallelized by breaking it into multiple newton methods
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
    delta.reset();
    G.reset();
    H.reset();
    //Pass gradient
    for(int i=0;i<n;++i){
        G(i) = g[i]; //

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
        cout<<"Badly Conditioned hessian matrix, all eigenvalues are near zero or negative"<<endl;//TODO CHANGE TO LOGGER
        return (2);
    }
    //make the deltas
    //cout<<H;

    matrixMultiply(H,G,delta);
    //cout<<"Deltas : "<<endl<<delta<<"--------------------"<<n<<endl;//TODO CHANGE TO LOGGER
    //Substract
    for(int i=0;i<n;++i){
    	//cout<<args[i]<<" ";
        args[i]= args[i]-delta(i);
        //cout<<args[i]<<" : "<<delta(i)<<endl;
        }

    return (inversionStatus);
}

#endif /* NEWTONOPTIMIZER_H_ */
