/*
 * BFGSOptimizer.h
 *
 *  Created on: May 29, 2014
 *      Author: mirt
 */

#ifndef BFGSOPTIMIZER_H_
#define BFGSOPTIMIZER_H_
#include <type/Constant.h>
#include <cmath>
#include <stdlib.h>
#include <iostream>
/**
 * Optimizer, takes a function to a pointer that specifies
 * the function to minimize, and the second parameter specifies a function pointer
 * to the gradient of the function, these pointers must be carefully crafted
 * to fit the arguments of the function
 * use a wrapper to call the functions if your function does not follow the form
 * of the input parameters
 * */
int static bfgs(double (*&fntomin)(double*, double*, int, int),
		void (*&gradient)(double*, double*, int, int, double*), double * args,
		double * pars, int nvars, int npars, int maxiter) {
	/*
	 //call of function
	 double result;
	 result = (*fntomin)(args, pars, n);
	 */
	//call of gradient
	/*
	 double result[n];
	 alloc memory in result
	 (*gradient) (args,pars,n,result);
	 */

	//Tolerances to end the algorithm are in constant.h
	//Allocate all needed arrays
	bool accpoint, enough;
	double *g, *t, *X, *c, *B;
	int count, funcount, gradcount;
	double f, gradproj;
	int i, j, ilast, iter = 0;
	double s, steplength;
	double D1, D2;
	int n, *l;
	double fmin = Constant::INFINITE;
	//carefull of negative max iterations
	if (maxiter <= 0) {
		return (1); //MAX_ITER_REACHED;
	}
	//allocate l
	l = (int *) malloc(nvars * sizeof(int));
	n = 0;

	for (i = 0; i < nvars; i++) {
		l[i] = i;
	}
	n = nvars;

	//Assign memory
	g = new double[n];
	t = new double[n];
	X = new double[n];
	c = new double[n];
	//Assign memory to the triangular lower matrix
	B = new double[n * (n + 1) / 2];
	//evaluate the function at the initial points
	f = (*fntomin)(args, pars, nvars, npars);
	if (!(f < Constant::INFINITE)) {
		return (2); //BAD_INITIAL_VALUES;
	}
	//the optimal point for now is f.
	fmin = f;
	funcount = gradcount = 1;
	//run the gradient
	(*gradient)(args, pars, nvars, npars, g);
	iter++;
	ilast = gradcount;

	//Main loop
	do {
		if (ilast == gradcount) {
			for (i = 0; i < n; i++) {
				for (j = 0; j < i; j++)
					B[(((i + 1) * (i + 2) / 2) - (i - j) - 1)] = 0.0;
				B[(((i + 1) * (i + 2) / 2) - (i - i) - 1)] = 1.0;
			}
		}
		for (i = 0; i < n; i++) {
			X[i] = args[l[i]];
			c[i] = g[l[i]];
		}
		gradproj = 0.0;
		for (i = 0; i < n; i++) {
			s = 0.0;
			for (j = 0; j <= i; j++)
				s -= B[(((i + 1) * (i + 2) / 2) - (i - j) - 1)] * g[l[j]];
			for (j = i + 1; j < n; j++)
				s -= B[(((j + 1) * (j + 2) / 2) - (j - i) - 1)] * g[l[j]];
			t[i] = s;
			gradproj += s * g[l[i]];
		}
		if (gradproj < 0.0) { /* search direction is downhill */
			steplength = 1.0;
			accpoint = false;
			do {
				count = 0;
				for (i = 0; i < n; i++) {
					args[l[i]] = X[i] + steplength * t[i];
					if (Constant::reltest + X[i]
							== Constant::reltest + args[l[i]]) { /* no change */
						count++;
					}
				}

//
				if (count < n) {
					f = (*fntomin)(args, pars, nvars, npars);
					funcount++;
					//
					accpoint = (f < Constant::INFINITE)
							&& (f
									<= fmin
											+ gradproj * steplength
													* Constant::acctol);
					if (!accpoint) {
						steplength *= Constant::stepredn;
					}
				}
			} while (!(count == n || accpoint));
			//
			enough = (f > Constant::abstol)
					&& fabs(f - fmin)
							> Constant::reltol
									* (fabs(fmin) + Constant::reltol);
			/* stop if value if small or if relative change is low */
			if (!enough) {
				count = n;
				fmin = f;
			}
			//
			if (count < n) {/* making progress */
				fmin = f;
				(*gradient)(args, pars, nvars, npars, g);
				gradcount++;
				iter++;
				D1 = 0.0;
				for (i = 0; i < n; i++) {
					t[i] = steplength * t[i];
					c[i] = g[l[i]] - c[i];
					D1 += t[i] * c[i];
				}
				if (D1 > 0) {
					D2 = 0.0;
					for (i = 0; i < n; i++) {
						s = 0.0;
						for (j = 0; j <= i; j++)
							s += B[(((i + 1) * (i + 2) / 2) - (i - j) - 1)]
									* c[j];
						for (j = i + 1; j < n; j++)
							s += B[(((j + 1) * (j + 2) / 2) - (j - i) - 1)]
									* c[j];
						X[i] = s;
						D2 += s * c[i];
					}
					D2 = 1.0 + D2 / D1;
					for (i = 0; i < n; i++) {
						for (j = 0; j <= i; j++)
							B[(((i + 1) * (i + 2) / 2) - (i - j) - 1)] += (D2
									* t[i] * t[j] - X[i] * t[j] - t[i] * X[j])
									/ D1;
					}
				} else { /* D1 < 0 */
					ilast = gradcount;
				}
			} else { /* no progress */
				if (ilast < gradcount) {
					count = 0;
					ilast = gradcount;
				}
			}
		} else { /* uphill search */
			count = 0;
			if (ilast == gradcount)
				count = n;
			else
				ilast = gradcount;
			/* Resets unless has just been reset */
		}
		//
		if (iter >= maxiter)
			break;
		if (gradcount - ilast > 2 * n)
			ilast = gradcount; /* periodic restart */
	} while (count != n || ilast != gradcount);
	if (iter < maxiter) {
		return (0); //SUCCESS;

	}
	return (3); //N_CONVERGENCE;
}

#endif /* BFGSOPTIMIZER_H_ */
