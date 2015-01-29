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
	//args = valores a optimizar (en 3pl a,b,c)
	//pars =variables necesarias por la función f pero no optimizables (en 3PL vector f,Matriz r, puntos de cuadrarura)

	//Tolerances to end the algorithm are in constant.h
	//Allocate all needed arrays
	bool accpoint, enough;
	double *g, *t, *X, *c, *B; //g = gradiente actual  , c almacena ultimo gradiente, B almacena matrix hessiana, 
								//t almacena el vector de trabajo de la busqueda lineal
	int count, funcount, gradcount; //conteo de parámetros que no cambian, llanmados a la logverosimilitud y llamados al gradiente
	double f, gradproj; // f almacena el lavor de el llamado a la logverosimilitud, gradproj es la proyección del gradiente
	int i, j, ilast, iter = 0; //ilast registra último paso en el que B fur inicializado a una matriz identidad
	double s, steplength; //steplength es el paso al que se avanza en el sentido de la proyección del gradiente, 
							// por defecto es 1, su redicción está dada por el valor de stepredn =	0.2
	double D1, D2; //D1 se utiliza para evaluar la dirección del gradiente, si es una dirección de ascenso se reinicializa 
					//B a la identidad
	int n, *l;
	double fmin = Constant::INFINITE; //Almacena valor actual de f
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
	B = new double[n * (n + 1) / 2]; //Debido a que B es simetrica solo se usa su parte triangular inferior
	//evaluate the function at the initial points
	/*Variables definidas en el archivo externo SICS/src/type/Constant.cpp	
	double Constant::INFINITE = 1e30;
	double Constant::stepredn =	0.2; constante a que se multimplica por steplength para obtener un nuevo steplength
	double Constant::acctol	=	0.0001; Tolerancia para aceptar un punto
	double Constant::reltest	=	10.0; Es para evaluar la igualdad de los parametros
	double Constant::abstol  =    0.00001;
	double Constant::reltol  =    1e-8;
	*/

	f = (*fntomin)(args, pars, nvars, npars);
	if (!(f < Constant::INFINITE)) {
		return (2); //BAD_INITIAL_VALUES; //Si f no es evaluable se retorna 2
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
		// Se inicializa B a la identidad
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

		//Se construye t_i y se calcula la proyección del gradiente
		//El calculo de t_i se encuentra el la linea siguiente a la ecuacion 2 del documento de Yuli
		for (i = 0; i < n; i++) {
			s = 0.0;
			//Debido a que la matriz B se definió triangular, el cálculo de t_i
			//se particiona en dos sumas, una para cada bloque (triangular inferior y triangular superior)
			for (j = 0; j <= i; j++)
				s -= B[(((i + 1) * (i + 2) / 2) - (i - j) - 1)] * g[l[j]];
			for (j = i + 1; j < n; j++)
				s -= B[(((j + 1) * (j + 2) / 2) - (j - i) - 1)] * g[l[j]];
			t[i] = s;
			gradproj += s * g[l[i]];
		}
		//Si la dirección es de descenso 
		if (gradproj < 0.0) { 
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

				//Si todavía quedan variables por optimizar
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
			
			//Lineas 144 a 152 no hace parte del algoritmo original
			//Si al cambio es muy pequeño se pone count = n e ilast = gradcount
			//ejemplificar caso hipotético con B diferente de la unidad.
			enough = (f > Constant::abstol)
					&& fabs(f - fmin)
							> Constant::reltol
									* (fabs(fmin) + Constant::reltol);
			/* stop if value if small or if relative change is low */
			if (!enough) {
				count = n;
				fmin = f;
			}
			
			//Si todavía quedan variables por optimizar se actualiza el gradiente y la matriz B
			if (count < n) {/* making progress */
				fmin = f;
				(*gradient)(args, pars, nvars, npars, g);
				gradcount++;
				iter++;

				//la actualización de B se da con los parámetros D1 y D2. Página 3 del documento de Yuli 
				D1 = 0.0;
				for (i = 0; i < n; i++) {
					t[i] = steplength * t[i];
					c[i] = g[l[i]] - c[i];
					D1 += t[i] * c[i];
				}
				//si se puede actualizar D
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
				} else { /* D1 < 0 */ //Si B no se puede actualizar se reinicia B if linea 171
					ilast = gradcount;
				}
			} else { /* no progress */ //Si ya no quedan variables por optimizar pero B no se acaba de inicializar if linea 157
				if (ilast < gradcount) { 
					count = 0;
					ilast = gradcount; //En caso de que count = n y ilast = gradcount se rompe el while principal
										//Y se retorna X
				}
			}
		} else { /* uphill search */ //if linea 114
			count = 0;
			if (ilast == gradcount) //termina el programa con el while principal
				count = n;
			else 
				ilast = gradcount; //si no solo reinicia B
			/* Resets unless has just been reset */
		}
		//
		if (iter >= maxiter) //Se rompe whike por iteraciones
			break;
		if (gradcount - ilast > 2 * n)
			ilast = gradcount; /* periodic restart */ //Reinicio periodico de B
	} while (count != n || ilast != gradcount); //En caso de que count = n y ilast = gradcount se rompe el while principal
										//Y se retorna X
	if (iter < maxiter) {
		return (0); //SUCCESS;
	}
	return (3); //N_CONVERGENCE;
}

#endif /* BFGSOPTIMIZER_H_ */
