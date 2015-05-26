#bfgs

Esta es la guia que muestra la implementacion de los criterios de wolfe en el bfgs.

##Criterios de wolfe

###Criterio de Wolfe f < f(x_k) + c1 * ak * pk * grad(f)
```
if (count < n) {
					f = (*fntomin)(args, pars, nvars, npars);
					funcount++;
					//Wolfe condition 1
					// f < f(x_k) + c1 * ak * pk * grad(f)
					// Where f(x_k) is fmin
					// c1 is acctol
					// ak is steplength
					// grad(f) is gradproj
					accpoint = (f < Constant::INFINITE)
							&& (f
									<= fmin
											+ gradproj * steplength
													* Constant::acctol);
					if (!accpoint) {
						// Reduce the steplength by the factor if wolfe condition is not met.
						steplength *= Constant::stepredn;
					}
```

### Actualizacion del B con d1 y d2

```
//la actualización de B se da con los parámetros D1 y D2. Página 3 del documento de Yuli 
				D1 = 0.0; //Paso 13 del algoritmo de nash, prepararse para la actualizacion de la matriz.
				for (i = 0; i < n; i++) {
					t[i] = steplength * t[i];
					c[i] = g[l[i]] - c[i];
					D1 += t[i] * c[i]; //Computa el producto interno t*y
				}
				//si se puede actualizar (Condicion)
				if (D1 > 0) {
					D2 = 0.0;// Calculo de D2 para actualizacion del B (Hessiana inv.)
					// d2 = (1+yT By/d1)d1
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
					//Actualizacion de la matriz B.
					for (i = 0; i < n; i++) {
						for (j = 0; j <= i; j++)
							B[(((i + 1) * (i + 2) / 2) - (i - j) - 1)] += (D2
									* t[i] * t[j] - X[i] * t[j] - t[i] * X[j])
									/ D1;
```
##Criterios de parada

### Caso 1 : Direccion hacia arriba
```
				
	//Case 1 for convergence criteria.
		} else { /* uphill search */ //if linea 114
			count = 0;
			if (ilast == gradcount) //termina el programa con el while principal
				count = n;
			else 
				ilast = gradcount; //si no solo reinicia B
			/* Resets unless has just been reset */
		}
		//
		if (iter >= maxiter) //Se rompe while por iteraciones
			break;
		if (gradcount - ilast > 2 * n)
			ilast = gradcount; /* periodic restart */ //Reinicio periodico de B
	} while (count != n || ilast != gradcount); //En caso de que count = n y ilast = gradcount se rompe el while principal
										//Y se retorna X
	if (iter < maxiter) {
		return (0); //SUCCESS; //Termino satisfactoriamente
	}
	return (3); //N_CONVERGENCE; //BFGS agoto iteraciones.


			
```
###Caso 2 :Cambio pequeño.
```

//Condicion de convergencia 2 :  Cambio pequeño

//Si al cambio es muy pequeño se pone count = n e ilast = gradcount
			//ejemplificar caso hipotético con B diferente de la unidad.
			enough = (f > Constant::abstol)
			&& fabs(f - fmin) > Constant::reltol
				* (fabs(fmin) + Constant::reltol);
			/* stop if value if small or if relative change is low */
			//Condicion de convergencia 2 :  Cambio pequeño
			if (!enough) {
				count = n;
				fmin = f;
			}
```
