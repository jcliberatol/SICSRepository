#Cotas y parametros usados en el algoritmo EM
Grupo de Investigacion SICS

##BFGS

###Numero de iteraciones del algoritmo 
El numero de iteraciones del BFGS.(Se espera implementacion adaptativa), actual = 15


```
//src/optimizer/Optimizer.cpp

	int r = bfgs(functionPtr,gradientPtr,args,pars,nargs,npars,15);

```
###Codigos de salida

Los codigos de salida son los siguientes :
*1 : Iteraciones Maximas Negativas
*2 : Valores Iniciales Malos, funcion no evaluable.
*3 : El algoritmo termino sin converger, Iteraciones Maximas
*0 : Convergencia

###Parametros del BFGS
```
//src/type/Constant.cpp

double Constant::stepredn =	0.020;
double Constant::acctol	=	0.0001;
double Constant::reltest	=	10.0;
double Constant::abstol  =    0.00001;
double Constant::reltol  =    1e-8;

```

##Estimacion de habilidades

###EAP y MAP

El optimizador brentfmin debe ser acotado asi :

```
//src/estimation/bayesian/LatentTraitEstimation.cpp

(*lt->traits)(counter, lt->dim - 1) = Brent_fmin(-5, 5, function,(void*) &temp, 0.0001220703);

```

###Modelos

Cota de success probability.

```
//src/model/parameter/_INSERTMODEL_.cpp

			tp = (ThreePLModel::successProbability ( theta[k], a,b,c));
			if (tp<1e-08)tp=1e-08;
			tq = 1-tp;

```
