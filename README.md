SICSRepository
==============

Sics Secret Repository


#### Using in C++.

##### Multidimensional Case : 

```
double fulloglik = em.getLoglik();
    std::cout<<"Ll : "<<fulloglik<<std::endl;
    double* returnpars;

    returnpars = new double[3*datSet->size];
    model->parameterModel->getParameters(returnpars);
    // For 2pl & 1pl
    if(e_model < 3)
        for (int i = 2*datSet->size;i < 3*datSet->size;i++)
            returnpars[i]=0;
            
        // Return in list
        for (int i = 0; i < 3*datSet->size; i++)
            pars[i] = returnpars[i];
            std::cout<<pars[i]<<" "<<std::endl;
```
