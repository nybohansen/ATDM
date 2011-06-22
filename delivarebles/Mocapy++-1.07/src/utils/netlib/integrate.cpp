// integrate.cpp --- Wrapper for fortran integration code
// Copyright (C) 2006-2008 Wouter Boomsma
//
// This file is part of BackboneDBN
//
// BackboneDBN is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// BackboneDBN is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.

// You should have received a copy of the GNU Lesser General Public License
// along with BackboneDBN.  If not, see <http://www.gnu.org/licenses/>.


#include <cassert>
#include "integrate.h"

// Computes a definite integral
// Integrate func from a to b (possibly infinite interval) using a technique
// from the Fortran library QUADPACK
double integrateQuad(double(*func)(double *), double a, double b, int *evaluations,
            double epsabs, double epsrel, int limit) {

     double res;
     double abserr;
     int ier;
  
     int *iwork = new int[limit];
     int lenw = limit*4;
     int last;
     double *work = new double[lenw];
     int neval;
  
     dqags_(func, &a, &b, &epsabs, &epsrel, &res, &abserr, &neval, &ier, &limit, &lenw, &last, iwork, work);

     // Check for error 
     assert(ier==0);
  
     delete[] work;
     delete[] iwork;

     if (evaluations) {
          *evaluations = neval;
     }
     return res;
}

// Global variables necessary to allow use to pass additional arguments with function
double (*quadFunc)(double *, double *);
void *quadExtraArgs = NULL;

// Single-argument funciton wrapper
// Calls the original function including additional arguments
double _quadSingleArgumentFunction(double *x) {
     return quadFunc(x, (double *)quadExtraArgs);
}

// Version of quad taking additional arguments of type double
double integrateQuad(double(*func)(double *, double *), double extraArguments[], double a, double b, int *evaluations,
            double epsabs, double epsrel, int limit) {
     quadFunc = func;                           // Original function pointer is saved
     quadExtraArgs = (void *)extraArguments;            // Additional argument values are saved

     // quad is called on single-argument-wrapper
     return integrateQuad(&_quadSingleArgumentFunction, a, b, evaluations, epsabs, epsrel, limit);
}



extern "C" int MAIN__() { return 1; }
