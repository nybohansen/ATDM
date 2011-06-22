/*
 * vonmises2d.cpp --- Bivariate von Mises distribution, cosine variant
 *
 *  Copyright (C) 2008, Wouter Boomsma, The Bioinformatics Centre, University of Copenhagen.
 *
 *  This file is part of Mocapy++.
 *
 *  Mocapy++ is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Mocapy++ is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with Mocapy++.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "vonmises2d.h"
#include <boost/math/special_functions/detail/bessel_i0.hpp>
#include "../utils/netlib/integrate.h"
#include "../utils/netlib/functions.h"
#include "vonmises2dess.h"

namespace mocapy {


// For values above this limit the calculations are scaled to avoid overflow
#define LIMITSCALING 300
// For values above this limit an asymptotic expression for the bessel function is used
#define LIMITAPPROXBESSEL 300

// Exact integrand for normalization
double logNormConstExactIntegrand(double *phi, double *extraArguments) {
     double k1 = extraArguments[0];
     double k2 = extraArguments[1];
     double k3 = extraArguments[2];
     double k13 = sqrt((k1*k1) + (k3*k3) - 2*k1*k3*cos(*phi));
     double res = (2*M_PI*exp(log(i0(k13)) + k2*cos(*phi)));
     return res;
}

// Integrand using approximation for Bessel function (works for large values of k13)
double logNormConstApproximateIntegrand(double *phi, double *extraArguments) {
     double k1 = extraArguments[0];
     double k2 = extraArguments[1];
     double k3 = extraArguments[2];
     double scaleTerm = extraArguments[3];
     double k13 = sqrt((k1*k1) + (k3*k3) - 2*k1*k3*cos(*phi));
     double e = exp(k13 + k2*cos(*phi) - scaleTerm);
     return sqrt((2*M_PI)/k13)*e;
}


// Calculate logarithm of normalize constant c(k1, k2, k3)
double computeLogNormConst(Vector_nD &k, double lower, double upper) {
     double k1 = k[0];
     double k2 = k[1];
     double k3 = k[2];
     double exactLimit = 400;
     if ((sqrt(k1*k1 + k2*k2 - 2*k1*k3) + k2) < exactLimit) {
          double extraArguments[3] = {k1, k2, k3};
          double resExact = integrateQuad(&logNormConstExactIntegrand, extraArguments, lower, upper);
          return (log(resExact));
     } else {
          double scaleTerm = sqrt(k1*k1 + k3*k3 - 2*k1*k3) + k2;
          double extraArguments[4] = {k1, k2, k3, scaleTerm};
          double resApproximate = integrateQuad(&logNormConstApproximateIntegrand, extraArguments, lower, upper);
          return (log(resApproximate) + scaleTerm);
     }
     return 0;
}

// Density function for a univariate von Mises distribution
double VM1CosDensity(double phi, double k, double mu, double scaleTerm=0.0) {
     double res;
     // For k>limitApproxBessel an asymptotic expression for the bessel function is used
     if (k<LIMITAPPROXBESSEL) {
          res = exp(k*cos(phi-mu) - k - log(i0e(k)) - log(2*M_PI));
     } else {
          res = exp(k*cos(phi-mu) - k - 0.5*log(4*M_PI*M_PI*k) - scaleTerm);
     }
     return res;
}

// Mixture of univariate von Mises. Used as comparison function in rejection sampling.
VM1CosMixtureDensity::VM1CosMixtureDensity(Vector_nD k, Vector_nD mu, Vector_nD proportion) {
          this->k = k;
          this->mu = mu;
          this->proportion = proportion / sum(proportion);
          this->cumProportion = cumSum(proportion);
     }

     double VM1CosMixtureDensity::density(double phi, double scaleTerm) {
          double res = 0.0;
          for (int i=0; i<k.size; i++) {
               res += proportion[i] * VM1CosDensity(phi, k[i], mu[i], scaleTerm);
          }
          return res;
     }

     double VM1CosMixtureDensity::logDensity(double phi) {
       Vector_nD tmp = fmap(&cos, -(mu-phi))-1;
       double scaleTerm = -(-0.5*(max(tmp)+min(tmp)))*max(k);

          double res;
          // For k>limitScaling the calculations are scaled to avoid overflow
          if (max(k) < LIMITSCALING) {
               res = log(density(phi));
          } else {
               res = log(density(phi, scaleTerm)) + scaleTerm;
          }
          return res;
     }

     Vector_nD VM1CosMixtureDensity::logDensity(Vector_nD phi) {
          Vector_nD res(phi.size);
          for (int i=0; i<phi.size; i++) {
               res[i] = logDensity(phi[i]);
          }
          return res;
     }







/*****************/
/**** DENSITY ****/
/*****************/

// Constructor
VM2cosDensity::VM2cosDensity(Vector_nD &k, Vector_nD &mu) {
     init(k, mu);
     this->logNormConst = computeLogNormConst(k);
}

// Constructor - logarithm of normalization constant specified
VM2cosDensity::VM2cosDensity(Vector_nD &k, Vector_nD &mu, double logNormConst) {
     init(k, mu);
     this->logNormConst = logNormConst;
}

// Destructor
VM2cosDensity::~VM2cosDensity() {}


// Initializer
void VM2cosDensity::init(Vector_nD &k, Vector_nD &mu) {
     this->k = k;
     this->mu = mu;
}

// Update given new k and mu values
void VM2cosDensity::update(Vector_nD &k, Vector_nD &mu) {
     this->k = k;
     this->mu = mu;
     this->logNormConst = computeLogNormConst(k);
}

// Update given new k and mu values
void VM2cosDensity::update(Vector_nD &k, Vector_nD &mu, double logNormConst) {
     this->k = k;
     this->mu = mu;
     this->logNormConst = logNormConst;
}

// Calculate log likelihood
double VM2cosDensity::get_log_likelihood(double *anglePair) {
     double k1 = k[0];
     double k2 = k[1];
     double k3 = k[2];
     double mu1 = mu[0];
     double mu2 = mu[1];
     double phi = anglePair[0];
     double psi = anglePair[1];
     double dens;
     // if (logSpace) {
     dens = k1*cos(phi-mu1) + k2*cos(psi-mu2) - k3*cos(phi-mu1-psi+mu2) - logNormConst;
     // } else {
     //      dens = exp(k1*cos(phi-mu1) + k2*cos(psi-mu2) - k3*cos(phi-mu1-psi+mu2))/exp(logNormConst);
     // }
     assert(std::isfinite(dens));
     return dens;
}

// Calculate likelihood
double VM2cosDensity::get_likelihood(double *anglePair) {
     return exp(get_log_likelihood(anglePair));
}




// Return a von Mises distribution pseudo-random variate on [-pi, +pi].
// The implementation is similar to the algorithm by Best and Fisher,
// 1979; see N.I. Fisher, "Statistical Analysis of Circular Data",
// Cambridge University Press, 1993, p. 49.
  double vonMisesSampler(double k, double mean, RandomGen* rg) {
     if (mean < -M_PI || mean > M_PI) {
          fprintf(stderr, "vonMises Error: mean must be in the interval (-pi,pi). Mean=%f\n", mean);
     }
     double res;
     double a = 1.0 + sqrt(1.0 + 4.0*k*k);
     double b = (a - sqrt(2*a)) / (2*k);
     double r = (1.0 + b*b)/(2*b);
     double f;
     while(true) {
          // double U1 = RandomLib::Random::Global.FixedN();
          double U1 = rg->get_rand();
          double z = cos(M_PI * U1);
          f = (1.0 + r*z)/(r + z);
          double c = k * (r - f);
          // double U2 = RandomLib::Random::Global.FixedN();
          double U2 = rg->get_rand();
          if (((c*(2.0-c) - U2) > 0.0) || ((log(c/U2) + 1.0 - c >= 0.0))){
               break;           // accept
          }
     }
     // double U3 = RandomLib::Random::Global.FixedN();
     double U3 = rg->get_rand();
     if (U3 > 0.5) {
          res = fmod(acos(f)+mean, 2*M_PI);
     } else {
          res = fmod(-acos(f)+mean, 2*M_PI);
     }
     return res;
}




     // Constructor
  VM1CosMixtureSampler::VM1CosMixtureSampler(Vector_nD &k, Vector_nD &mu, Vector_nD &proportion, RandomGen* rg) {
          this->k = k;
          this->mu = mu;
          this->proportion = proportion / sum(proportion);
          this->cumProportion = cumSum(proportion);
	  this->rg = rg;
     }

     // Return single sample
     double VM1CosMixtureSampler::operator()() {
          double k=0.0;
          double mu=0.0;
          // double u = RandomLib::Random::Global.FixedN();
          double u = rg->get_rand();
          for (int i=0; i<cumProportion.size; i++) {
               if (cumProportion[i] >= u) {
                    k = this->k[i];
                    mu = this->mu[i];
                    break;
               }
          }
          double res = vonMisesSampler(k, mu, rg);
          return res;
     }

     // Return n samples
     Vector_nD VM1CosMixtureSampler::operator()(int n) {
          Vector_nD res(n);
          for (int i=0; i<n; i++) {
               res[i] = (*this)();
          }
          return res;
     }


// Bivariate von Mises - cosine model. Marginal of phi. Density
class VM2CosMarginalDensity {
private:
     Vector_nD k;
     double logNormConst;
public:
     // Constructor
     VM2CosMarginalDensity(Vector_nD &k) {
          init(k);
          this->logNormConst = computeLogNormConst(k);
     }

     // Constructor - logarithm of normalization contant specified
     VM2CosMarginalDensity(Vector_nD &k, double logNormConst) {
          init(k);
          this->logNormConst = logNormConst;
     }

     // Initializer
     void init(Vector_nD &k) {
          this->k = k;
     }

     // Return log density for single phi
     double logDensity(double phi) {
          double k1 = k[0];
          double k2 = k[1];
          double k3 = k[2];
          double k13 = sqrt(k1*k1 + k3*k3 - 2*k1*k3*cos(phi));

          double res;
          // For values above limitApproxBessel an asymptotic expression for the bessel function is used
          if (k1*k1 + k2*k2 < LIMITAPPROXBESSEL) {
               res = log(2*M_PI) - logNormConst + log(i0(k13)) + k2*cos(phi);
                    } else {
               res = -logNormConst + log(sqrt((2*M_PI)/k13)) + k13 + k2*cos(phi);
          }
          return res;
     }

     // Return log densities for vector of phis
     Vector_nD logDensity(Vector_nD &phi) {
          Vector_nD res(phi.size);
          for (int i=0; i<phi.size; i++) {
               res[i] = logDensity(phi[i]);
          }
          return res;
     }

     // Retun density for single phi
     double density(double phi) {
          double k1 = k[0];
          double k2 = k[1];
          double k3 = k[2];
          double k13 = sqrt(k1*k1 + k3*k3 - 2*k1*k3*cos(phi));
          double res =  2*M_PI*i0(k13)*exp(k2*cos(phi));
          return res;
     }

     // When class is called as function it returns the negative density function
     double operator()(double phi) {
          return -logDensity(phi); // Note that this is the NEGATIVE density function
     }

};


// Bivariate von Mises - cosine model. Marginal of phi. Sampler
// Initializer
void VM2CosMarginalSampler::init(Vector_nD &k) {
          this->k = k;
          double k1 = k[0];
          double k2 = k[1];
          double k3 = k[2];

          // Determine whether marginal is bimodal and set roots accordingly
          // We can use the exponentially scaled Bessel functions
          // since the same scaling factor appears in the numerator
          // and the denominator
          root = 0.0;
          if (((i1e(fabs(k1 - k3)) / i0e(fabs(k1-k3))) > ((fabs(k1-k3) * k2)/(k1*k3))) && // Multimodal criterium Theorem 4
              ((derivativeFactorB(-M_PI) * derivativeFactorB(0.0)) < 0)) { // In case k1>k3>0 and k2>k3>0 is not fulfilled

               // Find root of b (theorem 4)
               VM2CosMarginalDensity densityObject(k, this->logNormConst);
               brentOptimizer.optimize(&densityObject, Bracket<double>(&densityObject, 0, -M_PI/2, -M_PI));

               root = brentOptimizer.xmin;
          }

          // find parameters for mixture of von Mises distributions that fits best
          envelopeSampler = NULL;
          findOptimalComparisonFunction();
     }

     // Constructor
VM2CosMarginalSampler::VM2CosMarginalSampler(Vector_nD &k, RandomGen* rg) {
	  this->rg = rg;
          this->logNormConst = computeLogNormConst(k);
          init(k);
     }

     // Constructor
VM2CosMarginalSampler::VM2CosMarginalSampler(Vector_nD &k, double logNormConst, RandomGen* rg) {
          this->logNormConst = logNormConst;
	  this->rg = rg;
          init(k);
     }

     // Copy constructor
VM2CosMarginalSampler::VM2CosMarginalSampler(const VM2CosMarginalSampler &sampler, RandomGen* rg) {
	  this->rg = rg;
          k = sampler.k;
          logNormConst = sampler.logNormConst;
          root = sampler.root;

          envelopeSampler = new VM1CosMixtureSampler(*sampler.envelopeSampler);
          envelopeDensity = new VM1CosMixtureDensity(*sampler.envelopeDensity);
          K = sampler.K;
     }

     // destructor
VM2CosMarginalSampler::~VM2CosMarginalSampler() {
          delete envelopeSampler;
          delete envelopeDensity;
     }

     // Generate a random sample from the marginal distribution.
     // Rejection sampling: 1) Sample random value r from comparison distribution
     //                     2) Sample random value u from uniform distribution
     //                     3) Accept r if u < f_target(r)/f_comparison(r) """
  double VM2CosMarginalSampler::operator()(bool *status) {


         VM2CosMarginalDensity density(k, logNormConst);
         int count = 0;
         while(1) {
              double r = (*this->envelopeSampler)();
              // double u = RandomLib::Random::Global.FixedN();
	      double u = rg->get_rand();
              double x = exp(density.logDensity(r) - (K + envelopeDensity->logDensity(r)));

              if (u<x) {
                   *status=true;
                   return(r);
              }

              if ((count > 0) && count%(10)==0) {
                   fprintf(stderr, "Large rejection count: %d\n", count);
                   if (count > 50) {
                        *status=false;
                        return 0.0;
                   }
              }
              count++;
         }
    }

    // Generate n samples
    Vector_nD VM2CosMarginalSampler::operator()(int n) {

         Vector_nD phi(n);
         int i=0;
         while (i<n) {
              bool success;
              phi[i] = (*this)(&success);

              if (success) {
                   i++;
              } else {
                   return Vector_nD();
              }
         }
         return phi;
     }

     // Factor b from derivative of density. Defined in proof of theorem 4
     double VM2CosMarginalSampler::derivativeFactorB(double phi){
          double k1 = k[0];
          double k2 = k[1];
          double k3 = k[2];
          double k13 = sqrt(k1*k1 + k3*k3 - 2*k1*k3*cos(phi));
          return (-k2 + k1*k3*(i1e(k13)/i0e(k13))/(k13));
     }


     // Find parameter k and scale factor K so the the target marginal density is at
     // no point above the comparison function.
     // Comparison function: Mixture of two univariate von Mises distributions
     VM2CosMarginalSampler::MaxRatio::MaxRatio(Vector_nD &k, double root, double logNormConst) {
               this->k = k;
               this->root = root;
               this->logNormConst = logNormConst;
          }

	  // Reparameterize k parameters
          void VM2CosMarginalSampler::MaxRatio::reparameterisation(double *k) {
               *k = fabs(*k);
          }

	  // Calculate parameters
          double VM2CosMarginalSampler::MaxRatio::operator()(Vector_nD x) {
               return (*this)(x[0]);
          }

	  // Calculate parameters
          double VM2CosMarginalSampler::MaxRatio::operator()(double x) {
               if (std::isnan(x)) {
                    fprintf(stderr, "k=nan proposed. Returning inf.");
                    assert(false);
                    return INFINITY;
               }

               reparameterisation(&x);

               Vector_nD phi = range(-M_PI, M_PI, 0.01);  // Hack: Only check function at 0.01 intervals
                                                         // instead of doing full minimization.

               VM1CosMixtureDensity envelopeDensity(Vector_nD(2,x,x), Vector_nD(2, root, -root), Vector_nD(2, 0.5, 0.5));
               VM2CosMarginalDensity density(k, logNormConst);
               double res = max(density.logDensity(phi) - envelopeDensity.logDensity(phi));
               return(res);
          }



     // Find parameter k and scale factor K so the the target marginal density is at
     // no point above the comparison function.
     // Comparison function: Mixture of four univariate von Mises distributions
     // This method is used for large values of k."""

VM2CosMarginalSampler::MaxRatioLargerMixture::MaxRatioLargerMixture(Vector_nD &k, double root, double logNormConst) {
               this->k = k;
               this->root = root;
               this->logNormConst = logNormConst;
          }

          // Reparameterisation of parameters.
          // k values are forced to be positive
          // proportion values are forced to be between 0 and 1"""
          void VM2CosMarginalSampler::MaxRatioLargerMixture::reparameterisation(Vector_nD *x, int kSize, bool inverse) {
               // first kSize elements of x are k values
               assert (kSize <= x->size);
               for (int i=0; i<kSize; i++) {
                    (*x)[i] = fabs((*x)[i]);
               }

               for (int i=kSize; i<x->size; i++) {
                    (*x)[i] = fabs((*x)[i]);

                    if (!inverse) {
                         (*x)[i] = 1-(1/exp((*x)[i]));
                    } else {
                         (*x)[i] = -log(1-(*x)[i]);
                    }
               }
          }

	  // Calculate parameters
          double VM2CosMarginalSampler::MaxRatioLargerMixture::operator()(Vector_nD x) {
               reparameterisation(&x, 2);

               Vector_nD phi = range(-M_PI, M_PI, 0.01);  // Hack: Only check function at 0.01 intervals
                                                          // instead of doing full minimization.
               VM1CosMixtureDensity envelopeDensity(Vector_nD(4, x[0], x[0], x[1], x[1]),
                                                    Vector_nD(4, root, -root, root, -root),
                                                    Vector_nD(4, 0.5, 0.5, x[2], x[2]));
               VM2CosMarginalDensity density(k, logNormConst);
               double res = max(density.logDensity(phi) - envelopeDensity.logDensity(phi));
               return(res);
          }


     // Minimize the necessary scale factor K
     void VM2CosMarginalSampler::findOptimalComparisonFunction() {

          MaxRatio function(k,root,logNormConst);

          Powell powell;
          powell.optimize(&function, Vector_nD(1, min(k)));

          function.reparameterisation(&powell.xmin[0]);

          Vector_nD k(2, powell.xmin[0], powell.xmin[0]);
          Vector_nD mu(2, root, -root);
          Vector_nD proportion(2, 0.5, 0.5);

          // Scale factor
          this->K = powell.fmin;

          if (K > 5) {
               // Retry with mixture model with four components
               MaxRatioLargerMixture function2(k,root,logNormConst);
               Vector_nD startValues(3, min(k), 0.1, 0.2);
               bool inverse = true;
               function2.reparameterisation(&startValues, 2, inverse);

               Powell powell2;
               powell2.optimize(&function2, startValues);

               function2.reparameterisation(&powell2.xmin, 2);
               k = Vector_nD(4, powell2.xmin[0], powell2.xmin[0], powell2.xmin[1], powell2.xmin[1]);
               mu = Vector_nD(4, root, -root, root, -root);
               proportion = Vector_nD(4, 0.5, 0.5, powell2.xmin[2], powell2.xmin[2]);

               this->K = powell2.fmin;
          }

          this->envelopeSampler = new VM1CosMixtureSampler(k, mu, proportion, rg);
          this->envelopeDensity = new VM1CosMixtureDensity(k, mu, proportion);
     }





/*****************/
/**** SAMPLER ****/
/*****************/

// Sample from a bivariate von Mises distritbution cosine model
// The distribution and the techniques used in this code are described in:
// 'Bivariate Densities for Angular Data' (draft), Mardia, Subramaniam, Taylor(2005)
// Any references found in comments are to this paper.
VM2cosSampler::VM2cosSampler(Vector_nD &k, Vector_nD &mu, double logNormConst, RandomGen* rg): marginalPhiSampler(NULL) {
     this->k = k;
     this->mu = mu;
     this->logNormConst = logNormConst;
     this->rg = rg;
}

VM2cosSampler::VM2cosSampler(Vector_nD &k, Vector_nD &mu, RandomGen* rg): marginalPhiSampler(NULL) {
     this->k = k;
     this->mu = mu;
     this->logNormConst = computeLogNormConst(k);
     this->rg = rg;
}

// Copy constructor
VM2cosSampler::VM2cosSampler(const VM2cosSampler &sampler, RandomGen* rg) {
     this->k = sampler.k;
     this->mu = sampler.mu;
     this->logNormConst = sampler.logNormConst;
     this->rg = rg;

     if (sampler.marginalPhiSampler) {
          marginalPhiSampler = new VM2CosMarginalSampler(*sampler.marginalPhiSampler);
     } else {
          marginalPhiSampler = NULL;
     }
}

// Destructor
VM2cosSampler::~VM2cosSampler() {
     delete marginalPhiSampler;
}

// Set random number generator
void VM2cosSampler::setRandomGen(RandomGen *rg) {
     this->rg = rg;
     if (marginalPhiSampler)
          marginalPhiSampler->setRandomGen(rg);
}

// Update given new k and mu values
void VM2cosSampler::update(Vector_nD &k, Vector_nD &mu) {
     this->k = k;
     this->mu = mu;
     this->logNormConst = computeLogNormConst(k);
     delete marginalPhiSampler;
     marginalPhiSampler = NULL;
}

// Update given new k and mu values
void VM2cosSampler::update(Vector_nD &k, Vector_nD &mu, double logNormConst) {
     this->k = k;
     this->mu = mu;
     this->logNormConst = logNormConst;
     delete marginalPhiSampler;
     marginalPhiSampler = NULL;
}

// Generate a single sample
void VM2cosSampler::sample(std::vector<double> &rVal) {
     // First sample the marginal of psi f(psi). Then use these
     // values to sample from the conditional f(phi|psi). See Section 2 in paper

     // Assign marginal sampler if not already present

  assert(rg!=NULL);

     if (!marginalPhiSampler) {
       assert(rg != NULL);
       marginalPhiSampler = new VM2CosMarginalSampler(k, this->logNormConst, rg);
     }


     bool success;

     double psi = (*this->marginalPhiSampler)(&success);

     if (!success)
          return;

     double k1 = k[0];
     double k3 = k[2];
     double mu1 = mu[0];
     double mu2 = mu[1];
     double k13 = sqrt(k1*k1 + k3*k3 - 2*k1*k3*cos(psi));
     double psi_mu = atan((-k3*sin(psi))/(k1 - k3*cos(psi)));
     double phi = vonMisesSampler(k13, psi_mu, rg);

     // Angles are in the interval [-pi, pi]. Add 3*pi to bring them
     // to: [2*pi, 4*pi] which with mu values in [-pi, pi] brings it
     // to [pi, 3pi], which can then readily be transformed back to [-pi, pi]
     psi = fmod(psi + (3*M_PI) + mu2, 2*M_PI) - M_PI;
     phi = fmod(phi + (3*M_PI) + mu1, 2*M_PI) - M_PI;

     rVal.resize(2);
     rVal[0] = phi;
     rVal[1] = psi;
}




/*******************/
/**** ESTIMATOR ****/
/*******************/

class logLikelihood_neg_opt {
private:
     Vector_nD mu;
     MDArray<double> &ess;
public:
     // Constructor
     logLikelihood_neg_opt(Vector_nD mu, MDArray<double> &ess): mu(mu), ess(ess) {}

     // Reparameterisation of parameters.
     // k1 and k2 values are forced to be positive
     // k3 is forced to be smaller than k1 and k2 and fabs(k3) < k1k2/(k1+k2)
     void reparameterization(Vector_nD &k, bool inverse=false) {

          // k1 and k2 must be positive (k3 can be negative)
          k[0] = fabs(k[0]);
          k[1] = fabs(k[1]);
               
	  if (!inverse) {
	       k[0] = k[0];
	       k[1] = k[1];
               k[2] = ((2/(1+exp(-k[2]/100)))-1) * min(Vector_nD(3, k[0], k[1], ((k[0]*k[1])/(k[0]+k[1]))));
          } else {
	       k[0] = k[0];
	       k[1] = k[1];
               double k12 = min(Vector_nD(3, k[0], k[1], ((k[0]*k[1])/(k[0]+k[1]))));
               k[2] = 100*(log(k12+k[2]) - log(k12-k[2]));
          }
     }


     // Return reasonable start values for optimization
     Vector_nD parameterStartValues_old() {

	  // Start values - assuming marginal distributions are von Mises
	  // see section 6

	  // Normalize
	  double C1 = ess[VonMises2dESS::C1]/ess[VonMises2dESS::COUNT];
	  double C2 = ess[VonMises2dESS::C2]/ess[VonMises2dESS::COUNT];;
	  double S1 = ess[VonMises2dESS::S1]/ess[VonMises2dESS::COUNT];;
	  double S2 = ess[VonMises2dESS::S2]/ess[VonMises2dESS::COUNT];;

	  // Approximations of start values from Dobson(1978)
	  double R1 = sqrt(C1*C1 + S1*S1);
	  double R2 = sqrt(C2*C2 + S2*S2);
	  double k1 = (1.28-0.53*R1)*tan(0.5*M_PI*R1);
	  double k2 = (1.28-0.53*R2)*tan(0.5*M_PI*R2);
	  double k3 = fmin(k1,k2)/10.0;

	  Vector_nD res(3, k1,k2,k3);

	  // Inverse reparameterisation
	  bool inverse = true;
	  reparameterization(res, inverse);

	  return res;
     }

     // Return reasonable start values for optimization
     // Use ME estimation
     Vector_nD parameterStartValues() {

	  Vector_nD estimate_ME = VonMises2dEstimate_ME(ess);

	  Vector_nD res(3, estimate_ME[0], estimate_ME[1], estimate_ME[2]);

          // Ensure that k3>k1, k3>k2 and k3<((k1*k2)/(k1+k2))
          double res12 = ((res[0]*res[1])/(res[0]+res[1]));
          if (res[2] > res[0] || res[2] > res[1]) {
               res[2] = std::min(res[0], res[1])-1e-9;
          } 

          if (res[2]>0 && res[2] > res12) {
               res[2] = res12-1e-9;
          } else if (res[2]<0 && res[2] < -res12) {
               res[2] = -res12+1e-9;
          }

          assert(std::isfinite(res[0]));
          assert(std::isfinite(res[1]));
          assert(std::isfinite(res[2]));

	  // Inverse reparameterisation
	  bool inverse = true;
	  reparameterization(res, inverse);

          assert(std::isfinite(res[0]));
          assert(std::isfinite(res[1]));
          assert(std::isfinite(res[2]));

	  return res;
     }

     // Evaluation - returns the negative log likelihood
     double operator()(Vector_nD k) {
	  reparameterization(k);

	  double logNormalizeContant = computeLogNormConst(k);
	  if (std::isinf(logNormalizeContant)) {
	       return std::numeric_limits<double>::infinity();
	  }

	  double logLikelihood = (k[0]*(cos(mu[0])*ess[VonMises2dESS::C1] + sin(mu[0])*ess[VonMises2dESS::S1]) +
				  k[1]*(cos(mu[1])*ess[VonMises2dESS::C2] + sin(mu[1])*ess[VonMises2dESS::S2]) -
				  k[2]*((cos(mu[0]-mu[1])*ess[VonMises2dESS::COS1MIN2]) +
					(sin(mu[0]-mu[1])*ess[VonMises2dESS::SIN1MIN2]))) - ess[VonMises2dESS::COUNT]*(logNormalizeContant);

	  // Return negative (for optimization purposes)
	  return -logLikelihood;
     }
};


std::vector<double> VonMises2dEstimate_ML(MDArray<double> &ess) {

     // No estimation if too few data points
     if (ess[VonMises2dESS::COUNT] < 10) {
	  return std::vector<double>();
     }

     // Estimation of mu1 and mu2 by the directional means
     double mu1 = atan2(ess[VonMises2dESS::S1], ess[VonMises2dESS::C1]);
     double mu2 = atan2(ess[VonMises2dESS::S2], ess[VonMises2dESS::C2]);

     logLikelihood_neg_opt func(Vector_nD(2, mu1, mu2), ess);

     // Vector_nD startValues = parameterStartValues_old(ess);
     Vector_nD startValues = func.parameterStartValues();

     // std::cout << "StartValues: " << startValues << " LL=" << func(startValues) << "\n";

     Powell powell;

     powell.optimize(&func, startValues);

     // Vector_nD optimizedValues(3, powell.xmin[0], powell.xmin[1], powell.xmin[2]);
     // std::cout << "optimizedValues: " << optimizedValues << " LL=" << func(optimizedValues) << "\n";


     func.reparameterization(powell.xmin);

     assert(std::isfinite(powell.xmin[0]));
     assert(std::isfinite(powell.xmin[1]));
     assert(std::isfinite(powell.xmin[2]));

     // return Vector_nD(5, powell.xmin[0], powell.xmin[1], powell.xmin[2], mu1, mu2);
     return vec(powell.xmin[0], powell.xmin[1], powell.xmin[2], mu1, mu2);
}


std::vector<double> VonMises2dEstimate_ME(MDArray<double> &ess) {

     // No estimation if too few data points
     if (ess[VonMises2dESS::COUNT] < 10) {
	  return std::vector<double>();
     }

     // Estimation of mu1 and mu2 by the directional means
     double mu1 = atan2(ess[VonMises2dESS::S1], ess[VonMises2dESS::C1]);
     double mu2 = atan2(ess[VonMises2dESS::S2], ess[VonMises2dESS::C2]);


     double S1 = 0.5 - ((1.0/(2*ess[VonMises2dESS::COUNT]))*cos(2*mu1)*ess[VonMises2dESS::_2COS1] +
			(1.0/(2*ess[VonMises2dESS::COUNT]))*sin(2*mu1)*ess[VonMises2dESS::_2SIN1]);
     double S2 = 0.5 - ((1.0/(2*ess[VonMises2dESS::COUNT]))*cos(2*mu2)*ess[VonMises2dESS::_2COS2] +
			(1.0/(2*ess[VonMises2dESS::COUNT]))*sin(2*mu2)*ess[VonMises2dESS::_2SIN2]);
     double S12 = ((1.0/(2*ess[VonMises2dESS::COUNT])) * cos(mu1-mu2) * ess[VonMises2dESS::COS1MIN2] +
		   (1.0/(2*ess[VonMises2dESS::COUNT])) * sin(mu1-mu2) * ess[VonMises2dESS::SIN1MIN2] -
		   (1.0/(2*ess[VonMises2dESS::COUNT])) * cos(mu1+mu2) * ess[VonMises2dESS::COS1PLUS2] -
		   (1.0/(2*ess[VonMises2dESS::COUNT])) * sin(mu1+mu2) * ess[VonMises2dESS::SIN1PLUS2]);

     double k1 = (S2-S12)/(S1*S2-(S12*S12));
     double k2 = (S1-S12)/(S1*S2-(S12*S12));
     double k3 = (-S12)/(S1*S2-(S12*S12));

     // Translation in case k3 is negative (this no longer occurs?)
     if (k1 < 0) {
	  mu1 =- M_PI;
	  k1 *= -1;
     }
     if (k2 < 0) {
	  mu2 =- M_PI;
	  k2 *= -1;
     }

     // return (Vector_nD(5, k1, k2, k3, mu1, mu2));
     return vec(k1, k2, k3, mu1, mu2);
}



}
