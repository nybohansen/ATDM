/*
 *  GDMdensity.cpp
 *
 *  Copyright (C) 2010, Simon H. Rasmussen, The Bioinformatics Centre,
 *  University of Copenhagen.
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

#include <boost/math/special_functions/gamma.hpp>
#include <vector>
#include <iostream>
#include <sstream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "GDMdensity.h"
#include "../utils/utils.h"
#include <boost/random/gamma_distribution.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <algorithm>
#include <utility>

namespace mocapy
{

  using namespace std;
  
  GDMdensity::GDMdensity(){
    dim=2;
    init_pars();
  }
  
  GDMdensity::GDMdensity(uint dim, ulong s,int nn){
    this->dim = dim;
    n = nn;
    ulong seed = s;
    cache_initialized = false;
    init_pars();
    rng1.seed(356);
    rng2.seed(2387);
  }
  
  void printv(vector<double> & v1, vector<double> & v2){
    int d = v1.size();
    cout << "[";
    for(int i = 0;i < d;i++){
      printf("%.2f ",v1[i]);
    }
    cout << ",\n ";
    for(int i = 0;i < d;i++){
      printf("%.2f ",v2[i]);
    }
    cout << "]" << endl;
  }
  double GDMdensity::sample_beta(double a,double b){
    // Method for sampling of the beta distribution.
    boost::gamma_distribution<double> gam1(a);     
    boost::variate_generator<boost::mt19937&, boost::gamma_distribution<double> > gamma1(rng1, gam1);
    boost::gamma_distribution<double> gam2(b);     
    boost::variate_generator<boost::mt19937&, boost::gamma_distribution<double> > gamma2(rng2, gam2);
    double X = gamma1();
    double Y = gamma2();
    double rn = X/(X+Y);
    if (a<0 || b<0){
      cout << "in sample_beta: rn alphas betas are " << rn << " " << a << " " << b << endl;
      return -1; 
    }
    return rn;
  }

  vector<double> GDMdensity::Sample_Gendir(vector<double> & as,vector<double> & bs){
    // Calculates a Generalized Dirichlet sample sampling from the beta distribution.
    double sum;
    vector<double> theta;
    for (int i = 0;i<dim;i++){
      theta.push_back(0);
    }
    theta[0] = sample_beta(as[0], bs[0]);
    sum = theta[0];
    
    for(int i = 1;i < dim-1;i++){
      theta[i] = sample_beta(as[i], bs[i]) * (1-sum);
      sum = sum + theta[i];
    }
    
    sum = 0;
    for (int i = 0;i<dim-1;i++){
      sum = sum + theta[i];
    }
    theta[dim-1] = (1-sum);
    return theta;
  }
  
  void GDMdensity::scanData(){
    // sets up cache for lnGammai()
    if (!cache_initialized){
      cache_initialized = true;
      int  max = 0;
      int sum;
      for(int i = 0;i<n;i++){
	sum = 0;
	for(int j = 0;j<dim;j++){
	  sum = sum + (int) (cvector[i][j] + .1);
	}
	if (sum > max){
	  max = sum;
	}
      }
      Nmax = max;
      for(int k = 0;k < Nmax-1;k++){
	lnGammaCache.push_back(-1);
      }
    }
  }

  vector<double> GDMdensity::Sample_mult(vector<double> & as, vector<double> & bs,int N){
    // A method for sampling multinomial vectors.
    CPD cumsum;
    double rnum;
    vector<uint> shape;
    shape.push_back(dim);
    cumsum.set_shape(shape);
    vector<double> theta = Sample_Gendir(as, bs);
    for(int i = 0;i<dim;i++){
      cumsum[i] =  theta[i];
    }
    cumsum.cumsum();
    vector<double> multsample;
    for(int i = 0;i < dim;i++){ 
      multsample.push_back(0);
    }
    int choice;
    for(int i = 0;i<N;i++){
      rnum = get_rand();
      choice = cumsum.bisect(rnum);
      multsample[choice] = multsample[choice] + 1; 
    }
    return multsample;
  }

  void GDMdensity::add_countv(vector<double> & c){
    // Add vector of counts.
    cvector.push_back(c);
    n = cvector.size();    
  }

vector< vector<double> > GDMdensity::reject( vector< vector<double> > & p){
  // In case estimates are nan of <= 0 we reject them in the current iteration
  // and use the ones from previous iteration.
  double init = 1.02;
  double lim = 0.1;
  for(int d = 0;d < dim-1;d++){
    if (p[0][d] < lim || isnan(p[0][d])){
      p[0][d] = init;
    }
    if (p[1][d] < lim || isnan(p[1][d])){
      p[1][d] = init;
    }
  }
  return p;
}

double GDMdensity::lngammai(int c){
  // Accessing the log(gamma(x))-function in gsl for integers, uses caching.
  if (c <= 2){
    return 0;
  }
  if (lnGammaCache[c-3] == -1){
    lnGammaCache[c-3] = boost::math::lgamma(c);
  }
  return lnGammaCache[c-3];
}

double GDMdensity::getMult_loglik(vector<double> & countv ,vector<double> & f){
  // Calculates the multinomial log-likelihood.
  double ll;
  int c;
  int s = 0;

  for(int i = 0;i < dim;i++){
    s = s + (int) (countv[i]+.1);
  }
  ll = lngammai(s+1);

  for(int i = 0;i < dim;i++){
    c = (int) (countv[i]+.1);
    ll = ll - lngammai(c+1);;
    ll = ll + c *  log(f[i]);
  }
  return ll;
}

double GDMdensity::loglik(vector<double> & countv){
  // Calculates the GDM log-likelihood
  int s = 0;
  int c;
  double ll = 0;
  double t2 = 0;
  double t1 = 0;
  double ap = 0;
  double bp = 0;
  for(int i = 0;i < dim;i++){
    s = s + (int) countv[i]+.1;
  }
  if(cache_initialized){
    ll = lngammai(s+1);
  }
  else{
    ll = lgamma(s+1);
  }
  for(int i = 0;i < dim;i++){
    c = (int) countv[i]+.1;
    if(cache_initialized){
      ll = ll - lngammai(c+1);
    }
    else{
      ll = ll - lgamma(c+1);
    } 
    
    if (i < dim - 1){
      t1 = t1 + boost::math::lgamma(alphas[i] + betas[i]) - boost::math::lgamma(alphas[i]) - boost::math::lgamma(betas[i]);
      ap = alphas[i] + c;
      s = s - c;
      bp = betas[i] + s;
      t2 = t2 - boost::math::lgamma(ap+bp) + boost::math::lgamma(ap) + boost::math::lgamma(bp);
    }
  } 
  ll = ll + t1 + t2;
  return ll;
}

double GDMdensity::sum_loglik(){
  // Calculates sum of loglikelihoods
  double sum=0;
  vector<double> vec;
  for(int i = 0;i < n;i++){
    sum = sum + loglik(cvector[i]);
  }
  return sum;
}

double GDMdensity::covsum(int d1,double *means,vector< vector<double> > & vec){
  // Calculates number of dimensions that each have negative covariance with.
  double cov;
  double dl1;
  double dl2;
  double cnum;

  cnum = 0;
  for (int d2 = 0;d2 < dim; d2++){
    cov = 0;
    for (int i = 0;i < n;i++){
      dl1 = (vec[i][d1] - means[d1]);
      dl2 = (vec[i][d2] - means[d2]);
      cov = cov + dl1*dl2;
      cov=cov/(double) (n - 1);
    }
    if (d2 != d1){
      if (cov < 0){
	cnum = cnum - 1;
      }
    }
  }
  return cnum; 
}

bool GDMdensity::sortdim(vector< vector<double> > & vec){
  // Sorts dimensions with regard to count of negative covariance with other dimensions.
  vector< pair<double, int> > covs; 
  double *means = new double[dim]; 
  for(int i = 0;i < dim;i++){
    means[i] = mean(i,vec);
  }
  for (int d1 = 0;d1 < dim; d1++){
    covs.push_back(pair<double, int>(covsum(d1,means,vec), d1));
  }
  stable_sort(covs.begin(),covs.end());

  vector<double> covvec;
  vector<pair<double, int> >::const_iterator itr;

  for(itr = covs.begin(); itr != covs.end(); ++itr){
    covvec.push_back((*itr).first);
    idxvec.push_back((*itr).second);
  }
  delete[] means;

  for(int i = 0;i < dim ;i++){
    if (idxvec[i] != i){
      return true;
    }
  }
  return false;
}
		
  void GDMdensity::swap_dims(vector< int > & idx,vector< vector<double> > & vec){
  // Swapps the dimensions of the dynamically allocated array given as argument
  // according to the order defined by the idx vector.
    vector<double> tmp;
    vector<double> tmp1;
  for (int i = 0;i < n; i++){
    tmp = vec[i];
    for (int j = 0;j < dim; j++){
      tmp1.push_back(tmp[idx[j]]);
    }
    vec[i] = tmp1;
    tmp1.clear();
  }
}

  void GDMdensity::set_parameters(vector<double> & as,vector<double> & bs){
    // Makes it posible to directly set parameters in a density.
    // Check cpd is right shape and not yet set.
    if ((int) as.size() == dim - 1 && (int) bs.size() == dim - 1 && alphas.empty() && betas.empty()){
      alphas=as;
      betas=bs;
    }
  }

  void GDMdensity::init_pars(){
    // Initializes the parameters to ones.
    int a = 6;
    if (alphas.size() < dim - 1 && betas.size() < dim - 1){
      alphas.clear();
      betas.clear();
      for(int i = 0;i<dim-1;i++){
	alphas.push_back(1);
	betas.push_back(1);
	//      alphas.push_back(a);
	//      betas.push_back(a*dim - a*i);
      }
    }
  }


  vector<int> GDMdensity::get_idxVec(){
    // get correct order of dimension in case they have been rearranged
    return idxvec;
  }

  vector< vector<double> > GDMdensity::metro_hast(int samples, bool init, bool sort){
  // Estimation of multinomial parameter set used for estimating GD parameters.
  // This is implemnted using Metropolis-Hastings sampling.
  double p1;
  double p2;
  double ran;
  int num = 0;
  pvector.clear();
  vector< vector<double> > p;
  vector< vector<double> > pars;
  vector< vector<double> > tmppars;
  if (init){
    scanData();
    init_pars();
    pars.push_back(alphas);
    pars.push_back(betas);

    if (sort && sortdim(cvector)){
      cout << "Dimension has been sorted, dimensions will be swapped" << endl;
      swap_dims(idxvec,cvector);
    }
  }

  for(int i = 0;i < n;i++){
    vector <double> f1 = Sample_Gendir(alphas, betas);
    p1 = getMult_loglik(cvector[i], f1);
    for(int j = 0;j < samples;j++){
      for(int j = 0;j<dim-1;j++){
	// Report if nan's appear in parameters. 
	//If this happens the user needs to rethink the model.
	if (isnan(alphas[i]) || isnan(betas[i])){
	  cout << "in Metro-Hastings: alphas betas i>"<< alphas[i] << " " << betas[i] << " " << i << endl;
	}
      }
      // Sample multinomial parameter vector.
      vector<double> f2 = Sample_Gendir(alphas, betas);
      p2 = getMult_loglik(cvector[i], f2);

      // if log-likelihood is better for f2, it is chosen.
      if (p2 > p1){
	f1 = f2;
	p1 = p2;
      }
      else{
	// if not we do a random choice.
	ran=get_rand();

	double rr = log((long double) ran);
	num++;
	if (rr < p2 - p1){
	  f1 = f2;
	  p1 = p2;
	}
      }
    }
    pvector.push_back(f1);
  }
  pars = estimate_pars(init);
  pars = reject(pars);
  alphas = pars[0];
  betas = pars[1];
  return pars;
}

  void GDMdensity::reset_cvec(){
    cvector.clear();
  }

double GDMdensity::mean(int d,vector< vector<double> > & vec){
  // Calculates mean of a vector.
  double mean = 0;
  int n = vec.size();
  for(int i = 0;i < n;i++){
    mean = mean + vec[i][d];
  }
  mean = mean/n;
  return mean;
}

double GDMdensity::var(double mean,int d,vector< vector<double> > & pvec){
  // Calculates variance of a vector, given the mean.
  double var=0;
  double delta;
  int n = pvec.size();
  for (int i = 0;i<n;i++){
    delta = (pvec[i][d] - mean);
    var = var + delta*delta;
  }
  var=var/(n-1);

  return var; 
}

vector< vector<double> > GDMdensity::estimate_pars(bool init){
  // Estimate GD parameters 
  vector< vector<double> > dpars;
  dpars.push_back(alphas);
  dpars.push_back(betas);
  double sumvar = 0;
  double t = 1;
  double u = 1;
  double e;
  double e2;
  double v;
  double den;
  vector< vector<double> > pars;

  for (int i = 0;i<dim-1;i++){
    e = mean(i,pvector);
    v = var(e,i,pvector);
    sumvar=sumvar + v;
    e2 = e*e;
    den = v*t+e2*t-e2*u;
    // If nan's appear here estimation has crashed.
    if (isnan(e) || isnan(v) || isnan(e2)){
      cout << "NAN e var e2: " << e << " " << v << " " << e2 << endl; 
    }   
    alphas[i] = (e*(-v+e*u-e2)/den);
    betas[i] = ((t-e)*(-v+e*u-e2)/den);
    // If nan's appear here estimation has crashed.
    if (isnan(alphas[i]) || isnan(betas[i])){
      cout << "in Estimate_pars: alphas betas i "<< alphas[i] << " " << betas[i] << " " << i << endl;
    }
    t = t * betas[i]/(alphas[i]+betas[i]);
    u = u * (betas[i]+1)/(alphas[i]+betas[i]+1);
  }
  
  pars.push_back(alphas);
  pars.push_back(betas);
  
  return pars;
}
}
