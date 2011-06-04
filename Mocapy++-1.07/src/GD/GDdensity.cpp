/*
 *  GDdensity.cpp
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
#include "GDdensity.h"
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
  
  GDdensity::GDdensity(){
    dim=2;
    GDinit_pars();
  }
  
  GDdensity::GDdensity(uint dim, ulong s,int nn){
    this->dim = dim;
    n = nn;
    cache_initialized = false;
    GDinit_pars();
    rng1.seed(356);
    rng2.seed(2387);
  }
  
  void GDprintv(vector<double> & v1, vector<double> & v2){
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

  double GDdensity::GDsample_beta(double a,double b){
    //draw a sample from beta distribution
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

  vector<double> GDdensity::GDSample_Gendir(vector<double> & as,vector<double> & bs){
    // Calculates a Generalized Dirichlet sample sampling from the beta distribution.
    double sum;
    vector<double> theta;
    for (int i = 0;i<dim;i++){
      theta.push_back(0);
    }
    theta[0] = GDsample_beta(as[0], bs[0]);
    sum = theta[0];
    
    for(int i = 1;i < dim-1;i++){
      theta[i] = GDsample_beta(as[i], bs[i]) * (1-sum);
      sum = sum + theta[i];
    }
    
    sum = 0;
    for (int i = 0;i<dim-1;i++){
      sum = sum + theta[i];
    }
    theta[dim-1] = (1-sum);
    return theta;
  }
  void GDdensity::GDadd_pvec(vector<double> & p){
    // Add vector of counts.
    pvector.push_back(p);
  }

vector< vector<double> > GDdensity::GDreject( vector< vector<double> > & p){
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

  double GDdensity::GDloglik(vector<double> & P,bool init){
    double p = 0;
    double pj = 0;
    double gami;
    
    if (init){
      GDinit_pars();
    }
    for(int i = 0;i < dim - 1;i++){
      if (i == dim - 2){
	gami = (betas[i] - 1.0);
      }
      else{
	gami = (betas[i] - alphas[i + 1] - betas[i + 1]);
      }
      
      pj = pj + P[i];
      p = p + boost::math::lgamma(alphas[i]+betas[i]) - boost::math::lgamma(alphas[i]) - boost::math::lgamma(betas[i]) + (alphas[i] - 1) * log(P[i]) + (gami * log(1.0 - pj));
    }
    return p;
  }

  double GDdensity::GDsum_loglik(){
    double sum=0;
    vector<double> vec;
    for(int i = 0;i < n;i++){
      sum = sum + GDloglik(cvector[i],false);
    }
    return sum;
  }
  
  double GDdensity::GDcovsum(int d1,double *means,vector< vector<double> > & vec){
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
  
  bool GDdensity::GDsortdim(vector< vector<double> > & vec){
    // Sorts dimensions with regard to count of negative covariance with other dimensions.
    vector< pair<double, int> > covs; 
    double *means = new double[dim]; 
    for(int i = 0;i < dim;i++){
      means[i] = mean(i,vec);
    }
    for (int d1 = 0;d1 < dim; d1++){
      covs.push_back(pair<double, int>(GDcovsum(d1,means,vec), d1));
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
    
  void GDdensity::GDswap_dims(vector< int > & idx,vector< vector<double> > & vec){
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
  void GDdensity::GDset_parameters(vector<double> & as,vector<double> & bs){
    // Makes it posible to directly set parameters in a density.
    // Check cpd is right shape and not yet set.
    if ((int) as.size() == dim - 1 && (int) bs.size() == dim - 1 && alphas.empty() && betas.empty()){
      alphas=as;
      betas=bs;
    }
  }

  void GDdensity::GDinit_pars(){
    // Initializes the parameters to Dirichlet parameters (not Generalized Dirichlet).
    //  cout << "init_pars" << endl;
    int a = 6;
    if (alphas.size() < dim - 1 && betas.size() < dim - 1){
      alphas.clear();
      betas.clear();    for(int i = 0;i<dim-1;i++){
	alphas.push_back(1);
	betas.push_back(1);
      }
    }
  }

  void GDdensity::reset_pvec(){
    // Calculates mean of a vector.
    pvector.clear();
  }
  
  double GDdensity::mean(int d,vector< vector<double> > & vec){
    // Calculates mean of a vector.
    double mean = 0;
    int n = vec.size();
    for(int i = 0;i < n;i++){
      mean = mean + vec[i][d];
    }
    mean = mean/n;
    return mean;
  }

  double GDdensity::var(double mean,int d,vector< vector<double> > & pvec){
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

  vector< vector<double> > GDdensity::GDestimate_pars(bool init){
    if (init){
      GDinit_pars();
    }
    
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
    
    return GDreject(pars);
  }

}
