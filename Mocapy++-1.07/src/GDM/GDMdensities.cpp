/*
 *  GDMdensities.cpp
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

#include <sstream>
#include <iostream>
#include "GDMdensities.h"

namespace mocapy
{

GDMdensities::GDMdensities()
{
    initialize();
}

  void GDMdensities::initialize()
  {
    type = GDM;
    output_size = dim;
    // setting default values which  can be changed withthe mutators
    metro_samples = 300;
    sumOfCounts=1000;
  }

  
  GDMdensities::GDMdensities(uint dim,uint sample_size)
  {
    this->dim = dim;
    this->sample_size = sample_size;
    initialize();
}

void GDMdensities::set_parameters(MDArray<double> & pars)
{
    // Check cpd is right shape
    vector<uint> shape = pars.get_shape();
    assert(shape[1] == 2 * (dim - 1));

    for(uint i = 0;i < node_size;i++){
      cpd.push_back(pars);
    }
    for(uint i = 0;i < node_size;i++){
      cpd[node_size - 1 - i] = pars.get_view(i);
    }
      assert(!cpd.empty());
      update_pars();
}

void GDMdensities::construct(vector<uint> & parent_sizes)
{
    // GDMNode always has a single parent
    assert(parent_sizes.size() == 1);
    node_size = parent_sizes[0];
    // vector of GDMDensity objects
    densities.resize(node_size);
    for(uint i = 0;i < parent_sizes.size();i++){
    }
    _create_densities();
}

void GDMdensities::_create_densities()
{
  first = true;
  ulong seed = 4450;
  MDArray<double> pars;
  vector<double> alphas;
  vector<double> betas;
  // Initialize the density objects and fill the density vector
  for(uint i=0; i<node_size; i++)
    {
      densities[i] = GDMdensity(dim, seed, sample_size);
      init.push_back(true);
    }
}

void GDMdensities::update_pars(){
  assert(!cpd.empty());
  MDArray<double> pars;
  vector<double> newparsalpha;
  vector<double> newparsbeta;
  for(uint j = 0; j < node_size; j++){
    pars = cpd[j];
    vector<uint> shape = pars.get_shape();
    assert(shape[0] == 2 * (dim - 1));
    for(uint i = 0; i < dim-1; i++){
      if (j == 0){
      newparsalpha.push_back(pars[i]);
      newparsbeta.push_back(pars[i+dim-1]);
      }
      else{
	newparsalpha[i] = pars[i];
	newparsbeta[i] = pars[i+dim-1];
      }
    }
     densities[j].set_parameters(newparsalpha,newparsbeta);
  }
}

void GDMdensities::estimate(vector<MDArray<double> > & ess)
{
  assert(!ess.empty());
  int n = ess.size();
  int p;
  for(uint j = 0;j < densities.size();j++){
    densities[j].reset_cvec();
  }
  MDArray<double> tmp;
  
  for(int i = 0; i < n; i++){
    tmp = ess[i];
    p = tmp[0];
    vector<double> cvec;
    for(uint j = 0; j < dim; j++){
      cvec.push_back(tmp[j + 1]);
    }
    densities[p].add_countv(cvec);
  }


  for(int j = 0;j < node_size;j++){
    
    vector< vector<double> > estpars = densities[j].metro_hast(metro_samples,init[j],false);
    MDArray<double> pars;
    
    pars.set_shape(2 * (dim - 1));
    
    for(uint i = 0; i < dim - 1; i++){
      pars[i] = estpars[0][i];
      pars[i + dim - 1] = estpars[1][i];
    }
    if (cpd.size() < node_size){
      for(int i = 0;i<node_size;i++){
	cpd.push_back(pars);
      }
    }
    cpd[j] = pars;
  }
  update_pars();
}
  
  void GDMdensities::set_SumOfCounts(int s){
    sumOfCounts = s;
  }

  void GDMdensities::set_metro_samples(int m){
    metro_samples = m;
  }

vector<double> GDMdensities::sample(vector<double> & pv){
  assert(!cpd.empty());
  uint p=(unsigned int) (pv[0]+0.1);

  MDArray<double> pars = cpd[p];
  vector<double> alphas;
  vector<double> betas;
  
  for(uint j = 0; j < dim-1; j++)
    {
      alphas.push_back(pars[j]);
      betas.push_back(pars[j + (dim - 1)]);
    }
  vector<double> sample;
  sample = densities[p].Sample_mult(alphas,betas,sumOfCounts);
  return sample;
}

double GDMdensities::get_lik(vector<double> & ptv, bool log_space)
{

  uint p=(unsigned int) (ptv[0]+0.1);

  // change container
  vector<double> newv;
  for(int i = 1;i<ptv.size();i++){
    newv.push_back(ptv[i]);
  }

  // get GDM loglik log(P(alpha|X))
  double ll = densities[p].loglik(newv);
  if (log_space)
    return ll;
  else
    return exp(ll);
}

vector<MDArray<double> > GDMdensities::get_parameters()
{
  return cpd;
}

ostream& operator<<(ostream& output, const GDMdensities& a)  {
    output << "Node: Generalized Dirichlet Multinomial, size: " << a.node_size << " " << a.dim << endl;

    if (a.cpd.empty()){
      output << "Uninitialized." << endl;
    }
    else{
      for(uint i = 0;i < a.cpd.size();i++){
        output << a.cpd[i].tostring();
      }
    }
    return output;
}

}
