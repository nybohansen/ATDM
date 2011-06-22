/*
 *  GDdensities.cpp
 *
 *  Copyright (C) 2008, Simon H. Rasmussen, The Bioinformatics Centre,
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
#include "GDdensities.h"

namespace mocapy
{

GDdensities::GDdensities()
{
    initialize();
}

  void GDdensities::initialize()
  {
    type = GD;
    output_size = dim;
  }

  
  GDdensities::GDdensities(uint dim,uint sample_size){
    this->dim = dim;
    this->sample_size = sample_size;
    initialize();
  }
  
void GDdensities::set_parameters(MDArray<double> & pars)
{
    // Check cpd is right shape
    vector<uint> shape = pars.get_shape();
    assert(shape[1] == 2 * (dim - 1));

    assert(shape[0] == node_size);
    for(uint i = 0;i < node_size;i++){
      cpd.push_back(pars.get_view(i));
    }
    if (!cpd.empty()){
      assert(!cpd.empty());
    }
    else{
      printf("cpd is empty \n");
    }
}

void GDdensities::construct(vector<uint> & parent_sizes){
    // GDNode always has a single parent
    assert(parent_sizes.size() == 1);
    node_size = parent_sizes[0];
    // vector of GDDensity objects
    densities.resize(node_size);
    _create_densities();
}

void GDdensities::_create_densities()
{
  first = true;
  ulong seed = 450;
  MDArray<double> pars;
  vector<double> alphas;
  vector<double> betas;
  // Initialize the density objects and fill the density vector
  for(uint i=0; i<node_size; i++)
    {
      densities[i] = GDdensity(dim, seed, sample_size);
      init.push_back(true);
    }
}

void GDdensities::update_pars(){
  //Update parameters in density class
  assert(!cpd.empty());
  MDArray<double> pars;
  vector<double> newparsalpha;
  vector<double> newparsbeta;
  for(uint j = 0; j < node_size; j++){
    pars = cpd[j];
    vector<uint> shape = pars.get_shape();
    assert(shape[0]==2 * (dim - 1));
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
     densities[j].GDset_parameters(newparsalpha,newparsbeta);
  }
}

void GDdensities::estimate(vector<MDArray<double> > & ess){
  assert(!ess.empty());
  int n = ess.size();

  if(!cpd.empty()){
    update_pars();
  }

  for(uint j = 0; j < node_size; j++){
    densities[j].reset_pvec();
  }
  int p;
  MDArray<double> tmp;
  vector<double> pvec;
  for(uint j = 1; j < dim + 1; j++){
    pvec.push_back(0);
  }
  
  for(int i = 0; i < n; i++){
    tmp = ess[i];
    p = tmp[0];
    for(uint j = 1; j < dim + 1; j++){
      pvec[j-1]=tmp[j];
    }
    double sum = 0;
    for(int k = 0;k < dim;k++){
      sum = sum + pvec[k];
    }
    assert(sum > .99 && sum < 1.01);
    densities[p].GDadd_pvec(pvec);
  }
  for(uint j = 0;j<node_size;j++){
    vector< vector<double> > estpars = densities[j].GDestimate_pars(init[j]);
  
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
    else{
      cpd[j] = pars;
    }
  }
  update_pars();
}

vector<double> GDdensities::sample(vector<double> & pv){
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
  sample = densities[p].GDSample_Gendir(alphas,betas);
  return sample;
}

double GDdensities::get_lik(vector<double> & ptv, bool log_space){

    uint p=(unsigned int) (ptv[0]+0.1);
    vector<double> newv;
    for(int i = 1;i<ptv.size();i++){
      newv.push_back(ptv[i]);
    }
  
    double ll = densities[p].GDloglik(newv,init[p]);
    if (init[p] == true){
      init[p] = false;
    } 
    if (log_space)
      return ll;
    else
      return exp(ll);
}

vector<MDArray<double> > GDdensities::get_parameters()
{
  return cpd;
}

ostream& operator<<(ostream& output, const GDdensities& a)  {
    output << "Node: Generalized Dirichlet, size: " << a.node_size << " " << a.dim << endl;

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
