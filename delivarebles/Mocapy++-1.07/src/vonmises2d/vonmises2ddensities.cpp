/*
 * vonmises2ddensities.cpp
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


#include "vonmises2ddensities.h"
#include "vonmises2dess.h"

using namespace std;

namespace mocapy {

     // Constructor
     VonMises2dDensities::VonMises2dDensities(MDArray<double> user_mus,
					      MDArray<double> user_kappas) {
	  this->user_mus = user_mus;
	  this->user_kappas = user_kappas;
	  this->kappa_size = 3;
	  initialize();
     }

     // Initializer
     void VonMises2dDensities::initialize() {
	  this->output_size = 2;
	  this->type = VONMISES2D;
     }

     // Override normal setRandomGen
     void VonMises2dDensities::setRandomGen(RandomGen* rg) {
          DensitiesBase::setRandomGen(rg);

          // Set random generator in all components
	  for (uint i = 0; i < node_size; i++) {
               dist_list[i].setRandomGen(rg);
          }
     }


     MDArray<double> VonMises2dDensities::make_uniform_kappas(uint node_size) {
	  std::vector<uint> shape(vec(node_size, kappa_size));
	  MDArray<double> kappas(shape);
	  for (uint i = 0; i < node_size; i++) {
	       for (uint j = 0; j < kappa_size; j++) {
		   kappas.set(i,j, 0);
	       }
	  }
	  return kappas;
     }


     MDArray<double> VonMises2dDensities::make_uniform_mus(uint node_size) {
	  std::vector<uint> shape(vec(node_size, output_size));
	  MDArray<double> mus(shape);
	  for (uint i = 0; i < node_size; i++) {
	       for (uint j = 0; j < output_size; j++) {
		    mus.set(i,j, 0);
	       }
	  }
	  return mus;
     }

     // Called in node.construct, and initializes the density arrays
     void VonMises2dDensities::construct(vector<uint> & parent_sizes) {

	  // Assume that node only has one parent
	  assert(parent_sizes.size() == 1);

	  // Set node size
	  this->node_size = parent_sizes[0];

	  parent_index = parent_sizes.front();

	  std::vector<uint> mus_shape(vec(node_size, output_size));
	  std::vector<uint> kappas_shape(vec(node_size, kappa_size));

	  // Check or construct mus
	  if (!user_mus.empty()) {
	       mus = user_mus;
	       assert(mus.get_shape() == mus_shape);
	       user_mus.clear();
	  } else {
	       mus = make_uniform_mus(node_size);
	  }

	  // Check or construct kappas
	  if (!user_kappas.empty()) {
	       kappas = user_kappas;
	       assert(kappas.get_shape() == kappas_shape);
	       user_kappas.clear();
	  } else {
	       kappas = make_uniform_kappas(node_size);
	  }

	  // Make a list of VonMises objects for sampling and calculations of data likelihood
	  dist_list = make_dist_list(node_size, mus, kappas);
     }

     // List of distribution objects from which you can sample and calculate densities
     vector<VonMises2d> VonMises2dDensities::make_dist_list(uint node_size, MDArray<double> &mus, MDArray<double> &kappas) {
	  assert(mus.get_shape() == vec(node_size, output_size));
	  assert(kappas.get_shape() == vec(node_size, kappa_size));

	  vector<VonMises2d> dist_list(node_size);
	  for (uint i = 0; i < node_size; i++) {
	    assert(randomGen != NULL);
	    dist_list[i] = VonMises2d(kappas.get_view(i), mus.get_view(i), randomGen);
	  }

	  return dist_list;
     }



     // Parameter estimation based on the ESS
     void VonMises2dDensities::estimate(vector<MDArray<double> > & ess) {
	for (uint i=0; i<node_size; i++) {
             // Use Moment Estimation only
	     // std::vector<double> parameters = VonMises2dEstimate_ME(ess[i]);

             // Use Maximum Likelihood (using Moment Estimation as start values)
	     std::vector<double> parameters = VonMises2dEstimate_ML(ess[i]);

	     if (parameters.size() == 5) {
		  kappas.set(i, 0, parameters[0]);
		  kappas.set(i, 1, parameters[1]);
		  kappas.set(i, 2, parameters[2]);
		  mus.set(i,0, parameters[3]);
		  mus.set(i,1, parameters[4]);
		  dist_list[i].update(kappas.get_view(i), mus.get_view(i));
	     } else {
                  kappas.set(i, 0, 0);                  
                  kappas.set(i, 1, 0);                  
                  kappas.set(i, 2, 0);                  
                  mus.set(i, 0, 0);
                  mus.set(i, 1, 0);
                  dist_list[i].update(kappas.get_view(i), mus.get_view(i));
             }
	}
     }

     // Return a sample, based on indicated parent values
     vector<double> VonMises2dDensities::sample(vector<double> & pv) {

       assert(randomGen != NULL);
	  uint parent_value = (uint) pv.front();

	  VonMises2d &component = dist_list[parent_value];
	 

	  std::vector<double> sampleVector;

	  assert(component.rg != NULL);
	  component.sampler->rg = randomGen;

	  component.sampler->sample(sampleVector);


	  return sampleVector;
     }

     // Return likelihood, that is: P(child|parents)
     double VonMises2dDensities::get_lik(vector<double> & pv, bool log) {
	  uint parent_value = (uint) pv.front();

	  VonMises2d &component = dist_list[parent_value];
	  double values[2] = {pv[1], pv[2]};
	  if (log)
	       return component.density->get_log_likelihood(values);
	  else
	       return component.density->get_likelihood(values);
     }

     // Return the distribution's parameters
     vector<MDArray<double> > VonMises2dDensities::get_parameters() {
	  vector<MDArray<double> > parameters;
	  parameters.push_back(kappas);
	  parameters.push_back(mus);
	  return parameters;
     }


     ostream& operator<<(ostream& output, const VonMises2dDensities& a)  {
    	 output << a.mus << endl;
    	 output << a.kappas << endl;
     	return output;
     }

}

