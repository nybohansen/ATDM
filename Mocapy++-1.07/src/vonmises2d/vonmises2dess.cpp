/*
 * vonmises2dess.cpp --- Bivariate von Mises distribution, cosine variant
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


#include "vonmises2dess.h"

using namespace std;

namespace mocapy {

// Initializes the ESS arrays.
void VonMises2dESS::construct(vector<uint> & parent_sizes, uint output_dim, uint node_size) {
	this->output_dim = output_dim;
	this->node_size = node_size;
	this->ess_size = ESSSIZE;

	ess_shape.resize(node_size);
	for (uint i=0; i<node_size; i++) {
	     ess_shape[i].set_shape(ess_size);
	}
	// ess_shape[_2COS1].set_shape(ess_size);
	// ess_shape[_2SIN1].set_shape(ess_size);
	// ess_shape[_2COS2].set_shape(ess_size);
	// ess_shape[_2SIN2].set_shape(ess_size);
	// ess_shape[COS1MIN2].set_shape(ess_size);
	// ess_shape[SIN1MIN2].set_shape(ess_size);
	// ess_shape[COS1PLUS2].set_shape(ess_size);
	// ess_shape[SIN1PLUS2].set_shape(ess_size);
	// ess_shape[C1].set_shape(ess_size);
	// ess_shape[S1].set_shape(ess_size);
	// ess_shape[C2].set_shape(ess_size);
	// ess_shape[S2].set_shape(ess_size);
	// ess_shape[COUNT].set_shape(ess_size);

	clear();
}

// Add another sample to the ESS
void VonMises2dESS::add_ptv(vector<double> ptv) {
	assert(ptv.size() == 3);

	double *node_value = &ptv[1];
	uint parent_value = (uint) ptv[0];

	ess[parent_value][_2COS1] += cos(2*node_value[0]);
	ess[parent_value][_2SIN1] += sin(2*node_value[0]);
	ess[parent_value][_2COS2] += cos(2*node_value[1]);
	ess[parent_value][_2SIN2] += sin(2*node_value[1]);
	ess[parent_value][COS1MIN2] += cos(node_value[0] - node_value[1]);
	ess[parent_value][SIN1MIN2] += sin(node_value[0] - node_value[1]);
	ess[parent_value][COS1PLUS2] += cos(node_value[0] + node_value[1]);
	ess[parent_value][SIN1PLUS2] += sin(node_value[0] + node_value[1]);
	ess[parent_value][C1] += cos(node_value[0]);
	ess[parent_value][S1] += sin(node_value[0]);
	ess[parent_value][C2] += cos(node_value[1]);
	ess[parent_value][S2] += sin(node_value[1]);
	ess[parent_value][COUNT] += 1;
}


void VonMises2dESS::save_ess(double weight) {

     uint node_size = ess.size();
     
     if (weight != 1.0) {
	  for (uint i=0; i<node_size; i++) {
	       for (uint j=0; j<ess_size; j++) {
		    ess[i][j] *= weight;
	       }
	  }
     }

     for (uint i=0; i<node_size; i++) {
	  for (uint j=0; j<ess_size; j++) {
	       saved_ess[i][j] += ess[i][j];
	  }
     }

     // Get ready for new ESS
     ess = ess_shape;
}
     
}
