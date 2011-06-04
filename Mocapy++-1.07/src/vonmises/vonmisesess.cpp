/*
 * vonmisesess.cpp
 *
 *  Copyright (C) 2008, Jes Frellsen, The Bioinformatics Centre, University of Copenhagen.
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


#include "vonmisesess.h"

using namespace std;

namespace mocapy {

// Initializes the ESS arrays.
void VonMisesESS::construct(vector<uint> & parent_sizes, uint new_output_dim, uint new_node_size) {
	output_dim = new_output_dim;
	node_size = new_node_size;
	ess_size = VM_ESSSIZE;

	ess_shape.resize(VM_ESSSIZE);
	ess_shape[VM_RX].set_shape(node_size);
	ess_shape[VM_RY].set_shape(node_size);
	ess_shape[VM_N].set_shape(node_size);

	clear();
}

// Add another sample to the ESS
void VonMisesESS::add_ptv(vector<double> ptv) {
	assert(ptv.size() == 2);

	double node_value = ptv[1];
	uint parent_value = TOINT(ptv[0]);

	assert(ess[VM_RX].size() > parent_value);
	ess[VM_RX][parent_value] += cos(node_value);
	assert(ess[VM_RY].size() > parent_value);
	ess[VM_RY][parent_value] += sin(node_value);
	assert(ess[VM_N].size() > parent_value);
	ess[VM_N][parent_value] += 1;
}

}
