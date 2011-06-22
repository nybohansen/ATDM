/*
 * discreteess.h
 *
 *  Copyright (C) 2008, Martin Paluszewski, The Bioinformatics Centre, University of Copenhagen.
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

#include "discreteess.h"

using namespace std;

namespace mocapy {

// Initializes the ESS arrays.
void DiscreteESS::construct(vector<uint> & parent_sizes, uint output_dim, uint node_size) {
	vector<uint> shape = vec_conc(parent_sizes, node_size);

	// Discrete ESS only has one MDArray
	ess_shape.resize(1);
	ess_shape[0].set_shape(shape);
	ess_size = 1;
	clear();
}

// Add another sample to the ESS
void DiscreteESS::add_ptv(vector<double> ptv) {
	vector<uint> int_indx;
	toint(ptv, int_indx);
	assert(!ess_shape.empty());

	ess[0][int_indx]++;
}

}
