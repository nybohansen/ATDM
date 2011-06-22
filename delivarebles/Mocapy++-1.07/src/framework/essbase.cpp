/*
 * essbase.h
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

#include "essbase.h"

using namespace std;

namespace mocapy {

// ESS has been used in M-step - get ready for new cycle
void ESSBase::clear() {
	assert(!ess_shape.empty());
	ess = ess_shape;
	saved_ess = ess_shape;
}

void ESSBase::save_ess(double weight) {
	if (weight != 1.0) {
		for (uint i=0; i<ess_size; i++) {
			ess[i].multiply_inplace(weight);
		}
	}

	for (uint i=0; i<ess_size; i++) {
		saved_ess[i].add_inplace(ess[i]);
	}

	// Get ready for new ESS
	ess = ess_shape;
}

// Return a vector of MDArrays containing the ESS
// This will be used in the density object for estimation
vector< MDArray<double> > ESSBase::get_array() {
	return saved_ess;
}

void ESSBase::combine(ESSBase & other_ess) {
	// For parallelization
	// TODO: Implement
	assert(false);
}

}
