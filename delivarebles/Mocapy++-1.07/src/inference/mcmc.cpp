/*
 * mcmc.cpp
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

#include "mcmc.h"

using namespace std;

namespace mocapy {

MCMC::MCMC(DBN * new_dbn) {
	dbn = new_dbn;
	total_output_size = dbn->nodes_0.size();
}

void MCMC::initialize(MDArray<eMISMASK> & mismask, uint start, uint end) {
	// Initialize unobserved nodes to random values.

	// mismask: array that is used to detect missing values in
	// the sequence. mismask[length, index] should return false for
	// missing values, and true for observed values.
	// @type mismask: numpy array, shape=(length, nr of nodes)

	// start: start position in sequence
	// end: end position in sequence


	for (uint i = start; i < end; i++) {
		vector<Node*> nodes;
		if (i == 0) {
			nodes = dbn->nodes_0;
		} else {
			nodes = dbn->nodes_1;
		}

		// Loop over nodes
		for (uint n = 0; n < dbn->nr_nodes; n++) {
			Node* node = nodes[n];
			if (mismask.get(i, n) == MOCAPY_HIDDEN) {
				// Initialize unobserved node value to random state
				node->sample(i);
			}
		}
	}
}

}
