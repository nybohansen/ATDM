/*
 * gibbsrandom.cpp
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

#include "gibbsrandom.h"

using namespace std;

namespace mocapy {

GibbsRandom::GibbsRandom(DBN* new_dbn) :
	MCMC(new_dbn) {
	first_sweep = true;
}

vector<pair<uint, uint> > GibbsRandom::traverse(uint nr_nodes, uint start,
		uint end) {
	//	Create a list of (slice index, node index) tuples
	// 	that covers all possible slice/node combinations
	//	in random order.

	//	nr_nodes: the number of nodes per slice
	//	start: start position in sequence
	//	end: end position in sequence

	vector<pair<uint, uint> > v;
	v.reserve((end - start) * dbn->nr_nodes);

	for (uint i = start; i < end; i++) {
		for (uint n = 0; n < dbn->nr_nodes; n++) {
			v.push_back(make_pair(i, n));
		}
	}

	//	randomize list
	if (first_sweep) {
		// First sweep - no shuffle
		first_sweep = false;
	} else {
		// After first sweep - shuffle
		// random_shuffle(l.begin(), l.end());
		random_shuffle(v.begin(), v.end());
	}
	return v;
}

pair<double, uint> GibbsRandom::sweep(uint mcmc_steps, bool burn_in_flag,
                                      Sequence &seq, MDArray<eMISMASK> & mismask, uint start, uint end) {
	//	Do one mcmc_steps of Gibbs sampling, and burn_in_steps
	//	of Gibbs sampling burn-in.

	//	mcmc_steps: number of Gibbs sampling sweeps through the data.

	//	burn_in_flag: if 1, we are doing burn in

	//	mismask: array that is used to detect missing values in
	//	the sequence. mismask[length, index] should return 0 for
	//	missing values, and 1 for observed values.

	//	start: start position in sequence
	//	end: end position in sequence

	//	Make some local copies for speed.
	vector<Node*> nodes_0 = dbn->nodes_0;
	vector<Node*> nodes_1 = dbn->nodes_1;
	uint nr_nodes = dbn->nr_nodes;

	//	Likelihood: sum[logP(H,O)]
	//	where sum runs over all samples
	//	P(H,O)

	double loglik = 0.0;
	uint slice_count = end - start;

	for (uint steps = 0; steps < mcmc_steps; steps++) {
		// Generate traversal
		vector<pair<uint, uint> > traversal = traverse(nr_nodes, start, end);

		for (uint i = 0; i < traversal.size(); i++) {
			pair<uint, uint> & trav = traversal[i];
			uint k = trav.first;
			uint n = trav.second;
			Node* node;
			if (k == 0) {
				node = nodes_0[n];
			} else {
				node = nodes_1[n];
			}

			// Loop over nodes
			uint index = node->node_index;
			eMISMASK mm;
			if (mismask.get_shape()[0] == 1) {
				// Mismask is only specified for a slice
				mm = mismask.get(0,index);
			} else {
				mm = mismask.get(k, index);
			}


			if (mm == MOCAPY_HIDDEN) {
//			if (mm != MOCAPY_OBSERVED) {
				// Unobserved node -- do Gibbs sampling
				node->blanket_sample(k, mismask);
			}
			if (!burn_in_flag && mm != MOCAPY_MISSING) {
//			if (!burn_in_flag ) {
				// Update likelihood
				pair<double, vector<double> > ll_fam =
						node->get_slice_log_likelihood_ptv(k);
				loglik += ll_fam.first;
				// update ESS
				if (!node->fixed) {
					node->update_ess(ll_fam.second);
				}
			}
		}
	}

	// Sequence is sampled - save the ESS
	if (!burn_in_flag) {
		for (uint i = 0; i < dbn->unique_nodes.size(); i++) {
			Node* node = dbn->unique_nodes[i];
			if (!node->fixed) {
				node->save_ess();
			}
		}
		loglik /= (double) mcmc_steps;
	} else {
		loglik = 0;
	}
	return make_pair(loglik, slice_count);
}

}
