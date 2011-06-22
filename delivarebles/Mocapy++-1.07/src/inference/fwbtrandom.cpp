/*
 * fwbtrandom.cpp
 *
 *  Copyright (C) 2009, Wouter Boomsma, Martin Paluszewski, The Bioinformatics Centre, University of Copenhagen.
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

#include "fwbtrandom.h"
#include "infenginehmm.h"

using namespace std;

namespace mocapy {

// Constructor
FwbtRandom::FwbtRandom(DBN* new_dbn) :
	MCMC(new_dbn) {
}

// Sample values for unobserved nodes 
pair<double, uint> FwbtRandom::sweep(uint mcmc_steps, bool burn_in_flag,
                                     Sequence &seq, MDArray<eMISMASK> & mismask, uint start, uint end) {
        //      Draw a sample using forward backtrack sampling

	//	mcmc_steps: Not applicable.

	//	burn_in_flag: Not applicable

	//	mismask: array that is used to detect missing values in
	//	the sequence. mismask[length, index] should return 0 for
	//	missing values, and 1 for observed values.

	//	start: start position in sequence
	//	end: end position in sequence


	//	Make some local copies for speed.
	vector<Node*> nodes_0 = dbn->nodes_0;
	vector<Node*> nodes_1 = dbn->nodes_1;

        // Create new sampler object
        SampleInfEngineHMM sampler = SampleInfEngineHMM(dbn, seq, mismask, 0, false);

        // Sample values
        Sequence &sampled_seq = sampler.sample_next();

        // Iterate over sequence
        for (uint i=0; i<sampled_seq.get_shape()[0]; i++) {

             // Iterate over nodes
             // k: node index
             // l: output-column index
             for (uint k = 0, l=0; k < dbn->nr_nodes; k++,l++) {

                  Node* node;
                  if (i == 0) {
                       node = nodes_0[k];
                  } else {
                       node = nodes_1[k];
                  }

                  eMISMASK mm;
                  if (mismask.get_shape()[0] == 1) {
                       // Mismask is only specified for a slice
                       mm = mismask.get(0, k);
                  } else {
                       mm = mismask.get(i, k);
                  }

                  // Insert sampled values in parent map
                  if (mm == MOCAPY_HIDDEN) {
                       std::vector<double> values;
                       for (uint j=0; j<node->get_output_size(); j++) {
                            values.push_back(sampled_seq.get(i,l+j));
                       }
                       node->parentmap.set(i, values);
                  } 

                  // Update ESS
                  if (mm != MOCAPY_MISSING) {
                       // update ESS
                       if (!node->fixed) {
                            vector<double> ptv;
                            node->parentmap.get(i, ptv);
                            node->update_ess(ptv);
                       }
                  }

                  // Adjust output column count l
                  for (uint j=0; j<node->get_output_size(); j++) {
                       if (j>0) {
                            l += 1;
                       }
                  }
             }
        }

	// Sequence is sampled - save the ESS
        for (uint i = 0; i < dbn->unique_nodes.size(); i++) {
             Node* node = dbn->unique_nodes[i];
             if (!node->fixed) {
                  node->save_ess();
             }
        }

        // Likelihood is stored internally in SampleInfEngineHMM object
        double loglik = sampler.fwbt.calc_ll();
        uint slice_count = end - start;

	return make_pair(loglik, slice_count);
}

}

