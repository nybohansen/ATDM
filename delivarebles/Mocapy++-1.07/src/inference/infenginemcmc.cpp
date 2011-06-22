/*
 * infenginemcmc.cpp
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

#include "infenginemcmc.h"

using namespace std;

namespace mocapy {

InfEngineMCMC::InfEngineMCMC(DBN* new_dbn, MCMC* new_sampler,
		Sequence * new_seq, MDArray<eMISMASK> & new_mismask, double new_weight) :
	AbstractInfEngine(new_dbn, new_seq, new_weight) {
	// dbn: DBN definition
	// sampler: sampler object that implements MCMC
	// seq: sequence data
	// mismask: mismask array
	sampler = new_sampler;
	mismask = new_mismask;

	// List of Family objects - one for each node
	// The Family objects tie the data and the nodes together
	parentmap_list = make_parentmap_list(seq, weight);
	is_initialized = false;
}

void InfEngineMCMC::initialize_sample_generator(uint burn_in_steps,
		bool init_random, MDArray<eMISMASK> new_mismask, uint new_start,
		int new_end) {
	if (!new_mismask.get_shape().empty())
		mismask = new_mismask;

	assert(!mismask.get_shape().empty());

	set_parentmap(parentmap_list);
	g_end = 0;
	if (new_end < 0)
		g_end = seq_len;

	g_start = new_start;

	assert(0 <= g_start && g_start < seq_len);
	assert(g_start < g_end && g_end <= seq_len);

	//  Initialize hidden nodes with random values?
	if (init_random) {
		sampler->initialize(mismask, g_start, g_end);
	}

	// Burn in?
	if (burn_in_steps > 0) {
		burn_in_flag = true;
		sampler->sweep(burn_in_steps, true, *seq, mismask, g_start, g_end);
	}
	is_initialized = true;
}

Sample InfEngineMCMC::sample_next() {
	// Put the original sequence back in the nodes
	assert(is_initialized);

	// (because calling viterbi can screw things up!)
	set_parentmap(parentmap_list);
	pair<double, uint> ll_sc =
             sampler->sweep(1, false, *seq, mismask, g_start, g_end);
	return Sample(*seq, ll_sc.first, ll_sc.second);
}

void InfEngineMCMC::initialize_viterbi_generator(uint new_mcmc_steps,
		uint burn_in_steps, bool init_random, MDArray<eMISMASK> mismask,
		uint start, int end, bool restart) {
	// The Viterbi path generator's 'next' method will return a sampled sequence.
	// The sequence will be sampled between 'start' and 'end'. Burn in
	// (burn_in_steps) and random value initialiazation of the hidden
	// node values is done only one.

	// When the 'next' method of the Viterbi path generator is called the
	// first time, 'mcmc_steps' sequences are sampled, and an initial
	// Viterbi path is calculated using these sequences. Subsequent calls
	// to 'next' sample again 'mcmc_steps' sequences, and add the Viterbi
	// path of the previous step to the list. Then, the next Viterbi path
	// is calculated. Subsequent calls to 'next' thus try to improve on
	// the previously generated Viterbi path. After a while, the returned
	// LogLik should become stationary, ie. the the procedure reached a
	// (possible local) minimum.
	// mcmc_steps: number of sequences used for the calculation of
	// the Viterbi path each time 'next' is called.
	// burn_in_steps: number of discarded burn-in steps
	// mismask: mismask array
	// start: start position in the sequence
	// end: end position in the sequence
	// init_random: initialize to random values flag
	// restart: if 1, re-initialize hidden nodes and do burn in for
	// each call to next

	mcmc_steps = new_mcmc_steps;
	prev_seq.clear();

	// The inferred sequences
	initialize_sample_generator(burn_in_steps, init_random, mismask, start, end);
}

Sample InfEngineMCMC::viterbi_next() {
	if (!prev_seq.empty()) {
		seq_list.clear();
		seq_list.push_back(prev_seq);
	} else
		seq_list.clear();

	for (uint i = 0; i < mcmc_steps; i++) {
		// Generate another batch of sampled sequences
		// Sampling starts where previous Viterbi call left off
		Sample s = sample_next();
		seq_list.push_back(s.seq);
	}

	nr_slices = seq_list.size();
	pair<Sequence, double> viterbi_ll = get_viterbi();
	prev_seq = viterbi_ll.first;
	return Sample(viterbi_ll.first, viterbi_ll.second, 0);
}

}
