/*
 * sampleinfengine.cpp
 *
 *  Copyright (C) 2009, Jes Frellsen, The Bioinformatics Centre, University of Copenhagen.
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

#include "sampleinfengine.h"

namespace mocapy {

SampleInfEngine::SampleInfEngine(DBN* new_dbn, Sequence & new_seq, MDArray<eMISMASK> & new_mismask, double new_weight) :
	MoreAbstractInfEngine(new_dbn), seq(new_seq), mismask(new_mismask), weight(new_weight) {

	randomGen = new_dbn->randomGen;

	// Check that seq and mismask are consistent
	check_seq_and_mismask(seq, mismask);

	// Make the parent map for the nodes
	parentmap_list = make_parentmap_list(&seq, weight);

	// Set the redefined parent map into the nodes
	set_parentmap(parentmap_list);

	// Set position specific variables
	seq_len = seq.get_shape().front();
	start = 0;
	end = seq_len;
}

void SampleInfEngine::set_seq_mismask(Sequence & new_seq, MDArray<eMISMASK> & new_mismask) {
	// Check that seq and mismask are consistent
	check_seq_and_mismask(new_seq, new_mismask);
	if (new_seq.get_shape() != seq.get_shape())         throw MocapyExceptions("SampleInfEngine: The new sequence must have same shape as the old sequence!");

	// Replace the this->seq and this->mismask
	seq = new_seq;
	mismask = new_mismask;

	// Disable undo
	prev_seq = Sequence();

	// Replace the current sequence with the previous in the parent map
	replace_seq_in_paremtmap_list(seq);

	// Set the redefined parent map into the nodes
	set_parentmap(parentmap_list);
}

void SampleInfEngine::undo() {
	// Check that there is a previous sequence
	if (prev_seq.empty())
		throw MocapyExceptions("SampleInfEngine: undo can only be called once for each call to sample!");

	// Replace the current sequence with the previous
	seq = prev_seq;
	prev_seq = Sequence();

	// Replace the current sequence with the previous in the parent map
	replace_seq_in_paremtmap_list(seq);

	// Set the redefined paremt map into the nodes
	set_parentmap(parentmap_list);
}

}
