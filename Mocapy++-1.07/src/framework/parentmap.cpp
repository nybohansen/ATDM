/*
 * ParentMap.cpp
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

// Fast parent node value/node value lookup for discrete nodes.

#include "parentmap.h"

using namespace std;

namespace mocapy {

ParentMap::ParentMap(Sequence * new_seq, vector<uint> & parents_0,
		vector<uint> & parents_1, uint new_this_index, double new_weight, uint new_dim) {
	// seq: the data sequence (2-dim)
	// parents_0: index of parent value in previous slice
	// parents_1: index of parent value in current slice
	// this_index: index of node value itself in current slice
	// weight: sequence weight

	assert(!new_seq->empty());
	lng = new_seq->get_shape().front();
	seq = new_seq->flat_iterator();
	dim = new_dim;

	// Step size
	vector<uint> shape = new_seq->get_shape();
	assert(shape.size() == 2);
	step = new_seq->get_shape()[1];

	// Self index
	// Cache the indices
	indices.clear();

	for (uint i = 0; i < lng; i++) {
		vector<uint> index;
		vector<uint*> index2;
		if (i > 0) {
			for (uint j = 0; j < parents_0.size(); j++) {
				uint pi = parents_0[j];
				index.push_back((i - 1) * step + pi);
			}
		}

		for (uint j = 0; j < parents_1.size(); j++) {
			uint pi = parents_1[j];
			index.push_back(i * step + pi);
		}
		index.push_back(i * step + new_this_index);
		indices.push_back(index);
	}
	this_index = new_this_index;
	weight = new_weight;
}

double ParentMap::get_weight() {
	return weight;
}

void ParentMap::set(uint h, MDArray<double> & v) {
	vector<uint> & v1 = indices[h];
	uint i = v1.back();
	assert(v.get_shape().size() == 1);
	uint dim = v.get_shape().front();
	uint k = 0;
	for (uint j = i; j < i + dim; j++, k++) {
		*(seq[j]) = v[k];
	}
}

void ParentMap::set(uint h, vector<double> & v) {
	vector<uint> & v1 = indices[h];
	uint i = v1.back();
	uint dim = v.size();
	uint k = 0;
	for (uint j = i; j < i + dim; j++, k++) {
		*(seq[j]) = v[k];
	}
}

void ParentMap::replace_seq(Sequence * new_seq) {
	// seq: the data sequence (2-dim)
	assert(!new_seq->empty());

	// Check that then new sequence is compatible with the old sequence
	vector<uint> new_shape = new_seq->get_shape();
	assert(new_shape.size() == 2);
	assert(new_shape[0] == lng);
	assert(new_shape[1] == step);

	// Replace the old sequence by the new
	seq = new_seq->flat_iterator();
}


}
