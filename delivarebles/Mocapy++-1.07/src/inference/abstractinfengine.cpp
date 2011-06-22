/*
 * abstractinfengine.cpp
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

#include "abstractinfengine.h"

using namespace std;

namespace mocapy {

AbstractInfEngine::AbstractInfEngine(DBN * new_dbn, Sequence * new_seq,
		double new_weight) :
	dbn(new_dbn), weight(new_weight), seq(new_seq) {

	assert(!seq->empty());
	seq_len = seq->get_shape().front();
	output_size = dbn->total_output_size;
}

vector<ParentMap> AbstractInfEngine::make_parentmap_list(Sequence * seq,
		double weight) {
	//	seq: sequence
	// 	weight: sequence weight

	vector<ParentMap> p;
	for (uint i = 0; i < dbn->unique_nodes.size(); i++) {
		Node* node = dbn->unique_nodes[i];
		ParentMap pm = node->get_parentmap(seq, weight);
		p.push_back(pm);
	}
	return p;
}

void AbstractInfEngine::set_parentmap(vector<ParentMap> & p) {
	// Put previously calculated L{ParentMap} objects (ie. for fast node value
	// lookup in sequences) into the nodes.
	assert(p.size() == dbn->unique_nodes.size());
	for (uint i = 0; i < dbn->unique_nodes.size(); i++) {
		Node* node = dbn->unique_nodes[i];
		ParentMap pm = p[i];
		node->set_parentmap(&pm);
	}
}

pair<Sequence, double> AbstractInfEngine::get_viterbi() {
	// Calculate the Viterbi path. This method can act on a set of samples
	// (resulting in a stochastic Viterbi path calculation) or on all possible
	// hidden node values for each slice (resulting in a deterministic Viterbi
	// path calculation).

	vector<uint> shs = seq->get_shape();
	Sequence core_seq(shs);
	vector<ParentMap> parentmap_list = make_parentmap_list(&core_seq, weight);
	set_parentmap(parentmap_list);

	Sequence viterbi;
	viterbi.set_shape(seq_len, output_size);

	MDArray<double> gain;
	gain.set_shape(seq_len, nr_slices);

	Sequence path;
	path.set_shape(seq_len-1, nr_slices);

	// Initialization
	// l=0
	for (uint i = 0; i < nr_slices; i++) {
		Sequence & slice = seq_list[i].get_view(0);
		core_seq.set(0, slice);
		double ll(0);
		for (uint n = 0; n < dbn->nodes_0.size(); n++) {
			Node* node = dbn->nodes_0[n];
			ll += node->get_slice_log_likelihood(0);
		}

		gain.set(0, i, ll);
	}
	// Propagation
	// l>0

	for (uint i = 1; i < seq_len; i++) {
		for (uint j = 0; j < nr_slices; j++) {
			Sequence & slice = seq_list[j].get_view(i);
			core_seq.set(i, slice);
			vector<double> m(nr_slices);
			for (uint k = 0; k < nr_slices; k++) {
				Sequence & cs = seq_list[k].get_view(i - 1);
				core_seq.set(i - 1, cs);
				double ll(0);
				for (uint n = 0; n < dbn->nodes_1.size(); n++) {
					Node* node = dbn->nodes_1[n];
					ll += node->get_slice_log_likelihood(i);
				}
				m[k] = gain.get(i - 1, k) + ll;
			}
			assert(!m.empty());
			uint index = argmax(m);
			gain.set(i, j, m[index]);
			path.set(i - 1, j, index);
		}
	}
	// Termination
	vector<double> & g = gain.get_view(seq_len - 1).get_values();
	uint index = argmax(g);
	Sequence & view = seq_list[index].get_view(seq_len - 1);
	viterbi.set(seq_len - 1, view);

	// Backtracking
	assert(seq_len>1);
	for (int i = seq_len - 2; i >= 0; i--) {
		index = (int) path.get((uint) i, index);
		Sequence & seq_view = seq_list[index].get_view((uint) i);
		viterbi.set(i, seq_view);
	}

	double ll = gain.get_view(seq_len - 1).get_max() / (double) seq_len;
	return make_pair(viterbi, ll);
}

}
