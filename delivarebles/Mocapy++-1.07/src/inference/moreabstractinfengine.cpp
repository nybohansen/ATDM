/*
 * moreabstractinfengine.cpp
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

#include "moreabstractinfengine.h"

using namespace std;

namespace mocapy {

MoreAbstractInfEngine::MoreAbstractInfEngine(DBN * new_dbn) :
	dbn(new_dbn) {}

vector<ParentMap> MoreAbstractInfEngine::make_parentmap_list(Sequence * seq, double weight) {
	//	seq: sequence
	// 	weight: sequence weight

	vector<ParentMap> p;
	for (vector<Node*>::iterator node = dbn->unique_nodes.begin(); node < dbn->unique_nodes.end(); node++) {
		ParentMap pm = (*node)->get_parentmap(seq, weight);
		p.push_back(pm);
	}
	return p;
}

void MoreAbstractInfEngine::set_parentmap(vector<ParentMap> & p) {
	// Put previously calculated L{ParentMap} objects (ie. for fast node value
	// lookup in sequences) into the nodes.
	assert(p.size() == dbn->unique_nodes.size());
	for (uint i = 0; i < dbn->unique_nodes.size(); i++) {
		Node* node = dbn->unique_nodes[i];
		ParentMap *pm = &p[i];
		node->set_parentmap(pm);
	}
}

void MoreAbstractInfEngine::set_parentmap(Sequence * seq, double weight) {
	for (vector<Node*>::iterator node = dbn->unique_nodes.begin(); node < dbn->unique_nodes.end(); node++) {
		ParentMap pm = (*node)->get_parentmap(seq, weight);
		(*node)->set_parentmap(&pm);
	}
}

void MoreAbstractInfEngine::check_seq_and_mismask(Sequence & seq, MDArray<eMISMASK> & mismask) {
	if (seq.empty() or mismask.empty()) throw MocapyExceptions("LikelihoodInfEngineHMM: sequence and mismask cannot be empty!");
	if (seq.get_shape().size()!=2)      throw MocapyExceptions("LikelihoodInfEngineHMM: The sequence must have exact 2 dimensions!");
	if (mismask.get_shape().size()!=2)  throw MocapyExceptions("LikelihoodInfEngineHMM: The mismask must have exact 2 dimensions!");
	if (seq.get_shape()[0] != mismask.get_shape()[0])  throw MocapyExceptions("LikelihoodInfEngineHMM: sequence and mismask must have same length!");
	

	if (seq.get_shape()[1] != dbn->total_output_size) {
	  std::cout << "seq shape: " << seq.get_shape()[1] << std::endl;
	  std::cout << "total output size: " << dbn->total_output_size << std::endl;
	  throw MocapyExceptions("LikelihoodInfEngineHMM: The sequence must have same height as the total output size!");
	}
	if (mismask.get_shape()[1] != dbn->nodes_0.size()) throw MocapyExceptions("LikelihoodInfEngineHMM: The mismask must have same height as number of nodes!");
}


}
