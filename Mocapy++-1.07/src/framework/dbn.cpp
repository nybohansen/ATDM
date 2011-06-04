/*
 * DBN.cpp
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
#include <assert.h>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include "dbn.h"

#include "mocapyexceptions.h"
#include "parentmap.h"

using namespace std;

namespace mocapy {
  
DBN::DBN(uint seed) {
  randomGen = new RandomGen(); 
  if (seed>0) {
    myseed=seed;
  }
  else {
    // global seed
    myseed=moc_seed;
  }
  
  randomGen->mocapy_seed(myseed);

};


DBN::DBN(vector<Node*> & new_nodes_0, vector<Node*> & new_nodes_1,
	 string new_name, uint seed) :
	name(new_name) {

	if (new_nodes_0.size() != new_nodes_1.size())
		throw MocapyExceptions("Slices need to have the same number of nodes.");

	set_slices(new_nodes_0, new_nodes_1);

	randomGen = new RandomGen(); 
	if (seed>0) {
	  myseed=seed;
	}
	else {
	  // global seed
	  myseed=moc_seed;
	}
	
	randomGen->mocapy_seed(myseed);

}

void DBN::seed(uint s) {
  myseed=s;
  randomGen->mocapy_seed(s);
}

void DBN::setRandomGen(RandomGen* rg) {
  randomGen = rg; 
  // Number of nodes/slice
  nr_nodes = nodes_0.size();
  
  for (uint i = 0; i < nr_nodes; i++) {
    Node* n0 = nodes_0[i];
    Node* n1 = nodes_1[i];
    

    assert(randomGen);
    n0->setRandomGen(randomGen);
    n1->setRandomGen(randomGen);
  }
}

void DBN::set_slices(vector<Node*> new_nodes_0, vector<Node*> new_nodes_1) {
	nodes_0 = new_nodes_0;
	nodes_1 = new_nodes_1;

	is_constructed = false;

	// Number of nodes/slice
	nr_nodes = nodes_0.size();

	// Set the node indices
	uint data_index = 0;
	for (uint i = 0; i < nr_nodes; i++) {
		Node* n0 = nodes_0[i];
		Node* n1 = nodes_1[i];

		// Identifier of the node
		n0->set_node_index(i);
		n1->set_node_index(i);

		// Where is the node data found?
		n0->set_data_index(data_index);
		n1->set_data_index(data_index);
		assert(n0->get_output_size() == n1->get_output_size());
		// Move to next free index in data
		data_index += n1->get_output_size();

		
	}
	make_index_map(nodes_0, nodes_1);
	total_output_size = data_index;
}


void DBN::make_index_map(vector<Node*> & nodes_0, vector<Node*> & nodes_1) {
	// Create a dicationary that maps the name of a node to its index.
	// nodes_0: nodes in slice 0.
	// nodes_1: nodes in slice 1.

	for (uint i = 0; i < nodes_0.size(); i++) {
		Node* node = nodes_0[i];
		if (node->get_name() != "")
			index_map[node->get_name()] = make_pair(i, 0);
		node = nodes_1[i];
		if (node->get_name() != "")
			index_map[node->get_name()] = make_pair(i, 1);
	}
}

pair<uint, uint> DBN::map_to_indices(NodeID & parent_i, NodeID & child_i) {
	// Return node indices associated with the two given node names 'parent'
	// and 'child'. If a name is an integer, it is returned as is.

	uint p_i = parent_i.integer;
	uint c_i = child_i.integer;

	if (!parent_i.isInteger) {
		IndexMap::iterator it = index_map.find(parent_i.name);
		if (it == index_map.end()) {
		  std::cout << parent_i.name << std::endl;
		  throw MocapyExceptions("Unknown node name");
		}
		else
			p_i = it->second.first;
	}
	if (!child_i.isInteger) {
		IndexMap::iterator it = index_map.find(child_i.name);
		if (it == index_map.end()) {
		  std::cout << child_i.name << std::endl;
		  throw MocapyExceptions("Unknown node name");
		}
		else
			c_i = it->second.first;
	}

	return make_pair(p_i, c_i);
}

pair<vector<Node*> , vector<Node*> > DBN::get_nodes() {
	return make_pair(nodes_0, nodes_1);
}

Node* DBN::get_node_by_name(string name) {
	map<string, pair<int, int> >::iterator it = index_map.find(name);
	assert(it != index_map.end());
	int index = it->second.first;
	int slice = it->second.second;

	if (slice == 0)
		return nodes_0[index];
	else
		return nodes_1[index];
}

pair<Sequence, double> DBN::sample_sequence(uint length) {
	// Return a sampled sequence sequence from the DBN with specified length.
	assert(is_constructed);

	Sequence seq;
	seq.set_shape(length, total_output_size);

	for (uint i = 0; i < unique_nodes.size(); i++) {
		Node* node = unique_nodes[i];
		ParentMap pm = node->get_parentmap(&seq, 1);
		node->set_parentmap(&pm);
	}
	vector<Node*> node_list;
	double ll(0);

	for (uint i = 0; i < length; i++) {
		if (i == 0) {
			node_list = nodes_0;
		} else {
			node_list = nodes_1;
		}

		for (uint n = 0; n < node_list.size(); n++) {
			Node* node = node_list[n];
			node->sample(i);
			double loglik = node->get_slice_log_likelihood(i);
			ll += loglik;
		}
	}
	return make_pair(seq, ll / (double) length);
}

// Syntactic sugar
void DBN::add_inter(int parent_i, int child_i) {
	NodeID pi(parent_i);
	NodeID ci(child_i);
	add_inter(pi, ci);
}

void DBN::add_inter(const char* parent_i, const char* child_i) {
	NodeID pi(parent_i);
	NodeID ci(child_i);
	add_inter(pi, ci);
}

// Add an edge between slices, from parent to child.
void DBN::add_inter(NodeID & parent_i, NodeID & child_i) {
	assert(!is_constructed);

	// Map names to indices if necessary
	pair<int, int> p_and_c = map_to_indices(parent_i, child_i);
	int p_i = p_and_c.first;
	int c_i = p_and_c.second;

	// slice 0
	Node* parent_0 = nodes_0[p_i];

	// slice 1
	Node* parent_1 = nodes_1[p_i];
	Node* child_1 = nodes_1[c_i];

	// Add children to parents
	parent_0->add_inter_child(child_1);

	// Note that child_0 cannot have parents
	child_1->add_inter_parent(parent_0->data_index, parent_0->get_node_size());

	// Parents tied?
	if (parent_0 != parent_1) {
		parent_1->add_inter_child(child_1);
	}
}

// Syntactic sugar
void DBN::add_intra(int parent_i, int child_i) {
	NodeID pi(parent_i);
	NodeID ci(child_i);
	add_intra(pi, ci);
}

void DBN::add_intra(const char* parent_i, const char* child_i) {
	NodeID pi(parent_i);
	NodeID ci(child_i);
	add_intra(pi, ci);
}

void DBN::add_intra(NodeID & parent_i, NodeID & child_i) {
	// Add an edge inside a slice, from parent to child.
	assert(!is_constructed);
	// Map names to indices if necessary
	pair<uint, uint> p_and_c = map_to_indices(parent_i, child_i);
	uint p_i = p_and_c.first;
	uint c_i = p_and_c.second;

	// slice 0
	assert(p_i < nodes_0.size());
	Node* parent_0 = nodes_0[p_i];
	assert(c_i < nodes_0.size());
	Node* child_0 = nodes_0[c_i];

	// slice 1
	Node* parent_1 = nodes_1[p_i];
	Node* child_1 = nodes_1[c_i];

	// Add children to parents
	parent_0->add_intra_child(child_0);

	if (parent_0 != parent_1) {
		parent_1->add_intra_child(child_1);
	}

	// Add parents to children
	child_0->add_intra_parent(parent_0->data_index, parent_0->get_node_size());

	if (child_0 != child_1) {
		child_1->add_intra_parent(parent_0->data_index, parent_0->get_node_size());
	}
}

void DBN::construct() {
	// Initialize the DBN data structures based on the added nodes and edges.
	// After calling this method no edges can be added.
  	assert(!is_constructed);
	// construct nodes in slice 0
	for (uint n = 0; n < nr_nodes; n++) {
		// Slice 0
		Node* node_0 = nodes_0[n];
		Node* node_1 = nodes_1[n];

		assert(randomGen);
		node_0->setRandomGen(randomGen);
		node_1->setRandomGen(randomGen);

		if (node_0 == node_1) {
			// Tied node
			node_0->construct(TIED);
			unique_nodes.push_back(node_0);
		} else {
			node_0->construct(START);
			node_1->construct(END);
			unique_nodes.push_back(node_0);
			unique_nodes.push_back(node_1);
		}
	}
	all_nodes = vec_concNode(nodes_0, nodes_1);
	is_constructed = true;
}


void DBN::save(const char* name) {
	std::ofstream ofs(name);
	boost::archive::text_oarchive oa(ofs);
	oa << *this;

}
void DBN::load(const char* name) {
  try {
    if (FileExists(name)) {
      std::ifstream ifs(name);
      boost::archive::text_iarchive ia(ifs);
      ia >> *this;
    }
    else {
      cout << "File: " << name << " does not exist" << endl;
      exit(0);
    }
  } catch(...) { 
    cout << "File: " << name << " could not be loaded." << endl;
    exit(0);
  }
}

uint DBN::get_total_output_size(){
        return total_output_size;
}





ostream& operator<<(ostream& output, const DBN& a)  {
	for (uint i=0; i<a.unique_nodes.size(); i++) {
		cout << "NODE: " << i << endl;
		cout << *(a.unique_nodes[i]) << endl;
	}
	return output;
}

}
