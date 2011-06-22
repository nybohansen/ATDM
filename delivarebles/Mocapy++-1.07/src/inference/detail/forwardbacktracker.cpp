/*
 * forwardbacktracker.cpp
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

#include "forwardbacktracker.h"

#include <math.h>

using namespace std;

namespace mocapy {

ForwardBacktracker::ForwardBacktracker(DBN* new_dbn, uint new_hidden_node, bool check) {
	// Set network and hidden node
	dbn = new_dbn;
	hidden_node = new_hidden_node;

	// Check that the sampler will work on the network
	if (check)
		check_dbn();

	// Store local variables
	nodes_0 = dbn->getNodes0();
	nodes_1 = dbn->getNodes1();
	hd_0 = dynamic_cast<DiscreteNode*> (dbn->getNodes0()[hidden_node]);
	hd_1 = dynamic_cast<DiscreteNode*> (dbn->getNodes1()[hidden_node]);
	hidden_node_size = hd_1->get_node_size();

	// Store CPDs
	cpd_0 = & dynamic_cast<DiscreteDensities*>(hd_0->get_densities())->getCPD();
	cpd_1 = & dynamic_cast<DiscreteDensities*>(hd_1->get_densities())->getCPD();

	cpd_1_swaped = dynamic_cast<DiscreteDensities*>(hd_1->get_densities())->getCPD().moveaxis(0, cpd_1->get_shape().size()-2);
}

void ForwardBacktracker::check_dbn() {
	// Get the node in both slices
	if (!(hidden_node < dbn->getNodes0().size()))
		throw MocapyExceptions("ForwardBacktrack: The hidden node index must be smaller than the number of hidden nodes in slice 0!");
	if (!(hidden_node < dbn->getNodes1().size()))
		throw MocapyExceptions("ForwardBacktrack: The hidden node index must be smaller than the number of hidden nodes in slice 1!");

	Node *node_0 = dbn->getNodes0()[hidden_node];
	Node *node_1 = dbn->getNodes1()[hidden_node];

	// Assert that the hidden node is discrete
	if (!(node_0->get_node_type() == DISCRETE && node_1->get_node_type() == DISCRETE))
		throw MocapyExceptions("ForwardBacktrack: The hidden node must be discrete!");

	// Check that the hidden node only depends on itself in previous slice
	if (node_1->parents_0  != vec(hidden_node) )
		throw MocapyExceptions("ForwardBacktrack: The hidden node can only (but must) depend on itself in previous slide!");

	// Check that the hidden node only has itself as a child in the next slice
	if (node_1->children_2 !=  vec(node_1) || node_0->children_2 != vec(node_1))
		throw MocapyExceptions("ForwardBacktrack: The hidden can only has itself as a child in the next slice!");

	// Check that the hidden node has same size in both slices
	if (node_0->get_node_size() != node_1->get_node_size())
		throw MocapyExceptions("ForwardBacktrack: The hidden node must have the same size in all slices!");

	// Check that the hidden node has the same children in all slices
	if (node_0->children_1 != node_1->children_1)
		throw MocapyExceptions("ForwardBacktrack: The hidden node must have the same children in all slices!");

	// Check that the children of the hidden node:
	// - only depends on the hidden node
	// - doesn't have any children
	for(vector<Node*>::iterator child = node_1->children_1.begin(); child < node_1->children_1.end(); child++ ){
		if ((*child)->parents_0.size() != 0)
			throw MocapyExceptions("ForwardBacktrack: The children of the hidden node cannot have parents in previous slices!");

		if ((*child)->parents_1 != vec(hidden_node))
			throw MocapyExceptions("ForwardBacktrack: The children of the hidden node can only have one parent (the hidden node)!");

		if ((*child)->children_1.size() != 0)
			throw MocapyExceptions("ForwardBacktrack: The children of the hidden node cannot have any children in the same slice!");

		if ((*child)->children_2.size() != 0)
			throw MocapyExceptions("ForwardBacktrack: The children of the hidden node cannot have any children in the next slice!");
	}
}

void ForwardBacktracker::forward_pass(uint start, uint end, uint seq_len, MDArray<eMISMASK> & mismask, bool multiply_by_parents) {
	// Set variables used in the forward algorithm
	uint slice_count = end - start;

	// Setup the forward array and the scales
	scales.clear();
	scales.reserve(slice_count);
	forward.set_shape(vec(slice_count, hidden_node_size));

	transition_matrices.clear();
	transition_matrices.reserve(slice_count);

	// Fill out the forward array
	for(uint i=0, l=start; i < slice_count; i++, l++) {
		DiscreteNode *node;
		MDArray<double> *cpd;
		vector<Node*> *nodes;

		// For the first slice in the network
		if(l==0 || i==0) {
			// Get the appropriate hidden node
			if(l==0) {
				node = hd_0;
				cpd = cpd_0;
			}
			else {
				node = hd_1;
				cpd = cpd_1;
			}

			// Get the values of the parents of the hidden node
			// (and remove the value of the hidden node)
			vector<double> dpv;
			vector<uint> ipv;

			node->parentmap.get(l, dpv);
			dpv.pop_back();
			toint(dpv, ipv);

			// Get the transition probability and store it in the forward array
			MDArray<double> *cpd_sliced = &cpd->get_view(ipv);
			forward.set(i, *cpd_sliced);

			// Store the transition matrix from a posible backtrack
			transition_matrices.push_back(cpd_sliced);
		}
		// For the remaining slices
		else {
			// Get the appropriate hidden node
			node = hd_1;

			// Get the values of the parents of the hidden node
			// (and remove the value of the hidden node)
			vector<double> dpv;
			vector<uint> ipv;

			node->parentmap.get(l, dpv);
			dpv.pop_back();
			dpv.erase(dpv.begin());
			toint(dpv, ipv);

			MDArray<double> *cpd_sliced = &cpd_1_swaped.get_view(ipv);

			// Store the transition matrix from a posible backtrack
			transition_matrices.push_back(cpd_sliced);

			// Multiply transitions probabilities P(HD_l = g | HD_{l-1} = h ) by the
			// forward values forward[i-1, g] and sum over node values of g

			// TODO: This can probably be done faster
			for(uint j=0; j < hidden_node_size; j++ ) {
				// Sum over the previous hidden value
				double sum = 0;
				for(uint k=0; k < hidden_node_size; k++ ) {
					sum += forward.get(i-1, k) * cpd_sliced->get(k, j);
				}
				forward.set(i, j, sum);
			}
		}

		// Multiply the forward value by the probability of the parents
		if(multiply_by_parents) {
			if(l==0)
				nodes = &nodes_0;
			else
				nodes = &nodes_1;

			for(vector<uint>::iterator parent=node->parents_1.begin(); parent < node->parents_1.end(); parent++) {
				double parent_likelihood = exp( (*nodes)[(*parent)]->get_slice_log_likelihood(l) );
				forward.get_view(i).multiply_inplace(parent_likelihood);
			}
		}

		// Multiply the forward values by the probability of the observed children
		// Note that the value of the hidden node is preserved in prev_node_value
		for(vector<Node*>::iterator child=node->children_1.begin(); child < node->children_1.end(); child++) {

			if (mismask.get(l, (*child)->node_index) == MOCAPY_OBSERVED) {
				// Get the original value of the hidden node
			 	vector<double> dpv;
				node->parentmap.get(l, dpv);

				// Multiply the forward vector with the likelihood values of the child
				for(uint j=0; j < hidden_node_size; j++ ) {
					node->parentmap.set(l, j);
					forward.get(i,j) *= exp( (*child)->get_slice_log_likelihood(l) );
				}

				// Restore the original value of the hidden node
				node->parentmap.set(l, dpv.back());
			}
		}

		// If the is the end slice but not the end of the DBN, take the next hidden value into account
		if(l==end-1 && l!=seq_len-1) {
			// The next node can only be hd_1
			vector<double> dpv;
			vector<uint> ipv;

			hd_1->parentmap.get(l+1, dpv);
			uint next_hidden_value = (uint)dpv.back();
			dpv.pop_back();
			dpv.erase(dpv.begin());
			toint(dpv, ipv);

			MDArray<double> * cpd_1_sliced = &cpd_1_swaped.get_view(ipv);

			// Multiply the forward vector with the probability of the transition to the next hidden state
			for(uint j=0; j < hidden_node_size; j++ ) {
				forward.get(i,j) *= cpd_1_sliced->get(j, next_hidden_value);
			}
		}

		// Scale forward to be a probability distribution (this is
		// allowed according to BSA sec. 3.6 (p. 78)
		vector<double> *column_i = (&(&forward.get_view(i))->get_values());
		double scale = 0;

		// Sum column i
		for(uint j=0; j < hidden_node_size; j++ ) {scale += (*column_i)[j];}

		// Normalise column i by it's sum
                bool observedNan=false;
		for(uint j=0; j < hidden_node_size; j++ ) {
                     (*column_i)[j] /= scale;
                     if (isnan((*column_i)[j])) 
                          observedNan = true;
                }

                if (observedNan) {
                     for(uint j=0; j < hidden_node_size; j++ ) {
                          (*column_i)[j] = 1.0/hidden_node_size;
                     }
                }

		// Save the scale
		scales.push_back(scale);
	}
}

  void ForwardBacktracker::backtrack_pass(uint start, uint end, uint seq_len, MDArray<eMISMASK> & mismask, RandomGen* rg) {
	uint choice = 0;
	MDArray<double> transition_matrix(vec(hidden_node_size));

	// The backtrack loop
	for(int l=end-1; l >= ((int) start); l--) {
		uint i = l - start;
		DiscreteNode *node;

		// Setup the node, cpd and nodes
		if(i==0) {
			node = hd_0;
		}
		else {
			node = hd_1;
		}

		// Choose whether to sample from the forward probability
		// distribution or take the previous choice into account
		if ( ((uint) l)==end-1 ) {
			transition_matrix = forward.get_view(i);
		}
		else {
			// TODO: Make this sum faster
			for (uint j=0; j<hidden_node_size; j++)
				transition_matrix[j] = forward.get(i,j) * transition_matrices[i+1]->get(j, choice);
			transition_matrix.normalize();
		}

		// Sample the hidden node from the forward probability distribution
		transition_matrix.cumsum();
		assert(rg);
		double r = rg->get_rand();
		choice = transition_matrix.bisect(r);
		node->parentmap.set(l, choice);

		// Sample the children
		for(vector<Node*>::iterator child=node->children_1.begin(); child < node->children_1.end(); child++) {
			if (mismask.get(l, (*child)->node_index) != MOCAPY_OBSERVED) {
				(*child)->sample(l);
			}
		}
	}
}

double ForwardBacktracker::calc_ll() {
	double sum = 0;
	for(vector<double>::iterator value = scales.begin(); value < scales.end(); value++) {
		sum += log(*value);
	}
	return sum;
}

}
