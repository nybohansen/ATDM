/*
 * generalinfenginemm.cpp
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

#include "generalinfenginemm.h"

#include <math.h>

using namespace std;

namespace mocapy {

GeneralInfEngineMM::GeneralInfEngineMM(DBN* new_dbn, uint new_hidden_node, bool check) {
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
}

void GeneralInfEngineMM::check_dbn() {
	// Get the node in both slices
	if (!(hidden_node < dbn->getNodes0().size()))
		throw MocapyExceptions("GeneralInfEngineMM: The hidden node index must be smaller than the number of hidden nodes in slice 0!");
	if (!(hidden_node < dbn->getNodes1().size()))
		throw MocapyExceptions("GeneralInfEngineMM: The hidden node index must be smaller than the number of hidden nodes in slice 1!");

	Node *node_0 = dbn->getNodes0()[hidden_node];
	Node *node_1 = dbn->getNodes1()[hidden_node];

	// Assert that the hidden node is discrete
	if (!(node_0->get_node_type() == DISCRETE && node_1->get_node_type() == DISCRETE))
		throw MocapyExceptions("GeneralInfEngineMM: The hidden node must be discrete!");

	// Check that the hidden node doesn't have any parents in previous slice
	if (node_1->parents_0.size() != 0 )
		throw MocapyExceptions("GeneralInfEngineMM: The hidden node cannot have any parents in previous slide!");

	// Check that the hidden node doesn't have any children in the next slice
	if (node_1->children_2.size() !=  0 || node_0->children_2.size() != 0)
		throw MocapyExceptions("GeneralInfEngineMM: The hidden cannot have any children in the next slice!");

	// Check that the hidden node has same size in both slices
	if (node_0->get_node_size() != node_1->get_node_size())
		throw MocapyExceptions("GeneralInfEngineMM: The hidden node must have the same size in all slices!");

	// Check that the hidden node has the same children in all slices
	if (node_0->children_1 != node_1->children_1)
		throw MocapyExceptions("GeneralInfEngineMM: The hidden node must have the same children in all slices!");

	// Check that the children of the hidden node:
	// - only depends on the hidden node
	// - doesn't have any children
	for(vector<Node*>::iterator child = node_1->children_1.begin(); child < node_1->children_1.end(); child++ ){
		if ((*child)->parents_0.size() != 0)
			throw MocapyExceptions("GeneralInfEngineMM: The children of the hidden node cannot have parents in previous slices!");

		if ((*child)->parents_1 != vec(hidden_node))
			throw MocapyExceptions("GeneralInfEngineMM: The children of the hidden node can only have one parent (the hidden node)!");

		if ((*child)->children_1.size() != 0)
			throw MocapyExceptions("GeneralInfEngineMM: The children of the hidden node cannot have any children in the same slice!");

		if ((*child)->children_2.size() != 0)
			throw MocapyExceptions("GeneralInfEngineMM: The children of the hidden node cannot have any children in the next slice!");
	}
}

double GeneralInfEngineMM::calc_ll(uint start, uint end, uint seq_len, MDArray<eMISMASK> & mismask, bool multiply_by_parents) {
	// Set variables used in the algorithm
	uint slice_count = end - start;
	double ll = 0;
	vector<Node*> *nodes;

	// Calculate the ll for each slice and sum it
	for(uint i=0, l=start; i < slice_count; i++, l++) {
		DiscreteNode *node;
		MDArray<double> *cpd;

		// Get the appropriate hidden node
		if(l==0) {
			node = hd_0;
			cpd = cpd_0;
		}
		else {
			node = hd_1;
			cpd = cpd_1;
		}

		// Save the original value of the hidden node
	 	vector<double> dpv;
		node->parentmap.get(l, dpv);

		// Sum out the hidden node by using the theorem of total probability:
		// ln(P(c1, c2, ...) = ln (sum_h exp[lnP(h) * lnP(c1) * lnP(c2) * ...] )
		//
		// slice_likelihood: P(c1, c2, ...)
		// h_likelihood: exp[lnP(h) * lnP(c1) * lnP(c2) * ...]

		double slice_likelihood = 0;

		for(uint j=0; j < hidden_node_size; j++ ) {
			// Calculate the likelihood for H=j
			double h_likelihood = 0;

			// Set H=j
			node->parentmap.set(l, j);

			// Add the likelihood of the hidden node
			h_likelihood += node->get_slice_log_likelihood(l);

			// Add the likelihood of the observed children
			for(vector<Node*>::iterator child=node->children_1.begin(); child < node->children_1.end(); child++) {
				if (mismask.get(l, (*child)->node_index) == MOCAPY_OBSERVED) {
					h_likelihood += (*child)->get_slice_log_likelihood(l);
				}
			}
			slice_likelihood += exp(h_likelihood);
		}

		// Transforming into log-space
		double slice_ll = log(slice_likelihood);

		// Restore the original value of the hidden node
		node->parentmap.set(l, dpv.back());

		// Add the ll of the parents to the slice ll
		if(multiply_by_parents) {
			if(l==0)
				nodes = &nodes_0;
			else
				nodes = &nodes_1;

			for(vector<uint>::iterator parent=node->parents_1.begin(); parent < node->parents_1.end(); parent++) {
				slice_ll += (*nodes)[(*parent)]->get_slice_log_likelihood(l);
			}
		}

		// Add the slice ll to the total ll
		ll += slice_ll;
	}

	return ll;
}


  void GeneralInfEngineMM::sample(uint start, uint end, uint seq_len, MDArray<eMISMASK> & mismask, RandomGen* rg) {

    assert(rg != NULL);
	// Set variables used in the algorithm
	uint slice_count = end - start;

	// Sample each slice independently
	for(uint i=0, l=start; i < slice_count; i++, l++) {
		DiscreteNode *node;
		MDArray<double> *cpd;

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
		// and slice the cpd of the hidden node
		vector<double> dpv;
		vector<uint> ipv;

		node->parentmap.get(l, dpv);
		dpv.pop_back();
		toint(dpv, ipv);

		MDArray<double> cpd_sliced = cpd->get_view(ipv);

		// Multiply the probability distribution of the hidden node
		// by the probability of the observed children
		// Node that we doen't need to store the value of the hidden node (we are going to sample it afterwards)
		for(vector<Node*>::iterator child=node->children_1.begin(); child < node->children_1.end(); child++) {
		  if (mismask.get(l, (*child)->node_index) == MOCAPY_OBSERVED) {
		    // Multiply the probability distribution over the hidden
		    // node by the probability of the child value
		    for(uint j=0; j < hidden_node_size; j++ ) {
		      node->parentmap.set(l, j);
		      cpd_sliced[j] *= exp((*child)->get_slice_log_likelihood(l));
		    }
		  }
		}

		// Scale the distribution over the hidden node to sum to one
		double scale = 0;
		for(uint j=0; j < hidden_node_size; j++ ) {scale += cpd_sliced[j];}
		for(uint j=0; j < hidden_node_size; j++ ) {cpd_sliced[j] /= scale;}


		// Sample the hidden node from the calculated distribution of the hidden node
		cpd_sliced.cumsum();
		double r = rg->get_rand();
		uint choice = cpd_sliced.bisect(r);
		node->parentmap.set(l, choice);

		// Sample the unobserved children
		for(vector<Node*>::iterator child=node->children_1.begin(); child < node->children_1.end(); child++) {
			if (mismask.get(l, (*child)->node_index) != MOCAPY_OBSERVED) {
			  (*child)->setRandomGen(rg);
			  (*child)->sample(l);
			}
		}
	}
}

}
