/*
 * discretedensities.h
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

#include "discretedensities.h"

using namespace std;

namespace mocapy {

DiscreteDensities::DiscreteDensities() {
	initialize();
}

DiscreteDensities::DiscreteDensities(uint new_node_size, Prior * new_prior, bool new_init_random) {
	// node_size is the number of states of the discrete node
	node_size = new_node_size;
	prior = new_prior;
	init_random = new_init_random;
	initialize();
}

DiscreteDensities::DiscreteDensities(uint new_node_size, CPD & new_user_cpd, Prior * new_prior) {
	// node_size is the number of states of the discrete node
	node_size = new_node_size;
	prior = new_prior;
	user_cpd = new_user_cpd;
	init_random = false;
	initialize();
}


void DiscreteDensities::initialize() {
	type = DISCRETE;
	output_size = 1;
}

// Normalize CPD and make sure that CPD is 'well'
void DiscreteDensities::set_cpd(CPD & new_cpd) {
	assert(new_cpd.size()> 0);
	cpd = new_cpd;
	cpd.add_inplace(_MIN_TRANSITION);
	cpd.normalize();
	cpd.clip(_MIN_TRANSITION, 1000);

	log_cpd = cpd;
	log_cpd.log_all();

	cum_cpd = cpd;
	cum_cpd.cumsum();
}

CPD DiscreteDensities::make_random_cpd(vector<uint> & shape, bool no_zeroes) {
	// Return a random CPD.
	CPD c(shape);
	c.randomize(randomGen);

	// avoid zero entries
	if (no_zeroes) {
		c.clip(0.01, 2);
	}
	c.normalize();
	return c;
}

CPD DiscreteDensities::make_uniform_cpd(const vector<uint> & shape) {
	// Return a uniform CPD.
	CPD c(shape);

	// avoid zero entries
	c.clip(0.1, 2);
	c.normalize();
	return c;
}


// Called in node.construct, and initializes the density arrays
void DiscreteDensities::construct(vector<uint> & parent_sizes) {
	CPD_shape = vec_conc(parent_sizes, node_size);

	if(user_cpd.empty()) {
		CPD cpd;
		if(!init_random) {
			cpd = make_uniform_cpd(CPD_shape);
		}
		else {
			cpd = make_random_cpd(CPD_shape, true);
		}
		set_cpd(cpd);
	}
	else {
		assert(user_cpd.get_shape() == CPD_shape);
		set_cpd(user_cpd);
	}
}

// Parameter estimation based on the ESS
void DiscreteDensities::estimate(vector<MDArray<double> > & ess) {
	assert(!ess.empty());
	if (prior) {
		prior->apply_prior(ess.front());
	}

	set_cpd(ess.front());
}

// Return a sample, based on indicated parent values
vector<double> DiscreteDensities::sample(vector<double> & pv) {
	MDArray<double>* cumulative;
	if (pv.empty())
		cumulative = &cum_cpd;
	else {
		vector<uint> ipv;
		toint(pv, ipv);
		cumulative = &(cum_cpd.get_view(ipv));
	}

	double r = randomGen->get_rand();

	assert(r<=1 && r>=0);
	vector<double> choice;
	choice.push_back( cumulative->bisect(r) );
	return choice;
}

// Return likelihood, that is: P(child|parents)
double DiscreteDensities::get_lik(vector<double> & ptv, bool log) {
	// The following size==2 special case is not strictly necessary.
	// However, often a node has only one parent, and avoiding the
	// toint casting function gives a huge speedup.
	if (ptv.size() == 2) {
		if (log) {
			return log_cpd.get((uint)ptv[0], (uint)ptv[1]);
		}
		else {
			return cpd.get((uint)ptv[0], (uint)ptv[1]);
		}
	}
	else {
		vector<uint> ipv;
		toint(ptv, ipv);
		if (log)
			return log_cpd[ipv];
		else
			return cpd[ipv];
	}
}

// Return the distribtion's parameters
vector<MDArray<double> > DiscreteDensities::get_parameters() {
   vector<MDArray<double> > ret;
   ret.push_back(cpd);
   return ret;
}


void DiscreteDensities::set_user_cpd(CPD & new_user_cpd) {
	user_cpd = new_user_cpd;
}

void DiscreteDensities::set_prior(Prior * new_prior) {
	prior = new_prior;
}


ostream& operator<<(ostream& output, const DiscreteDensities& a)  {
	output << "Node: Discrete, size: ";

	if (a.cpd.empty()) {
		output << "Empty" << endl;
	}
	else {
		vector<uint> shape = a.cpd.get_shape();
		uint len = shape.size();

		for(uint i=0; i<len; i++)
			output << shape[i] << " ";

        output << endl;

	// a.cpd.tostring() does not work properly when the shape size is more than 2
	//        output << a.cpd.tostring();

	output << a.cpd << endl;
	}

	return output;
}


}
