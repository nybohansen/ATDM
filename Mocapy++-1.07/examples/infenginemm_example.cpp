/*
 * infenginemm_example.cpp
 *
 *  Copyright (C) 2008, Jes Frellsen, The Bioinformatics Centre, University of Copenhagen.
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

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include "mocapy.h"
#include <time.h>

using namespace mocapy;
using namespace std;

template<typename T>
vector<T> vect(T t1, T t2, T t3) {
	vector<T> v;
	v.push_back(t1);
	v.push_back(t2);
	v.push_back(t3);
	return v;
}

void copy_flat(const vector<double*> & flat_iterator, double array[]) {
	for(uint i=0; i<flat_iterator.size(); i++) {
		*flat_iterator[i] = array[i];
	}
}


int main(void) {
	time_t t = time(NULL);
	mocapy_seed(t);
	cout << "Mocapy seed = " << t << endl;

	cout << "******************" << endl;
	cout << "Setting up network" << endl;
	cout << "******************" << endl;

	// Set cpd for in-node
	double in_cpd_array[] = {0.5, 0.5};
	CPD in_cpd;
	in_cpd.set_shape(2);
	copy_flat(in_cpd.flat_iterator(), in_cpd_array);


	// Set cpd for hidde node
	double hd_cpd_array[] = {1.0, 0.0, 0.0, 0.0,
			              0.25, 0.25, 0.25, 0.25};
	CPD hd_cpd;
	hd_cpd.set_shape(2, 4);
	copy_flat(hd_cpd.flat_iterator(), hd_cpd_array);

	// Set cpd for out node
	double out_cpd_array[] = {0.7, 0.0, 0.0, 0.0,
				              0.0, 0.7, 0.0, 0.0,
		                      0.0, 0.0, 0.7, 0.0,
		                      0.0, 0.0, 0.0, 0.7};


	CPD out_cpd;
	out_cpd.set_shape(4, 4);
	copy_flat(out_cpd.flat_iterator(), out_cpd_array);

	// Setup nodes
	DiscreteNode* in = new DiscreteNode();
	in->set_densities(DiscreteDensities(2, in_cpd));

	DiscreteNode* hd = new DiscreteNode();
	hd->set_densities(DiscreteDensities(4, hd_cpd));

	DiscreteNode* out = new DiscreteNode();
	out->set_densities(DiscreteDensities(4, out_cpd));

	// Setup dbn
	vector<Node*> nodes = vect(in->base(), hd->base(), out->base());

	DBN dbn = DBN(nodes, nodes);

	dbn.add_intra(0, 1);
	dbn.add_intra(1, 2);
	dbn.construct();

	// Output the CPDs
	cout << "in:  \n" << in->get_densities()->getCPD() << endl;
	cout << "hd:  \n" << hd->get_densities()->getCPD() << endl;
	cout << "out: \n" << out->get_densities()->getCPD() << endl;
	cout << endl;

	cout << "******************" << endl;
	cout << "Setting up dataset" << endl;
	cout << "******************" << endl;

	Sequence data;
	data.set_shape(5, 3);
	double seq_array[] = {1,0,1,
	                      1,0,1,
	                      1,0,1,
	                      1,0,1,
	                      1,0,1};
	copy_flat(data.flat_iterator(), seq_array);


	MDArray<eMISMASK> mism;
	mism.set_shape(5, 3);
	mism.set_wildcard(vec(-1, 0), MOCAPY_OBSERVED);
	mism.set_wildcard(vec(-1, 1), MOCAPY_HIDDEN);
	mism.set_wildcard(vec(-1, 2), MOCAPY_OBSERVED);

	MDArray<eMISMASK> mism_sample;
	mism_sample.set_shape(5, 3);
	mism_sample.set_wildcard(vec(-1, 0), MOCAPY_OBSERVED);
	mism_sample.set_wildcard(vec(-1, 1), MOCAPY_HIDDEN);
	mism_sample.set_wildcard(vec(-1, 2), MOCAPY_HIDDEN);

	cout << "Data: \n" << data << endl;
	cout << "Mism: \n" << mism << endl;
	cout << "Mism_sample: \n" << mism_sample << endl;

	cout << "*************************" << endl;
	cout << "Example of SampleInfEngineMM" << endl;
	cout << "*************************" << endl;

	// Setup the sampler
	SampleInfEngineMM sampler(&dbn, data, mism_sample, 1);

	MDArray<double> sample = sampler.sample_next();
	cout << "Sample:\n" << sample;
	cout << "ll of sample = " << sampler.calc_ll(mism) << endl << endl;

	cout << "undo()" << endl;
	sampler.undo();
	cout << "ll of initial values =" << sampler.calc_ll(mism) << endl << endl;

	cout << "Setting start=0 and end=1" << endl << endl;
	sampler.set_start_end(0, 1);

	sample = sampler.sample_next();
	cout << "Sample:\n" << sample;
	cout << "ll of sample = " << sampler.calc_ll(mism) << endl << endl;

	cout << "undo()" << endl;
	sampler.undo();


	// Use the other mismask, so that we can just sample the hidden states
	sampler.set_seq_mismask(data, mism);
	sampler.set_start_end();

	sample = sampler.sample_next();
	cout << "Sample only hidden states:\n" << sample;
	cout << "ll of sample = " << sampler.calc_ll(mism) << endl << endl;


	cout << "********************************" << endl;
	cout << "Example of LikelihoodInfEngineMM" << endl;
	cout << "********************************" << endl;

	double ll;
	LikelihoodInfEngineMM infengine(&dbn, 1);

	ll = infengine.calc_ll(data, mism, false);

	cout << "ll of initail values = " << ll << endl;

	ll = infengine.calc_ll(data, mism, 0, -1, true);

	cout << "ll of initail values (including parents) = " << ll << endl;

	// Clean up!
	delete in;
	delete hd;
	delete out;

	return EXIT_SUCCESS;
}
