/*
 * test_vonmises.cpp
 *
 *  Created on: Nov 24, 2008
 *      Author: frellsen
 */

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include "mocapy.h"
#include <time.h>

using namespace mocapy;
using namespace std;

template<typename T>
vector<T> vect(T t1, T t2) {
	vector<T> v;
	v.push_back(t1);
	v.push_back(t2);
	return v;
}



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

	cout << "**************" << endl;
	cout << "Setting up DBN" << endl;
	cout << "**************" << endl;

	// Set cpd for hd0
	double hd0_cpd_array[] = {0.25, 0.25, 0.25, 0.25};
	CPD hd0_cpd;
	hd0_cpd.set_shape(4);
	copy_flat(hd0_cpd.flat_iterator(), hd0_cpd_array);

	// Set cpd for hd1
	double hd1_cpd_array[] = {0.90, 0.10, 0.00, 0.00,
		                      0.00, 0.80, 0.10, 0.10,
		                      0.10, 0.10, 0.70, 0.10,
		                      0.20, 0.10, 0.10, 0.60};

	CPD hd1_cpd;
	hd1_cpd.set_shape(4, 4);
	copy_flat(hd1_cpd.flat_iterator(), hd1_cpd_array);

	// Set cpd for out node
	double out_cpd_array[] = {0.7, 0.1, 0.1, 0.1,
				              0.1, 0.7, 0.1, 0.1,
		                      0.1, 0.1, 0.7, 0.1,
		                      0.1, 0.1, 0.1, 0.7};


	CPD out_cpd;
	out_cpd.set_shape(4, 4);
	copy_flat(out_cpd.flat_iterator(), out_cpd_array);

	// Setup nodes
	DiscreteNode* hd0 = new DiscreteNode();
	hd0->set_densities(DiscreteDensities(4, hd0_cpd));

	DiscreteNode* hd1 = new DiscreteNode();
	hd1->set_densities(DiscreteDensities(4, hd1_cpd));

	DiscreteNode* out = new DiscreteNode();
	out->set_densities(DiscreteDensities(4, out_cpd));

	// Setup dbn
	vector<Node*> start_nodes = vect(hd0->base(), out->base());
	vector<Node*> end_nodes = vect(hd1->base(), out->base());

	DBN dbn = DBN(start_nodes, end_nodes);

	dbn.add_intra(0, 1);
	dbn.add_inter(0, 0);
	dbn.construct();

	// Output the CPDs
	cout << "hd0: \n" << hd0->get_densities()->getCPD() << endl;
	cout << "hd1: \n" << hd1->get_densities()->getCPD() << endl;
	cout << "out:  \n" << out->get_densities()->getCPD() << endl;
	cout << endl;

	cout << "******************" << endl;
	cout << "Setting up dataset" << endl;
	cout << "******************" << endl;

	Sequence data;
	data.set_shape(3, 2);
	double seq_array[] = {0,0,
	                      0,0,
	                      0,0};
	copy_flat(data.flat_iterator(), seq_array);


	MDArray<eMISMASK> mism;
	mism.set_shape(3, 2);
	mism.set_wildcard(vec(-1, 0), MOCAPY_HIDDEN);
	mism.set_wildcard(vec(-1, 1), MOCAPY_OBSERVED);

	MDArray<eMISMASK> mism_sample;
	mism_sample.set_shape(3, 2);
	mism_sample.set_wildcard(vec(-1, 0), MOCAPY_HIDDEN);
	mism_sample.set_wildcard(vec(-1, 1), MOCAPY_HIDDEN);

	cout << "Data: \n" << data << endl;
	cout << "Mism: \n" << mism << endl;
	cout << "Mism_sample: \n" << mism_sample << endl;


	cout << "*********************************************************" << endl;
	cout << "Testing sampler and likelihood inf engine for consistency" << endl;
	cout << "*********************************************************" << endl;

	//
	// Calculate likelihood of data
	//
	printf("Calculating likelihood by the forward algorithm ...\n");

	double ll;
	LikelihoodInfEngineHMM infengine(&dbn, 0);
	ll = infengine.calc_ll(data, mism);

	printf("likelihood of data = %.8f\n\n", exp(ll));

	//
	// Sample likelihood by forward algorithm
	//
	printf("Sampling likelihood by the forward-backtrack algorithm ...\n");

	// Setup the sampler
	SampleInfEngineHMM sampler(&dbn, data, mism_sample, 0);

	const int total = 1000000;
	int count = 0;

	for(int i=0; i<total; i++) {
		MDArray<double> sample = sampler.sample_next();
		if(TOINT(sample.get(0,1))==TOINT(data.get(0,1)) &&
		   TOINT(sample.get(1,1))==TOINT(data.get(1,1)) &&
		   TOINT(sample.get(2,1))==TOINT(data.get(2,1))) {
			count++;
		}
	}

	printf("estimated likelihood = %.8f\n\n", 1.0*count/total);

	//
	// Sample likelihood by sampling the network from one end to the other
	//
	cout << "Sampling likelihood by sampling nodes from one end to the other ..." << endl;

	count = 0;

	for(int i=0; i<total; i++) {
		MDArray<double> sample = dbn.sample_sequence(3).first;

		if(TOINT(sample.get(0,1))==TOINT(data.get(0,1)) &&
		   TOINT(sample.get(1,1))==TOINT(data.get(1,1)) &&
		   TOINT(sample.get(2,1))==TOINT(data.get(2,1))) {
			count++;
		}
	}

	printf("estimated likelihood = %.8f\n", 1.0*count/total);

	// Clean up!
	delete hd0;
	delete hd1;
	delete out;

	return EXIT_SUCCESS;
 }
