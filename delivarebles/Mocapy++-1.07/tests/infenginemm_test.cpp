	/*
 * infenginemm_test.cpp
 *
 *  Created on: Mar 24, 2009
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

	// Set cpd for hidde node
	double hd_cpd_array[] = {0.5, 0.5};
	CPD hd_cpd;
	hd_cpd.set_shape(2);
	copy_flat(hd_cpd.flat_iterator(), hd_cpd_array);

	// Set cpd for out1 node
	double out1_cpd_array[] = {0.8, 0.2,
				               0.2, 0.8};

	CPD out1_cpd;
	out1_cpd.set_shape(2, 2);
	copy_flat(out1_cpd.flat_iterator(), out1_cpd_array);

	// Set cpd for out2 node
	double out2_cpd_array[] = {0.8, 0.2,
				               0.2, 0.8};

	CPD out2_cpd;
	out2_cpd.set_shape(2, 2);
	copy_flat(out2_cpd.flat_iterator(), out2_cpd_array);

	// Setup nodes
	DiscreteNode* hd = new DiscreteNode();
	hd->set_densities(DiscreteDensities(2, hd_cpd));

	DiscreteNode* out1 = new DiscreteNode();
	out1->set_densities(DiscreteDensities(2, out1_cpd));

	DiscreteNode* out2 = new DiscreteNode();
	out2->set_densities(DiscreteDensities(2, out2_cpd));

	// Setup dbn
	vector<Node*> nodes = vect(hd->base(), out1->base(), out2->base());

	DBN dbn = DBN(nodes, nodes);

	dbn.add_intra(0, 1);
	dbn.add_intra(0, 2);
	dbn.construct();

	// Output the CPDs
	cout << "hd:  \n" << hd->get_densities()->getCPD() << endl;
	cout << "out1: \n" << out1->get_densities()->getCPD() << endl;
	cout << "out2: \n" << out2->get_densities()->getCPD() << endl;
	cout << endl;

	cout << "******************" << endl;
	cout << "Setting up dataset" << endl;
	cout << "******************" << endl;

	Sequence data;
	data.set_shape(1, 3);
	double seq_array[] = {1,0,0};
	copy_flat(data.flat_iterator(), seq_array);


	MDArray<eMISMASK> mism;
	mism.set_shape(1, 3);
	mism.set_wildcard(vec(-1, 0), MOCAPY_HIDDEN);
	mism.set_wildcard(vec(-1, 1), MOCAPY_OBSERVED);
	mism.set_wildcard(vec(-1, 2), MOCAPY_OBSERVED);

	MDArray<eMISMASK> mism_marginal;
	mism_marginal.set_shape(1, 3);
	mism_marginal.set_wildcard(vec(-1, 0), MOCAPY_HIDDEN);
	mism_marginal.set_wildcard(vec(-1, 1), MOCAPY_OBSERVED);
	mism_marginal.set_wildcard(vec(-1, 2), MOCAPY_HIDDEN);

	MDArray<eMISMASK> mism_sample;
	mism_sample.set_shape(1, 3);
	mism_sample.set_wildcard(vec(-1, 0), MOCAPY_HIDDEN);
	mism_sample.set_wildcard(vec(-1, 1), MOCAPY_HIDDEN);
	mism_sample.set_wildcard(vec(-1, 2), MOCAPY_OBSERVED);


	cout << "Data: \n" << data << endl;
	cout << "Mism: \n" << mism << endl;
	cout << "Mism_cond: \n" << mism_marginal << endl;
	cout << "Mism_sample: \n" << mism_sample << endl;

	cout << "*********************************************************" << endl;
	cout << "Testing sampler and likelihood inf engine for consistency" << endl;
	cout << "*********************************************************" << endl;

	//
	// Calculate likelihoods by using the LikelihoodInfEngineMM
	//

	printf("Calculating likelihood by summing over the hidden node using the LikelihoodInfEngineMM ...\n");

	LikelihoodInfEngineMM infengine(&dbn, 0);

	double ll;
	ll = infengine.calc_ll(data, mism);
	printf("P(out1, out2) = %.8f\n", exp(ll));

	double ll_marginal;
	ll_marginal = infengine.calc_ll(data, mism_marginal);
	printf("P(out1) = %.8f\n", exp(ll_marginal));

	double ll_cond = ll - ll_marginal;
	printf("P(out1|out2) = %.8f\n\n", exp(ll_cond));

	//
	// Calculate likelihoods by sampling directly from P(out1|out2) using SampleInfEngineMM
	//

	cout << "Estimating P(out1|out2) by sampling directly from P(out1|out2) using SampleInfEngineMM" << endl;

	SampleInfEngineMM sampler(&dbn, data, mism_sample, 0);

	vector<int> counter2(2,0);

	for(int i=0; i<100000; i++) {
		MDArray<double> sample = sampler.sample_next();
		counter2[TOINT(sample.get(0,1))]++;
	}

	double estimate = 1.0 * counter2[data.get(0,1)] / (counter2[0] + counter2[1]);
	printf("Estimate of P(out1|out2) = %.8f\n\n", estimate);

	//
	// Calculate likelihoods by sampling from P(out1,out2) using the standard sample
	//

	cout << "Estimating P(out1|out2) by sampling from P(out1,out2) using the standard DBN sampler" << endl;

	vector<int> counter4(4,0);

	for(int i=0; i<200000; i++) {
		MDArray<double> sample = dbn.sample_sequence(1).first;
		counter4[TOINT(sample.get(0,1) + 2*TOINT(sample.get(0,2)))]++;
	}

	estimate = 1.0 * counter4[data.get(0,1)] / (counter4[0] + counter4[1]);
	printf("P(out1|out2) = %f\n", estimate);


	// Clean up!
	delete hd;
	delete out1;
	delete out2;

	return EXIT_SUCCESS;
}
