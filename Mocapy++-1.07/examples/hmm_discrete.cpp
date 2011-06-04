/*
 * discrete_hmm.cpp
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
 *
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
using namespace mocapy;
using namespace std;

int main(void) {
   	mocapy_seed((uint)5556574);

	// Number of trainining sequences
	int N = 100;

	// Sequence lengths
	int T = 100;

	// Gibbs sampling parameters
	int MCMC_BURN_IN = 10;

	// HMM hidden and observed node sizes
	uint H_SIZE=2;
	uint O_SIZE=2;
	bool init_random=false;

	CPD th0_cpd;
	th0_cpd.set_shape(2); th0_cpd.set_values(vec(0.1, 0.9));

	CPD th1_cpd;
	th1_cpd.set_shape(2,2); th1_cpd.set_values( vec(0.95, 0.05, 0.1, 0.9));

	CPD to0_cpd;
	to0_cpd.set_shape(2,2); to0_cpd.set_values( vec(0.1, 0.9, 0.8, 0.2));

	// The target DBN (This DBN generates the data)
	Node* th0 = NodeFactory::new_discrete_node(H_SIZE, "th0", init_random, th0_cpd);
	Node* th1 = NodeFactory::new_discrete_node(H_SIZE, "th1", init_random, th1_cpd);
	Node* to0 = NodeFactory::new_discrete_node(O_SIZE, "to0", init_random, to0_cpd);

	DBN tdbn;
	tdbn.set_slices(vec(th0, to0), vec(th1, to0));

	tdbn.add_intra("th0", "to0");
	tdbn.add_inter("th0", "th1");
	tdbn.construct();

	// The model DBN (this DBN will be trained)
	// For mh0, get the CPD from th0 and fix parameters
	Node* mh0 = NodeFactory::new_discrete_node(H_SIZE, "mh0", init_random, CPD(), th0, true );
	Node* mh1 = NodeFactory::new_discrete_node(H_SIZE, "mh1", init_random);
	Node* mo0 = NodeFactory::new_discrete_node(O_SIZE, "mo0", init_random);

	DBN mdbn;
	mdbn.set_slices(vec(mh0, mo0), vec(mh1, mo0));

	mdbn.add_intra("mh0", "mo0");
	mdbn.add_inter("mh0", "mh1");
	mdbn.construct();

	cout << "*** TARGET ***" << endl;
	cout << *th0 << endl;
	cout << *th1 << endl;
	cout << *to0 << endl;

	cout << "*** MODEL ***" << endl;
	cout << *mh0 << endl;
	cout << *mh1 << endl;
	cout << *mo0 << endl;

	vector<Sequence> seq_list;
	vector< MDArray<eMISMASK> > mismask_list;

	cout << "Generating data" << endl;

	MDArray<eMISMASK> mismask;
	mismask.repeat(T, vec(MOCAPY_HIDDEN, MOCAPY_OBSERVED));

	// Generate the data
	double sum_LL(0);
	for (int i=0; i<N; i++) {
		pair<Sequence, double>  seq_ll = tdbn.sample_sequence(T);
		sum_LL += seq_ll.second;
 		seq_list.push_back(seq_ll.first);
 		mismask_list.push_back(mismask);
	}
	cout << "Average LL: " << sum_LL/N << endl;

	GibbsRandom mcmc = GibbsRandom(&mdbn);
	EMEngine em = EMEngine(&mdbn, &mcmc, &seq_list, &mismask_list);

	InfEngineMCMC inf = InfEngineMCMC(&mdbn, &mcmc, &(seq_list[0]), mismask_list[0]);
	inf.initialize_viterbi_generator(5, 5, true);

	for (uint i=0; i<3; i++) {
		Sample s = inf.viterbi_next();
		cout << "Viterbi LL= " << s.ll << endl;
	}
	cout << "Ending Viterbi" << endl;

	cout << "Starting EM loop" << endl;
	double bestLL=-1000;
	uint it_no_improvement(0);
	uint i(0);
	// Start EM loop
	while (it_no_improvement<100) {
		em.do_E_step(1, MCMC_BURN_IN, true);

		double ll = em.get_loglik();

		cout << "LL= " << ll;

		if (ll > bestLL) {
			cout << " * saving model *" << endl;
			mdbn.save("discrete_hmm.dbn");
			bestLL = ll;
			it_no_improvement=0;
		}
		else { it_no_improvement++; cout << endl; }

		i++;
		em.do_M_step();
	}

	cout << "DONE" << endl;

	mdbn.load("discrete_hmm.dbn");

	cout << "*** TARGET ***" << endl;
	cout << *th0 << endl;
	cout << *th1 << endl;
	cout << *to0 << endl;

	cout << "*** MODEL ***" << endl;
	cout << *mh0 << endl;
	cout << *mh1 << endl;
	cout << *mo0 << endl;


	delete th0;
	delete th1;
	delete to0;

	delete mh0;
	delete mh1;
	delete mo0;
	return EXIT_SUCCESS;
}
