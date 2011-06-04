/*
 * test example
 *
 *  Copyright (C) 2010, Simon H. Rasmussen, The Bioinformatics Centre, University of Copenhagen.
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

int main(void) {
   	mocapy_seed((uint)924573);

	// Number of trainining sequences
	int N = 2000;

	// Allocate storage in each node.
	int max_num_dat_points = N;

	// Sequence lengths
	int T = 1;

	// Gibbs sampling parameters
	int MCMC_BURN_IN = 10;

	// HMM hidden and observed node sizes
	uint H_SIZE=2;
	uint O_SIZE=20;
	bool init_random=true;

	// The target DBN (This DBN generates the data)
	Node* th0 = NodeFactory::new_discrete_node(H_SIZE, "th0", init_random);
	Node* th1 = NodeFactory::new_discrete_node(H_SIZE, "th1", init_random);

	Node* to0 = NodeFactory::new_GD_node("to0", O_SIZE , max_num_dat_points);
	cout<< "Done making GD node" << endl;

	DBN tdbn;
	tdbn.set_slices(vec(th0, to0), vec(th1, to0));

	tdbn.add_intra("th0", "to0");
	tdbn.add_inter("th0", "th1");


	cout<< "model has been defined" << endl;
	tdbn.construct();
	
	// The model DBN (this DBN will be trained)
	// For mh0, get the CPD from th0 and fix parameters
	Node* mh0 = NodeFactory::new_discrete_node(H_SIZE, "mh0", init_random, CPD(), th0, true );
	Node* mh1 = NodeFactory::new_discrete_node(H_SIZE, "mh1", init_random);
	Node* mo0 = NodeFactory::new_GD_node("mo0",O_SIZE, max_num_dat_points);

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

	// Setting parameters
	//CPD pars2;
	CPD pars;
	int t = 4;
	pars.set_shape(H_SIZE,2*(O_SIZE-1));
	//pars2.set_shape(H_SIZE,2*(O_SIZE-1));
	vector<double> parvec;
	for(uint i = 0;i < O_SIZE-1;i++){
	  parvec.push_back(t* 1.0);
	}
	for(double i = O_SIZE-1;i >= 1;i--){
	  parvec.push_back(t*i);
	}
	//	pars.set_values(parvec);
	double s = 2;
	for(uint i = 0;i < O_SIZE-1;i++){
	  parvec.push_back(s);
	}
	for(double i = s*(O_SIZE-1);i >= s*1;i=i-s){
	  parvec.push_back(i);
	}
	pars.set_values(parvec);

	((GDNode*)to0)->get_densities()->set_parameters(pars);
	pair<Sequence, double>  seq_ll;	
	// Generate the data
	double sum_LL(0);
	for (int i=0; i<N; i++) {
	  seq_ll = tdbn.sample_sequence(T);
	  sum_LL += seq_ll.second;
	  seq_list.push_back(seq_ll.first);
	  mismask_list.push_back(mismask);
	}
	cout << "Average LL: " << sum_LL/N << endl;

	GibbsRandom mcmc = GibbsRandom(&mdbn);
	EMEngine em = EMEngine(&mdbn, &mcmc, &seq_list, &mismask_list);

	cout << "Starting EM loop" << endl;
	double bestLL=-100000;
	uint it_no_improvement(0);
	uint i(0);
	bool l;
	// Start EM loop

	while (it_no_improvement<100) {
		em.do_E_step(1, MCMC_BURN_IN, true);
		double ll = em.get_loglik();
		cout << "LL= " << ll;

		if (ll > bestLL) {
		  cout << " * saving model *" << endl;
		  cout << *mo0 << endl;	
		  l = true;
		  //mdbn.save("discrete_hmm.dbn");
		  bestLL = ll;
		  it_no_improvement=0;
		}
		else { 
		  it_no_improvement++; cout << endl;
		}
		i++;
		em.do_M_step();
		if (l){
		  l=false;
		}
	}

	cout << "DONE" << endl;

	//	mdbn.load("discrete_hmm.dbn");

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
