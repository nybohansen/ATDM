/*
 * hmm_mixed4.cpp
 *
 *  Created on: 16. jun 2011
 *      Author: Kasper Nybo Hansen
 */

#include <iostream>
#include <stdio.h>
#include <stdlib.h>

#include "mocapy.h"
using namespace mocapy;
using namespace std;

int main(void) {
	// The dynamic Bayesian network
	DBN dbn;

	uint H_SIZE=5;
	uint O_SIZE=2;

	// Nodes in slice 1
	Node* h1 = NodeFactory::new_discrete_node(H_SIZE, "h1");
	Node* o1 = NodeFactory::new_mixed_node(O_SIZE, "o1");

	// Nodes in slice 2
	Node* h2 = NodeFactory::new_discrete_node(H_SIZE, "h2");

	// Set architecture
	dbn.set_slices(vec(h1, o1), vec(h2, o1));

	dbn.add_intra("h1", "o1");
	dbn.add_inter("h1", "h2");
	dbn.construct();

	cout << "Loading traindata" << endl;
	GibbsRandom mcmc = GibbsRandom(&dbn);

	EMEngine em = EMEngine(&dbn, &mcmc);

    //Test data for mixed
    em.load_mismask("data/energy_NH.mismask");
    em.load_sequences("data/energy_NH.data");


	double bestLL=-1000;
	uint it_no_improvement(0);

	cout << "Starting EM loop" << endl;
	while (it_no_improvement<30) {
        //While we don't converge
		em.do_E_step(1, 10, true);
		double ll = em.get_loglik();

        cout.precision(15);
		cout << "LL= " << ll ;
		if (ll > bestLL) {
			cout << " * saving model *" << endl;
			dbn.save("mixed_hmm.dbn");
			bestLL = ll;
			it_no_improvement=0;
		} else { 
            cout << endl; 
            it_no_improvement++;
        }
        cout.precision(5);
        em.do_M_step();
    }

	    
	cout << "h1: " << *h1 << endl;
	cout << "o1: " << *o1 << endl;
	cout << "h2: " << *h2 << endl;


	return EXIT_SUCCESS;
}
