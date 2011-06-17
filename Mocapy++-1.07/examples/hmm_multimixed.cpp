                                     /*
 * hmm_mixed.cpp
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


    //Gaussian used in multimixed_node, this should be made so any node could be used here!
	Node* to0 = NodeFactory::new_gaussian_node(1, "to0");

	// Nodes in slice 1
	Node* h1 = NodeFactory::new_discrete_node(5, "h1");
	Node* o1 = NodeFactory::new_multimixed_node(2, "o1", false, to0);

	// Nodes in slice 2
	Node* h2 = NodeFactory::new_discrete_node(5, "h2");

	// Set architecture
	dbn.set_slices(vec(h1, o1), vec(h2, o1));

	dbn.add_intra("h1", "o1");
	dbn.add_inter("h1", "h2");
	dbn.construct();

	cout << "Loading traindata" << endl;
	GibbsRandom mcmc = GibbsRandom(&dbn);

	EMEngine em = EMEngine(&dbn, &mcmc);

    //Test data for mixed
    em.load_mismask("data/energy_CO_test.mismask");
    em.load_sequences("data/energy_CO_test.data");

    // em.load_mismask("data/mismask.dat");
    // em.load_weights("data/weights.dat");
    // em.load_sequences("data/traindata.dat");

	cout << "Starting EM loop" << endl;
	for (uint i=0; i<100; i++) {
		em.do_E_step(20, 10);
		double ll = em.get_loglik();
		cout << "LL= " << ll << endl;
        em.do_M_step();
	}

	cout << "h1: " << *h1 << endl;
	cout << "o1: " << *o1 << endl;
	cout << "h2: " << *h2 << endl;

	return EXIT_SUCCESS;
}
