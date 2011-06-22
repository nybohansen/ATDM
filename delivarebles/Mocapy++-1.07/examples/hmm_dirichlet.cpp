#include <iostream>
#include <stdio.h>
#include <stdlib.h>

#include "mocapy.h"
using namespace mocapy;
using namespace std;

int main1(void) {
     mocapy_seed((uint)94574);
    
     /****
       Estimate the parameters of the sequence 

       ****/
    MDArray<double> a;
    a.set_shape(5);
    a[0] = 0.28;
    a[1] = 5;
    a[2] = 0.32;
    a[3] = 0.55;
    a[4] = .07;
        
    cout << "hello" << endl;
    vector<uint> parent_sizes(1);
    parent_sizes[0] = 1;
    DirichletDensity density = DirichletDensity(a);

    DirichletESS dess = DirichletESS();
    dess.construct(parent_sizes, (uint)5, (uint)1);

    vector<double> ps = density.sample(2);
    cout << ps << endl;
    cout << density.get_parameters() << endl;
    cout << "lik: " << density.get_lik(ps, true, 0) << endl;

    vector<double> ps1(5);
    for (int i = 0; i < 5; i++) 
        ps1[i] = .2;

    cout << "lik1: " << density.get_lik(ps1, true, 0) << endl;
    int num_samples = 500;

    for (int i = 0; i < num_samples; i++) {
        vector<double> sampled = density.sample(1);

        vector<double> b(sampled.size() + 1);
        b[0] = 0;
        for (uint j = 0; j < sampled.size(); j++)
            b[j+1] = sampled[j];
        dess.add_ptv(b);
    }

    dess.save_ess(1.0);
    DirichletDensities densities = DirichletDensities(5, 0, num_samples, true);
    cout << "0" << endl;
    densities.construct(parent_sizes);

    cout << "1" << endl;
    vector <MDArray<double> > local_ess = dess.get_array();

    cout << "2" << endl;
    cout << "local_ess:" << local_ess << endl;
    densities.estimate(local_ess);

    cout << "Estimated:" << densities.get_parameters() << endl;
    exit(0);
    
     ////////////////////////////////////////
     // End estimation
     ///////////////////////////////////////


    
}

int main(void) {
     mocapy_seed((uint)94574);
	// Number of trainining sequences
	int N = 300;

	// Sequence lengths
	int T = 100;

	// Gibbs sampling parameters
	int MCMC_BURN_IN = 10;

	// HMM hidden and observed node sizes
    // Node size indicates how many values a node can take on
	uint H_SIZE=3;
	uint O_SIZE=5;
	bool init_random=false;

    //Conditional probability density
	CPD th0_cpd;
	th0_cpd.set_shape(3); th0_cpd.set_values(vec(0.1, 0.2, 0.7));

    vector<double> vals;
    vals.push_back(0.85);
    vals.push_back(0.10);
    vals.push_back(0.05);
    vals.push_back(0.1);
    vals.push_back(0.2);
    vals.push_back(0.7);
    vals.push_back(0.33);
    vals.push_back(0.33);
    vals.push_back(0.33);

    cout << "hello" << endl;
	CPD th1_cpd;
	th1_cpd.set_shape(3,3); th1_cpd.set_values( vals ); 

    vals.clear();
    vals.push_back(1);
    vals.push_back(2);
    vals.push_back(0.7);
    vals.push_back(5);
    vals.push_back(3);
    vals.push_back(0.15);
    vals.push_back(6);
    vals.push_back(7);
    vals.push_back(8);
    vals.push_back(.1);
    vals.push_back(0.2);
    vals.push_back(10);
    vals.push_back(1);
    vals.push_back(4.5);
    vals.push_back(9);

    cout << "hello1" << endl;
	CPD to0_cpd;
	to0_cpd.set_shape(3,5); to0_cpd.set_values( vals );

	// The target DBN (This DBN generates the data)
	Node* th0 = NodeFactory::new_discrete_node(H_SIZE, "th0", init_random, th0_cpd);
	Node* th1 = NodeFactory::new_discrete_node(H_SIZE, "th1", init_random, th1_cpd);

	Node* to0 = NodeFactory::new_dirichlet_node(O_SIZE, "to0", true);
	((DirichletNode*)to0)->get_densities()->set_parameters(to0_cpd);

	DBN tdbn;
	tdbn.set_slices(vec(th0, to0), vec(th1, to0));

	tdbn.add_intra("th0", "to0");
	tdbn.add_inter("th0", "th1");

	tdbn.construct();

	// The model DBN (this DBN will be trained)
	// For mh0, get the CPD from th0 and fix parameters
	Node* mh0 = NodeFactory::new_discrete_node(H_SIZE, "mh0", init_random, CPD(), th0, true );
	Node* mh1 = NodeFactory::new_discrete_node(H_SIZE, "mh1", init_random);
	Node* mo0 = NodeFactory::new_dirichlet_node(O_SIZE, "mo0", init_random);

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
        //cout << "seq:" << seq_ll.first << endl;
 		seq_list.push_back(seq_ll.first);
        //cout << "mismask:" << mismask << endl;
 		mismask_list.push_back(mismask);
	}
	cout << "Average LL: " << sum_LL/N << endl;

	GibbsRandom mcmc = GibbsRandom(&mdbn);
	EMEngine em = EMEngine(&mdbn, &mcmc, &seq_list, &mismask_list);

	cout << "Starting EM loop" << endl;
	double bestLL=-1000;
	uint it_no_improvement(0);
	uint i(0);
	// Start EM loop
	while (it_no_improvement<80) {
		em.do_E_step(1, MCMC_BURN_IN, true);

		double ll = em.get_loglik();

		cout << "LL= " << ll;

		if (ll > bestLL) {
			cout << " * saving model *" << endl;
			mdbn.save("dirichlet_hmm.dbn");
			bestLL = ll;
			it_no_improvement=0;
		}
		else { it_no_improvement++; cout << endl; }

		i++;
		em.do_M_step();
	}

	cout << "DONE" << endl;

	mdbn.load("dirichlet_hmm.dbn");

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
   cout << "hello" << endl;
}

