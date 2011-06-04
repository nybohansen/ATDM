/*
 * nodefactory.cpp
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

#include "nodefactory.h"

using namespace std;

namespace mocapy {

NodeFactory::NodeFactory() {
}

NodeFactory::~NodeFactory() {
}

DiscreteNode* NodeFactory::new_discrete_node(uint node_size, const char* name, bool init_random, CPD new_cpd, Node* discrete_node, bool fixed) {
	DiscreteNode* n = new DiscreteNode();
	n->set_densities( DiscreteDensities(node_size, NULL, init_random ) );

	assert(!discrete_node || new_cpd.empty());

	if (discrete_node) {
		CPD new_CPD = ((DiscreteNode*)discrete_node)->get_densities()->getCPD();
		n->get_densities()->set_user_cpd(new_CPD);
	}
	if (!new_cpd.empty()) {
		n->get_densities()->set_user_cpd(new_cpd);
	}
	n->fixed=fixed;
	n->set_name(name);
	return n;
}

  GaussianNode* NodeFactory::new_gaussian_node(uint dimension, const char* name, bool init_random, bool shrinkage, MDArray<double> means, MDArray<double> covs, eCovType new_cov_type) {

  bool useShrinkage = shrinkage;

	GaussianNode* n = new GaussianNode();
	GaussianDensities gd = GaussianDensities(dimension, init_random, new_cov_type, false, means, covs);
	gd.useShrinkage=useShrinkage;
	  
	n->set_densities( gd );
	n->set_name(name);
	n->get_ess_ref().useShrinkage = useShrinkage;
	return n;
}

MultinomialNode* NodeFactory::new_multinomial_node(uint dimension, const char* name, bool init_random, double pseudocount) {
	MultinomialNode* n = new MultinomialNode();
	n->set_densities( MultinomialDensities(dimension, pseudocount, 10, init_random) );
	n->set_name(name);
	return n;
}

DirichletNode* NodeFactory::new_dirichlet_node(uint dimension, const char* name, bool init_random, double pseudocount) {
	DirichletNode* n = new DirichletNode();
	n->set_densities( DirichletDensities(dimension, pseudocount, 1, init_random) );
	n->set_name(name);
	return n;
}


  /*
AllCountsMultinomialNode* NodeFactory::new_multinomial_ac_node(uint dimension, const char* name, bool init_random, double pseudocount) {
	AllCountsMultinomialNode* n = new AllCountsMultinomialNode();
	n->set_densities( AllCountsMultinomialDensities(dimension, pseudocount, 10, init_random) );
	n->set_name(name);
	return n;
}
  */

VonMises2dNode* NodeFactory::new_vonmises2d_node(const char* name, MDArray<double> user_mus, MDArray<double> user_kappas) {
	VonMises2dNode* n = new VonMises2dNode();
	n->set_densities(VonMises2dDensities(user_mus, user_kappas) );
	n->set_name(name);
	return n;
}

VonMisesNode* NodeFactory::new_vonmises_node(const char* name, vector<double> user_mus, vector<double> user_kappas) {
	VonMisesNode* n = new VonMisesNode();
	n->set_densities(VonMisesDensities(user_mus, user_kappas) );
	n->set_name(name);
	return n;
}

KentNode* NodeFactory::new_kent_node(const char* name, vector<double> kappas, vector<double> betas, vector<MDArray<double> > es, double kappa_max, double beta_max) {
     KentNode* n = new KentNode();
     n->set_densities( KentDensities(kappas, betas, es, kappa_max, beta_max) );
     n->set_name(name);
     return n;
}

PoissonNode* NodeFactory::new_poisson_node(const char* name, vector<double> new_user_means) {
	PoissonNode* n = new PoissonNode();
	n->set_densities( PoissonDensities(new_user_means) );
	n->set_name(name);
	return n;
}

  BippoNode* NodeFactory::new_bippo_node(const char* name, vector<double> new_lambdas, vector<double> new_thetas,  vector<double> new_ns ) {
	BippoNode* n = new BippoNode();
	n->set_densities( BippoDensities(new_lambdas, new_thetas, new_ns) );
	n->set_name(name);
	return n;
}


}
