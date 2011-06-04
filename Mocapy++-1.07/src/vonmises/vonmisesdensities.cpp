/*
 * vonmisesdensities.h
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

#include "vonmisesdensities.h"
#include "vonmisesess.h"

using namespace std;

namespace mocapy {

#define MAX_KAPPA_VALUE 600.0

VonMisesDensities::VonMisesDensities(vector<double> new_user_mus,
		vector<double> new_user_kappas) {

	// node_size: number of mixture components.
	// user_mus: node_size array of mus
	// user_kappas: node_size array kappas
	output_size = 1;
	user_mus = new_user_mus;
	user_kappas = new_user_kappas;

	initialize();
}

void VonMisesDensities::initialize() {
	type = VONMISES;
}


vector<double> VonMisesDensities::make_uniform_kappas(uint node_size) {
	vector<double> kappas(node_size);
	for (uint i = 0; i < node_size; i++) {
		kappas[i] = 0;
	}
	return kappas;
}


vector<double> VonMisesDensities::make_uniform_mus(uint node_size) {
	vector<double> mus(node_size);
	for (uint i = 0; i < node_size; i++) {
		mus[i] = 0;
	}
	return mus;
}


// Called in node.construct, and initializes the density arrays
void VonMisesDensities::construct(vector<uint> & parent_sizes) {
	parent_index = parent_sizes.front();
	node_size = parent_sizes.front();

	// Check or construct mus
	if (!user_mus.empty()) {
		assert(user_mus.size() == node_size);
		mus = user_mus;
	}
	else {
		mus = make_uniform_mus(node_size);
	}

	// Check or construct kappas
	if (!user_kappas.empty()) {
		assert(user_kappas.size() == node_size);
		kappas = user_kappas;
	}
	else {
		kappas = make_uniform_kappas(node_size);
	}

	// Make a list of VonMises objects for sampling and calculations of data likelihood
	vm_list = make_vm_list(node_size, mus, kappas);
}


vector<VonMises> VonMisesDensities::make_vm_list(uint node_size, vector<double> mus, vector<double> kappas) {
	assert(mus.size() == node_size);
	assert(kappas.size() == node_size);
	vector<VonMises> vm_list(node_size);

	for (uint i = 0; i < node_size; i++) {
		vm_list[i] = VonMises(mus[i], kappas[i]);
	}

	return vm_list;
}


// Parameter estimation based on the ESS
void VonMisesDensities::estimate(vector<MDArray<double> > & ess) {
	// Make the estimation
	uint size = ess[VM_N].size();

	vector<double> new_mus(size);
	vector<double> new_kappas(size);
	vector<bool> valid_estimate(size);

	for(uint i=0; i<size; i++) {
		pair<pair<double, double>, bool> result = estimate_mu_kappa(ess[VM_RX][i], ess[VM_RY][i], (uint)ess[VM_N][i]);
		new_mus[i] = result.first.first;
		new_kappas[i] = result.first.second;
		valid_estimate[i] =result.second;
	}

	// Set the kappas and mus with non-valid estimates to previous estimate
	for (uint i=0; i<node_size; i++) {
		if (!valid_estimate[i]) {
			new_mus[i] = mus[i];
			new_kappas[i] = kappas[i];
		}
	}

	// Update the VM list and values
	vm_list = make_vm_list(node_size, new_mus, new_kappas);
	mus = new_mus;
	kappas = new_kappas;
}

// Return a sample, based on indicated parent values
vector<double> VonMisesDensities::sample(vector<double> & pv) {
	assert(pv.size()==1);
	uint parent_value = (uint) pv.front();

	// Get the Von Mises object and sample a value
	VonMises &vm = vm_list[parent_value];
	double sample = vm.sample(randomGen);
	return vec(sample);
}

// Return likelihood, that is: P(child|parents)
double VonMisesDensities::get_lik(vector<double> & ptv, bool log) {
	uint ipv = (uint) ptv.front();

	VonMises &vm  = vm_list[ipv];

	ptv.erase(ptv.begin());
	double ll = vm.get_lik(ptv[0], true);
	return ll;
}

// Return the distribtion's parameters
vector< MDArray<double> > VonMisesDensities::get_parameters() {
	vector<MDArray<double> > parameters;

	MDArray<double> MD_mus(vec<uint>(mus.size()));
	MD_mus.set_values(mus);
	parameters.push_back(MD_mus);

	MDArray<double> MD_kappas(vec<uint>(kappas.size()));
	MD_kappas.set_values(kappas);
	parameters.push_back(MD_kappas);

	return parameters;
}


pair<pair<double, double>, bool> VonMisesDensities::estimate_mu_kappa(double rx, double ry, int n) {
	// Check if a valid estimation can be done
	if(n==0 || (rx==0.0 && ry==0.0)) {
		return make_pair(make_pair(0,0), false);
	}
	double mu, kappa;
	double c = rx / n;
	double s = ry / n;
	double rho = sqrt(c*c + s*s);

	// Estimate mu
	if(s > 0) {
		mu = acos(c / rho);
	}
	else {
		mu = 2*M_PI - acos(c / rho);
	}

	// Estimate kappa (using approximation from P.M. Lee. Bayesian Statistics: An Introduction. Hodder Arnold, 2004.)
	if(rho < 2.0/3.0) {
		kappa = rho * ((2.0 - rho*rho) / (1.0 - rho*rho));
	}
	else {
		kappa = (rho + 1.0) / (4.0 * rho * (1 - rho));
	}

	// Check if the estimated kappa will be greater than MAX_KAPPA_VALUE
	if (kappa > MAX_KAPPA_VALUE) {
		std::cout << "WARNING: Capping kappa from " << kappa << " to " << MAX_KAPPA_VALUE << " in the Von Mises node!" << std::endl;
		kappa = MAX_KAPPA_VALUE;
	}

	return make_pair(make_pair(mu, kappa), true);
}

ostream& operator<<(ostream& output, const VonMisesDensities& a)  {
	output << a.mus << endl;
	output << a.kappas << endl;
	return output;
}


}
