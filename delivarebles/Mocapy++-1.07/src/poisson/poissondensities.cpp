/*
 * poissondensities.h
 *
 *  Copyright (C) 2008, Mikael Borg, The Bioinformatics Centre, University of Copenhagen.
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

#include <sstream>
#include "poissondensities.h"
#include <math.h>

using namespace std;

namespace mocapy {
PoissonDensities::PoissonDensities(vector<double> new_user_means) {

    // node_size is the number of states of the discrete node
    user_means = new_user_means;
	initialize();
}

void PoissonDensities::initialize() {
	type = POISSON;
	output_size = 1;
}

// Called in node.construct, and initializes the density arrays
void PoissonDensities::construct(vector<uint> & parent_sizes) {
	// Set node size
	node_size = parent_sizes[0];

	// Check or construct means
	if (!user_means.empty()) {
		assert(user_means.size() == node_size);
		means = user_means;
	}
	else {
		means.resize(node_size);
		for(uint i = 0; i < means.size(); i++)
		  means[i]= 1.0;

		  //	means[i]= randomGen->get_rand()*10;
	}
}

// Parameter estimation based on the ESS
void PoissonDensities::estimate(vector<MDArray<double> > & ess) {
	assert(!ess.empty());
        for(uint p = 0; p < node_size; p++){
                assert(ess[P_S][p] > 0);
                means[p]=ess[P_SY][p]/ess[P_S][p];
        }
}

// Return a sample, based on indicated parent values
vector<double> PoissonDensities::sample(vector<double> & pv) {
        uint parent_value = (uint) pv.front();
        assert(parent_value < means.size());
	vector<double> s;
	s.push_back( poissonsample(means.at(parent_value), randomGen));
	return s;
}

// Return likelihood, that is: P(child|parents)
double PoissonDensities::get_lik(vector<double> & ptv, bool logflag) {
        uint p = (uint) ptv.at(0);
        uint k = (uint) ptv.at(1);
        double lambda = means.at(p);
        double loglambda = log(lambda);
        //double log_dens=-self.lam+k*self.log_lam-gammaln(k+1);
        double log_dens=-lambda+k*loglambda-lgamma(k+1);
	if (logflag)
		return(log_dens);
	else
		return(exp(log_dens));
}

// Return the distribtion's parameters
vector< MDArray<double> > PoissonDensities::get_parameters() {
        assert(false);
        return vector< MDArray<double> >();
}


ostream& operator<<(ostream& output, const PoissonDensities& a)  {

	for(uint i=0; i<a.means.size(); i++){
		output << a.means[i];
	    output << " ";
	}
	return output;
}



}
