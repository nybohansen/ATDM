/*
 * bippodensities.h
 *
 *  Copyright (C) 2009, Thomas Hamelryck, The Bioinformatics Centre, 
 *  University of Copenhagen.
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
#include <math.h>
#include "bippodensities.h"

using namespace std;

namespace mocapy {
BippoDensities::BippoDensities(vector<double> new_user_lambdas, 
    vector<double> new_user_thetas, vector<double> new_user_ns) {

    // node_size is the number of states of the discrete node
    user_lambdas = new_user_lambdas;
    user_thetas = new_user_thetas;
    user_ns = new_user_ns;
	initialize();
}

void BippoDensities::initialize() {
	type = BIPPO;
	output_size = 1;
}

// Called in node.construct, and initializes the density arrays
void BippoDensities::construct(vector<uint> & parent_sizes) {
	// Set node size
	node_size = parent_sizes[0];

	// Check or construct means
	if (!user_lambdas.empty() && !user_thetas.empty() && !user_ns.empty()) {
		assert(user_lambdas.size() == node_size);
		assert(user_thetas.size() == node_size);
		assert(user_ns.size() == node_size);
		lambdas = user_lambdas;
		thetas = user_thetas;
		ns = user_ns;
	}
	else {
		lambdas.resize(node_size);
		thetas.resize(node_size);
		ns.resize(node_size);
		for(uint i = 0; i < node_size; i++) {
			lambdas[i]= 0.1+5*randomGen->get_rand(); 
			thetas[i]= 0.5;  
			ns[i]= 5.0;
        }
	}

    // Now make density objects (one for each parent value)
    densities.resize(node_size);
    for(uint i = 0; i < node_size; i++)
        densities[i]=Bippo(lambdas[i], thetas[i], (uint) ns[i]);
}

// Parameter estimation based on the ESS
void BippoDensities::estimate(vector<MDArray<double> > & ess) {
	assert(!ess.empty());
    for(uint p = 0; p < node_size; p++) {
            double m1, m2, m3;
            double s1, s2, s3, count;
            s1 = ess[BP_S][p];
            s2 = ess[BP_S2][p];
            s3 = ess[BP_S3][p];
            count = ess[BP_N][p];

            if (count>20)
            {
                //mean
                m1 = s1/count;
                // second moment
                m2 = s2/count-m1*m1; 
                // third moment
                m3 = (s3-3*m1*s2)/count+2*m1*m1*m1;

                // Update density objects
                Bippo &dens = densities[p];
                dens.estimate(m1, m2, m3);
                vector<double> parameters = dens.get_parameters();
                // update parameters
                lambdas[p]=parameters[0];
                thetas[p]=parameters[1];
                ns[p]=parameters[2];
            }
    }
}

// Return a sample, based on indicated parent values
vector<double> BippoDensities::sample(vector<double> & pv) {
    uint parent_value = (uint) pv.front();
    assert(parent_value < node_size);
	vector<double> s;
    Bippo &dens = densities[parent_value];
	s.push_back( dens.sample(randomGen) );
	return s;
}

// Return likelihood, that is: P(child|parents)
double BippoDensities::get_lik(vector<double> & ptv, bool logflag) {
    uint p = (uint) ptv.at(0);
    uint k = (uint) ptv.at(1);
    Bippo &dens = densities[p];
    return dens.get_lik(k, logflag);
}

// Return the distribution's parameters
vector< MDArray<double> > BippoDensities::get_parameters() {
    vector< MDArray<double> > pars(node_size);
    MDArray<double> ar;
    ar.set_shape(3);
	for(uint i=0; i<node_size; i++){
        ar[0]=lambdas[i]; 
        ar[1]=thetas[i]; 
        ar[2]=(double) ns[i]; 
        pars.push_back(ar);
    }
    return pars;
}

ostream& operator<<(ostream& output, const BippoDensities& a)  {
	for(uint i=0; i<a.node_size; i++){
		output << a.lambdas[i] << " " << a.thetas[i] << " " << a.ns[i] << endl;
	}
	return output;
}

}
