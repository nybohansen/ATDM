/*
 * mixeddensities.h
 *
 *  Copyright (C) 2011, Kasper Nybo Hansen, Department of Computer Science, University of Copenhagen.
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

#include "mixeddensities.h"

using namespace std;

namespace mocapy {

MixedDensities::MixedDensities() {
    // cout << "MixedDensities::MixedDensities()" << endl;
}

MixedDensities::MixedDensities(uint new_node_size, Prior * new_prior, bool new_init_random) {
    // cout << "MixedDensities(uint new_node_size, Prior * new_prior, bool new_init_random)" << endl;
    // node_size is the number of states of the discrete node
    node_size = new_node_size;
    prior = new_prior;
    init_random = new_init_random;
}

MixedDensities::MixedDensities(uint new_node_size, CPD & new_user_cpd, Prior * new_prior) {
	// node_size is the number of states of the discrete node
	node_size = new_node_size;
	prior = new_prior;
	user_cpd = new_user_cpd;
	init_random = false;
}



// Normalize CPD and make sure that CPD is 'well'
void MixedDensities::set_cpd(CPD & new_cpd) {
    cout << "MixedDensities::set_cpd(CPD & new_cpd)" << endl;
    assert(new_cpd.size()> 0);
    cpd = new_cpd;
    cpd.add_inplace(_MIN_TRANSITION);
    cpd.normalize();
    cpd.clip(_MIN_TRANSITION, 1000);
    
    cout << "cpd = " << cpd << endl;
    
    log_cpd = cpd;
    log_cpd.log_all();
    
    cum_cpd = cpd;
    cum_cpd.cumsum();
    
    
}

CPD MixedDensities::make_random_cpd(vector<uint> & shape, bool no_zeroes) {
    cout << "MixedDensities::make_random_cpd(vector<uint> & shape, bool no_zeroes) " << endl;
    // Return a random CPD.
    CPD c(shape);
    c.randomize(randomGen);

    // avoid zero entries
    if (no_zeroes) {
        c.clip(0.01, 2);
    }
    c.normalize();
    return c;
}

CPD MixedDensities::make_uniform_cpd(const vector<uint> & shape) {
    cout << "MixedDensities::make_uniform_cpd(const vector<uint> & shape)" << endl;
    // Return a uniform CPD.
    CPD c(shape);

    // avoid zero entries
    c.clip(0.1, 2);
    c.normalize();
    return c;
}


// Called in node.construct, and initializes the density arrays
void MixedDensities::construct(vector<uint> & parent_sizes) {
    cout << "MixedDensities::construct(vector<uint> & parent_sizes)" << endl;

    //Save the parrent size, we need it later
    parent_size = parent_sizes[0];

    //Initialize the mean and variance arrays
    means.set_shape(parent_size);
    variance.set_shape(parent_size);
        
    type = MIXED;
    output_size = 2;     
    
    CPD_shape = vec_conc(parent_sizes, output_size); 
    cout << "CPD_shape = " << CPD_shape << endl;

    if(user_cpd.empty()) {
        CPD cpd;
        if(!init_random) {
            cpd = make_uniform_cpd(CPD_shape);
        }
        else {
            cpd = make_random_cpd(CPD_shape, true);
        }
        set_cpd(cpd);
    } else {
        assert(user_cpd.get_shape() == CPD_shape);
        set_cpd(user_cpd);
    }
}

// Parameter estimation based on the ESS
void MixedDensities::estimate(vector<MDArray<double> > & ess) {
    cout << "MixedDensities::estimate(vector<MDArray<double> > & ess)" << endl;
    assert(!ess.empty());

    for(uint i = 0; i < parent_size; i++){
        //For each parrent value, calculate mean and variance
        if( ess[M_V].get_view(i)[SUM]!=0 ){
            //If the for value i is >0 we know that there must be at least one observation of value i    
            //The total number of observations made with a value i
            uint total = ess[M_D].get_view(i)[DISCRETE_TYPE]+ess[M_D].get_view(i)[CONTINUOUS_TYPE];             
            //Calculate the mean for value i
            means[i] = ess[M_V].get_view(i)[SUM]/total;
            //Calculate the variance for value i
            variance[i] = ess[M_V].get_view(i)[SUM_SQUARED]/total - pow(means[i],2); 
        }
    }

    // cout << "ess[M_D]" << endl << ess[M_D] << endl;
    // cout << "ess[M_V]" << endl << ess[M_V] << endl;
    // cout << "means" << endl << means << endl;
    // cout << "variance" << endl << variance << endl;
        
    set_cpd(ess[M_D]);
}

// Return a sample, based on indicated parent values
vector<double> MixedDensities::sample(vector<double> & pv) {
     cout << "MixedDensities::sample(vector<double> & pv)" << endl;
    MDArray<double>* cumulative;
    if (pv.empty())
     cumulative = &cum_cpd;
    else {
     vector<uint> ipv;
     toint(pv, ipv);
     cumulative = &(cum_cpd.get_view(ipv));
    }
    
    double r = randomGen->get_rand();
    
    assert(r<=1 && r>=0);
    vector<double> choice;
    choice.push_back( cumulative->bisect(r) );
    return choice;
}

// Return likelihood, that is: P(child|parents)
double MixedDensities::get_lik(vector<double> & ptv, bool log) {  
    cout << "MixedDensities::get_lik called with: " << ptv << " and log = " << log << endl;
    if(ptv[INDICATOR]){
        return 0;
    }else{
        if (log) {
            cout << "log_cpd" << endl << log_cpd << endl;
            cout << "log_cpd " << log_cpd.get((uint)ptv[PV], DISCRETE_TYPE) << endl;
            return log_cpd.get((uint)ptv[PV], DISCRETE_TYPE);
        } else {
            cout << "cpd" << endl << cpd << endl;
            cout << "cpd " << cpd.get((uint)ptv[PV], DISCRETE_TYPE) << endl;
            return cpd.get((uint)ptv[PV], DISCRETE_TYPE);
        }        
    }

}

// Return the distribtion's parameters
vector<MDArray<double> > MixedDensities::get_parameters() {
    cout << "MixedDensities::get_parameters()" << endl;
   vector<MDArray<double> > ret;
   ret.push_back(cpd);
   return ret;
}


void MixedDensities::set_user_cpd(CPD & new_user_cpd) {
    cout << "MixedDensities::set_user_cpd(CPD & new_user_cpd)" << endl;
    user_cpd = new_user_cpd;
}

void MixedDensities::set_prior(Prior * new_prior) {
    cout << "MixedDensities::set_prior(Prior * new_prior)" << endl;
    prior = new_prior;
}

//Override operator <<
ostream& operator<<(ostream& output, const MixedDensities& a)  {
    cout << "Override operator << called";
	output << "Node: Discrete, size: ";

    if (a.cpd.empty()) {
     output << "Empty" << endl;
    }
    else {
     vector<uint> shape = a.cpd.get_shape();
     uint len = shape.size();
    
     for(uint i=0; i<len; i++)
         output << shape[i] << " ";
    
            output << endl;
    
    // a.cpd.tostring() does not work properly when the shape size is more than 2
    //        output << a.cpd.tostring();
    
    output << a.cpd << endl;
    }

	return output;
}


}
