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
#include "utils/random_data.h"

using namespace std;

namespace mocapy {

MixedDensities::MixedDensities() {
    // cout << "MixedDensities::MixedDensities()" << endl;
    initialize();
}

MixedDensities::MixedDensities(uint new_node_size, Prior * new_prior, bool new_init_random) {
    // cout << "MixedDensities(uint new_node_size, Prior * new_prior, bool new_init_random)" << endl;
    // node_size is the number of states of the discrete node
    node_size = new_node_size;
    prior = new_prior;
    init_random = new_init_random;
    initialize();
}

MixedDensities::MixedDensities(uint new_node_size, CPD & new_user_cpd, Prior * new_prior) {
	// node_size is the number of states of the discrete node
	node_size = new_node_size;
	prior = new_prior;
	user_cpd = new_user_cpd;
	init_random = false;
    initialize();
}

void MixedDensities::initialize() {
    type = MIXED;
    output_size = 2;     
}

// Normalize CPD and make sure that CPD is 'well'
void MixedDensities::set_cpd(CPD & new_cpd) {
    // cout << "MixedDensities::set_cpd(CPD & new_cpd)" << endl;
    assert(new_cpd.size()> 0);
    cpd = new_cpd;
    cpd.add_inplace(_MIN_TRANSITION);
    cpd.normalize();
    cpd.clip(_MIN_TRANSITION, 1000);
    
    // cout << "cpd = " << cpd << endl;
    
    log_cpd = cpd;
    log_cpd.log_all();
    
    cum_cpd = cpd;
    cum_cpd.cumsum();
    
    
}

CPD MixedDensities::make_random_cpd(vector<uint> & shape, bool no_zeroes) {
    // cout << "MixedDensities::make_random_cpd(vector<uint> & shape, bool no_zeroes) " << endl;
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
    // cout << "MixedDensities::make_uniform_cpd(const vector<uint> & shape)" << endl;
    // Return a uniform CPD.
    CPD c(shape);

    // avoid zero entries
    c.clip(0.1, 2);
    c.normalize();
    return c;
}


// Called in node.construct, and initializes the density arrays
void MixedDensities::construct(vector<uint> & parent_sizes) {
    // cout << "MixedDensities::construct(vector<uint> & parent_sizes)" << endl;

    //Save the parrent size, we need it later
    parent_size = parent_sizes[0];

    //Initialize the mean and variance arrays
    means.set_shape(parent_size);
    variance.set_shape(parent_size);


    means.randomize(randomGen);
    variance.randomize(randomGen);    
            
    
    CPD_shape = vec_conc(parent_sizes, output_size); 
    // cout << "CPD_shape = " << CPD_shape << endl;

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
    // cout << "MixedDensities::estimate(vector<MDArray<double> > & ess)" << endl;
    assert(!ess.empty());

    for(uint i = 0; i < parent_size; i++){
        //For each parrent value, calculate mean and variance
        if( ess[M_CV].get_view(i)[SUM]!=0 ){
            //If the for value i is >0 we know that there must be at least one observation of value i    
            //The total number of observations made with a value i
            uint total = ess[M_D].get_view(i)[CONTINUOUS_TYPE];             
            //Calculate the mean for value i
            means[i] = ess[M_CV].get_view(i)[SUM]/total;
            //Calculate the variance for value i
            variance[i] = ess[M_CV].get_view(i)[SUM_SQUARED]/total - pow(means[i],2); 
        }
    }

    // cout << "ess[M_D]" << endl << ess[M_D] << endl;
    // cout << "ess[M_CV]" << endl << ess[M_CV] << endl;
    // cout << "means" << endl << means << endl;
    // cout << "variance" << endl << variance << endl;
        
    set_cpd(ess[M_D]);
    // cout << "CPD" << endl;
    // cout << cpd << endl;    
}


double MixedDensities::sample_1d_gauss(vector<double> & pv){
    //Draw sample from the gaussian distribution with estimated mean and variance
    double* mean = new double[1];
    double *var = new double[1];
    mean[0] = means[(uint)pv[PV]];
    var[0] = variance[(uint)pv[PV]];    
    double* s = normal_multivariate(1, 1, var, mean, &(randomGen->moc_seed2));
    return s[0];
}

double MixedDensities::sample_discrete(vector<double> & pv){
    return 0;
}


// Return a sample, based on indicated parent values
vector<double> MixedDensities::sample(vector<double> & pv) {
    // cout << "MixedDensities::sample(vector<double> & pv) called with pv = " << pv << endl;

    //Draw random number in [0,1], see if above threshold. If it is, then sample from gauss, otherwise samplediscrete          
    vector<double> choice;

    double r = randomGen->get_rand();
    
    if(cpd.get((uint)pv[PV], DISCRETE_TYPE)<r){
        //Sample from discrete distribution
        choice.push_back( DISCRETE_TYPE );
        choice.push_back( sample_discrete( pv ) );        
    }else{
        //Sample from continouos distribution
        choice.push_back( CONTINUOUS_TYPE );
        choice.push_back( sample_1d_gauss( pv ) );                
    }
    return choice;
}

// Return likelihood, that is: P(child|parents)
double MixedDensities::get_lik(vector<double> & ptv, bool log_space) {  
    // cout << "double MixedDensities::get_lik(vector<double> & ptv, bool log_space)" << endl;
    if(ptv[INDICATOR]){
        
        double a = ptv[ENERGY];
        
        double mu = means[(uint)ptv[PV]]; //Mean
        double v = variance[(uint)ptv[PV]]; //Variance
        double s = sqrt(v); // the standard deviation

        double x = 1/(sqrt(2*M_PI*v))*exp(-pow(a-mu,2)/(2*v));

        x = x*cpd.get((uint)ptv[PV], CONTINUOUS_TYPE);

        if(x== 0){
            x = _MIN_TRANSITION;
            // cout << "x = _MIN_TRANSITION" << x << endl;
        } 

		if (log_space){
            // cout << "Returning likelihood " << log(x) << " for ptv = " << ptv << endl;
		    return log(x);
		}else{
            // cout << "Returning likelihood " << x << " for ptv = " << ptv << endl;
		    return x;
		}
		
    }else{
        if (log_space) {
            // cout << "log_cpd" << endl << log_cpd << endl;
            // cout << "log_cpd " << log_cpd.get((uint)ptv[PV], DISCRETE_TYPE) << endl;
            // cout << "Returning likelihood " << log_cpd.get((uint)ptv[PV], DISCRETE_TYPE) << " for ptv = " << ptv << endl;
            return log_cpd.get((uint)ptv[PV], DISCRETE_TYPE);
        } else {
            // cout << "cpd" << endl << cpd << endl;
            // cout << "cpd " << cpd.get((uint)ptv[PV], DISCRETE_TYPE) << endl;
            // cout << "Returning likelihood " << cpd.get((uint)ptv[PV], DISCRETE_TYPE) << " for ptv = " << ptv << endl;
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
    output << "Node: Mixed, size: ";

    if (a.cpd.empty()) {
        output << "Empty" << endl;
    } else {
        vector<uint> shape = a.cpd.get_shape();
        uint len = shape.size();

        for(uint i=0; i<len; i++){
            output << shape[i] << " ";         
        }
        
        output << endl;

        output <<  "Mixed CPD: " << endl;
        output << a.cpd << endl;
    
        output <<  "Gaussian means: " << endl;
        output << a.means << endl;
    
        output <<  "Gaussian variances: " << endl;
        output << a.variance << endl;
    
    
    }
    return output;
}


}
