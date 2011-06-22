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
    initialize();
}

MixedDensities::MixedDensities(uint new_node_size, bool new_init_random, MDArray<double> user_means, MDArray<double> user_variance){
    node_size = new_node_size;
    init_random = new_init_random;
    means = user_means;
    variance = user_variance;
    initialize();
}

MixedDensities::MixedDensities(uint new_node_size, CPD & new_user_cpd, MDArray<double> user_means, MDArray<double> user_variance){
	node_size = new_node_size;
	user_cpd = new_user_cpd;
	init_random = false;
    means = user_means;
    variance = user_variance;	
    initialize();
}

void MixedDensities::initialize() {
    type = MIXED;
    output_size = 2;     
}

// Normalize CPD and make sure that CPD is 'well'
void MixedDensities::set_cpd(CPD & new_cpd) {
    assert(new_cpd.size()> 0);
    cpd = new_cpd;
    cpd.add_inplace(_MIN_TRANSITION);
    cpd.normalize();
    cpd.clip(_MIN_TRANSITION, 1000);
    
    log_cpd = cpd;
    log_cpd.log_all();
    
    cum_cpd = cpd;
    cum_cpd.cumsum();    
}

CPD MixedDensities::make_random_cpd(vector<uint> & shape, bool no_zeroes) {
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
    // Return a uniform CPD.
    CPD c(shape);

    // avoid zero entries
    c.clip(0.1, 2);
    c.normalize();
    return c;
}


// Called in node.construct, and initializes the density arrays
void MixedDensities::construct(vector<uint> & parent_sizes) {
    //Save the parrent size, we need it later
    parent_size = parent_sizes[0];
    //Initialize the mean and variance arrays
    means.set_shape(parent_size);
    variance.set_shape(parent_size);
    //Randomize the arrays
    means.randomize(randomGen);
    variance.randomize(randomGen);    
    
    //Calculate the CPD shape            
    CPD_shape = vec_conc(parent_sizes, output_size); 

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
    
    set_cpd(ess[M_D]);
}


double MixedDensities::sample_1d_gauss(vector<double> & pv){
    //Draw sample from the gaussian distribution with estimated mean and variance
    double* mean = new double[1];
    double *sd = new double[1];
    mean[0] = means[(uint)pv[PV]];
    sd[0] = sqrt(variance[(uint)pv[PV]]);    
    double* s = normal_multivariate(1, 1, sd, mean, &(randomGen->moc_seed2));
    return s[0];
}

double MixedDensities::sample_discrete(vector<double> & pv){
    //Dummy function, just return 0 since in the discrete case the only energy is 0.
    return 0;
}


// Return a sample, based on indicated parent values
vector<double> MixedDensities::sample(vector<double> & pv) {
    //Draw random number in [0,1], see if above threshold. If it is, then sample from gauss, otherwise samplediscrete          
    vector<double> choice;

    double r = randomGen->get_rand();
    if(cpd.get((uint)pv[PV], DISCRETE_TYPE)>r){
        //Sample from discrete distribution
        choice.push_back( DISCRETE_TYPE );
        choice.push_back( sample_discrete( pv ) );        
    }else{
        //Sample from continouos distribution
        choice.push_back( CONTINUOUS_TYPE );
        choice.push_back( sample_1d_gauss( pv ) );                
    }             
    cout << choice << endl;
    return choice;
}

// Return likelihood, that is: P(child|parents)
double MixedDensities::get_lik(vector<double> & ptv, bool log_space) {  
    if(ptv[INDICATOR]){
        
        double a = ptv[ENERGY]; //Sample point we would like to test the likelihood of
        
        double mu = means[(uint)ptv[PV]]; //Mean
        double v = variance[(uint)ptv[PV]]; //Variance
        double s = sqrt(v); // the standard deviation

        //Calculate the gaussian likelihood
        double x = 1/(sqrt(2*M_PI*v))*exp(-pow(a-mu,2)/(2*v));

        //Multiply with the probability that we see a continous node
        x = x*cpd.get((uint)ptv[PV], CONTINUOUS_TYPE);

        if(x == 0){
            //We are in trouble if x==0 since log(0) is illegal
            //We still want to "Punish" the system so set x to a very low number
            x = _MIN_TRANSITION;
        } 

		if (log_space){
		    return log(x);
		}else{
		    return x;
		}
		
    }else{
        if (log_space) {
            return log_cpd.get((uint)ptv[PV], DISCRETE_TYPE);
        } else {
            return cpd.get((uint)ptv[PV], DISCRETE_TYPE);
        }        
    }

}

// Return the distribtion's parameters
vector<MDArray<double> > MixedDensities::get_parameters() {
   cout << "MixedDensities::get_parameters()" << endl;
   vector<MDArray<double> > ret;
   //Return the CPD
   ret.push_back(cpd);
   //Return the means
   ret.push_back(means);   
   //Return the variance
   ret.push_back(variance);
   return ret;
}


void MixedDensities::set_user_cpd(CPD & new_user_cpd) {
    user_cpd = new_user_cpd;
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
