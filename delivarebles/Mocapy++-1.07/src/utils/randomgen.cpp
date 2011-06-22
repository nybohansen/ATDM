/*
 * randomgen.cpp
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


#include "randomgen.h"

namespace mocapy {

// Initialize the random number generators that are used by Mocapy.
void RandomGen::mocapy_seed(uint s) {
    // srand(s);
    moc_seed1 = s;
    moc_seed2 = s;
    rng.seed(moc_seed1);
    rng_lagged_fibonacci607.seed(moc_seed1);
}

double RandomGen::get_rand(const bool use_lagged_fibonacci607_flag) {
    double result=0;
    if (use_lagged_fibonacci607_flag==true) {
        // lagged_fibonacci607
        assert(r_lagged_fibonacci607);
        result = (*r_lagged_fibonacci607)();
    }
    else {
        // Mersenne Twister
        assert(r);
        result = (*r)();
    }
    return result;

    // Using C standard rand()
    //        return ((double) rand()) / ((double) (RAND_MAX + 1.0));
}

double RandomGen::get_rand_normal(
        double mean, double sigma,
        const bool use_lagged_fibonacci607_flag
        ) {
    double result=0;
    if (use_lagged_fibonacci607_flag==true) {
        // select Gaussian probability distribution
        boost::normal_distribution<double> norm_dist(mean, sigma);
        // bind random number generator to distribution, forming a function
        boost::variate_generator<boost::lagged_fibonacci607&, boost::normal_distribution<double> >  normal_sampler(rng_lagged_fibonacci607, norm_dist);
        // sample from the distribution
        result = normal_sampler();
    }
    else {
        // class default generator?
        if (mean==0.0 and sigma==1.0) {
            result = (*norm)();  // class default generator
        }
        // not default: build a new generator
        else {  
            // select Gaussian probability distribution
            boost::normal_distribution<double> norm_dist(mean, sigma);
            // bind random number generator to distribution, forming a function
            boost::variate_generator<boost::mt19937&, boost::normal_distribution<double> >  normal_sampler(rng, norm_dist);
            // sample from the distribution
            result = normal_sampler();
        }
    }
    return result;
}

}
