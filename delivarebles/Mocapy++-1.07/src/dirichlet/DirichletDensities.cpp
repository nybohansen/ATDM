/*
 *  DirichletDensities.cpp
 *
 *  Copyright (C) 2010, Peter Kerpedjiev, The Bioinformatics Centre,
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
#include <iostream>
#include "DirichletDensities.h"
#include <boost/math/special_functions/digamma.hpp>
#include <boost/math/special_functions/zeta.hpp>
#include "../utils/netlib/functions.h"

using namespace std;

namespace mocapy
{

DirichletDensities::DirichletDensities()
{
    initialize();
}

void DirichletDensities::initialize()
{
    type = DIRICHLET;
    output_size = dim;
}


DirichletDensities::DirichletDensities(uint dim,
    double pseudo_count,
    uint sample_size,
    bool init_random)
{
    this->dim = dim;
    this->pseudo_count = pseudo_count;
    this->sample_size = sample_size;
    this->init_random = init_random;
    initialize();
}

void DirichletDensities::set_parameters(MDArray<double> & cpd)
{
    // Check cpd is right shape
    vector<uint> shape = cpd.get_shape();
    assert(shape[1] == dim);

    this->cpd = cpd;
    //this->cpd.normalize();
    // make new DirichletDensity objects
    _create_densities();
}

void DirichletDensities::construct(vector<uint> & parent_sizes)
{
    // DirichletNode always has a single parent
    assert(parent_sizes.size()==1);
	node_size = parent_sizes[0];
    // vector of DirichletDensity objects
    densities.resize(node_size);


    if(cpd.empty())
    {
        vector<uint> shape(2);
        shape[0] = node_size;
        shape[1] = dim;
        cpd.set_shape(shape);

        if(!init_random)
            // Uniform probabilities
            cpd.add_inplace(1.0);
        else
        {
            // Random probabilities
            cpd.randomize(randomGen);
            cpd.add_inplace(0.05);
        }
    }

    _create_densities();
}

void DirichletDensities::_create_densities()
{
    // Initialize the density objects and fill the density vector
    for(uint i=0; i<node_size; i++)
    {
        densities[i] = DirichletDensity(cpd.get_view(i), sample_size);
    }
}

int factorial(int n) {
    int total = 1;

    for (int i = 1; i <= n; i++)
        total *= i;

    return total;
}

/* 
   The polygamma function for n >= 0 
   Only tested for n == 1
   :-(
*/
double polygamma(int n, double x) {
    if (n == 0)
        return boost::math::digamma(x);
    if (n < 0)
        return -1;

    double nn = (double)n;
    double nx = x;
   
    if (x > 0)
      return pow((double)-1, (double)(n+1)) * boost::math::tgamma(n+1.0) * zeta(nn+1, x);
    return -1;

}


double inv_digamma(double y) {
    double x;

    if (y >= -2.22)
        x = exp(y) + 0.5;
    else
        x = -1 / (y - boost::math::digamma(1));

    for (int i = 0; i < 5; i++) 
        x = x - (boost::math::digamma(x) - y)/polygamma(1, x);

   return x;
}

void DirichletDensities::estimate(vector<MDArray<double> > & ess)
{
	assert(!ess.empty());
    uint iter = 3000;
    double tolerance = 0.000001;

    MDArray<double> &ess_sumofdata = ess[G_SUMOFDATA];
    MDArray<double> &ess_sumoflogs = ess[G_SUMOFLOG];
    MDArray<double> &ess_datapoints = ess[G_DATAPOINTS];
    MDArray<double> a;
    MDArray<double> prev_a;
    
    a.set_shape(ess_sumofdata.get_shape());
    prev_a.set_shape(ess_sumofdata.get_shape());


    for (uint i = 0; i < node_size; i++) {
        vector<double> logp(ess_sumoflogs.get_shape()[1]);
        double sum_a = 0, max_diff=0;

        for (uint j = 0; j < dim; j++) {
            a.set(i, j, ess_sumofdata.get(i, j) / ess_datapoints.get(i, j));
            logp[j] = ess_sumoflogs.get(i, j) / ess_datapoints.get(i, j);
            sum_a += a.get(i, j); 
        }

        for (uint j = 0; j < iter; j++) {

            for (uint k = 0; k < dim; k++) {
                double digamma_a;

                digamma_a = boost::math::digamma(sum_a) + logp[k];
                a.set(i,k, inv_digamma(digamma_a));
            }

            sum_a = 0.0;
            max_diff = 0.0;

            for (uint k = 0; k < dim; k++) {
                double curr_a = a.get(i,k);
                double p_a = prev_a.get(i,k);

                double diff = abs(curr_a - p_a); 
                max_diff = diff > max_diff ? diff : max_diff;

                sum_a += a.get(i, k);
                prev_a.set(i,k, curr_a);

            }

            if (max_diff < tolerance) 
                break;      //close enough
        }
    }

	set_parameters(a);
}

vector<double> DirichletDensities::sample(vector<double> & pv)
{
    uint p=(unsigned int) (pv[0]+0.1);
    return densities[p].sample();
}

double DirichletDensities::get_lik(vector<double> & ptv, bool log_space)
{
    uint p=(unsigned int) (ptv[0]+0.1);
    // Start at 1 (last argument): first value is parent
    double ll = densities[p].get_lik(ptv, true, 1);

    if (log_space)
        return ll;
    else
        return exp(ll);
}

vector<MDArray<double> > DirichletDensities::get_parameters()
{
    vector<MDArray<double> > v(1);
    v[0] = cpd;
    return v;
}

ostream& operator<<(ostream& output, const DirichletDensities& a)  {
    output << "Node: Dirichlet, size: " << a.node_size << " " << a.dim << endl;

    if (a.cpd.empty())
        output << "Uninitialized." << endl;
    else
        output << a.cpd.tostring();

	return output;
}



}
