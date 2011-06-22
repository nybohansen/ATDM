/*
 *  multinomialdensities.cpp
 *
 *  Copyright (C) 2008, Thomas Hamelryck, The Bioinformatics Centre,
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
#include "MultinomialDensities.h"

using namespace std;

namespace mocapy
{

MultinomialDensities::MultinomialDensities()
{
    initialize();
}

void MultinomialDensities::initialize()
{
    type = MULTINOMIAL;
    output_size = dim;
}


MultinomialDensities::MultinomialDensities(uint dim,
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

void MultinomialDensities::set_parameters(MDArray<double> & cpd)
{
    // Check cpd is right shape
    vector<uint> shape = cpd.get_shape();
    assert(shape[1] == dim);

    this->cpd = cpd;
    this->cpd.normalize();
    // make new MultinomialDensity objects
    _create_densities();
}

void MultinomialDensities::construct(vector<uint> & parent_sizes)
{
    // MultinomialNode always has a single parent
    assert(parent_sizes.size()==1);
	node_size = parent_sizes[0];
    // vector of MultinomialDensity objects
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

    cpd.normalize();
    _create_densities();
}

void MultinomialDensities::_create_densities()
{
    // Initialize the density objects and fill the density vector
    for(uint i=0; i<node_size; i++)
    {
        densities[i] = MultinomialDensity(cpd.get_view(i), sample_size);
    }
}

void MultinomialDensities::estimate(vector<MDArray<double> > & ess)
{
	assert(!ess.empty());

    MDArray<double> e = ess[0];
    e.add_inplace(pseudo_count);
	set_parameters(e);
}

vector<double> MultinomialDensities::sample(vector<double> & pv)
{
    uint p=(unsigned int) (pv[0]+0.1);
    return densities[p].sample(randomGen);
}

double MultinomialDensities::get_lik(vector<double> & ptv, bool log_space)
{
    uint p=(unsigned int) (ptv[0]+0.1);
    // Start at 1 (last argument): first value is parent
    double ll = densities[p].get_lik(ptv, true, 1);

    if (log_space)
        return ll;
    else
        return exp(ll);
}

vector<MDArray<double> > MultinomialDensities::get_parameters()
{
    vector<MDArray<double> > v(1);
    v[0] = cpd;
    return v;
}

ostream& operator<<(ostream& output, const MultinomialDensities& a)  {
    output << "Node: Multinomial, size: " << a.node_size << " " << a.dim << endl;

    if (a.cpd.empty())
        output << "Uninitialized." << endl;
    else
        output << a.cpd.tostring();

	return output;
}



}
