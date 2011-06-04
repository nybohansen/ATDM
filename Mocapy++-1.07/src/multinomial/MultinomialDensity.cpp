/*
 *  multinomialdensity.cpp
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


#include "MultinomialDensity.h"


namespace mocapy
{

MultinomialDensity::MultinomialDensity(MDArray<double> & p, uint sample_size)
{
    if (sample_size > 0)
        this->sample_size = sample_size;
    else
        this->sample_size = 1;

    // Make copy of p and normalize
    // Dimension of p and count vector
    dim = p.get_shape()[0];

    set_parameters(p);

    // Initialize sample container
    sampled.resize(dim);
    sampled.assign(dim, 0);
}

void MultinomialDensity::set_parameters(MDArray<double> &a)
{
    assert(a.get_shape().size() == 1);
    assert(a.get_shape()[0] == dim);
    assert(dim > 1);

    p = a;
    p.normalize();

    // Make cumulative p for sampling
    cum_p = p;
    cum_p.cumsum();

    // Make log(p) for loglik calculation
    log_p = p;
    log_p.log_all();
}

}
