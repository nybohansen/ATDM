/*
 *  dirichletdensity.cpp
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
#include "DirichletDensity.h"
#include <boost/math/special_functions/gamma.hpp>

namespace mocapy
{

DirichletDensity::DirichletDensity(MDArray<double> & a, uint sample_size)
{
    if (sample_size > 0)
        this->sample_size = sample_size;
    else
        this->sample_size = 1;

    dim = a.get_shape()[0];

    set_parameters(a);

    // Initialize sample container
    sampled.resize(dim);
    sampled.assign(dim, 0);
}

void DirichletDensity::set_parameters(MDArray<double> &a_new)
{
    assert(a_new.get_shape().size() == 1);
    assert(a_new.get_shape()[0] == dim);
    assert(dim > 1);

    a = a_new;
    am.set_shape(a.get_shape());

    double y = 0;
    double a_sum = 0;

    for (uint i = 0; i < dim; i++)  {
        y += boost::math::lgamma(a[i]);
        a_sum += a[i];
        am[i] = a[i] - 1;
    }
   
    lnc = boost::math::lgamma(a_sum) - y;
}

}
