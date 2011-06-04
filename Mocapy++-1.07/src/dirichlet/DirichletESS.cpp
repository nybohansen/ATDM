/*
 *  dirichletess.cpp
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

#include "DirichletESS.h"

using namespace std;

namespace mocapy
{

void DirichletESS::construct(vector<uint> & parent_sizes, uint output_dim,
    uint node_size)
{
    dim = output_dim;

    // DirichletNode always has a single parent
    assert(parent_sizes.size()==1);
    // Make sure parent size and child size are compatible
    assert(node_size==parent_sizes[0]);

    // Specify shape of counts array
    vector<uint> shape(2);
    shape[0] = node_size;
    shape[1] = dim;

    // Only one count array
    ess_shape.resize(G_DIRESSSIZE);

    ess_shape[G_SUMOFDATA].set_shape(shape);
    ess_shape[G_SUMOFLOG].set_shape(shape);
    ess_shape[G_DATAPOINTS].set_shape(shape);

    ess_size = G_DIRESSSIZE;
    clear();
}

void DirichletESS::add_ptv(vector<double> ptv)
{
    // Get parent value
    uint p=(unsigned int) (ptv[0]+0.1);

    // Reference to ess array to avoid multiple vector access
    MDArray<double> &ess_sumofdata = ess[G_SUMOFDATA];
    MDArray<double> &ess_sumoflogs = ess[G_SUMOFLOG];
    MDArray<double> &ess_datapoints = ess[G_DATAPOINTS];

    // Put counts from child value into ess
    for (uint i=0; i<dim; i++)
    {
        // Note i+1 because of first value in the vector, which is the
        // parent
        /*
        double new_value = a.get(p, i) + ptv[i+1];
        a.set(p, i, new_value);
        */

        double new_datapoints = ess_datapoints.get(p, i) + 1;
        ess_datapoints.set(p, i, new_datapoints);

        double new_sumoflogs = ess_sumoflogs.get(p, i) + log(ptv[i+1]);
        ess_sumoflogs.set(p, i, new_sumoflogs);
        
        double new_sumofdata = ess_sumofdata.get(p, i) + ptv[i+1];
        ess_sumofdata.set(p, i, new_sumofdata);

    }

}

}
