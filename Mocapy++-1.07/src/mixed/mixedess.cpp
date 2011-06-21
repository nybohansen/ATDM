/*
 * discreteess.h
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

#include "mixedess.h"

using namespace std;

namespace mocapy {

// Initializes the ESS arrays.
void MixedESS::construct(vector<uint> & parent_sizes, uint output_dim, uint node_size) {

    //Sorry, only nodes with sizes two are supported.
    assert(node_size == 2);

    //Calculate shape of node
	vector<uint> shape = vec_conc(parent_sizes, node_size);
		
	ess_shape.resize(M_ESSSIZE);
    //M_D = indicator based on value
    ess_shape[M_D].set_shape(shape);
    //M_CV = Gaussian Energy, E and E^2 based on value
    ess_shape[M_CV].set_shape(shape);
          
	ess_size = M_ESSSIZE;
	clear();
}

// Add another sample to the ESS
void MixedESS::add_ptv(vector<double> ptv) {
    //Update the indicator table    
    ess[M_D].get_view( (uint)ptv[PV] )[ (uint) ptv[INDICATOR] ]++;
    if(ptv[INDICATOR]){
        //Update the energy table for the continous case
        ess[M_CV].get_view( (uint)ptv[PV] )[ SUM ] += ptv[ENERGY];        
        ess[M_CV].get_view( (uint)ptv[PV] )[ SUM_SQUARED ] += pow(ptv[ENERGY],2);
    }    
}

}
