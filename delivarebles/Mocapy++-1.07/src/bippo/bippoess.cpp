/*
 * bippoess.h
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

#include "bippoess.h"

using namespace std;

namespace mocapy{

// Initializes the ESS arrays.
void BippoESS::construct(vector<uint> & parent_sizes, uint output_dim, 
    uint node_size) {

	ess_shape.resize(BP_ESSSIZE);
	ess_shape[BP_N].set_shape(node_size);  // Nr of data points
	ess_shape[BP_S].set_shape(node_size);  // Sum for mean
	ess_shape[BP_S2].set_shape(node_size); // Sum of squares for second moment 
	ess_shape[BP_S3].set_shape(node_size); // Sum of cubes for third moment
	ess_size = BP_ESSSIZE;
	clear();
}

// Add another sample to the ESS
void BippoESS::add_ptv(vector<double> ptv) {
        assert(ptv.size() == 2);
        // parent hidden node value
        uint p = (uint) ptv.at(0);
        // bippo node value
        double y = ptv.at(1);
        // save count
        ess[BP_N][p]++;
        double y2=y*y; // for efficiency
        ess[BP_S][p] += y;
        ess[BP_S2][p] += y2;
        ess[BP_S3][p] += y2*y;
}
}
