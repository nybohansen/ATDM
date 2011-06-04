/*
 * poissoness.h
 *
 *  Copyright (C) 2008, Mikael Borg, The Bioinformatics Centre, University of Copenhagen.
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

#include "poissoness.h"

using namespace std;

namespace mocapy{

// Initializes the ESS arrays.
void PoissonESS::construct(vector<uint> & parent_sizes, uint output_dim, uint node_size) {

	ess_shape.resize(P_ESSSIZE);
	ess_shape[P_S].set_shape(node_size);
	ess_shape[P_SY].set_shape(node_size);
        ess_shape[P_DATALENGTH].set_shape(1);
	ess_size = P_ESSSIZE;
	clear();
}

// Add another sample to the ESS
void PoissonESS::add_ptv(vector<double> ptv) {
        assert(ptv.size() == 2);
        uint p = (uint) ptv.at(0);
        uint y = (uint) ptv.at(1);
        ess[P_DATALENGTH][0]++;
        ess[P_S][p]++;
        ess[P_SY][p]+=y;
}
}
