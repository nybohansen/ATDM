/*
 *  kentess.cpp
 *
 *  Copyright (C) 2008, Kasper Stovgaard, The Bioinformatics Centre, University of Copenhagen.
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

#include "kentess.h"

using namespace std;

namespace mocapy {

void KentESS::construct(vector<uint> & parent_sizes, uint new_output_dim, uint new_node_size) {

     // Check parent vs child size
     assert(new_node_size==parent_sizes[0]);

     node_size = new_node_size;
     output_dim = new_output_dim;
     ess_size = K_ESSSIZE;

//     cout<<"Kent ess constructed with output_dim, node_size: "<<output_dim<<", "<<node_size<<endl;

     ess_shape.resize(K_ESSSIZE);
     ess_shape[K_N].set_shape(node_size);
     ess_shape[K_YBAR].set_shape(node_size, output_dim);
     ess_shape[K_S].set_shape(node_size, output_dim, output_dim);
     ess_shape[K_DATALENGTH].set_shape(1);

     clear();
}

void KentESS::add_ptv(vector<double> ptv) {

    assert(ptv.size() == 4);

    //cout<<"Kent ess add_ptv "<<ptv<<endl;

    uint p = TOINT(ptv[0]);
    ptv.erase(ptv.begin());

    //cout<<"Kent ess with x, y,z  "<<ptv[0]<<", "<<ptv[1]<<", "<<ptv[2]<<endl;

    ess[K_N][p]++; //count for parent p

    vector<double> ybar;
    ybar.push_back(ptv[0]); ybar.push_back(ptv[1]); ybar.push_back(ptv[2]);
    ess[K_YBAR].get_view(p).add_inplace(ybar);

    vector<vector<double> > S;
    dot2(ybar, ybar, S);
    ess[K_S].get_view(p).add_inplace(S);

    ess[K_DATALENGTH][0]++; //Adding a data point

}

}


