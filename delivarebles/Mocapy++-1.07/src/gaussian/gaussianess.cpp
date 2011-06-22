/*
 * GaussianESS.cpp
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

#include "gaussianess.h"

using namespace std;

namespace mocapy {

// Initializes the ESS arrays.
void GaussianESS::construct(vector<uint> & parent_sizes, uint new_output_dim, uint new_node_size) {
	output_dim = new_output_dim;
	node_size = new_node_size;

	ess_shape.resize(G_ESSSIZE);
	ess_shape[G_S].set_shape(node_size);
	ess_shape[G_SY].set_shape(node_size, output_dim);
	ess_shape[G_SYTY].set_shape(node_size, output_dim, output_dim);
	ess_shape[G_STYY].set_shape(node_size);
	ess_shape[G_DATALENGTH].set_shape(1);

	ess_size = G_ESSSIZE;

	clear();
}


void GaussianESS::save_ess(double weight) {
  if (weight != 1.0) {
    for (uint i=0; i<ess_size; i++) {
      ess[i].multiply_inplace(weight);
    }
  }

  for (uint i=0; i<ess_size; i++) {
    saved_ess[i].add_inplace(ess[i]);
  }


  for (uint i=G_ESSSIZE; i<ess.size(); i++) {
    saved_ess.push_back(ess[i]);
  }
  ess = ess_shape;
}


void GaussianESS::clear() {
  ess = ess_shape;
  saved_ess = ess_shape;

  /*
  if (saved_ess.empty())
    saved_ess = ess_shape;

  for (uint i=0; i<G_ESSSIZE; i++) {
    saved_ess[i] = ess_shape[i];
  }
  */
  //  saved_ess.clear();

}


// Add another sample to the ESS
void GaussianESS::add_ptv(vector<double> ptv) {
	assert(ptv.size() >= 2);

	if (useShrinkage) {
	  // Push back the actual data point (for shrinkage estimation)
	  MDArray<double> mdPtv;
	  mdPtv.set_shape(ptv.size());
	  mdPtv.set_values(ptv);
	  ess.push_back(mdPtv);
	}

	uint p = (uint) ptv.front();

	// We don't need the parent value
	ptv.erase(ptv.begin());

	// One more data point
	ess[G_DATALENGTH][0]++;

	// update s for obs m with parent p
	ess[G_S][p]++;

	// update styy for obs m with parent p
	ess[G_STYY][p] += dot(ptv, ptv);

	// update syty for obs m with parent p
	vector<vector<double> > mda;
	dot2(ptv, ptv, mda);
	ess[G_SYTY].get_view(p).add_inplace(mda);

	// update sy for obs m with parent p
	ess[G_SY].get_view(p).add_inplace(ptv);
}

}
