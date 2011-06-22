/*
 *  GDMess.cpp
 *
 *  Copyright (C) 2010, Simon H. Rasmussen, The Bioinformatics Centre,
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

#include "GDMess.h"

namespace mocapy
{


// ESS has been used in M-step - get ready for new cycle
void GDMess::clear() {
  ess.clear();
  saved_ess.clear();
}


void GDMess::construct(vector<uint> & parent_sizes, uint output_dim,
    uint node_size)
{

  // Each element in the ESS is a data point. 
  // A data point is a 1-dim MDArray with ptvs
  clear();
}

void GDMess::save_ess(double weight) {
  for (uint i=0; i<ess.size(); i++) {
    saved_ess.push_back(ess[i]);
  }

  ess.clear();
}


void GDMess::add_ptv(vector<double> ptv)
{
  MDArray<double> mdPtv;
  mdPtv.set_shape(ptv.size());
  mdPtv.set_values(ptv);
  ess.push_back(mdPtv);
}


}
