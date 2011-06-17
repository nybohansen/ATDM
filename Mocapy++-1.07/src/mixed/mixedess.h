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
#ifndef MIXEDESS_H_
#define MIXEDESS_H_

#include "../framework/essbase.h"

// Expected Sufficient Statistics

namespace mocapy {

// Indices for the ess std::vector
//M_D = Table of indicator values
//M_V = Table of values of the
enum MIXED_ESS_INDEX {M_D, M_V, M_DATALENGTH, M_ESSSIZE};
//Enum used to define how the PTV is indexed
//PV = Parent Value, Indicator = indicator, ENERGY = energy 
enum PTV_INDEX {PV, INDICATOR, ENERGY};
enum M_V_INDEX {SUM, SUM_SQUARED};
enum M_D_INDEX {DISCRETE_TYPE, CONTINUOUS_TYPE};

class MixedESS: public ESSBase {
public:
	MixedESS() {};
	virtual ~MixedESS() {};

	// Initializes the ESS arrays.
	void construct(std::vector<uint> & parent_sizes, uint output_dim, uint node_size);

	// Add another sample to the ESS
	void add_ptv(std::vector<double> ptv);
};

}

#endif /* MIXEDESS_H_ */
