/*
 *  vonmises.cpp
 *
 *  Copyright (C) 2008, Jes Frellsen, The Bioinformatics Centre,
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

#include "vonmises.h"
#include <boost/math/special_functions/detail/bessel_i0.hpp>

namespace mocapy {

VonMises::VonMises(double new_mu, double new_kappa) : mu(new_mu), kappa(new_kappa) {
	// Variables for the density function
	log_denominator = log(2.0 * M_PI * boost::math::detail::bessel_i0(kappa));

	// Variables for the sampler
	a = 1.0 + sqrt(1.0 + 4.0 * kappa * kappa);
	b = (a - sqrt(2.0 * a)) / (2.0 * kappa);
	r = (1.0 + b * b) / (2.0 * b);
}

}
