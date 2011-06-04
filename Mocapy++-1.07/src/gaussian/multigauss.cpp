/*
 * MultiGauss.cpp
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

#include "multigauss.h"

namespace mocapy {

MultiGauss::MultiGauss(MDArray<double> & new_mean, MDArray<double> & new_sigma) {
	// mean: mean of Gaussian density
	// sigma: sigma of Gaussian density
	sigma = new_sigma;
	mean = new_mean;
	c=0;
	c2=0;
	b=0;
	sq_det=0;
	logb=0;
	
	if (mean.get_values().size() == 1) {
		// 1D Gaussian
        assert(sigma.size()==1);
        double v = sigma.get_values_flat().front(); // the variance
        double s = sqrt(v); // the standard deviation

		c = -log(s * sqrt(2 * M_PI));
		c2 = 2.0 * s * s;
		one_dim = true;
	} else {
		// Multidimensional Gaussian
	  //	assert(sigma.det()> 0);
		dim = mean.get_shape()[0];
		assert(sigma.get_shape().front() == dim);
		assert(sigma.get_shape().back() == dim);
		assert(sigma.get_shape().size() == 2);
		double d = sigma.det();
		sq_det = sqrt(d);
		b = 1.0 / (pow(2.0 * M_PI, (double) dim / 2.0) * sq_det);
		logb = log(b);
		inv_sigma.copy(sigma);
		inv_sigma.inv();
		one_dim = false;
	}
}

}
