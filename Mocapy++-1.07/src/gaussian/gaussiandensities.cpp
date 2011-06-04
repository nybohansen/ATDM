/*
 * gaussiandensities.cpp
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

#include "gaussiandensities.h"
#include "gaussianess.h"
#include "../utils/random_data.h"
#include "shrinkage/Shrinkage.h"

using namespace std;

namespace mocapy {

  GaussianDensities::GaussianDensities(uint dim, bool new_init_random,
				       eCovType new_cov_type, bool new_cov_tied,
				       MDArray<double> new_user_means, MDArray<double> new_user_covs) {
    // Dimension of output vector
    output_size = dim;
    // Covariance tied?
    cov_tied = new_cov_tied;
    // Covariance SPHE, DIAG or FULL
    cov_type = new_cov_type;
    user_means = new_user_means;
    user_covs = new_user_covs;
    if (cov_type >= WRONG_COV_TYPE) {
      throw MocapyExceptions("Unknown covariance type.");
    }

    init_random = new_init_random;
    initialize();
    useShrinkage=false;
  }

  void GaussianDensities::initialize() {
    type = GAUSSIAN;
  }

  vector<MultiGauss> GaussianDensities::make_mgauss_list(MDArray<double> & means, MDArray<double> & covs) {
    vector<MultiGauss> mgauss_list;
    for (uint i = 0; i < node_size; i++) {
      MDArray<double> & m = means.get_view(i);
      MDArray<double> & c = covs.get_view(i);
      MultiGauss mgauss = MultiGauss(m, c);
      mgauss_list.push_back(mgauss);
    }
    return mgauss_list;
  }

  MDArray<double> GaussianDensities::make_rnd_covs(uint node_size, uint dim,
						   vector<uint> & cov_shape, bool tied) {
    // Return random covariance matrices.

    // node_size: number of mixture components
    // dim: dimension of distribution
    // cov_shape: shape of array of covs (node_size, dim, dim)
    // tied: flags tied or untied covariance matrices

    MDArray<double> covs(cov_shape);
    MDArray<double> c(cov_shape);
    for (uint i = 0; i < node_size; i++) {
      if (i == 0 || !tied) {
	if (cov_type == SPHE) {
	  // Make sure entries in c>0
	  double r = randomGen->get_rand();
	  double r2 = (r + 0.1) / 1.1;
	  c.eye(dim, r2);
	} else {
	  c.rnd_cov(dim, randomGen);
	  if (cov_type == DIAG)
	    c.keep_eye();
	}
      }
      covs.set(i, c);
    }
    return covs;
  }

  MDArray<double> GaussianDensities::make_uniform_covs(uint node_size, uint dim,
						       vector<uint> & cov_shape) {
    // Return uniform covariance matrices.
    // node_size: number of mixture components
    // dim: dimension of distribution
    // cov_shape: shape of array of covs (node_size, dim, dim)

    MDArray<double> covs(cov_shape);
    MDArray<double> c(cov_shape);
    for (uint i = 0; i < node_size; i++) {
      c.eye(dim, 1);
      covs.set(i, c);
    }
    return covs;
  }

  // Initializes the Density arrays
  void GaussianDensities::construct(vector<uint> & parent_sizes) {
    // Set node size
    node_size = parent_sizes[0];
    assert(parent_sizes.size() == 1);
    parent_index = parent_sizes.front();

    mean_shape = vec(node_size, output_size);
    cov_shape = vec(node_size, output_size, output_size);

    // Check or construct means
    if (!user_means.empty()) {
      means = user_means;
      assert(means.get_shape() == mean_shape);
      user_means.clear();
    } else {
      means.set_shape(node_size, output_size);
      if (init_random) {
	means.randomize(randomGen);
	means.multiply_inplace(50.0);
      }
    }

    // Check or construct covariances
    if (!user_covs.empty()) {
      covs = user_covs;
      assert(covs.get_shape() == cov_shape);
    } else {
      if (init_random) {
	//  Random covarience matrices
	covs = make_rnd_covs(node_size, output_size, cov_shape, cov_tied);
      }
      else {
	//  Uniform covarience matrices
	covs = make_uniform_covs(node_size, output_size, cov_shape);
      }
    }
    // Make a list of MultiGauss objects for calculations of P(obs)
    mgauss_list = make_mgauss_list(means, covs);
  }

  // Parameter estimation based on the ESS
  void GaussianDensities::estimate(vector<MDArray<double> > & ess) {
    // Backup old data - new covs might be illegal
    MDArray<double> old_means = means;
    MDArray<double> old_covs = covs;

    means = calc_means(ess);
    if (cov_tied && cov_type==FULL) {
      covs = calc_cov_full_tied(ess);
    }
    else if (!cov_tied && cov_type==FULL) {
      covs = calc_cov_full(ess);
    }
    else {
      // Implement other estimation functions
      assert(false);
      exit(0);
    }

    for (uint i=0; i<node_size; i++) {
      MDArray<double> & c = covs.get_view(i);
      bool ok;

      if (useShrinkage)
	ok = c.det() > 0 && ess[G_S][i] >= 1;
      else
	ok = c.det() > 0 && ess[G_S][i] > output_size+1;

      if (!ok) {
	covs.get_view(i) = old_covs.get_view(i);
	means.get_view(i) = old_means.get_view(i);
      }
    }

    // update MultiGauss list
    mgauss_list = make_mgauss_list(means, covs);
  }

  MDArray<double> GaussianDensities::multivariate_normal(MDArray<double> & m, MDArray<double> & c) {
    assert(!m.get_shape().empty());

    uint dim = m.get_shape().front();
    double* mean = new double[dim];
    double *cov = new double[dim * dim];

    for (uint i = 0; i < dim; i++) {
      mean[i] = m[i];
      for (uint j = 0; j < dim; j++) {
	double v = c.get(i, j);
	cov[i * dim + j] = v;
      }
    }
    dpofa(cov, dim, dim);
    // Must point to an int
    double* s = normal_multivariate(dim, 1, cov, mean, &(randomGen->moc_seed2));
    vector<uint> sh = m.get_shape();
    MDArray<double> sa(sh);
    for (uint i = 0; i < dim; i++) {
      sa.set(i, s[i]);
    }
    delete[] s;
    delete[] mean;
    delete[] cov;
    return sa;
  }


  // Return a sample, based on indicated parent values
  vector<double> GaussianDensities::sample(vector<double> & pv) {
    uint parent_value = (uint) pv.front();
    MDArray<double> & m = means.get_view(parent_value);
    MDArray<double> & c = covs.get_view(parent_value);

    MDArray<double> s = multivariate_normal(m, c);
    return s.get_values();
  }

  // Return likelihood, that is: P(child|parents)
  double GaussianDensities::get_lik(vector<double> & ptv, bool log) {
    uint ipv = (uint) ptv.front();
    MultiGauss & mgauss = mgauss_list[ipv];

    //ptv.erase(ptv.begin());
    vector<double> tv;
    for (uint i=1; i<ptv.size(); i++) {
      tv.push_back(ptv[i]);
    }

    double ll = mgauss.get_lik(tv, log);
    return ll;
  }

  // Return the distribtion's parameters
  vector< MDArray<double> > GaussianDensities::get_parameters() {
    vector<MDArray<double> > ret;
    ret.push_back(means);
    ret.push_back(covs);
    return ret;
  }

  MDArray<double> GaussianDensities::calc_means(vector<MDArray<double> > & ess) {
    //Calculate means from ESS.
    vector<uint> sh = ess[G_SY].get_shape();
    assert(sh.size() == 2);

    MDArray<double> m(sh);
    for (uint i = 0; i < sh[0]; i++) {
      for (uint j = 0; j < sh[1]; j++) {
	double v = ess[G_SY].get(i, j) / ess[G_S][i];
	m.set(i, j, v);
      }
    }
    return m;
  }

  MDArray<double> GaussianDensities::calc_cov_full(vector<MDArray<double> > & ess) {
    // Calculate full covariance matrices from ESS.
    MDArray<double> covs;	
    covs.set_shape(node_size, output_size, output_size);
      
    // Standard covariance implementation
    for (uint i=0; i<node_size; i++) {
      MDArray<double> d;
      MDArray<double> m = means.get_view(i);
      dot2(m, m, d);
      MDArray<double> v = ess[G_SYTY].get_view(i)/ess[G_S][i];
      MDArray<double> a = v-d;
      
      if (!useShrinkage) {
	covs.set(i, a);
      }
      else {
	vector< vector<double> > datapoints;
	for (uint j=G_ESSSIZE; j<ess.size(); j++) {
	  MDArray<double> point = ess[j];
	  if (point[0] == i) {
	    // Data point for this component
	    vector<double> data = point.get_slice(1, point.size());
	    datapoints.push_back(data);
	  }
	}

	// With shrinkage
	MDArray<double> md_data;
	if (!datapoints.empty()) {

	  //	  std::cout << "DATA: " << datapoints << std::endl;

	  md_data = toMDArray(datapoints);
	  
	  Shrinkage SHR(md_data,a);
	  double lambda=SHR.get_lambda_cor();  

	  MDArray<double> target(SHR.get_target_D());
	  MDArray<double> aux=SHR.combine(target,a,lambda);
	  covs.set(i,aux);
	  //	  cout << "old: " << a << endl;
	  //	  cout << "new: " << aux << endl;
	}
      }
    }
    return covs;
  }



  MDArray<double> GaussianDensities::calc_cov_full_tied(vector<MDArray<double> > & ess) {
    assert(cov_tied);
    assert(cov_type==FULL);

    uint dim = output_size;
    // Calculate full, tied covariance matrix from ESS.
    MDArray<double> a;
    a.set_shape(dim, dim);
    MDArray<double> b;
    b.set_shape(dim, dim);

    for (uint i = 0; i < node_size; i++) {
      MDArray<double> & m = means.get_view(i);
      MDArray<double> d;
      dot2(m.get_values(), m.get_values(), d);
      d.multiply_inplace(ess[G_S][i]);
      b.add_inplace(d);
      a.add_inplace(ess[G_SYTY].get_view(i));
    }

    MDArray<double> cov;
    cov.set_shape(dim, dim);
    a.sub_inplace(b);
    double sumv = sumvect(ess[G_S].get_values());
    a.multiply_inplace(1.0 / sumv);
    cov = a;

    MDArray<double> covs;
    covs.set_shape(node_size, dim, dim);

    for (uint i = 0; i < node_size; i++) {
      covs.set(i, cov);
    }

    return covs;
  }

  ostream& operator<<(ostream& output, const GaussianDensities& a)  {
    output << "Means: " << endl << a.means << endl;
    output << "Covariances: " << endl << a.covs << endl;
    return output;
  }


}
