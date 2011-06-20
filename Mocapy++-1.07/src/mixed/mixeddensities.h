/*
 * mixeddensities.h
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

#ifndef MIXEDDENSITIES_H_
#define MIXEDDENSITIES_H_

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <fstream>

#include "../framework/densitiesbase.h"
#include "mixedess.h"

namespace mocapy {

class MixedDensities;
std::ostream& operator<<(std::ostream&, const MixedDensities&);

class MixedDensities : public DensitiesBase {
public:
	MixedDensities();                                                                       

    MixedDensities(uint new_node_size, 
                   bool new_init_random = false,
                   MDArray<double> user_means = MDArray<double> (), 
                   MDArray<double> user_variance = MDArray<double> ());
                   
    MixedDensities(uint new_node_size, 
                   CPD & new_user_cpd, 
                   MDArray<double> user_means = MDArray<double> (), 
                   MDArray<double> user_variance = MDArray<double> ()); 
	virtual ~MixedDensities() {};

	// Initializes the Density arrays
	void construct(std::vector<uint> & parent_sizes);

	// Parameter estimation based on the ESS
	void estimate(std::vector<MDArray<double> > & ess);

	// Return a sample, based on indicated parent values
	std::vector<double> sample(std::vector<double> & ptv);

	// Return likelihood, that is: P(child|parents)
	double get_lik(std::vector<double> & ptv, bool log=false);

	// Return the distribution's parameters
	std::vector< MDArray<double> > get_parameters();

	void set_user_cpd(CPD & userCPD);
	void set_prior(Prior * new_prior);
	CPD & getCPD() {return cpd;}

	// Persistence
	friend class boost::serialization::access;
	template<class Archive>
	void serialize(Archive & ar, const unsigned int version);

	friend std::ostream& operator<< (std::ostream& output, const MixedDensities& a);

private:
	void test_cpd(CPD & cpd, std::vector<uint> & shape);
	void set_cpd(CPD & cpd);
	void initialize();

	CPD make_uniform_cpd(const std::vector<uint> & shape);
	CPD make_random_cpd(std::vector<uint> & shape, bool no_zeroes = true);

    double sample_1d_gauss(std::vector<double> & pv);
    double sample_discrete(std::vector<double> & pv);
    
	Prior* prior;
	bool init_random;

    // MDArray<double> user_means;
    // MDArray<double> user_variance;
	MDArray<double> means;
    MDArray<double> variance;
    
    uint parent_size;
 
	// CPDs
	CPD cpd;
	CPD user_cpd;
	CPD cum_cpd;
	CPD log_cpd;
	std::vector<uint> CPD_shape;
};

template<class Archive>
void MixedDensities::serialize(Archive & ar, const unsigned int version) {
    ar & boost::serialization::base_object<DensitiesBase>(*this);
    ar & init_random;

    ar & means;                
    ar & variance;

    ar & cpd;
    ar & user_cpd;
    ar & cum_cpd;
    ar & log_cpd;
    ar & CPD_shape;
}

}

#endif /* MIXEDDENSITIES_H_ */
