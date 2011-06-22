/*
 * EMEngine.cpp
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

#include "emengine.h"
#include "mocapyexceptions.h"

using namespace std;

namespace mocapy {

EMEngine::EMEngine(DBN* new_dbn, MCMC* new_sampler,
		vector<Sequence> * new_seq_list,
		vector<MDArray<eMISMASK> > * new_mismask_list, vector<double> * new_weight_list) {

  	vector<Sequence*> seq_list2;
	if (new_seq_list != NULL) {
		for (uint i=0; i<new_seq_list->size(); i++)  {
			Sequence* s = new Sequence;
			s->copy( (*new_seq_list)[i] );
			seq_list2.push_back(s);
		}
	}

	initialize(new_dbn, new_sampler, &seq_list2, new_mismask_list, new_weight_list);
}

void EMEngine::initialize(DBN* new_dbn, MCMC* new_sampler,
		vector<Sequence*> * new_seq_list,
		vector<MDArray<eMISMASK> > * new_mismask_list, vector<double> * new_weight_list) {

	if (new_weight_list != NULL)
		weight_list = *new_weight_list;

	if (weight_list.empty() && new_seq_list != NULL) {
		weight_list.resize(new_seq_list->size(), 1.0);
	}

	dbn = new_dbn;
	sampler = new_sampler;

	if (new_seq_list != NULL)
		seq_list = *new_seq_list;

	if (new_mismask_list != NULL)
		mismask_list = *new_mismask_list;

	if (new_mismask_list != NULL && new_seq_list != NULL)
		inf_list = make_inf_list();

	nr_seqs = seq_list.size();
	// This flags if an E-step has been done (true) or not (false)
	E_done = false;
	keep_sequences = false;
	dataIsVerified=false;
}

bool EMEngine::load_sequences(const char* filename, char sep, vector<uint> columns) {
	vector<Sequence> data = data_loader(filename, sep, columns);
	seq_list.clear();
	for (uint i=0; i<data.size(); i++) {
		MDArray<double>* mda = new MDArray<double>;
		*mda = data[i];
		seq_list.push_back(mda);
	}

	nr_seqs = seq_list.size();

	if (weight_list.empty()) {
		weight_list.resize(seq_list.size(), 1.0);
	}

	inf_list = make_inf_list();
	cout << "Read " << nr_seqs << " sequences." << endl;
	return true;
}

bool EMEngine::load_mismask(const char* filename, char sep, vector<uint> columns) {
	if (!seq_list.empty()) {
		cout << "Load mismask and optional weights before loading sequences" << endl;
		return false;
	}

	vector<MDArray<double> > data = data_loader(filename, sep, columns);
	mismask_list = toMismask(data);
	return true;
}

bool EMEngine::load_weights(const char* filename) {
  // Loading weights of sequences. Each line in the file corresponds to the weight of a sequence.
	if (!seq_list.empty()) {
		cout << "Load mismask and optional weights before loading sequences" << endl;
		return false;
	}

	vector<MDArray<double> > data = data_loader(filename);
	assert(data.size()==1);
	vector<double*> pWeights = data[0].flat_iterator();
	weight_list.clear();
	for (uint i=0; i<pWeights.size(); i++) {
		weight_list.push_back(*pWeights[i]);
	}
	return true;
}


vector<InfEngineMCMC> EMEngine::make_inf_list() {
	vector<InfEngineMCMC> k;

	cout << "Sequences: " << seq_list.size() << endl;
	cout << "Mismasks: " << mismask_list.size() << endl;

	assert(seq_list.size() == mismask_list.size() || mismask_list.size() == 1);
	for (uint i = 0; i < seq_list.size(); i++) {
		Sequence * seq = seq_list[i];
		MDArray<eMISMASK> mismask;

		if (mismask_list.size() == 1) {
			mismask = mismask_list[0];
		} else {
			mismask = mismask_list[i];
		}


		double weight = weight_list[i];
		InfEngineMCMC inf_engine = InfEngineMCMC(dbn, sampler, seq, mismask,
				weight);
		k.push_back(inf_engine);
	}
	return k;
}

void EMEngine::do_E_step(uint mcmc_steps, uint burn_in_steps, bool init_random) {
	//	Perform the E step.
	//	mcmc_steps: number of Gibbs sampling steps
	//	burn_in_steps: number of Gibbs sampling burn-in steps
	//	init_random: if true, initialize hidden nodes to random values


  if (!dataIsVerified) {
    verifyData();
    dataIsVerified=true;
  }
    

	//	LogLik sum
	loglik = 0;

	//	Total number of sequence slices sampled
	slice_count = 0;
	for (uint i = 0; i < nr_seqs; i++) {
		// Pick a sequence, ie. a family object and a mismask array
		InfEngineMCMC & inf_engine = inf_list[i];
		// Do one sweep of Gibbs sampling over all nodes and whole seq
		inf_engine.initialize_sample_generator(burn_in_steps, init_random); //inf_engine.get_sample_generator(burn_in_steps, init_random);
		for (uint j = 0; j < mcmc_steps; j++) {
			Sample s = inf_engine.sample_next();
			loglik += s.ll;
			slice_count += s.slice_count;
		}
	}
	E_done = true;
}

  void EMEngine::verifyData() {
    // Check that each sequence has a datapoint for each node in a slice
    if (seq_list.empty()) {
      cout << "No data has been loaded" << endl;
      exit(1);
    }

    if (mismask_list.empty()) {
      cout << "No mismask has been loaded" << endl;
      exit(1);
    }
    
    for (uint i=0; i<seq_list.size(); i++) {
      Sequence* seq = seq_list[i];
      vector<uint> shape = seq->get_shape();
      assert(shape.size() == 2);
      
      int slice0=0;
      for (uint j=0; j<dbn->nodes_0.size(); j++) {
	slice0 += dbn->nodes_0[j]->get_output_size();
      }
      
      int sliceN=0;
      for (uint j=0; j<dbn->nodes_1.size(); j++) {
	sliceN += dbn->nodes_1[j]->get_output_size();
      }
      
      // Run through all slices
      for (uint s=0; s<shape[0]; s++) {
	MDArray<double> slice = seq->get_view(s);
	if ((i==0 && slice.size() != slice0) || (i>0 && slice.size() != sliceN)) {
	  cout << "Wrong number of data values in training sequence " << i+1 << " slice " << s+1 << endl;

	  cout << "Expected ";
	  if (i==0)
	    cout << slice0;
	  else
	    cout << sliceN;

	  cout << " data values, got " << slice.size() << " data values" << endl;
	  
	  exit(1);
	}
      }
    } 
  
    // Check that each slice has a valid mismask
    for (uint i=0; i<mismask_list.size(); i++) {
      MDArray<eMISMASK> mm = mismask_list[i];
      vector<uint> shape = mm.get_shape();
      assert(shape.size() == 2);
      
      int slice0 = dbn->nodes_0.size();
      int sliceN = dbn->nodes_1.size();
      
      // Run through all slices
      for (uint s=0; s<shape[0]; s++) {
	MDArray<eMISMASK> slice = mm.get_view(s);
	if ((i==0 && slice.size() != slice0) || (i>0 && slice.size() != sliceN)) {
	  cout << "Wrong number of mismask values for training sequence " << i+1 << " slice " << s+1 << endl;

	  cout << "Expected ";
	  if (i==0)
	    cout << slice0;
	  else
	    cout << sliceN;

	  cout << " mismask values, got " << slice.size() << " mismask values" << endl;
	  
	  exit(1);
	}
      }
    }
  }


double EMEngine::get_loglik() {
	//	Return LogLik per slice during E-step.
	if (!E_done) {
		throw MocapyExceptions("E step was not performed.");
	}

	// Return LogLik/slice count
	return loglik / (double) slice_count;
}

void EMEngine::do_M_step() {
	//	Perform M-step.
	if (!E_done) {
		throw MocapyExceptions("ESS step was not performed.");
	}

	for (uint n = 0; n < dbn->nodes_1.size(); n++) {
		Node* node_0 = dbn->nodes_0[n];
		Node* node_1 = dbn->nodes_1[n];
		vector<Node*> node_list;

		if (node_0 == node_1) {
			// Tied nodes
			assert(node_0->slice == TIED);
			node_list.push_back(node_0);
		} else {
			// Untied nodes
			node_list.push_back(node_0);
			node_list.push_back(node_1);
		}
		for (uint i = 0; i < node_list.size(); i++) {
			Node* node = node_list[i];
			if (!node->fixed) {
				// Collect ESS for one node
				vector<MDArray<double> > local_ess = node->get_ess();
				node->do_M_step(local_ess);
			}
		}
	}
	E_done = false;
}

}
