/*
 * mlr-uni.cpp
 *
 *  Copyright (C) 2010, Christian Andreetta, Thomas Hamelryck,
 *    The Bioinformatics Centre, University of Copenhagen.
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
//*----------------------------------------------------------------
//* GENERAL OVERVIEW
//* 
//* Implements multiniomial logistic model.
//* 
//* This is an empirical Bayes model. The parameters are estimated using S-EM.
//* The elements in the Beta matrix are modelled by Gaussians. The means/variances of
//* these are estimated via S-EM. E-step is done using Metropolis-Hastings sampling. 
//* 
//* The mean can be used as a emirical Bayes MAP estimate.
//*----------------------------------------------------------------
//* Python version: Thomas Hamelryck <thamelry -at- binf.ku.dk>
//* C++ version: Christian Andreetta <chrandr -at- binf.ku.dk>
//*----------------------------------------------------------------
//#define BOOST_DISABLE_ASSERTS    // after debug, for performance increase using arrays
#include <boost/algorithm/string.hpp>
#include <boost/program_options.hpp>
#include <boost/format.hpp>
#include <boost/thread.hpp>
#include <boost/bind.hpp>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cstdlib>
#include <mocapy.h>
#define FLOAT_TYPE float

using namespace std;

//#define FLOAT_TYPE double
//*----------------------------------------------------------------
//* GENERAL OVERVIEW
//* 
//* Implements multiniomial logistic model.
//* 
//* This is an empirical Bayes model. The parameters are estimated using S-EM.
//* The elements in the Beta matrix are modelled by Gaussians. The means/variances of
//* these are estimated via S-EM. E-step is done using Metropolis-Hastings sampling. 
//* 
//* The mean can be used as a emirical Bayes MAP estimate.
//* 
//* Seems to work just fine.
//*----------------------------------------------------------------
unsigned int _DEBUG=0;
std::string data_separator_field=",\t ";
std::string data_separator_record=";\n";
//*----------------------------------------------------------------
//* PROGRAM OPTIONS
//*----------------------------------------------------------------
namespace po = boost::program_options;
// parser
po::variables_map getCmdOptions(int ac, char* av[]) {
    // COMMAND LINE PARSER
    po::variables_map options;
    // definition
    try {
        po::options_description desc("Allowed options");
        desc.add_options()
            // HELP
            ("help",  "help message")
            // GENERAL
            ("debug",    po::value<unsigned int>()->default_value(0), "debug level (150: full detail)")
            ("threads",  po::value<unsigned int>()->default_value(1), "number of threads used in the proposal phase")
            // RANDOM NUMBER GENERATOR
            ("seed",  po::value<unsigned int>()->default_value(12345), "random generator: seed")
            ("generator-lagged-fibonacci607",  po::value<bool>()->default_value(true), "random generator: use Lagged Fibonacci607")
            // INPUT
            ("data-input-filename",    po::value<std::string>()->default_value("./mlr.dat"), "data I/O: input file. format= 'state;input_vector'")
            ("data-separator-field",   po::value<std::string>()->default_value(data_separator_field), "data I/O: separators between record fields")
            ("data-separator-record",  po::value<std::string>()->default_value(data_separator_record), "data I/O: separators between records")
            // ML
            ("ml-proposal-init-means",     po::value<FLOAT_TYPE>()->default_value(0.0), "input: ML proposal: initial means value")
            ("ml-proposal-init-std_devs",  po::value<FLOAT_TYPE>()->default_value(1.0), "input: ML proposal: initial std_devs value")
            ("ml-proposal-filename-input",  po::value<std::string>()->default_value(""), "input: ML proposal file with format: 'mean_vector[i];std_dev_vector[i]'")
            ("ml-proposal-filename-output", po::value<std::string>()->default_value(""), "output: ML proposal file with format: 'mean_vector[i];std_dev_vector[i]'")
            // EM DETAILS
            ("em-iterations-max",                  po::value<unsigned int>()->default_value(1e5), "EM: max number of iterations")
            ("em-iterations-without-improvements", po::value<unsigned int>()->default_value(1e3), "EM: max number of iterations without LL improvements")
            ("em-progress-steps",                  po::value<unsigned int>()->default_value(20), "EM: show progress and dump output files every N iterations")
            ("em-step-gibbs-iterations", po::value<unsigned int>()->default_value(7), "EM step: Gibbs iterations")
            ("em-step-resample-perc",    po::value<FLOAT_TYPE>()->default_value(0.3), "EM step: probability of resampling an element in the beta matrices")
        ;
        po::store(po::command_line_parser(ac, av).options(desc).run(), options);
        po::notify(options);
        // parse options
        if (options.count("help") or not(options.count("data-input-filename") ) ) {
            printf("# USAGE: <program> [options] --data-input-filename <input.dat>\n");
            std::cout << desc << "\n";
            exit(0);
        }
        // fill global variables
        _DEBUG=options["debug"].as<unsigned int>();
        data_separator_field=options["data-separator-field"].as<std::string>();
        data_separator_record=options["data-separator-record"].as<std::string>();
    } catch(std::exception& e) {
        std::cerr << "error: " << e.what() << "\n";
        exit(1);
    } catch(...) {
        std::cerr << "Exception of unknown type!\n";
        exit(1);
    }
    return options;
}
void printOptions(po::variables_map &options) {
    printf("# PARAMETERS values:\n");
    for ( po::variables_map::iterator vm_it=options.begin() ; vm_it != options.end(); vm_it++ ) {
        std::cout << "# PARAM: " << (*vm_it).first << " => ";
        const po::variable_value& v = (*vm_it).second;
        if ( ! v.empty() ) {
            const std::type_info& type = v.value().type() ;
            if ( type == typeid( std::string ) ) { std::cout <<"'"<<v.as<std::string>()<<"'"; }
            else if ( type == typeid( bool ) ) { std::cout << v.as<bool>(); }
            else if ( type == typeid( int ) ) { std::cout << v.as<int>(); }
            else if ( type == typeid( unsigned int ) ) { std::cout << v.as<unsigned int>(); }
            else if ( type == typeid( float ) ) { std::cout << v.as<float>(); }
            else if ( type == typeid( double ) ) { std::cout << v.as<double>(); }
            else if ( type == typeid( FLOAT_TYPE ) ) { std::cout << v.as<FLOAT_TYPE>(); }
            else {
                printf("unmanaged container...");
            }
        }
        printf("\n");
    }
}
//*----------------------------------------------------------------
//* MDArray initializations
//*----------------------------------------------------------------
template <class T>
mocapy::MDArray<T> vec_2_MDArray( const std::vector< T > &vec ) {
    mocapy::MDArray<T> mda;
    mda.set_shape( vec.size() );
    mda.set_values(vec);
    return mda;
}
//*----------------------------------------------------------------
template <class T>
mocapy::MDArray<T> vecVec_2_MDArray( const std::vector< std::vector< T > > &vec, bool add_one_flag=false ) {
    mocapy::MDArray<T> mda;
    // shape
    unsigned int add_one_flag_offset=0;
    if (add_one_flag==true) { add_one_flag_offset=1; }
    //
    std::vector<unsigned int> shape;
    shape.push_back( vec.size() );
    shape.push_back( add_one_flag_offset+vec[0].size() );
    mda.set_shape( shape );
    // fill values
    for (unsigned int j=0; j<vec.size(); j++) {
        if (add_one_flag==true) { mda.set(j,0,1.0); }
        for (unsigned int k=0; k<vec[0].size(); k++) {
            mda.set( j, add_one_flag_offset+k, vec[j][k] );
        }
    }
    //
    return mda;
}
//*----------------------------------------------------------------
//* Returns array of probabilities over all output states.
//*
//* size: number of output states (e. the multinomial dimension)
//* dim: dimension of input vector
//* beta: (size-1) x (dim+1) matrix, the model's parameters 
//* add_one: flag to add a "1" to the start of the input vector
void mldens(
        mocapy::MDArray<FLOAT_TYPE> &p_states_mda,
        mocapy::MDArray<FLOAT_TYPE> &invec_mda_elem,
        mocapy::MDArray<FLOAT_TYPE> &betas_mda_elem2D,
        const unsigned int &dim,
        const unsigned int &size,
        const FLOAT_TYPE tiny_probability=1e-5) {
    // loops over all states
    p_states_mda[0]=1;
    for (unsigned int j=1; j<size; j++) {
        p_states_mda[j]=betas_mda_elem2D.get_view(j-1).dot( invec_mda_elem );
        if (_DEBUG>=100) {
            printf("# mldens: dot product: p_states_mda[%d]: %.3f\n# betas_elem: %s# invec_elem: %s",
                j,p_states_mda[j],
                betas_mda_elem2D.get_view(j-1).tostring().c_str(),
                invec_mda_elem.tostring().c_str()
            );
        }
    }
    if (_DEBUG>=90) { printf("# mldens pre-norm(): p_states_mda: %s",p_states_mda.tostring().c_str()); }
    // normalization of probabilities
    // avoids overflow, subtracts largest value in log space
    p_states_mda.add_inplace(-1*p_states_mda.get_max());
    p_states_mda.exp_all();
    p_states_mda.normalize();
    // avoids underflow
    p_states_mda.add_inplace(tiny_probability);
    // store normalized probabilities
    p_states_mda.normalize();
    if (_DEBUG>=90) { printf("# mldens post-norm(): p_states_mda: %s",p_states_mda.tostring().c_str()); }
}
//*----------------------------------------------------------------
//* DATAFILES
//*----------------------------------------------------------------
std::vector< std::vector<std::string> > datafile_parse(
        const std::string &data_input_filename,
        const std::string sep_record=";"
        ) {
    if (_DEBUG>=10) { printf("# DEBUG: entering function: datafile_parse\n"); }
    std::vector< std::vector<std::string> > parse_vec;
    std::string line;
    unsigned int line_count=0;
    // open file
    std::ifstream fileH(data_input_filename.c_str());
    if (!fileH.is_open()) {
        printf("# ERROR: Unable to open file: '%s'\n",data_input_filename.c_str());
        exit(1);
    }
    // parse file
    while (! fileH.eof() ) {
        getline (fileH,line);
        std::vector<std::string> remove_comment;
        // remove comments
        boost::split( remove_comment, line, boost::is_any_of("#") );
        // process line
        if (remove_comment[0].length()>0) {
            if (_DEBUG>=100) {
                printf("# line %3d: '%s' -> '%s'\n",
                    line_count,line.c_str(),remove_comment[0].c_str());
            }
            std::vector<std::string> split_records;
            boost::split( split_records, remove_comment[0], boost::is_any_of(sep_record.c_str()) );
            parse_vec.push_back(split_records);
        }
        line_count++;
    }
    fileH.close();
    //
    return parse_vec;
}
void input_datafile_load(
        const std::string &data_input_filename,
        std::vector<unsigned int> &state_vec,
        std::vector< std::vector<FLOAT_TYPE> > &invec_vec,
        const std::string sep_record=data_separator_record,
        const std::string sep_fields=data_separator_field
        ) {
    if (_DEBUG>=10) { printf("# DEBUG: entering function: input_datafile_load\n"); }
    // parse input file
    std::vector< std::vector<std::string> > parse_vec = datafile_parse(
        data_input_filename, sep_record );
    // parse outer records vector
    unsigned int state=0;
    double val=0;
    for (unsigned int j=0; j<parse_vec.size(); j++) {
        // state
        state=atoi(parse_vec[j][0].c_str());
        state_vec.push_back( state );
        // input vector
        std::vector<std::string> fields_rec;
        std::vector<FLOAT_TYPE> val_vec;
        boost::split( fields_rec, parse_vec[j][1], boost::is_any_of(sep_fields.c_str()) );
        for (unsigned int j=0; j<fields_rec.size(); j++) {
            val=atof(fields_rec[j].c_str());
            val_vec.push_back(val);
            if (_DEBUG>=100) { printf("%.1f, ", val ); }
        }
        invec_vec.push_back(val_vec);
        if (_DEBUG>=100) { printf("\n"); }
    }
}
//*----------------------------------------------------------------
template <class T>
void print_vec( const std::vector< T > &vec, bool print_index_flag=false) {
    if (print_index_flag) { printf("# "); }
    for (unsigned int j=0; j< vec.size(); j++) {
        if (print_index_flag) { printf("%d: ",j); }
        //printf("%.2f",vec[j]);
        std::cout << vec.at(j);
        if (j<vec.size()-1) { printf(", "); }
    }
    printf("\n");
}
template <class T>
void print_vec_vec( const std::vector< std::vector<T> > &vec) {
    for (unsigned int j=0; j< vec.size(); j++) {
        printf("# %d: ",j);
        print_vec( vec[j] );
    }
}
//*----------------------------------------------------------------
//* Multinomial Logistics - univariate
//*----------------------------------------------------------------
/*
 *  Estimate the parameters of a multinomial logistic model using EM
 *  The elements of the beta matrix are described by univariate Gaussians distributions
 *  Point estimate: the matrix with the means of the Gaussians
 */
class MLR_uni {
    private:
        // data base size
        unsigned int dim;  // size of a single input vector
        unsigned int size; // size of output
        unsigned int N;    // number of input data
        unsigned int threads_number; // number of threads used in the proposal phase
        std::vector< pair<unsigned int,unsigned int> > threads_input_intervals;
        std::vector<FLOAT_TYPE> ll_vec;
        // random generator: use Lagged Fibonacci 607?
        bool use_lagged_fibonacci607;
        // dimensions of means, std_dev, betas matrices
        unsigned int means_shape0;
        unsigned int means_shape1;
        // default values
        FLOAT_TYPE means_default;
        FLOAT_TYPE std_devs_default;
        // input
        mocapy::MDArray<unsigned int> states_mda;
        mocapy::MDArray<FLOAT_TYPE> states_stats_mda;
        mocapy::MDArray<FLOAT_TYPE> invec_mda;
        // point estimate: mean and standard deviation
        mocapy::MDArray<FLOAT_TYPE> means_mda;
        mocapy::MDArray<FLOAT_TYPE> std_devs_mda;
        // betas matrix: keeps the sampled beta matrices: N x (size-1) x (dim-1)
        mocapy::MDArray<FLOAT_TYPE> betas_mda;
        // beta_tmp matrix: keeps the proposal beta matrix: threads_number x (size-1) x (dim-1)
        // cache_tmp matrix: as beta_tmp. used to compute the threaded standard deviation
        mocapy::MDArray<FLOAT_TYPE> beta_tmp_mda;
        mocapy::MDArray<FLOAT_TYPE> cache_tmp_mda;
        // probabilities P(state|input_vector,beta)
        mocapy::MDArray<FLOAT_TYPE> p_states_accepted_mda;
        mocapy::MDArray<FLOAT_TYPE> p_states_stats_mda;
        // probabilities over the state: used as a temporary storage. specific for each thread
        mocapy::MDArray<FLOAT_TYPE> p_states_mda;
        // random generator
        std::vector< mocapy::RandomGen* > randomGen_vec;
        // functions
        void matrices_reset(void);
        void init(void);
        FLOAT_TYPE get_rand_uniform(
            const unsigned int thread=0
        );
        FLOAT_TYPE get_rand_normal(
            const FLOAT_TYPE mu=0.0,
            const FLOAT_TYPE sigma=1.0,
            const unsigned int thread=0
        );
        void em_step_do_parallel(
            const unsigned int input_idx_start,
            const unsigned int input_idx_end,
            const unsigned int thread,
            const unsigned int gibbs_iter_num,
            const FLOAT_TYPE resample_perc,
            const bool em_init_flag
        );
        void resample_beta_slice(
            const unsigned int i,
            const unsigned int gibbs_iter_num,
            const FLOAT_TYPE resample_perc,
            const unsigned int thread=0
        );
        FLOAT_TYPE sample_b_elem(
            const unsigned int j, const unsigned int k,
            const unsigned int thread=0
        );
        void update_parameters_parallel_mean(
            const unsigned int input_idx_start,
            const unsigned int input_idx_end,
            mocapy::MDArray<FLOAT_TYPE> *beta_tmp_mda_elem
        );
        void update_parameters_parallel_std_dev(
            const unsigned int input_idx_start,
            const unsigned int input_idx_end,
            mocapy::MDArray<FLOAT_TYPE> *beta_tmp_mda_elem,
            mocapy::MDArray<FLOAT_TYPE> *cache_tmp_mda_elem
        );
        void update_parameters(void);
        void get_log_lik_parallel(
            FLOAT_TYPE *ll, 
            const unsigned int input_idx_start,
            const unsigned int input_idx_end,
            mocapy::MDArray<FLOAT_TYPE> &p_states_mda_elem
        );
    public:
        unsigned int get_N(void)    {return N;}
        unsigned int get_dim(void)  {return dim;}
        unsigned int get_size(void) {return size;}
        // constructor, destructor
        MLR_uni(
            const unsigned int threads_number=1,
            const unsigned int seed=0,
            const bool use_lagged_fibonacci607=true,
            const FLOAT_TYPE means_default=0.0,
            const FLOAT_TYPE std_devs_default=1.0
        );
        ~MLR_uni();
        // functions
        void set_input_from_datafile(const std::string &filename);
        void set_input_from_vectors(
            const std::vector<unsigned int> &state_vec,
            const std::vector< std::vector<FLOAT_TYPE> > &invec_vec,
            bool add_one_flag=false
        );
        void set_ML_proposal(
            const std::string &filename,
            const std::string sep_record=data_separator_record,
            const std::string sep_fields=data_separator_field
        );
        void save_ML_proposal(
            const std::string &filename,
            const unsigned int precision=8,
            const std::string sep_record=data_separator_record.substr(0,1),
            const std::string sep_fields=data_separator_field.substr(0,1)
        );
        void em_do(
            const unsigned int em_iterations_max=1e3,
            const unsigned int em_iterations_without_improvements=500,
            const unsigned int em_progress_steps=20,
            const unsigned int em_step_gibbs_iterations=5,
            const FLOAT_TYPE em_step_resample_perc=0.2,
            const std::string ML_filename_out="",
            bool init_flag=true
        );
        void em_step_do(
            const unsigned int gibbs_iter_num,
            const FLOAT_TYPE resample_perc,
            const bool em_init_flag=false
        );
        FLOAT_TYPE calc_p(
            const unsigned int &state,
            mocapy::MDArray<FLOAT_TYPE> &p_states_mda_elem,
            mocapy::MDArray<FLOAT_TYPE> &invec_mda_elem,
            mocapy::MDArray<FLOAT_TYPE> &betas_mda_elem2D
        );
        std::string evaluate_classification(bool betas_efficiency_check=false);
        FLOAT_TYPE get_log_lik(void);
        std::string get_ML(void);
};
//*----------------------------------------------------------------
MLR_uni::MLR_uni(
        const unsigned int threads_number,
        const unsigned int seed,
        const bool use_lagged_fibonacci607,
        const FLOAT_TYPE means_default,
        const FLOAT_TYPE std_devs_default
        ) {
    if (_DEBUG>=10) { printf("# DEBUG: entering function: MLR_uni constructor\n"); }
    // default values
    this->means_default=means_default;
    this->std_devs_default=std_devs_default;
    this->threads_number=threads_number;
    this->use_lagged_fibonacci607=use_lagged_fibonacci607;
    // random number generator
    randomGen_vec.resize(threads_number);
    for (unsigned int thread=0; thread<threads_number; thread++) {
        randomGen_vec[thread] = new mocapy::RandomGen();
        if (seed>0) { randomGen_vec[thread]->mocapy_seed(seed+thread); }
    }
    if (_DEBUG>=100) {
        printf("# DEBUG: TEST: random number generator: normal distribution: mu=100, sigma=20\n");
        for (unsigned int j=0; j<100;j++) { printf( "%.2f ", get_rand_normal( 100,20) ); }
        printf("\n");
    }
}
MLR_uni::~MLR_uni() {
    for (unsigned int thread=0; thread<threads_number; thread++) {
        delete randomGen_vec[thread];
    }
}
//*----------------------------------------------------------------
void MLR_uni::matrices_reset(void) {
    if (_DEBUG>=10) { printf("# DEBUG: entering function: matrices_reset\n"); }
    // constants initialization
    N=states_mda.get_shape()[0];
    size=states_mda.get_max()+1;
    dim=invec_mda.get_shape()[1]-1;
    // sanity checks
    assert(N>0);
    assert(size>0);
    assert(size==dim);
    assert(N==invec_mda.get_shape()[0]);
    assert(threads_number<N);
    // threads: intervals definition
    unsigned int idx_start=0;
    unsigned int idx_end=0;
    unsigned int input_portion_per_thread=N/threads_number;
    for (unsigned int thread=0; thread<threads_number; thread++) {
        idx_start=input_portion_per_thread * thread;
        if (thread<threads_number-1) {
            idx_end=input_portion_per_thread * (thread+1);
        } else { idx_end=N; }
        threads_input_intervals.push_back( std::pair<unsigned int,unsigned int>(idx_start,idx_end) );
        if (_DEBUG>=50) { printf("# thread %d: start at input: %u, end: %u\n",thread,idx_start,idx_end); }
    }
    
    // ML data base matrices
    means_shape0=size-1;
    means_shape1=dim+1;
    // means: start equal to 0
    means_mda.set_shape( means_shape0, means_shape1 );
    means_mda.add_inplace(means_default);
    // standard deviations: start equal to 1
    std_devs_mda.set_shape( means_shape0, means_shape1 );
    std_devs_mda.add_inplace(std_devs_default);
    // probabilities over the states
    p_states_mda.set_shape( threads_number,size );
    p_states_stats_mda.set_shape( size );
    // beta matrix: keeps the sampled beta matrices: N x (size-1) x (dim-1)
    betas_mda.set_shape( N, means_shape0, means_shape1 );
    beta_tmp_mda.set_shape( threads_number, means_shape0, means_shape1 );
    cache_tmp_mda.set_shape( threads_number, means_shape0, means_shape1 );
    // probabilities P(state|input_vector,beta)
    p_states_accepted_mda.set_shape( N );
    // log likelihood: computed thread-wise
    ll_vec.resize(threads_number);
    // output
    if (_DEBUG>=50) {
        printf("# states_mda: shape: (%d). dimensions: %d\n",states_mda.get_shape()[0],(int)states_mda.get_shape().size());
        if (_DEBUG>=100) { printf("# states_mda: %s",states_mda.tostring().c_str()); }
        printf("# invec_mda: shape: (%d,%d). dimensions: %d\n",invec_mda.get_shape()[0],invec_mda.get_shape()[1],(int)invec_mda.get_shape().size());
        if (_DEBUG>=100) { printf("# invec_mda:\n%s",invec_mda.tostring().c_str()); }
        printf("# means_mda: shape: (%d,%d). dimensions: %d\n",means_mda.get_shape()[0],means_mda.get_shape()[1],(int)means_mda.get_shape().size());
        if (_DEBUG>=100) { printf("# means_mda:\n%s",means_mda.tostring().c_str()); }
        printf("# std_devs_mda: shape: (%d,%d). dimensions: %d\n",std_devs_mda.get_shape()[0],std_devs_mda.get_shape()[1],(int)std_devs_mda.get_shape().size());
        if (_DEBUG>=100) { printf("# std_devs_mda:\n%s",std_devs_mda.tostring().c_str()); }
        printf("# betas_mda: shape: (%d,%d,%d). dimensions: %d\n",
            betas_mda.get_shape()[0],betas_mda.get_shape()[1],betas_mda.get_shape()[2],
            (int)betas_mda.get_shape().size()
            );
        printf("# p_states_accepted_mda: shape: (%d). dimensions: %d\n",
            p_states_accepted_mda.get_shape()[0],(int)p_states_accepted_mda.get_shape().size()
            );
        if (_DEBUG>=100) { printf("# p_states_accepted_mda:\n%s",p_states_accepted_mda.tostring().c_str()); }
    }
}
void MLR_uni::init(void) {
    if (_DEBUG>=10) { printf("# DEBUG: entering function: init\n"); }
    // initialization of beta matrices
    // loops over the input vectors
    for (unsigned int i=0; i<N ; i++) {
        // loops over the beta matrix
        for (unsigned int j=0; j<means_shape0; j++) {
            for (unsigned int k=0; k<means_shape1; k++) {
                betas_mda.set( i,j,k, sample_b_elem(j,k) );
            }
        }
        // kickstart P(state|input_vector,beta_matrix_slice)
        p_states_accepted_mda[i]=calc_p(
            states_mda[i],
            p_states_mda.get_view(0),
            invec_mda.get_view(i),
            betas_mda.get_view(i)
        );
    }
    if (_DEBUG>=120) { printf("# betas_mda:\n%s",betas_mda.tostring().c_str()); }
    if (_DEBUG>=100) {
        printf("# init: p_states_accepted_mda:\n%s",p_states_accepted_mda.tostring().c_str());
        printf("# init: states_reference_mda:\n%s",states_mda.tostring().c_str());
    }
}
//*----------------------------------------------------------------
void MLR_uni::set_input_from_datafile(const std::string &filename) {
    if (_DEBUG>=10) { printf("# DEBUG: entering function: set_input\n"); }
    std::vector<unsigned int> state_vec;
    std::vector< std::vector<FLOAT_TYPE> > invec_vec;
    // load input
    input_datafile_load( filename, state_vec, invec_vec );
    // convert input to MDArrays
    set_input_from_vectors( state_vec,invec_vec );
    if (_DEBUG>=100) {
        printf("# states: "); print_vec( state_vec );
        printf("# inputs:\n"); print_vec_vec( invec_vec );
    }
}
//*----------------------------------------------------------------
void MLR_uni::set_input_from_vectors(
        const std::vector<unsigned int> &state_vec,
        const std::vector< std::vector<FLOAT_TYPE> > &invec_vec,
        bool add_one_flag
        ) {
    if (_DEBUG>=10) { printf("# DEBUG: entering function: MLR_uni constructor\n"); }
    // MDArray initialization
    states_mda=vec_2_MDArray( state_vec );
    // global variable size
    size=states_mda.get_max()+1;
    // input data: statistics
    states_stats_mda.set_shape(size);
    for (unsigned int j=0; j<state_vec.size(); j++) {
        states_stats_mda[ state_vec[j] ]++;
    }
    states_stats_mda.normalize();
    // do input size and number of states match?
    if ( ( add_one_flag==false) and ( size==invec_vec[0].size() ) ) {
        add_one_flag=true;
    }
    invec_mda=vecVec_2_MDArray( invec_vec, add_one_flag );
    // initialize ML storage matrices
    matrices_reset();
    //
    printf("# input data loaded: model is: N: %d, size: %d, dim: %d\n",
            get_N(),get_size(),get_dim() );
}
//*----------------------------------------------------------------
void MLR_uni::set_ML_proposal(
        const std::string &filename,
        const std::string sep_record,
        const std::string sep_fields
        ) {
    // TODO: vec<vec<str>> -> vec<vec<float>>
    if (_DEBUG>=10) { printf("# DEBUG: entering function: set_ML_proposal\n"); }
    // parse file
    std::vector< std::vector<std::string> > parse_vec =
        datafile_parse(filename);
    // parse records
    double val=0;
    // means
    std::vector< std::vector<FLOAT_TYPE> > means_vecVec;
    for (unsigned int j=0; j<parse_vec.size(); j++) {
        // input vector
        std::vector<std::string> fields_rec;
        std::vector<FLOAT_TYPE> val_vec;
        boost::split( fields_rec, parse_vec[j][0], boost::is_any_of(sep_fields.c_str()) );
        for (unsigned int j=0; j<fields_rec.size(); j++) {
            val=atof(fields_rec[j].c_str());
            val_vec.push_back(val);
            if (_DEBUG>=100) { printf("%.1f, ", val ); }
        }
        means_vecVec.push_back(val_vec);
        if (_DEBUG>=100) { printf("\n"); }
    }
    // means: mda
    means_mda=vecVec_2_MDArray( means_vecVec );
    // standard deviations
    std::vector< std::vector<FLOAT_TYPE> > std_devs_vecVec;
    for (unsigned int j=0; j<parse_vec.size(); j++) {
        // input vector
        std::vector<std::string> fields_rec;
        std::vector<FLOAT_TYPE> val_vec;
        boost::split( fields_rec, parse_vec[j][1], boost::is_any_of(sep_fields.c_str()) );
        for (unsigned int j=0; j<fields_rec.size(); j++) {
            val=atof(fields_rec[j].c_str());
            val_vec.push_back(val);
            if (_DEBUG>=100) { printf("%.1f, ", val ); }
        }
        std_devs_vecVec.push_back(val_vec);
        if (_DEBUG>=100) { printf("\n"); }
    }
    // std_devs: mda
    std_devs_mda=vecVec_2_MDArray( std_devs_vecVec );
    // initialize
    init();
    std::string eval_init=evaluate_classification();
    printf("# init: set_ML_proposal: classification:\n%s\n",eval_init.c_str()); 
}
//*----------------------------------------------------------------
// returns P(state|input_vector,beta_matrix_slice)
//     note: get_view() in mldens() expect NON-const arguments
FLOAT_TYPE MLR_uni::calc_p(
        const unsigned int &state,
        mocapy::MDArray<FLOAT_TYPE> &p_states_mda_elem,
        mocapy::MDArray<FLOAT_TYPE> &invec_mda_elem,
        mocapy::MDArray<FLOAT_TYPE> &betas_mda_elem2D ) {
    // compute probabilities over all states, stores into p_states_mda_elem
    if (_DEBUG>=20) { printf("# DEBUG: entering function: calc_p\n"); }
    mldens(p_states_mda_elem, invec_mda_elem, betas_mda_elem2D, this->dim, this->size);
    if (_DEBUG>=50) {
        printf("# calc_p: state %d -> %.3f\n", state, p_states_mda[state]);
        if (_DEBUG>=100) {
            printf("# calc_p: p_states_mda_elem: shape: (%d). dimensions: %d\n",
                p_states_mda_elem.get_shape()[0],(int)p_states_mda_elem.get_shape().size() );
            printf("# calc_p: p_states_mda_elem: %s",p_states_mda_elem.tostring().c_str());
        }
    }
    // return relevant probability
    return p_states_mda_elem[state];
}
//*----------------------------------------------------------------
std::string MLR_uni::evaluate_classification(bool betas_efficiency_check) {    
    if (_DEBUG>=10) { printf("# DEBUG: entering function: evaluate_classification\n"); }
    // clear statistics
    mocapy::MDArray<FLOAT_TYPE> predictions_count;
    predictions_count.set_shape( 2 );
    p_states_stats_mda.multiply_inplace(0.0);
    // loops over input data, evaluate ML: beta=means
    mocapy::MDArray<FLOAT_TYPE> p_states_mda_elem=p_states_mda.get_view(0);
    for (unsigned int i=0; i<N; i++) {
        if (betas_efficiency_check==true) {
            mldens(p_states_mda_elem, invec_mda.get_view(i), betas_mda.get_view(i), dim, size);
        } else {
            mldens(p_states_mda_elem, invec_mda.get_view(i), means_mda, dim, size);
        }
        // get most probable predicted state
        unsigned int argmax_idx=0;
        FLOAT_TYPE argmax_val=0;
        for (unsigned int j=0; j<size; j++) {
            if (argmax_val<p_states_mda_elem[j]) {
                argmax_val=p_states_mda_elem[j];
                argmax_idx=j;
            }
        }
        p_states_stats_mda[argmax_idx]++;
        // correct?
        if (argmax_idx == states_mda[i]) { predictions_count[0]++; }
        else { predictions_count[1]++; }
        if (_DEBUG>=80) {
            printf("# evaluate_classification: input %3d: state_input: %2d, proposed: %d\n",
                i,states_mda[i],argmax_idx
            );
            printf("# states distrib. prediction: %s",p_states_mda_elem.tostring().c_str());
            printf("# states distrib. from input: %s",states_stats_mda.tostring().c_str());
            printf("# predictions_count: %s",predictions_count.tostring().c_str());
        }
    }
    predictions_count.normalize();
    p_states_stats_mda.normalize();
    // output
    std::string result_str="";
    if (betas_efficiency_check==true) {
        result_str+=(boost::format("# betas_mda:  statistics: correct: %.3f, wrong: %.3f\n") 
            % predictions_count[0] % predictions_count[1]).str();
    } else {
        result_str+=(boost::format("# prediction: statistics: correct: %.3f, wrong: %.3f\n") 
            % predictions_count[0] % predictions_count[1]).str();
        result_str+="# states distrib. prediction: "+p_states_stats_mda.tostring();
        result_str+="# states distrib. from input: "+states_stats_mda.tostring();
    }
    return result_str;
}
//*----------------------------------------------------------------
void MLR_uni::get_log_lik_parallel(
        FLOAT_TYPE *ll, 
        const unsigned int input_idx_start,
        const unsigned int input_idx_end,
        mocapy::MDArray<FLOAT_TYPE> &p_states_mda_elem
        ) {
    if (_DEBUG>=10) { printf("# DEBUG: entering function: get_log_lik_parallel\n"); }
    //
    *ll=0.0;
    for (unsigned int i=input_idx_start; i<input_idx_end;i++) {
        *ll+=log(
            calc_p(
                states_mda[i],
                p_states_mda_elem,
                invec_mda.get_view(i),
                means_mda
            )
        );
    }
    if (_DEBUG>=100) { 
        printf("# start: %d, end: %d. LL partial: %.3f\n",input_idx_start,input_idx_end,*ll);
    }
}
//*----------------------------------------------------------------
FLOAT_TYPE MLR_uni::get_log_lik(void) {
    if (_DEBUG>=10) { printf("# DEBUG: entering function: get_log_lik\n"); }
    // loops over input data, uses means_mda as ML matrix
    boost::thread_group Tg;
    // splits the input data: Gibbs sampling is done independently for each slice
    for (unsigned int thread=0; thread<threads_number; thread++) {
        Tg.create_thread(
            boost::bind(
                &MLR_uni::get_log_lik_parallel,
                this,
                &ll_vec[thread],
                threads_input_intervals[thread].first,
                threads_input_intervals[thread].second,
                p_states_mda.get_view(thread)
            )
        );
    }
    Tg.join_all();
    FLOAT_TYPE ll=0;
    if (_DEBUG>=50) { printf("# LL partials: "); }
    for (unsigned int thread=0; thread<threads_number; thread++) {
        ll+=ll_vec[thread];
        if (_DEBUG>=50) { printf("%.3f ",ll_vec[thread]); }
    }
    if (_DEBUG>=50) { printf("-> LL tot: %.3f\n",ll); }
    return ll;
}
//*----------------------------------------------------------------
FLOAT_TYPE MLR_uni::sample_b_elem(
        const unsigned int j, const unsigned int k,
        const unsigned int thread
        ) {
    // used to resample a single element in the beta matrices
    if (_DEBUG>=25) { printf("# DEBUG: entering function: sample_b_elem\n"); }
    return get_rand_normal( means_mda.get(j,k),std_devs_mda.get(j,k), thread );
}
//*----------------------------------------------------------------
FLOAT_TYPE MLR_uni::get_rand_uniform(
        const unsigned int thread
        ) {
    if (_DEBUG>=30) { printf("# DEBUG: entering function: get_rand_uniform\n"); }
    return randomGen_vec[thread]->get_rand(use_lagged_fibonacci607);
}
FLOAT_TYPE MLR_uni::get_rand_normal(
        const FLOAT_TYPE mu, const FLOAT_TYPE sigma,
        const unsigned int thread
        ) {
    // gets a random sample from a normal(mu,sigma) distribution
    if (_DEBUG>=30) { printf("# DEBUG: entering function: get_rand_normal\n"); }
    return randomGen_vec[thread]->get_rand_normal(mu,sigma,use_lagged_fibonacci607);
}
//*----------------------------------------------------------------
void MLR_uni::resample_beta_slice(
        const unsigned int i,           // input slice ID
        const unsigned int gibbs_iter_num,
        const FLOAT_TYPE resample_perc,
        const unsigned int thread
        ) {
    // gibbs sampling over each single slice
    if (_DEBUG>=10) {
        printf("# DEBUG: entering function: resample_beta_slice\n");
        if (_DEBUG>=100) {
            printf("#    params: i: %d, gibbs_iter: %d, resample_perc: %.2f, thread: %d\n",
                i, gibbs_iter_num,
                resample_perc,
                thread
            );
        }
    }
    for (unsigned int gibbs_j=0; gibbs_j<gibbs_iter_num; gibbs_j++) {
        FLOAT_TYPE p1=0;
        FLOAT_TYPE p2=0;
        // clear cache, set to a slice of the complete betas matrix
        beta_tmp_mda.get_view(thread).multiply_inplace(0.0);
        beta_tmp_mda.get_view(thread).add_inplace( betas_mda.get_view(i) );
        // loops over beta_tmp elements
        for (unsigned int j=0; j<means_shape0; j++) {
            for (unsigned int k=0; k<means_shape1; k++) {
                // shall we resample element?
                if ( get_rand_uniform(thread) < resample_perc ) {
                    beta_tmp_mda.get_view(thread).set( j,k, sample_b_elem(j,k,thread) );
                }
            }
        }
        // P(state|input_vector,previous_beta_matrix)
        p1=p_states_accepted_mda[i];
        // P(state|input_vector,beta_matrix)
        p2=calc_p(
            states_mda[i],
            p_states_mda.get_view(thread),
            invec_mda.get_view(i),
            beta_tmp_mda.get_view(thread)
            );
        // Standard Metropolis-Hastings
        if ( (p2>p1) or (get_rand_uniform(thread)<(p2/p1)) ) {
            // accepted
            p_states_accepted_mda[i]=p2;
            betas_mda.get_view(i).multiply_inplace(0.0);
            betas_mda.get_view(i).add_inplace(beta_tmp_mda.get_view(thread));
        }
    }
}
//*----------------------------------------------------------------
void MLR_uni::em_step_do_parallel(
        const unsigned int input_idx_start,
        const unsigned int input_idx_end,
        const unsigned int thread,
        const unsigned int gibbs_iter_num,
        const FLOAT_TYPE resample_perc,
        const bool em_init_flag
        ) {
    // iterate over the splitted input
    if (_DEBUG>=10) {
        printf("# DEBUG: entering function: em_step_do_parallel\n");
        if (_DEBUG>=100) {
            printf("#    params: %d, %d, %d, %d, %f, %d\n",
                input_idx_start,input_idx_end,
                thread,
                gibbs_iter_num,resample_perc,
                em_init_flag
            );
        }
    }
    for (unsigned int i=input_idx_start; i<input_idx_end; i++) {
        resample_beta_slice( i, gibbs_iter_num, resample_perc, thread );
    }
}
//*----------------------------------------------------------------
void MLR_uni::em_step_do(
        const unsigned int gibbs_iter_num,
        const FLOAT_TYPE resample_perc,
        const bool em_init_flag
        ) {
    if (_DEBUG>=10) { printf("# DEBUG: entering function: do_em\n"); }
    // initialize beta matrices?
    if (em_init_flag==true) { init(); }
    // threaded?
    boost::thread_group Tg;
    // splits the input data: Gibbs sampling is done independently for each slice
    for (unsigned int thread=0; thread<threads_number; thread++) {
        Tg.create_thread(
            boost::bind(
                &MLR_uni::em_step_do_parallel,
                this,
                threads_input_intervals[thread].first,
                threads_input_intervals[thread].second,
                thread,
                gibbs_iter_num,
                resample_perc,
                em_init_flag
            )
        );
    }
    Tg.join_all();
    // update parameters for ML
    if ( not(em_init_flag) ) {
        update_parameters();
    }
}
void MLR_uni::em_do(
        const unsigned int em_iterations_max,
        const unsigned int em_iterations_without_improvements,
        const unsigned int em_progress_steps,
        const unsigned int em_step_gibbs_iterations,
        const FLOAT_TYPE em_step_resample_perc,
        const std::string ML_filename_out,
        bool init_flag
        ) {
    // EM loop
    // exit from EM if there have been no improvements in LL for the last N iterations
    FLOAT_TYPE ll=-1e30, ll_best=-1e30;
    unsigned int ll_best_step=0;
    for (unsigned int em_step_j=0; em_step_j<em_iterations_max; em_step_j++) {
        em_step_do(
            em_step_gibbs_iterations,
            em_step_resample_perc,
            init_flag
            );
        init_flag=false;
        ll=get_log_lik();
        printf( "step: %6d. LL: %.3f\n",em_step_j,ll );
        // store best LL
        if ( ll_best<ll) {
            ll_best=ll;
            ll_best_step=em_step_j;
            if ( ML_filename_out.size() > 0 ) {
                save_ML_proposal(ML_filename_out+".best_ll");
            }
        }
        if ( em_step_j % em_progress_steps ==0) {
            if ( em_step_j % (5*em_progress_steps) ==0) {
                printf("%s",evaluate_classification(true).c_str());
            }
            printf("%s",evaluate_classification().c_str());
            printf("# Beta matrix (ML):\n%s",get_ML().c_str());
            // should we save the ML matrices?
            if ( ML_filename_out.size() > 0 ) {
                save_ML_proposal(ML_filename_out);
            }
        }
        std::cout.flush();
        // shall we exit because no improvements?
        if ( (em_step_j-ll_best_step)>=em_iterations_without_improvements ) {
            printf("# step: %d. Exiting: %d iterations without LL improvements\n",em_step_j,em_step_j-ll_best_step);
            break;
        }
    }
}
//*----------------------------------------------------------------
void MLR_uni::update_parameters_parallel_mean(
        const unsigned int input_idx_start,
        const unsigned int input_idx_end,
        mocapy::MDArray<FLOAT_TYPE> *beta_tmp_mda_elem
        ) {
    if (_DEBUG>=10) { printf("# DEBUG: entering function: update_parameters_parallel_mean\n"); }
    // clear cache
    beta_tmp_mda_elem->multiply_inplace(0.0);
    for (unsigned int i=input_idx_start; i<input_idx_end; i++) {
        beta_tmp_mda_elem->add_inplace( betas_mda.get_view(i) );
    }
}
void MLR_uni::update_parameters_parallel_std_dev(
        const unsigned int input_idx_start,
        const unsigned int input_idx_end,
        mocapy::MDArray<FLOAT_TYPE> *beta_tmp_mda_elem,
        mocapy::MDArray<FLOAT_TYPE> *cache_tmp_mda_elem
        ) {
    if (_DEBUG>=10) { printf("# DEBUG: entering function: update_parameters_parallel_std_dev\n"); }
    // clear cache
    beta_tmp_mda_elem->multiply_inplace(0.0);
    for (unsigned int i=input_idx_start; i<input_idx_end; i++) {
        cache_tmp_mda_elem->multiply_inplace(0.0);
        cache_tmp_mda_elem->add_inplace( betas_mda.get_view(i) );
        cache_tmp_mda_elem->sub_inplace( means_mda );
        cache_tmp_mda_elem->sqr_all();
        // add to thread summation matrix
        beta_tmp_mda_elem->add_inplace(*cache_tmp_mda_elem);
    }
}
//*----------------------------------------------------------------
void MLR_uni::update_parameters(void) {
    // update ML parameters means and std_devs
    // loops over betas matrices elements
    // loops over input data
    if (_DEBUG>=10) { printf("# DEBUG: entering function: update_parameters\n"); }
    boost::thread_group Tg;
    // two pass algorithm: mean
    // splits the input data: uses beta_tmp_mda as temp matrix
    for (unsigned int thread=0; thread<threads_number; thread++) {
        Tg.create_thread(
            boost::bind(
                &MLR_uni::update_parameters_parallel_mean,
                this,
                threads_input_intervals[thread].first,
                threads_input_intervals[thread].second,
                &beta_tmp_mda.get_view(thread)
            )
        );
    }
    Tg.join_all();
    // assemble means_mda
    means_mda.multiply_inplace(0.0);
    for (unsigned int thread=0; thread<threads_number; thread++) {
        means_mda.add_inplace(beta_tmp_mda.get_view(thread));
    }
    means_mda.div_inplace(N);
    // two pass algorithm: std_dev
    for (unsigned int thread=0; thread<threads_number; thread++) {
        Tg.create_thread(
            boost::bind(
                &MLR_uni::update_parameters_parallel_std_dev,
                this,
                threads_input_intervals[thread].first,
                threads_input_intervals[thread].second,
                &beta_tmp_mda.get_view(thread),
                &cache_tmp_mda.get_view(thread)
            )
        );
    }
    Tg.join_all();
    // assemble std_devs_mda
    std_devs_mda.multiply_inplace(0.0);
    for (unsigned int thread=0; thread<threads_number; thread++) {
        std_devs_mda.add_inplace(beta_tmp_mda.get_view(thread));
    }
    std_devs_mda.div_inplace(N-1);
    //std_devs_mda.sqrt_all();
    for (unsigned int j=0; j<means_shape0; j++) {
        for (unsigned int k=0; k<means_shape1; k++) {
            std_devs_mda.set( j,k, sqrt(std_devs_mda.get(j,k)) );
        }
    }
}
//*----------------------------------------------------------------
std::string MLR_uni::get_ML(void) {
    if (_DEBUG>=10) { printf("# DEBUG: entering function: get_ML\n"); }
    return means_mda.tostring();
}
//*----------------------------------------------------------------
void MLR_uni::save_ML_proposal(
        const std::string &filename,
        const unsigned int precision,
        const std::string sep_record,
        const std::string sep_fields
        ) {
    if (_DEBUG>=10) { printf("# DEBUG: entering function: save_ML_proposal\n"); }
    // open output file
    std::ofstream fH(filename.c_str());
    if (!fH) { 
        printf("# ERROR: Cannot open file: '%s'\n",filename.c_str());
        exit(1);
    }
    fH << "# format: means[i,0],means[i,1],...;std_dev[i,0],...\n";
    fH << "# LL: " << get_log_lik() << "\n";
    fH << "# classification:\n" << evaluate_classification();
    // iterate over outer dimension
    string out_line="";
    for (unsigned int j=0; j<means_shape0; j++) {
        out_line=means_mda.get_view(j).tostring(precision,0,sep_fields,"") +
            sep_record +
            std_devs_mda.get_view(j).tostring(precision,0,sep_fields,"");
        if (_DEBUG>=100) { printf( "%s\n",out_line.c_str() ); }
        fH << out_line << "\n";
    }
    fH.close();
}
//*----------------------------------------------------------------
//* MAIN
//*----------------------------------------------------------------
int main(int argc, char *argv[]) {
    po::variables_map options=getCmdOptions(argc, argv);
    printOptions(options);
    // MLR_uni construction
    MLR_uni mlr_uni(
        options["threads"].as<unsigned int>(),
        options["seed"].as<unsigned int>(),
        options["generator-lagged-fibonacci607"].as<bool>(),
        options["ml-proposal-init-means"].as<FLOAT_TYPE>(),
        options["ml-proposal-init-std_devs"].as<FLOAT_TYPE>()
        );
    // input datafile loading
    mlr_uni.set_input_from_datafile( options["data-input-filename"].as<std::string>() );
    // ML proposal loading?
    if ( options["ml-proposal-filename-input"].as<std::string>().size() > 0 ) {
        mlr_uni.set_ML_proposal(options["ml-proposal-filename-input"].as<std::string>());
    }
    // REMOVE: printf("# DEVELOP: exiting\n");exit(0);
    // EM loop
    mlr_uni.em_do(
        options["em-iterations-max"].as<unsigned int>(),
        options["em-iterations-without-improvements"].as<unsigned int>(),
        options["em-progress-steps"].as<unsigned int>(),
        options["em-step-gibbs-iterations"].as<unsigned int>(),
        options["em-step-resample-perc"].as<FLOAT_TYPE>(),
        options["ml-proposal-filename-output"].as<std::string>()
    );
    printf( "# current LL: %.3f\n",mlr_uni.get_log_lik() );
    printf( "%s",mlr_uni.evaluate_classification(true).c_str() );
    printf( "%s",mlr_uni.evaluate_classification().c_str() );
    printf("# Beta matrix (ML):\n%s",mlr_uni.get_ML().c_str());
    // save ML matrices?
    if ( options["ml-proposal-filename-output"].as<std::string>().size() > 0 ) {
        mlr_uni.save_ML_proposal(options["ml-proposal-filename-output"].as<std::string>());
    }
    // exit
    return 0;
}
