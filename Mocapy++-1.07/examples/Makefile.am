AM_CXXFLAGS = -ggdb -Wall -I$(top_srcdir)/src $(BOOST_CPPFLAGS)
AM_LDFLAGS =  $(BOOST_LDFLAGS) $(LAPACK_LIBS) $(BLAS_LIBS) $(LIBS) 
LDADD = $(top_srcdir)/src/libMocapy.la $(BOOST_SERIALIZATION_LIB)

noinst_PROGRAMS = discrete_hmm_with_prior hmm_bippo hmm_discrete hmm_factorial hmm_gauss hmm_gauss_1d hmm_kent hmm_multinomial hmm_poisson hmm_simple hmm_torus hmm_vonmises infenginehmm_example infenginemm_example pythoninterface

discrete_hmm_with_prior_SOURCES = discrete_hmm_with_prior.cpp 

hmm_bippo_SOURCES = hmm_bippo.cpp 

hmm_discrete_SOURCES = hmm_discrete.cpp 

hmm_factorial_SOURCES = hmm_factorial.cpp 

hmm_gauss_SOURCES = hmm_gauss.cpp 

hmm_gauss_1d_SOURCES = hmm_gauss_1d.cpp 

hmm_kent_SOURCES = hmm_kent.cpp 

hmm_multinomial_SOURCES = hmm_multinomial.cpp 

hmm_poisson_SOURCES = hmm_poisson.cpp 

hmm_simple_SOURCES = hmm_simple.cpp 

hmm_torus_SOURCES = hmm_torus.cpp 

hmm_vonmises_SOURCES = hmm_vonmises.cpp 

infenginehmm_example_SOURCES = infenginehmm_example.cpp 

infenginemm_example_SOURCES = infenginemm_example.cpp 

pythoninterface_SOURCES = pythoninterface.cpp

