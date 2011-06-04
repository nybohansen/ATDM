SET(BUILD_RELEASE ON)
SET(BUILD_DEBUG OFF)
SET(BUILD_PROFILE OFF)
SET(BUILD_DUMMYRANDOM OFF)

SET(Boost_USE_MULTITHREAD OFF)
       
FIND_PACKAGE(Boost 1.39)

INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS})
LINK_DIRECTORIES(${Boost_LIBRARY_DIRS})

SET(Mocapy_src
utils/LogFactorial.cpp
multinomial/MultinomialDensities.cpp
multinomial/MultinomialDensity.cpp
multinomial/MultinomialESS.cpp
dirichlet/DirichletDensities.cpp
dirichlet/DirichletDensity.cpp
dirichlet/DirichletESS.cpp
inference/abstractinfengine.cpp
framework/dbn.cpp
discrete/discretedensities.cpp
discrete/discreteess.cpp
discrete/discretepriors.cpp
framework/emengine.cpp
framework/essbase.cpp
inference/detail/forwardbacktracker.cpp
inference/detail/generalinfenginemm.cpp
gaussian/gaussiandensities.cpp
gaussian/gaussianess.cpp
gaussian/shrinkage/Shrinkage.cpp
gaussian/shrinkage/ShrinkageHelpers.cpp
inference/gibbsrandom.cpp
inference/fwbtrandom.cpp
inference/infenginemm.cpp
inference/infenginehmm.cpp
inference/infenginemcmc.cpp
inference/mcmc.cpp
inference/moreabstractinfengine.cpp
gaussian/multigauss.cpp
framework/node.cpp
framework/nodefactory.cpp
framework/parentmap.cpp
poisson/poissondensities.cpp
poisson/poissoness.cpp
poisson/poissonsample.cpp
bippo/bippodensities.cpp
bippo/bippoess.cpp
bippo/bippo.cpp
utils/random_data.cpp
utils/randomgen.cpp
inference/sampleinfengine.cpp
utils/utils.cpp
vonmises/vonmises.cpp
vonmises2d/vonmises2d.cpp
vonmises2d/vonmises2ddensities.cpp
vonmises2d/vonmises2dess.cpp
vonmises/vonmisesdensities.cpp
vonmises/vonmisesess.cpp
kent/kentdensities.cpp
kent/kentess.cpp
)

SET(Integrate_src
utils/netlib/integrate.cpp
utils/netlib/integrate.h
utils/netlib/slatec/d1mach.f
utils/netlib/slatec/d9lgmc.f
utils/netlib/slatec/dcsevl.f
utils/netlib/slatec/dgamlm.f
utils/netlib/slatec/dgamma.f
utils/netlib/slatec/dqags.f
utils/netlib/slatec/dqagse.f
utils/netlib/slatec/dqelg.f
utils/netlib/slatec/dqk21.f
utils/netlib/slatec/dqpsrt.f
utils/netlib/slatec/fdump.f
utils/netlib/slatec/i1mach.f
utils/netlib/slatec/initds.f
utils/netlib/slatec/j4save.f
utils/netlib/slatec/xercnt.f
utils/netlib/slatec/xerhlt.f
utils/netlib/slatec/xermsg.f
utils/netlib/slatec/xerprn.f
utils/netlib/slatec/xersve.f
utils/netlib/slatec/xgetua.f
)

SET(Functions_src
utils/vector_nD.h
utils/optimize.h
utils/netlib/cephes/chbevl.c
utils/netlib/cephes/const.c
utils/netlib/cephes/gamma.c
utils/netlib/cephes/zeta.c
utils/netlib/cephes/i0.c
utils/netlib/cephes/i1.c
utils/netlib/cephes/isnan.c
utils/netlib/cephes/mconf.h
utils/netlib/cephes/mtherr.c
utils/netlib/cephes/polevl.c
)

SET(LIBRARY_OUTPUT_PATH ${MOCAPY_SOURCE_DIR}/libs)

# build release library
IF(BUILD_RELEASE)
add_library (Mocapy STATIC
        ${Mocapy_src}
        ${Integrate_src}
        ${Functions_src}
        )

SET_TARGET_PROPERTIES(Mocapy PROPERTIES COMPILE_FLAGS " -O3 -fPIC -c -fmessage-length=0 -I${MOCAPY_SOURCE_DIR}/src")
ENDIF(BUILD_RELEASE)

#IF(BUILD_RELEASE)
#add_library (Mocapy SHARED
#        ${Mocapy_src}
#        ${Integrate_src}
#        ${Functions_src}
#        )

#SET_TARGET_PROPERTIES(Mocapy PROPERTIES COMPILE_FLAGS " -O3 -fPIC -c -fmessage-length=0 -I${MOCAPY_SOURCE_DIR}/src")
#ENDIF(BUILD_RELEASE)


IF(BUILD_DEBUG)
add_library (MocapyD 
        ${Mocapy_src}
        ${Integrate_src}
        ${Functions_src}
        )
SET_TARGET_PROPERTIES(MocapyD PROPERTIES COMPILE_FLAGS "-O0 -fPIC -g -Wall -c -fmessage-length=0")
ENDIF(BUILD_DEBUG)

IF(BUILD_PROFILE)
add_library (MocapyP 
        ${Mocapy_src}
        ${Integrate_src}
        ${Functions_src}
        )
SET_TARGET_PROPERTIES(MocapyP PROPERTIES COMPILE_FLAGS "-pg -fPIC -O3 -Wall -c -fmessage-length=0")
ENDIF(BUILD_PROFILE)

IF(BUILD_DUMMYRANDOM)
add_library (MocapyDR 
        ${Mocapy_src}
        ${Integrate_src}
        ${Functions_src}
        )
SET_TARGET_PROPERTIES(MocapyDR PROPERTIES COMPILE_FLAGS "-DRANDDUMMY  -O3 -fPIC -Wall -c -fmessage-length=0")
ENDIF(BUILD_DUMMYRANDOM)


# install commands
#INSTALL(TARGETS Mocapy 
#        LIBRARY DESTINATION lib
#        ARCHIVE DESTINATION lib
#        )
