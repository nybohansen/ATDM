# The name of our project is "MOCAPY". CMakeLists files in this project can
# refer to the root source directory of the project as ${MOCAPY_SOURCE_DIR} and
# to the root binary directory of the project as ${MOCAPY_BINARY_DIR}.
cmake_minimum_required (VERSION 2.6.2)
project (MOCAPY C CXX Fortran)
# version number
set(MOCAPY_VERSION 0.1)

# Recurse into the "src" and "examples" subdirectories. This does not actually
# cause another cmake executable to run. The same process will walk through
# the project's entire directory structure.
add_subdirectory (src)
add_subdirectory (examples)
add_subdirectory (tests)

# settings for making distribution package
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "Mocapy is a Dynamic Bayesian Network toolkit")
set(CPACK_PACKAGE_VENDOR "Thomas Hamelryck and Martin Paluszewski")
#set(CPACK_PACKAGE_DESCRIPTION_FILE ${CMAKE_CURRENT_SOURCE_DIR}/README)
set(CPACK_GENERATOR TGZ)
set(
        CPACK_SOURCE_PACKAGE_FILE_NAME
        "Mocapy-${MOCAPY_VERSION}"
        CACHE INTERNAL "tarball basename"
        )
set(CPACK_SOURCE_GENERATOR TGZ)
include(CPack)

