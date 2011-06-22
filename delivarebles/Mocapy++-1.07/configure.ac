AC_INIT([Mocapy], [0.01], [palu@binf.ku.dk], [mocapy])
AM_INIT_AUTOMAKE([-Wall -Werror foreign])
AC_PROG_CC
AC_PROG_CXX
AC_PROG_F77([gfortran g77])
AC_F77_LIBRARY_LDFLAGS
AC_F77_MAIN
AC_FC_LIBRARY_LDFLAGS

AC_PROG_LIBTOOL

AC_LAPACK

AX_BOOST_BASE
AX_BOOST_SERIALIZATION

AM_CONDITIONAL([COND_G77], [test $F77 = g77])

AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES([
 Makefile
 src/utils/netlib/Makefile
 src/Makefile
 examples/Makefile
])

AC_OUTPUT