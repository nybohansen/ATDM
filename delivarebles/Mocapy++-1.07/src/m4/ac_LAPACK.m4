# ac_LAPACK.m4 --- LAPACK m4 macros
# Copyright (C) 2008 Wouter Boomsma
#
# This file is part of Phaistos 
#
# Phaistos is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Phaistos is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with Phaistos.  If not, see <http://www.gnu.org/licenses/>.


# Check whether framework contains function
AC_DEFUN([AC_CHECK_FUNC_IN_FRAMEWORK],
[
  AC_CACHE_CHECK([for $1 in framework $2], 
  		 [ac_check_func_in_fw_$2],
  ac_check_func_in_fw_tmp=$LIBS
  ac_check_func_in_fw_$2=[no]
  LIBS="-framework $2 $LIBS"
  AC_TRY_LINK_FUNC([$1], [ac_check_func_in_fw_$2="-framework $2"], [])
  LIBS=$ac_check_func_in_fw_tmp)
])


# Wrapper for the ACX_LAPACK function. This inludes a test
# for the vecLib library available on OS X
AC_DEFUN([AC_LAPACK], 
[

  AC_REQUIRE([ACX_LAPACK])
  AC_SUBST([LAPACK_LIBS], ["$LAPACK_LIBS $BLAS_LIBS"])
  
  AC_CHECK_FUNC_IN_FRAMEWORK([sgemm_], [vecLib])

  if test "$ac_check_func_in_fw_vecLib" != no; then
    AC_SUBST([LAPACK_LIBS], [$ac_check_func_in_fw_vecLib])
  fi
])

