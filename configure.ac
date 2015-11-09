#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.
AC_PREREQ(2.59) ## specificy version of autoconf
AC_INIT(swap, 0.0, jwala@broadinstitute.org)
AM_INIT_AUTOMAKE(foreign)
AC_CONFIG_SRCDIR([src/swap.cpp]) ## just checks that this file existss
AM_MAINTAINER_MODE([disable])
##AC_CONFIG_HEADER([config.h])
##m4_include([m4/m4_ax_openmp.m4])

# Checks for programs.
AC_PROG_CXX ## test for cpp compiler
AC_PROG_CC ## test for C compiler
AC_PROG_RANLIB ## required if libraries are built in package

# Check for headers
AC_LANG([C++])
AC_CHECK_HEADER([zlib.h])

# Check for libraries
##AC_SEARCH_LIBS([library],[function], [action-if-found], [action if not found])
AC_SEARCH_LIBS([gzopen],[z],,[AC_MSG_ERROR([libz not found, please install zlib (http://www.zlib.net/)])])
AC_SEARCH_LIBS([clock_gettime], [rt], [AC_DEFINE([HAVE_CLOCK_GETTIME], [1], [clock_getttime found])], )
AC_SEARCH_LIBS([pthread_create], [pthread], [AC_DEFINE([HAVE_PTHREAD], [1], [pthread found])], )

# Only fail on warnings when the --enable-development flag is passed into configure
AC_ARG_ENABLE(development, AS_HELP_STRING([--enable-development],
	[Turn on development options, like failing compilation on warnings]))
if test "$enable_development"; then
    fail_on_warning="-Werror"
fi

# Set compiler flags.
AC_SUBST(AM_CXXFLAGS, "-Wall -Wextra $fail_on_warning -Wno-unknown-pragmas -std=c++11")
AC_SUBST(CXXFLAGS, "-Ofast")
AC_SUBST(CFLAGS, "-Ofast")
AC_SUBST(CPPFLAGS, "$CPPFLAGS")
AC_SUBST(LDFLAGS, "$LDFLAGS")

# We always need to specify to link in bamtools
AC_SUBST(LIBS, "$LIBS")

# Make sure the bamtools headers can be found
##AC_CHECK_HEADERS([plinkio/bed.h],,[AC_MSG_ERROR([The plinkio library must be installed (https://github.com/mfranberg/libplinkio). You can specify its path with the --with-plinkio=PATH option])])

AC_CONFIG_FILES([Makefile
	         SnowTools/src/Makefile
		 src/Makefile])

AC_OUTPUT
