AC_PREREQ(2.69) ## specificy version of autoconf
AC_INIT(ginseng, 0.0, jwala@broadinstitute.org)
AM_INIT_AUTOMAKE(foreign)
AC_CONFIG_SRCDIR([src/ginseng.cpp]) ## just checks that this file existss
AM_MAINTAINER_MODE([disable])

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

AC_ARG_WITH(boost, AS_HELP_STRING([--with-boost=PATH],
	[specify directory containing the boost library)]))
if test "$with_boost" -a -d "$with_boost"; then
    boost_include="-I$with_boost"
    AC_CHECK_FILE("$with_boost/stage/lib", STAGE_PATH=1, STAGE_PATH=0)
    if test ${STAGE_PATH} = 1; then
        boost_lib="$with_boost/stage/lib"  
    else
        boost_lib="$with_boost/lib"        
    fi
    boost_lib_flag="-L/$boost_lib"
    AC_SUBST(boost_lib)
fi

AC_ARG_WITH(apophenia, AS_HELP_STRING([--with-apophenia=PATH],
	[specify directory containing the Apophenia library)]))
if test "$with_apophenia" -a -d "$with_apophenia"; then
    apophenia_include="-I$with_apophenia/include"
    apophenia_lib="$with_apophenia/lib/libapophenia.a"        
    AC_SUBST(apophenia_lib)
    AC_SUBST(apophenia_include)
fi

AC_ARG_WITH(gsl, AS_HELP_STRING([--with-gsl=PATH],
	[specify directory containing the Gsl library)]))
if test "$with_gsl" -a -d "$with_gsl"; then
    gsl_include="-I$with_gsl/include"
    gsl_lib="-L$with_gsl/lib"
    AC_SUBST(gsl_lib)
fi

# Only fail on warnings when the --enable-development flag is passed into configure
AC_ARG_ENABLE(development, AS_HELP_STRING([--enable-development],
	[Turn on development options, like failing compilation on warnings]))
if test "$enable_development"; then
    fail_on_warning="-Werror"
fi

# Set compiler flags.
AC_SUBST(AM_CXXFLAGS, "-Wall -Wextra $fail_on_warning -Wno-unknown-pragmas -std=c++11 -g")
AC_SUBST(CXXFLAGS, "")
AC_SUBST(CFLAGS, "")
AC_SUBST(CPPFLAGS, "$CPPFLAGS $boost_include $apophenia_include $gsl_include")
AC_SUBST(LDFLAGS, "-fopenmp $LDFLAGS")

AC_SUBST(LIBS, " $apophenia_lib $gsl_lib  -lgsl -lgslcblas $LIBS")

# Ensure the GSL is available
AC_CHECK_HEADERS([gsl/gsl_rng.h],,[AC_MSG_ERROR([GSL library is required: See http://apophenia.info/setup.html])])

# Ensure the Apophenia is available
AC_CHECK_HEADERS([apop.h],,[AC_MSG_ERROR([Apophenia library is required: https://github.com/b-k/apophenia])])

AC_CONFIG_FILES([Makefile
		 SeqLib/src/Makefile
		 src/Makefile])

AC_OUTPUT
