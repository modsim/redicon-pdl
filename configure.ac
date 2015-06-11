#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.61)
AC_INIT(pdl, 0.0.0, valiska@gmail.com)
AC_CONFIG_SRCDIR([src/system.hh])
AC_CONFIG_HEADER([config.h])

PACKAGE=mplib
VERSION=0.0.0

dnl IFACE_NUMBER (current) : REV : AGE (< IFACE_NUMBER)
MPLIB_LT_VERSION="0:0:0"
AC_SUBST(MPLIB_LT_VERSION)

AC_SUBST(PACKAGE)
AC_SUBST(PACKAGE_STRING)
AC_SUBST(VERSION)

#AM_INIT_AUTOMAKE($PACKAGE, $VERSION)
AM_INIT_AUTOMAKE()


dnl maintainer
AC_SUBST(PDL_MAINTAINER)
PDL_MAINTAINER="valiska@gmail.com"

# Checks for programs.
AC_PROG_CC
AC_PROG_CXX
#AC_PROG_RANLIB
AC_PROG_LIBTOOL

# Checks for libraries.
#
# GNU Scientific Library
#
#AM_PATH_GSL(1.3, AC_DEFINE_UNQUOTED([GSL_VERSION],
#                    "`$ac_cv_path_GSL_CONFIG --version`" , [GSL version we use to compile this package]),
#  AC_MSG_ERROR([gsl is missing. You need gsl to compile this package.
#      See http://www.gnu.org/software/gsl/])
#)


# Checks for header files.

# Checks for typedefs, structures, and compiler characteristics.

# Checks for library functions.

# 'Optimizization' must be beore all the other flags!
AC_ARG_ENABLE([optimize],
    [AC_HELP_STRING([--enable-optimize],
            [disable debugging (removes -g flags) and set optimization level O3 (disabled)])],
    [CFLAGS="-O3"], [CFLAGS="$CFLAGS_save"])

# Enable/disable assetions
AC_ARG_ENABLE(assertions,
	      AC_HELP_STRING([--disable-assertions],
			     [disable all assertions (enabled)]), 
			     [disable_assertions=yes AC_DEFINE([MPLIB_NO_ASSERT],,[no assertions])],
			      disable_assertions=no)

if test "x$disable_assertions" = "xyes"; then
	CFLAGS="-DMPLIB_NO_ASSERT $CFLAGS"
#	CFLAGS="-DNDEBUG $CFLAGS"
fi

PDLIB_CFLAGS="-I$includedir"
PDLIB_LIBS="-L$libdir -lpdlib"

AC_SUBST(PDLIB_CFLAGS)
AC_SUBST(PDLIB_LIBS)

AC_CONFIG_FILES([Makefile
		pdlib.pc pdlib-config
		src/Makefile
		src/tests/Makefile
		])
AC_OUTPUT
