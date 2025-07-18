#
#				configure.ac
#
# Process this file with autoconf to generate a configure script.
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#	This file part of:	SkyMaker
#
#	Copyright:		(C) 1998-2021 IAP/CNRS/SorbonneU
#	          		(C) 2021-2023 CFHT/CNRS
#	          		(C) 2023-2025 CEA/AIM/UParisSaclay
#
#	License:		GNU General Public License
#
#	SkyMaker is free software: you can redistribute it and/or modify
#	it under the terms of the GNU General Public License as published by
#	the Free Software Foundation, either version 3 of the License, or
# 	(at your option) any later version.
#	SkyMaker is distributed in the hope that it will be useful,
#	but WITHOUT ANY WARRANTY; without even the implied warranty of
#	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#	GNU General Public License for more details.
#	You should have received a copy of the GNU General Public License
#	along with SkyMaker.  If not, see <http://www.gnu.org/licenses/>.
#
#	Last modified:		02/04/2025
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# First, disable the annoying config.cache
define([AC_CACHE_LOAD],)
define([AC_CACHE_SAVE],)

# This is your standard Bertin source code...
AC_INIT([SkyMaker],[4.3.0], [https://github.com/astromatic/skymaker],
	skymaker, [https://astromatic.net/software/skymaker])
AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_AUX_DIR(autoconf)
AC_CONFIG_SRCDIR(src/makeit.c)
AC_CONFIG_HEADERS(config.h)
AM_INIT_AUTOMAKE
date=`date +%Y-%m-%d`
date1=`date -R`
date2=`date +"%a %b %d %Y"`
date3=`date +"%B %Y"`
AC_DEFINE_UNQUOTED(DATE, "$date", [Archive creation date])
AC_SUBST(PACKAGER, "AstrOmatic")
AC_SUBST(DATE1, "$date1")
AC_SUBST(DATE2, "$date2")
AC_SUBST(DATE3, "$date3")

# Checks for programs.
# AC_LANG(C)
AC_PROG_CC
AC_PROG_CXX

# Provide a special option for setting the package release number
AC_ARG_WITH(release,
	[AS_HELP_STRING([--with-release=<release number>],
	[set the package release number (default = 1)])],
	[],
	[with_release="no"])
if test "$with_release" = "no"; then
  with_release=1
fi
AC_SUBST(PACKAGE_RELEASE, "$with_release")

# Display pakage and version number
AC_MSG_RESULT([******** Configuring:  $PACKAGE_NAME $PACKAGE_VERSION - $PACKAGE_RELEASE ($date) ********])

# Provide special option for choosing automatically the compilation flags.
AC_MSG_CHECKING([if compilation flags are set automatically])
AC_ARG_ENABLE(auto-flags,
	[AS_HELP_STRING([--enable-auto-flags],
	[Let the configure script choose the compilation flags (default = no)])],
	AC_MSG_RESULT([yes]),
	AC_MSG_RESULT([no]))

# Provide special option for the Linux Intel C "classic" compiler
AC_MSG_CHECKING([whether the classic INTEL compiler is enabled])
AC_ARG_ENABLE(icc,
	[AS_HELP_STRING([--enable-icc],
	[Use the Intel compiler (default = no)])],
        CC="icc"
        enable_iccx="yes"
	AC_MSG_RESULT([yes]),
	AC_MSG_RESULT([no]))

# Provide special option for the Linux Intel OneAPI C compiler
AC_MSG_CHECKING([whether the INTEL compiler is enabled])
AC_ARG_ENABLE(icx,
	[AS_HELP_STRING([--enable-icx],
	[Use the Intel compiler (default = no)])],
        CC="icx"
        enable_iccx="yes"
	AC_MSG_RESULT([yes]),
	AC_MSG_RESULT([no]))

# Provide special options for INTEL MKL
# We use icc if available, if not Intel OneAPI's icx
# (we don't use AC_PROG_CC as it does not play nice
# in a conditional block)
AC_MSG_CHECKING([whether INTEL's MKL is enabled])
AC_ARG_ENABLE(mkl,
	[AS_HELP_STRING([--enable-mkl],
	[Use INTEL's MKL for FFTs (default = no)])],
	    AC_MSG_RESULT([yes])
        CC=""
        [AC_CHECK_PROGS(CC, [icc icx cc])]
        enable_iccx="yes",
	AC_MSG_RESULT([no]))

AC_SEARCH_LIBS([strerror],[cposix])
if test "$enable_auto_flags" = "yes"; then
  CFLAGS=""
  LDFLAGS=""
  ACX_PROG_CC_OPTIM
fi
LT_INIT
AC_PROG_INSTALL

# Checks for libraries.
AC_CHECK_LIB(m, sin)
# the Calabretta wcslib
AC_CHECK_LIB([wcs], [wcsinit],[],[],[-lwcs])

AC_CHECK_HEADERS([wcslib/wcs.h],[],[
AC_MSG_ERROR([Error: Could  not find header wcslib/wcs.h ... aborting...])
])



# Checks for header files.
AC_CHECK_HEADERS([limits.h malloc.h stdlib.h string.h sys/mman.h \
		sys/time.h  sys/types.h unistd.h])
# Checks for INTEL math header files.
if test "$enable_iccx" = "yes"; then
  AC_CHECK_HEADERS(mathimf.h)
fi

# Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST
AC_TYPE_OFF_T
AC_TYPE_SIZE_T
AC_TYPE_LONG_LONG_INT
AC_TYPE_UNSIGNED_LONG_LONG_INT
AC_STRUCT_TM
AC_TYPE_UID_T

# Checks for library functions.
AC_FUNC_ERROR_AT_LINE
AC_FUNC_MMAP
AC_FUNC_STAT
AC_FUNC_STRFTIME
AC_CHECK_FUNCS([atexit getenv gettimeofday memcpy memmove memset munmap \
		posix_memalign sincos sincosf strstr sysconf])
 
# Check support for large files
AC_SYS_LARGEFILE
AC_FUNC_FSEEKO

# Provide special options for FFTW
AC_ARG_WITH(fftw-libdir,
	[AS_HELP_STRING([--with-fftw-libdir=<FFTW library path>],
	[Provide an alternative path to the FFTW library])])
AC_ARG_WITH(fftw-incdir,
	[AS_HELP_STRING([--with-fftw-incdir=<FFTW include dir>],
	[Provide an alternative path to the FFTW include directory])])

# Provide special options for the MKL library
AC_ARG_WITH(mkl-dir,
	[AS_HELP_STRING([--with-mkl-dir=<MKL path>],
	[Provide an alternative path to the MKL library])])

# Provide special option for CFITSIO
AC_MSG_CHECKING([whether CFITSIO support should be disabled (default=enabled)])
AC_ARG_ENABLE(cfitsio,
	[AS_HELP_STRING([--disable-cfitsio],
	[Disable support for compressed FITS files through the CFITSIO library (default=enabled)])],
	if test "$enable_cfitsio" = "no"; then
	  AC_MSG_RESULT([yes])
	else
	  AC_MSG_RESULT([no])
	fi,
	AC_MSG_RESULT([no]))
AC_ARG_WITH(cfitsio-libdir,
	[AS_HELP_STRING([--with-cfitsio-libdir=<CFITSIO library path>],
	[Provide an alternative path to the CFITSIO library])])
AC_ARG_WITH(cfitsio-incdir,
	[AS_HELP_STRING([--with-cfitsio-incdir=<CFITSIO include dir>],
	[Provide an alternative path to the CFITSIO include directory])])


# Set flags for multithreading
n_pthreads=1024
AC_ARG_ENABLE(threads,
	[AS_HELP_STRING([--enable-threads@<:@=<max_number_of_threads>@:>@],
	[Enable multhreading (on with up to 1024 threads by default)])],
    if test "$enableval" = "no"; then
      use_pthreads="no"
    else
      use_pthreads="yes"
      if test "$enableval" != "yes"; then
        n_pthreads=$enableval
      fi
    fi,
    use_pthreads=yes
    )

# Provide special option for profiling
AC_MSG_CHECKING([for profiler mode])
AC_ARG_ENABLE(profiling,
	[AS_HELP_STRING([--enable-profiling],
	[Enable special mode for profiling (default = no)])],
	AC_MSG_RESULT([yes]),
	AC_MSG_RESULT([no]))

# Enable linking options for making the executable as portable as possible.
AC_MSG_CHECKING([best linking option])
AC_ARG_ENABLE(best-link,
	[AS_HELP_STRING([--enable-best-link],
	[Choose the right combination of static and dynamic linking to make \
the executable as portable as possible (default = no)])],
	AC_MSG_RESULT([yes]),
	AC_MSG_RESULT([no]))

 AX_OPENMP

################# Actions to complete in case of multhreading ################
AC_DEFINE_UNQUOTED(THREADS_NMAX, $n_pthreads,[Maximum number of POSIX threads])
if test "$use_pthreads" = "yes"; then
  AC_MSG_CHECKING([for multithreading])
  AC_MSG_RESULT([maximum of $n_pthreads thread(s)])
  AC_DEFINE(USE_THREADS, 1, [Triggers multhreading])
# CC, CFLAGS and LIBS are system and compiler-dependent
  ACX_PTHREAD
  CC="$PTHREAD_CC"
  [AM_CFLAGS="$AM_CFLAGS $OPENMP_CFLAGS $PTHREAD_CFLAGS -D_REENTRANT"]
  [AM_CXXFLAGS="$AM_CXXFLAGS $OPENMP_CFLAGS $PTHREAD_CFLAGS"]
  LIBS="$LIBS $PTHREAD_LIBS"
fi
AM_CONDITIONAL(USE_THREADS, test $use_pthreads = "yes")

#################### handle the INTEL MKL library (FFTW) ####################
  if test "$enable_mkl" = "yes"; then
    convlibs="${srcdir}/wcs/libwcs_c.a"
    ACX_MKL($with_mkl_dir,,$enable_best_link,$convlibs)
    if test "$MKL_WARN" != ""; then
      AC_MSG_WARN([$MKL_WARN])
    fi
    AM_CFLAGS="$AM_CFLAGS $MKL_CFLAGS "
    AM_LDFLAGS="$AM_LDFLAGS $MKL_LDFLAGS "
    LIBS="$LIBS $MKL_LIBS"
  else
################ handle the FFTW library (Fourier transforms) ################
    ACX_FFTW($with_fftw_libdir, $with_fftw_incdir, $use_pthreads,yes,
      [
        [LIBS="$FFTW_LIBS $LIBS"]
        if test "$FFTW_WARN" != ""; then
          AC_MSG_WARN([$FFTW_WARN])
        fi
      ],
      AC_MSG_ERROR([$FFTW_ERROR Exiting.])
    )
fi

########################## handle the CFITSIO library ########################
if test "$enable_cfitsio" != "no"; then
  ACX_CFITSIO($with_cfitsio_libdir, $with_cfitsio_incdir,
	[LIBS="$CFITSIO_LIBS $LIBS"]
	if test "$CFITSIO_WARN" != ""; then
	  AC_MSG_WARN([$CFITSIO_WARN])
	fi,
	AC_MSG_ERROR([$CFITSIO_ERROR Exiting.]))
fi

AM_CONDITIONAL(USE_MODEL, [test "$enable_model_fitting" != "no"])

# Compile with profiling option
if test "$enable_profiling" = "yes"; then
  if test "$enable_iccx" = "yes"; then
    AM_CFLAGS="$AM_CFLAGS -pq"
  else
    AM_CFLAGS="$AM_CFLAGS -pg"
  fi
  enable_best_link="no"
fi

# "Best" linking option
if test "$enable_best_link" = "yes"; then
  if test "$enable_iccx" = "yes"; then
    AM_LDFLAGS="-static-intel -qopenmp-link=static -shared-libgcc \
	-static-libtool-libs -avoid-version $AM_LDFLAGS"
  else
    AM_LDFLAGS="-shared-libgcc -static-libtool-libs -avoid-version $AM_LDFLAGS"
  fi
fi

AC_SUBST(AM_CFLAGS)
AC_SUBST(AM_CXXFLAGS)
AC_SUBST(AM_CPPFLAGS)
AC_SUBST(AM_LDFLAGS)

# Display compiler and linker flags
AC_MSG_RESULT([***************************************************************])
AC_MSG_RESULT([Compile cmdline:  $CC $AM_CPPFLAGS $CPPFLAGS $AM_CFLAGS $CFLAGS])
AC_MSG_RESULT([Link    cmdline:  $CC $AM_LDFLAGS $LDFLAGS $LIBS])
AC_MSG_RESULT([Default XSLT URL: $xsl_url])
AC_MSG_RESULT([***************************************************************])

AC_CONFIG_FILES([
Makefile
skymaker.spec
debian/Makefile
debian/changelog
doc/Makefile
doc/src/conf.py
man/Makefile
man/sky.1
src/Makefile
src/fits/Makefile
])
AC_OUTPUT
