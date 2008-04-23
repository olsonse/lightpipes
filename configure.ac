#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.59)
AC_INIT(lightpipes, rev-1, olsonse@umich.edu)
AM_INIT_AUTOMAKE(lightpipes,rev-1)
AC_CONFIG_SRCDIR([src/lib/Field.h])
AC_CONFIG_HEADER([src/lib/config.h])
AX_PREFIX_CONFIG_H([src/lib/$PACKAGE-config.h])
AC_PREFIX_DEFAULT(`pwd`)



# Checks for programs.
AC_PROG_CXX
AC_PROG_CC
AC_PROG_MAKE_SET
AC_PROG_F77(gfortran g77)
AC_PROG_RANLIB
AC_PROG_MAKE_SET
AC_PROG_INSTALL
AM_PROG_MKDIR_P

# Check for pkg-config
AC_PATH_PROG(PKGCONFIG, pkg-config, false)



#                      BEGIN OCTAVE CHECK
probe_for_mkoctfile(){
    # Check for mkoctfile
    AC_PATH_PROG(MKOCTFILE, mkoctfile, false)
    AM_CONDITIONAL([OCTAVE_ENABLED], [test "$MKOCTFILE" != "false"])

    if test "$MKOCTFILE" != "false" ; then
        FPICFLAG=`$MKOCTFILE -p FPICFLAG`
    fi
}

AC_ARG_ENABLE(octave,
    AC_HELP_STRING([--disable-octave], [Disable octave components.]),
[ if test "$enableval" = "yes" ; then
    probe_for_mkoctfile
  else
    echo "octave extensions disabled"
  fi
],[
    probe_for_mkoctfile
])
#                        END OCTAVE CHECK

CFLAGS="$CFLAGS $FPICFLAG"
CXXFLAGS="$CXXFLAGS $FPICFLAG"



# find out the exe extension for this platform.
AC_EXEEXT
AC_OBJEXT

AC_LANG(C++)

# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([fcntl.h stdlib.h string.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_CONST
AC_C_INLINE
AC_TYPE_MODE_T
AC_TYPE_SIZE_T

# Checks for library functions.
AC_FUNC_ERROR_AT_LINE
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_CHECK_FUNCS([floor memset pow sqrt strstr])


AC_CONFIG_FILES([Makefile
                 lightpipes.pc
                 src/Makefile
                 src/lib/Makefile
                 src/progs/Makefile
                 octave/Makefile
                 ])
AC_OUTPUT
