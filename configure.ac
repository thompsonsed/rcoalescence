AC_INIT(all, 1.0)
AC_CONFIG_MACRO_DIR([inst/extdata/m4])
dnl AC_CONFIG_MACRO_DIR([extdata/m4])
CPPFLAGS="$CPPFLAGS"
dnl Switch to a C++ compiler, and check if it works.
m4_include(inst/extdata/m4/boost.m4)
# find R home and set correct compiler + flags
: ${R_HOME=`R RHOME`}
if test -z "${R_HOME}"; then
  AC_MSG_ERROR([cannot determine R_HOME. Make sure you use R CMD INSTALL!])
fi
AC_MSG_NOTICE([R_HOME: ${R_HOME}])
RBIN="${R_HOME}/bin/R"
# pick all flags for testing from R
: ${CC=`"${RBIN}" CMD config CC`}
: ${CXX=`"${RBIN}" CMD config CXX`}
: ${CPP=`"${RBIN}" CMD config CPP`}
: ${CFLAGS=`"${RBIN}" CMD config CFLAGS`}
: ${CPPFLAGS=`"${RBIN}" CMD config CPPFLAGS`}
: ${LDFLAGS=`"${RBIN}" CMD config LDFLAGS`}

PKG_CXXFLAGS="-I. -std=c++14"

AC_LANG(C++)
AC_PROG_CXX
PKG_LIBS=${PKGLIBS}
dnl Check sqlite headers exist
AC_CHECK_HEADERS(sqlite3.h, [has_sqlite3=1], AC_ERROR([A working sqlite3 library is required. Make sure it is correctly installed on your system.]))
with_sqlite3="yes"
LIBS="$LIBS -lsqlite3"

dnl Check boost headers
AC_LANG(C++)
BOOST_REQUIRE([1.54.0], AC_ERROR([A working boost library v1.54 or greater is required. Make sure it is correctly installed on your system.]))

ifdef(BOOST_CPPFLAGS,[boost_found="yes"],[boost_found="no"])
substring="-I/"
if [[ "$BOOST_CPPFLAGS" == "${BOOST_CPPFLAGS%$substring*}" ]];
then
	LIBS="$LIBS -lboost_filesystem -lboost_system -lpthread"
	AC_MSG_WARN([Library locations not detecting: using default locations.])
	boost_found="install detected: location unknown"
else
	dnl goofy workaround for boost libraries
	lib_boost="${BOOST_CPPFLAGS}"
	lib_boost="${lib_boost/-I/-L}"
	lib_boost="${lib_boost/include/lib}"
	CPPFLAGS="$CPPFLAGS $BOOST_CPPFLAGS"
	LDFLAGS="$LDFLAGS $lib_boost"
	LIBS="$LIBS -lboost_filesystem -lboost_system -lpthread"
	boost_found="yes"
fi
BOOST_FIND_HEADER([boost/filesystem.hpp])

dnl Check for gdal library
GDAL_CONFIG="gdal-config"

GDAL_CONFIG_SET="no"

AC_ARG_WITH([gdal-config],
    AC_HELP_STRING([--with-gdal-config=GDAL_CONFIG],
           [the location of gdal-config]),
           [gdal_config=$withval])
if test [ -n "$gdal_config" ] ; then
    GDAL_CONFIG_SET="yes"
    AC_SUBST([GDAL_CONFIG],["${gdal_config}"])
    AC_MSG_NOTICE(gdal-config set to $GDAL_CONFIG)
fi

if test ["$GDAL_CONFIG_SET" = "no"] ; then
  AC_PATH_PROG([GDAL_CONFIG], ["$GDAL_CONFIG"],["no"])
  if test ["$GDAL_CONFIG" = "no"] ; then
    AC_MSG_RESULT(no)
    AC_MSG_ERROR([gdal-config not found or not executable.])
  fi
else
  AC_MSG_CHECKING(gdal-config exists)
  if test -r "${GDAL_CONFIG}"; then
    AC_MSG_RESULT(yes)
  else
    AC_MSG_RESULT(no)
    AC_MSG_ERROR([gdal-config not found - configure argument error.])
  fi
  AC_MSG_CHECKING(gdal-config executable)
  if test -x "${GDAL_CONFIG}"; then
    AC_MSG_RESULT(yes)
  else
    AC_MSG_RESULT(no)
    AC_MSG_ERROR([gdal-config not executable.])
  fi
fi

AC_MSG_CHECKING(gdal-config usability)
if test `${GDAL_CONFIG} --version`;
then

	GDAL_CPPFLAGS=`${GDAL_CONFIG} --cflags`
	GDAL_LIBS=`${GDAL_CONFIG} --libs`
  GDAL_VERSION=`${GDAL_CONFIG} --version`
  GDAL_DEP_LIBS=`${GDAL_CONFIG} --dep-libs`
  GDAL_DATADIR=`${GDAL_CONFIG} --datadir`
  gdalok=yes
  AC_MSG_RESULT(yes)
else
  AC_MSG_RESULT(no)
	echo "Error: gdal-config not found"
	echo "The gdal-config script distributed with GDAL could not be found."
	echo "If you have not installed the GDAL libraries, you can"
	echo "download the source from  http://www.gdal.org/"
	echo "If you have installed the GDAL libraries, then make sure that"
	echo "gdal-config is in your path. Try typing gdal-config at a"
	echo "shell prompt and see if it runs. If not, use:"
	echo " --configure-args='--with-gdal-config=/usr/local/bin/gdal-config'"
	echo "with appropriate values for your installation."
	echo ""

	exit 1

fi

AC_MSG_NOTICE([GDAL: ${GDAL_VERSION}])
AC_MSG_CHECKING([GDAL version >= 2.1.0])
GDAL_VER_DOT=`echo $GDAL_VERSION | tr -d "."`
if test ${GDAL_VER_DOT} -lt 163 ; then
  AC_MSG_RESULT(no)
  AC_MSG_ERROR([upgrade GDAL to 2.1.0 or later])
else
  AC_MSG_RESULT(yes)
fi
use_gdal=0


CPPFLAGS="${PKG_CXXFLAGS} ${CPPFLAGS} ${GDAL_CPPFLAGS}"
LIBS="${PKG_LIBS} ${LIBS} ${GDAL_LIBS}"
AC_MSG_NOTICE([CPP flags: ${CPPFLAGS}])

PKG_CPPFLAGS="${PKG_CPPFLAGS} ${GDAL_CPPFLAGS}"
CPPFLAGS="${CPPFLAGS} ${GDAL_CPPFLAGS}"
PKG_LIBS="${PKG_LIBS} ${GDAL_LIBS}"
LIBS="${LIBS} ${GDAL_LIBS}"

GDAL_DATA_TEST_FILE="${GDAL_DATADIR}/pcs.csv"
AC_MSG_CHECKING(GDAL: ${GDAL_DATADIR}/pcs.csv readable)
if test -r "${GDAL_DATA_TEST_FILE}" ; then
  AC_MSG_RESULT(yes)
else
  AC_MSG_RESULT(no)
  AC_MSG_ERROR([pcs.csv not found in GDAL data directory.])
fi
AC_CHECK_HEADERS([gdal_priv.h], [use_gdal=1], AC_MSG_ERROR([Ensure gdal is properly installed - cannot find gdal_priv.h]))
AC_CHECK_HEADERS([gdal.h], [use_gdal=1], AC_MSG_ERROR([Ensure gdal is properly installed - cannot find gdal_priv.h]))
defines="$defines -Dwith_gdal -Dsql_ram -Dwith_rcoalescence"
AC_SUBST(defines)
dnl Create the output variables
AC_SUBST(PKG_CXXFLAGS)
AC_SUBST(PKG_CPPFLAGS)
AC_SUBST(PKG_LIBS)
dnl Process Makevars.in to create Makevars
AC_CONFIG_FILES([src/Makevars])
AC_OUTPUT
