AC_INIT([promising], [0.1], [kelsey.anderson@ucdenver.edu])
AM_INIT_AUTOMAKE([-Wall -Werror foreign])
AC_PROG_CXX
##AC_CONFIG_HEADERS([config.h])
AC_CHECK_HEADERS([iostream])
AC_SEARCH_LIBS([sgetrf_], [lapack], [], [AC_MSG_ERROR([unable to find lapack])])
AC_SEARCH_LIBS([sgetri_], [lapack], [], [AC_MSG_ERROR([unable to find lapack])])
AC_CONFIG_FILES([
 Makefile
 src/Makefile
])
AC_OUTPUT
