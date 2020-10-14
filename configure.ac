AC_INIT([dtf-bench], [1.0])

AC_ARG_VAR(MPICC, MPI C Compiler)
AC_ARG_VAR(CC, Overwritten by MPICC if MPICC is set)
if test "x$MPICC"  = x && test "x$CC"  != x ; then MPICC=$CC; fi
if test "x$MPICC" = x && test "x$CC" = x; then MPICC=mpiicc; fi

CC=${MPICC}
AC_SUBST(MPICC)
AC_PROG_CC

AC_ARG_WITH([pnetcdf], AC_HELP_STRING([--with-pnetcdf=<PNETCDF-INSTALL-PATH>], [To Specify PnetCDF Install Path]),,)
if test "x${with_pnetcdf}" != "x" && test -d ${with_pnetcdf}; then
	pnetcdf_ins=${with_pnetcdf}
else
	AC_MSG_ERROR([PnetCDF Library Missing])
fi

AC_ARG_WITH([dtf], AC_HELP_STRING([--with-dtf=<DTF-INSTALL-PATH>], [To Specify DTF Install Path]),,)
if test "x${with_dtf}" != "x" && test -d ${with_dtf}; then
	dtf_ins=${with_dtf}
else
	AC_MSG_ERROR([DTF Library Missing])
fi

AC_ARG_WITH([split], AC_HELP_STRING([--with-split=<SPLITWRAPPER-PATH>], [To Specify SPLIT-WRAPPER Path]),,)
if test "x${with_split}" != "x" && test -e ${with_split}; then
	split_wrapper=${with_split}
else
	split_wrapper=
fi

AC_SUBST(split_wrapper)
AC_SUBST(pnetcdf_ins)
AC_SUBST(dtf_ins)

AC_CONFIG_FILES([Makefile])
AC_CONFIG_FILES([src/Makefile])
AC_CONFIG_FILES([script/exec.sh])
AC_CONFIG_FILES([src/include/bpconf.h])
AC_OUTPUT
