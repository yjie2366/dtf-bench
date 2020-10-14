#ifndef __UTIL_H_
#define __UTIL_H_

#include <stdio.h>
#include <mpi.h>
#include <errno.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <ctype.h>
#include <pnetcdf.h>
#include "bpconf.h"
#include "iobench.h"
#include "dtf.h"

#define check_error(exp, func) do { \
	if (!(exp)) { \
		fprintf(stderr, #func" at line %d in file %s [ERROR]: %s\n", \
			__LINE__, __FILE__, (char *)strerror(errno)); \
		MPI_Abort(MPI_COMM_WORLD, -1);				\
	}								\
} while (0)

#define check_mpi(errcode, func) do {						\
	if (errcode != MPI_SUCCESS) {					\
		char error_msg[MPI_MAX_ERROR_STRING];			\
		int str_len;						\
		MPI_Error_string((errcode), error_msg, &str_len);		\
		fprintf(stderr, #func" at line %d in file %s [MPI ERROR]: %s\n", \
			__LINE__, __FILE__, error_msg); \
		MPI_Abort(MPI_COMM_WORLD, (errcode));			\
	}								\
} while (0)

#define check_io(errcode, func) do {						\
	if (errcode != NC_NOERR) {					\
		fprintf(stderr, #func" at line %d in file %s [IO ERROR]: %s\n", \
			__LINE__, __FILE__, ncmpi_strerror(errcode)); \
		MPI_Abort(MPI_COMM_WORLD, (errcode));			\
	}								\
} while (0)

#define TOLOWER(str) ({							\
	char *tmp = str;						\
	for (; *tmp != '\0'; tmp++) *tmp = tolower(*tmp);		\
	str;})

#endif // __UTIL_H_
