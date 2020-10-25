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
#include "constants.h"
#include "dtf.h"

#ifdef __GNUC__
#define likely(x)       __builtin_expect(!!(x), 1)
#define unlikely(x)     __builtin_expect(!!(x), 0)
#else
#define likely(x)       (x)
#define unlikely(x)     (x)
#endif

#define check_error(exp, func) do { \
	if (unlikely(!(exp))) { \
		fprintf(stderr, #func" at line %d in file %s [ERROR]: %s\n", \
			__LINE__, __FILE__, (char *)strerror(errno)); \
		MPI_Abort(MPI_COMM_WORLD, -1);				\
	}								\
} while (0)

#define check_mpi(errcode, func) do {						\
	if (unlikely(errcode != MPI_SUCCESS)) {					\
		char error_msg[MPI_MAX_ERROR_STRING];			\
		int str_len;						\
		MPI_Error_string((errcode), error_msg, &str_len);		\
		fprintf(stderr, #func" at line %d in file %s [MPI ERROR]: %s\n", \
			__LINE__, __FILE__, error_msg); \
		MPI_Abort(MPI_COMM_WORLD, (errcode));			\
	}								\
} while (0)

#define check_io(errcode, func) do {						\
	if (unlikely(errcode != NC_NOERR)) {					\
		fprintf(stderr, #func" at line %d in file %s [IO ERROR]: %s\n", \
			__LINE__, __FILE__, ncmpi_strerror(errcode)); \
		MPI_Abort(MPI_COMM_WORLD, (errcode));			\
	}								\
} while (0)

#define TOLOWER(str) ({							\
	char *tmp = str;						\
	for (; *tmp != '\0'; tmp++) *tmp = tolower(*tmp);		\
	str;})

void init_pd(int argc, char **argv, PD *pd);
int finalize_pd(PD *pd);

#endif // __UTIL_H_
