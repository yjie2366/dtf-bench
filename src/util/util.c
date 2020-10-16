#include <string.h>
#include <stdlib.h>
#include <sys/stat.h>

#include "util.h"

extern struct dim_pair anal_dims[];
extern struct dim_pair hist_dims[];

int create_dirs(char *path)
{
	int offset = 0;
	char tmp[MAX_PATH_LEN] = { 0 };
	char *pt = NULL;
	
	if (strlen(path) >= MAX_PATH_LEN) {
		fprintf(stderr, "[ERROR] %s Path %s is too long to process.\n", __func__, path);
		return EINVAL;
	}
	
	memcpy(tmp, path, strlen(path)+1);

	// Ignore leading slashes
	while (tmp[offset] == '/') offset++;
	for (pt = tmp + offset; *pt != '\0'; pt++) {
		if (*pt != '/') continue;
		
		*pt = '\0';
		if (mkdir(tmp, 0750)) {
			if (errno != EEXIST) {
				fprintf(stderr, "[ERROR] Failed to create folder %s\n", tmp);
				return errno;
			}
		}
		*pt = '/';
	}
	
	return 0;
}

void fmt_filename(int cycle, int id, int total_chrs, char *prefix, char *suffix, char *name )
{
	char tmp[MAX_NAME_LEN] = { 0 };
	char digits[MAX_NAME_LEN] = { 0 };
	int num_digits = 0;
	int i, off; 
	
	if (prefix) {
		int l = strlen(prefix);
		memcpy(tmp, prefix, l);
		off = l;
	}

	int size = sprintf(tmp+off, "%d/", cycle);
	off += size;

	do {
		digits[num_digits++] = '0' + id % 10;
	} while (id /= 10);

	for (i = 0; i < (total_chrs - num_digits); i++) tmp[off++] = '0';
	for (i = num_digits-1; i >= 0; i--) tmp[off++] = digits[i];

	if (suffix) {
		int l = strlen(suffix);
		memcpy(tmp+off, suffix, l);
		off += l;
	}

	memcpy(name, tmp, strlen(tmp));
}

MPI_Offset get_databuf_size(PD *pd, int file_idx)
{
	struct file_info *file = NULL;
	MPI_Offset total_size = 0;
	int idx;

	if (!pd->files) {
		fprintf(stderr, "[ERROR]: files info are not filled yet\n");
		return -1;
	}
	
	file = &pd->files[file_idx];
	if (!file->dims || !file->vars) {
		fprintf(stderr, "[ERROR]: dims and vars are not filled yet\n");
		return -1;
	}

	idx = (file_idx == ANAL)? ANAL_DATA_VARS_OFFSET : HIST_DATA_VARS_OFFSET;
	for (; idx < file->nvars; idx++) {
		int j;
		struct var_pair *var = &file->vars[idx];
		MPI_Offset total_dim_size = 1;

		for (j = 0; j < var->ndims; j++) {
			MPI_Offset dim_size = 0;
			char *name = var->dim_name[j];

			if (!(strcmp(name, "x")) || !(strcmp(name, "xh"))) {
				dim_size = (file_idx == ANAL)? IA(pd) : IMAX(pd);
			}
			else if (!(strcmp(name, "y")) || !(strcmp(name, "yh"))) {
				dim_size = (file_idx == ANAL)? JA(pd) : JMAX(pd);
			}
			else if (!(strcmp(name, "z"))) {
				dim_size = (file_idx == ANAL)? KA : KMAX;
			}
			else if (!(strcmp(name, "zh"))) {
				dim_size = KA + 1;
			}
			else if (!(strcmp(name, "oz"))) {
				dim_size = OKMAX;
			}
			else if (!(strcmp(name, "lz"))) {
				dim_size = LKMAX;
			}
			else if (!(strcmp(name, "uz"))) {
				dim_size = UKMAX;
			}
			else if (!(strcmp(name, "time"))) {
				dim_size = 1;
			}
			total_dim_size *= dim_size;
		}
		total_size += total_dim_size;
	}

	return total_size * sizeof(float);
}

int prepare_file(struct file_info *file, MPI_Comm comm, char *file_path, int flag, int *ncid)
{
	int i, ret, fid;
	
	if (flag & FILE_CREATE) {
		ret = ncmpi_create(comm, file_path, NC_CLOBBER | NC_64BIT_OFFSET,
					MPI_INFO_NULL, &fid);
		check_io(ret, ncmpi_create);
	}
	else {
		int mode = (flag & FILE_OPEN_W) ? NC_WRITE : NC_NOWRITE;

		ret = ncmpi_open(comm, file_path, mode, MPI_INFO_NULL, &fid);
		check_io(ret, ncmpi_open);
	}
	
	for (i = 0; i < file->ndims; i++) {
		struct dim_pair *dim = &file->dims[i];

		if (flag & FILE_CREATE) {
			ret = ncmpi_def_dim(fid, dim->name, dim->length, &dim->dimid);
			check_io(ret, ncmpi_def_dim);
		}
		else {
			ret = ncmpi_inq_dimid(fid, dim->name, &dim->dimid);
			check_io(ret, ncmpi_inq_dimid);
		}
	}

	for (i = 0; i < file->nvars; i++) {
		int j;
		struct var_pair *var = &file->vars[i];

		for (j = 0; j < var->ndims; j++) {
			ret = ncmpi_inq_dimid(fid, var->dim_name[j], &var->dims[j]);
			check_io(ret, ncmpi_inq_dimid);
		}

		if (flag & FILE_CREATE) {
			ret = ncmpi_def_var(fid, var->name, var->type, var->ndims,
					var->dims, &var->varid);
			check_io(ret, ncmpi_def_var);
			int rank = -1;
			MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		}
		else {
			ret = ncmpi_inq_varid(fid, var->name, &var->varid);
			check_io(ret, ncmpi_inq_varid);
		}
	}

	*ncid = fid;

	return 0;
}
