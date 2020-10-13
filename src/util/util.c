#include <string.h>
#include <stdlib.h>
#include <sys/stat.h>

#include "util.h"

extern struct dim_pair anal_dims[];
extern struct dim_pair hist_dims[];

inline char *get_type(nc_type type) {
	switch(type) {
		case NC_BYTE:
			return "BYTE";
		case NC_CHAR:
			return "CHAR";
		case NC_SHORT:
			return "SHORT";
		case NC_INT:
			return "INT";
		case NC_FLOAT:
			return "FLOAT";
		case NC_DOUBLE:
			return "DOUBLE";
		case NC_NAT:
		default:
			return NULL;
	}
}

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

void print_structs(PD *pd)
{
	int i, j;
	int num_files = pd->nfiles;

	printf("Process info: \n\n"
		"IMAX: %ld  JMAX: %ld  KMAX: %ld\n"
		"Number of Ensembles: %d\n"
		"Number of Cycles: %d\n"
		"Universal Process Rank: %d\n"
		"Ensemble Process Rank: %d\n"
		"Size: %d(ens) / %d(world)\n"
		"Ensemble ID: %d\n"
		"Number of Processes on X COORD: %d\n"
		"Numbre of Processes on Y COORD: %d\n"
		"\n",
		pd->imax, pd->jmax, pd->kmax,
		pd->num_ens, pd->cycles, pd->world_rank, pd->ens_rank,
		pd->ens_size, pd->world_size,
		pd->ens_id, pd->proc_num_x, pd->proc_num_y);

	printf("Dimension of anal: \n\n");
	for (i = 0; i < NUM_ANALDIMS; i++) {
		printf("name: %s length: %ld\n", anal_dims[i].name, anal_dims[i].length);
	}
	
	printf("\nDimension of hist: \n\n");
	for (i = 0; i < NUM_HISTDIMS; i++) {
		printf("name: %s length: %ld\n", hist_dims[i].name, hist_dims[i].length);
	}

	for (i = 0; i < num_files; i++) {
		int num_vars = (i == ANAL) ? NUM_ANALVARS : NUM_HISTVARS;
		printf("Variables of File %s:\n", (i == ANAL)? "anal" : "hist");

		for (j = 0; j < num_vars; j++) {
			int d;
			struct var_pair *var = &pd->files[i].vars[j];
			
			printf("name: %s type: %s ndims: %d dims: ",
				var->name, get_type(var->type), var->ndims);
			for (d = 0; d < var->ndims; d++) {
				printf("%s ", var->dim_name[d]);
			}
			printf("\n");
		}
	}
}

MPI_Offset get_databuf_size(PD *pd, int file_idx)
{
	MPI_Offset ret = 0;
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
