#include "util.h"

#define LVARS_ANAL 	11
#define LVARS_HIST 	20
#define LETKF_WEIGHT 	10

static char *anal_vars[LVARS_ANAL] = {
	"DENS", "MOMX", "MOMY", "MOMZ", "RHOT", "QV",
	"QC", "QR", "QI", "QS", "QG"
};

static char *hist_vars[LVARS_HIST] = {
	"U", "V", "W", "T", "PRES", "QV", "QC",
	"QR", "QI", "QS", "QG", "RH", "height",
	"topo", "SFC_PRES", "PREC", "U10", "V10", "T2", "Q2"
};

static inline int find_var(PD *pd, int file_idx, char *var_name)
{
	int idx = -1, j;
	int offset = NUM_AXIS_VARS;
	struct file_info *file = &pd->files[file_idx];

	for (j = offset; j < file->nvars; j++) {
		if (strcmp(file->vars[j].name, var_name)) continue;
		idx = j; 
		break;
	}

	return idx;
}

int read_hist(PD *pd, char *dir_path, int cycle)
{
	int i, ncid;
	int ret = 0;
	char file_path[MAX_PATH_LEN] = { 0 };
	struct file_info *file = &pd->files[HIST];
	float **var_arrays = file->var_read_buffers;
	int num_var = LVARS_HIST;
	
	fmt_filename(cycle, pd->ens_id, 6, dir_path, ".hist.nc", file_path);

	prepare_file(file, pd->ens_comm, file_path, FILE_OPEN_R, &ncid);

	if (!var_arrays) {
		var_arrays = (float **)malloc(sizeof(float *) * num_var);
		check_error(var_arrays, malloc);
		memset(var_arrays, 0, sizeof(float *) * num_var);
		file->nvar_read_buf = num_var;
	}

	if (cycle) dtf_time_start();

	for (i = 0; i < num_var; i++) {
		int j;
		int idx = -1;
		MPI_Offset total_count = 1;
		struct var_pair *var = NULL;
		MPI_Offset *start = NULL;
		MPI_Offset *count = NULL;
		float *array = NULL;

		ret = ncmpi_inq_varid(ncid, hist_vars[i], &idx);
		check_io(ret, ncmpi_inq_varid);
		var = &file->vars[idx];

		if (strcmp(var->name, hist_vars[i])) {
			fprintf(stderr, "[WARNING]: Unmatched var %s ID and idx\n", hist_vars[i]);
			idx = find_var(pd, HIST, hist_vars[i]);
			check_error(idx >= 0, find_var);
			var = &file->vars[idx];
		}

		start = (MPI_Offset *)malloc(sizeof(MPI_Offset) * var->ndims * 2);
		check_error(start, malloc);
		count = start + var->ndims;

		/* History files have no halo */
		for (j = 0; j < var->ndims; j++) {
			if (strchr(var->dim_name[j], 'z')) {
				start[j] = 0;
				count[j] = KMAX;
			}
			else if (strchr(var->dim_name[j], 'y')) {
				start[j] = pd->proc_rank_y * JMAX(pd);
				count[j] = JMAX(pd);
			}
			else if (strchr(var->dim_name[j], 'x')) {
				start[j] = pd->proc_rank_x * IMAX(pd);
				count[j] = IMAX(pd);
			}
			else if (strstr(var->dim_name[j], "time")) {
				start[j] = 0;
				count[j] = 1;
			}
			total_count *= count[j];
		}
		
		array = var_arrays[i];
		if (!array) {
			array = (float *)malloc(sizeof(float) * total_count);
			check_error(array, malloc);
			memset(array, 0, total_count * sizeof(float));
		}

		ret = ncmpi_iget_vara_float(ncid, var->varid, start, count, array, NULL);
		check_io(ret, ncmpi_iget_vara_float);

		var_arrays[i] = array;
		free(start);
	
	}

	ret = ncmpi_wait_all(ncid, NC_REQ_ALL, NULL, NULL);
	check_io(ret, ncmpi_wait_all);

	ret = dtf_transfer(file_path, ncid);
	check_error(!ret, dtf_transfer);
	
	file->var_read_buffers = var_arrays;

	ret = ncmpi_close(ncid);
	check_io(ret, ncmpi_close);

	if (cycle) dtf_time_end();
	
	MPI_Barrier(pd->ens_comm);

	return ret;
}

int read_anal(PD *pd, char *dir_path, int cycle)
{
	int i, ncid;
	int ret = 0;
	char file_path[MAX_PATH_LEN] = { 0 };
	struct file_info *file = &pd->files[ANAL];
	float **var_arrays = file->var_read_buffers;

	fmt_filename(cycle, pd->ens_id, 6, dir_path, ".anal.nc", file_path);
	
	prepare_file(file, pd->ens_comm, file_path, FILE_OPEN_R, &ncid);

	if (!var_arrays) {
		var_arrays = (float **)malloc(sizeof(float *) * LVARS_ANAL);
		check_error(var_arrays, malloc);
		memset(var_arrays, 0, sizeof(float *) * LVARS_ANAL);
		file->nvar_read_buf = LVARS_ANAL;
	}

	if (cycle) dtf_time_start();

	for (i = 0; i < LVARS_ANAL; i++) {
		int j;
		int idx = -1;
		MPI_Offset total_count = 1;
		struct var_pair *var = NULL;
		MPI_Offset *start = NULL;
		MPI_Offset *count = NULL;
		float *array = NULL;

		ret = ncmpi_inq_varid(ncid, anal_vars[i], &idx);
		check_io(ret, ncmpi_inq_varid);
		var = &file->vars[idx];

		if (strcmp(var->name, anal_vars[i])) {
			fprintf(stderr, "[WARNING]: Unmatched var %s ID and idx\n", anal_vars[i]);
			idx = find_var(pd, ANAL, anal_vars[i]);
			check_error(idx >= 0, find_var);
			var = &file->vars[idx];
		}

		start = (MPI_Offset *)malloc(sizeof(MPI_Offset) * var->ndims * 2);
		check_error(start, malloc);
		count = start + var->ndims;

		for (j = 0; j < var->ndims; j++) {
			if (strchr(var->dim_name[j], 'z')) {
				start[j] = 0;
				count[j] = KMAX;
			}
			else if (strchr(var->dim_name[j], 'y')) {
				start[j] = pd->proc_rank_y * JMAX(pd) + JHALO;
				count[j] = JMAX(pd);
			}
			else if (strchr(var->dim_name[j], 'x')) {
				start[j] = pd->proc_rank_x * IMAX(pd) + IHALO;
				count[j] = IMAX(pd);
			}
			total_count *= count[j];
		}

		array = var_arrays[i];
		if (!array) {
			array = (float *)malloc(sizeof(float) * total_count);
			check_error(array, malloc);
		}
		
		ret = ncmpi_iget_vara_float(ncid, var->varid, start, count, array, NULL);
		check_io(ret, ncmpi_iget_vara_float);

		var_arrays[i] = array;
		free(start);

	}
	ret = ncmpi_wait_all(ncid, NC_REQ_ALL, NULL, NULL);
	check_io(ret, ncmpi_wait_all);

	ret = dtf_transfer(file_path, ncid);
	check_error(!ret, dtf_transfer);

	file->var_read_buffers = var_arrays;

	ret = ncmpi_close(ncid);
	check_io(ret, ncmpi_close);

	if (cycle) dtf_time_end();

	MPI_Barrier(pd->ens_comm);

	return ret;
}

int write_anal(PD *pd, char *dir_path, int cycle)
{
	int ret = 0, i;
	int ncid = -1;
	char file_path[MAX_PATH_LEN] = { 0 };
	struct file_info *file = &pd->files[ANAL];
	float **var_arrays = file->var_write_buffers;

	fmt_filename(cycle, pd->ens_id, 6, dir_path, ".anal.nc", file_path);

	prepare_file(file, pd->ens_comm, file_path, FILE_OPEN_W, &ncid);

	if (!var_arrays) {
		var_arrays = (float **)malloc(sizeof(float *) * LVARS_ANAL);
		check_error(var_arrays, malloc);
		memset(var_arrays, 0, sizeof(float *) * LVARS_ANAL);
		file->nvar_write_buf = LVARS_ANAL;
	}

	if (cycle) dtf_time_start();

	for (i = 0; i < LVARS_ANAL; i++) {
		int j;
		int idx = -1;
		MPI_Offset total_count = 1;
		struct var_pair *var = NULL;
		MPI_Offset *start = NULL;
		MPI_Offset *count = NULL;
		float *array = NULL;

		ret = ncmpi_inq_varid(ncid, anal_vars[i], &idx);
		check_io(ret, ncmpi_inq_varid);
		var = &file->vars[idx];

		if (strcmp(var->name, anal_vars[i])) {
			fprintf(stderr, "[WARNING]: Unmatched var %s ID and idx\n", anal_vars[i]);
			idx = find_var(pd, ANAL, anal_vars[i]);
			check_error(idx >= 0, find_var);
			var = &file->vars[idx];
		}

		start = (MPI_Offset *)malloc(sizeof(MPI_Offset) * var->ndims * 2);
		check_error(start, malloc);
		count = start + var->ndims;

		for (j = 0; j < var->ndims; j++) {
			if (strchr(var->dim_name[j], 'z')) {
				start[j] = 0;
				count[j] = KMAX;
			}
			else if (strchr(var->dim_name[j], 'y')) {
				start[j] = pd->proc_rank_y * JMAX(pd) + JHALO;
				count[j] = JMAX(pd);
			}
			else if (strchr(var->dim_name[j], 'x')) {
				start[j] = pd->proc_rank_x * IMAX(pd) + IHALO;
				count[j] = IMAX(pd);
			}
			total_count *= count[j];
		}

		array = var_arrays[i];
		if (!array) {
			array = (float *)malloc(sizeof(float) * total_count);
			check_error(array, malloc);
		}
		
		for (j = 0; j < total_count; j++) {
			array[j] = (float)((pd->ens_id - 1 + j) * (pd->world_rank + 1)) * LETKF_WEIGHT;
		}

		ret = ncmpi_iput_vara_float(ncid, var->varid, start, count, array, NULL);
		check_io(ret, ncmpi_iput_vara_float);

		var_arrays[i] = array;
		free(start);
	}

	ret = ncmpi_wait_all(ncid, NC_REQ_ALL, NULL, NULL);
	check_io(ret, ncmpi_wait_all);
	
	ret = dtf_transfer(file_path, ncid);
	check_error(!ret, dtf_transfer);

	file->var_write_buffers = var_arrays;

	ret = ncmpi_close(ncid);
	check_io(ret, ncmpi_close);

	if (cycle) dtf_time_end();

	MPI_Barrier(pd->ens_comm);

	return ret;
}
