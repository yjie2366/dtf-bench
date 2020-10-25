#include "util.h"

#define LVARS_ANAL 	11
#define LVARS_HIST 	20

static char *anal_vars[LVARS_ANAL] = {
	"DENS", "MOMX", "MOMY", "MOMZ", "RHOT", "QV",
	"QC", "QR", "QI", "QS", "QG"
};

static char *hist_vars[LVARS_HIST] = {
	"U", "V", "W", "T", "PRES", "QV", "QC",
	"QR", "QI", "QS", "QG", "RH", "height",
	"topo", "SFC_PRES", "PREC", "U10", "V10", "T2", "Q2"
};

int read_hist(PD *pd, char *dir_path, int cycle)
{
	int i, ncid, first_run = 0;
	int ret = 0;
	char file_path[MAX_PATH_LEN] = { 0 };
	struct file_info *file = &pd->files[HIST];
	struct data_buf *arrays = file->var_read_buffers;
	int num_var = LVARS_HIST;
	
	if (!arrays) {
		arrays = (struct data_buf *)malloc(sizeof(struct data_buf) * num_var);
		check_error(arrays, malloc);
		memset(arrays, 0, sizeof(struct data_buf) * num_var);
		file->nvar_read_buf = num_var;
		first_run = 1;
	}

	fmt_filename(cycle, pd->ens_id, 6, dir_path, ".hist.nc", file_path);
	prepare_file(file, pd->ens_comm, file_path, FILE_OPEN_R, &ncid);

	if (cycle) dtf_time_start();

	for (i = 0; i < num_var; i++) {
		struct data_buf *array = &arrays[i];
		int varid = array->varid, ndims = array->ndims;
		struct var_pair *var = NULL;
		MPI_Offset *start, *count;
		float *data = array->data;

		if (first_run) {
			ret = ncmpi_inq_varid(ncid, hist_vars[i], &varid);
			check_io(ret, ncmpi_inq_varid);
		}

		var = &file->vars[varid];
		if (unlikely(strcmp(var->name, hist_vars[i]))) {
			int idx;
			fprintf(stderr, "[WARNING]: Unmatched var %s ID and idx\n", hist_vars[i]);
			idx = find_var(file, hist_vars[i]);
			check_error(idx >= NUM_AXIS_VARS, find_var);

			var = &file->vars[idx];
			array->varid = var->varid;
		}
		
		if (first_run) {
			int j;
			MPI_Offset total_count = 1;
			MPI_Offset *s_idx, *e_idx;

			ndims = var->ndims;
			
			start = (MPI_Offset *)malloc(sizeof(MPI_Offset) * ndims * 4);
			check_error(start, malloc);

			count = start + ndims;
			s_idx = count + ndims;
			e_idx = s_idx + ndims;

			/* History files have no halo */
			for (j = 0; j < ndims; j++) {
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
				s_idx[j] = 0;
				e_idx[j] = count[j];
				total_count *= count[j];
			}
			
			data = (float *)malloc(sizeof(float) * total_count);
			check_error(data, malloc);
			memset(data, 0, total_count * sizeof(float));
		
			array->data = data;
			array->shape = start;
			array->idxes = s_idx;
			array->ndims = ndims;
			array->varid = varid;
			array->nelems = total_count;
		}
		else {
			start = array->shape;
			count = start + ndims;
		}

		ret = ncmpi_iget_vara_float(ncid, varid, start, count, data, NULL);
		check_io(ret, ncmpi_iget_vara_float);
	}

	ret = ncmpi_wait_all(ncid, NC_REQ_ALL, NULL, NULL);
	check_io(ret, ncmpi_wait_all);

	ret = dtf_transfer(file_path, ncid);
	check_error(!ret, dtf_transfer);
	
	ret = ncmpi_close(ncid);
	check_io(ret, ncmpi_close);

	if (cycle) dtf_readtime_end();
	
//	if (!pd->world_rank) {
//		int idx;
//		for (idx = 0; idx < num_var; idx++) {
//			int j;
//			fprintf(stderr, "Variable: %s\n", hist_vars[idx]);
//			for (j = 0; j < arrays[idx].nelems; j++) {
//				fprintf(stderr, "%f ", arrays[idx].data[j]);
//			}
//			fprintf(stderr, "\n");
//		}
//	}

	file->var_read_buffers = arrays;

	/* Check validity of read values */
	for (i = 0; i < num_var; i++) {
		struct data_buf *read_buf = &file->var_read_buffers[i];

	//	ret = compare_buffer(pd, read_buf, cycle, SCALE_WEIGHT);
	//	check_error(!ret, compare_buffer);
	}

	return ret;
}

int read_anal(PD *pd, char *dir_path, int cycle)
{
	int i, ncid, first_run = 0;
	int ret = 0;
	char file_path[MAX_PATH_LEN] = { 0 };
	struct file_info *file = &pd->files[ANAL];
	struct data_buf *arrays = file->var_read_buffers;
	int num_var = LVARS_ANAL;

	fmt_filename(cycle, pd->ens_id, 6, dir_path, ".anal.nc", file_path);
	prepare_file(file, pd->ens_comm, file_path, FILE_OPEN_R, &ncid);

	if (!arrays) {
		arrays = (struct data_buf *)malloc(sizeof(struct data_buf) * num_var);
		check_error(arrays, malloc);
		memset(arrays, 0, sizeof(struct data_buf) * num_var);
		file->nvar_read_buf = num_var;
		first_run = 1;
	}

	if (cycle) dtf_time_start();

	for (i = 0; i < num_var; i++) {
		int j;
		struct data_buf *array = &arrays[i];
		int varid = array->varid;
		int ndims = array->ndims;
		float *data = array->data;

		struct var_pair *var = NULL;
		MPI_Offset *start, *count;

		if (first_run) {
			ret = ncmpi_inq_varid(ncid, anal_vars[i], &varid);
			check_io(ret, ncmpi_inq_varid);
		}

		var = &file->vars[varid];
		if (unlikely(strcmp(var->name, anal_vars[i]))) {
			int idx;
			fprintf(stderr, "[WARNING]: Unmatched var %s ID and idx\n", anal_vars[i]);
			idx = find_var(file, anal_vars[i]);
			check_error(idx >= NUM_AXIS_VARS, find_var);
			
			var = &file->vars[idx];
			array->varid = var->varid;
		}

		if (first_run) {
			MPI_Offset *s_idx, *e_idx;
			MPI_Offset total_count = 1;

			ndims = var->ndims;
			start = (MPI_Offset *)malloc(sizeof(MPI_Offset) * ndims * 4);
			check_error(start, malloc);
			
			count = start + ndims;
			s_idx = count + ndims;
			e_idx = s_idx + ndims;

			for (j = 0; j < ndims; j++) {
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
				s_idx[j] = 0;
				e_idx[j] = count[j];
			}

			data = (float *)malloc(sizeof(float) * total_count);
			check_error(data, malloc);
			memset(data, -1, sizeof(float) * total_count);

			array->ndims = ndims;
			array->shape = start;
			array->idxes = s_idx;
			array->varid = varid;
			array->data = data;
			array->nelems = total_count;
		}
		else {
			start = array->shape;
			count = start + ndims;
		}

		ret = ncmpi_iget_vara_float(ncid, var->varid, start, count, data, NULL);
		check_io(ret, ncmpi_iget_vara_float);
	}

	ret = ncmpi_wait_all(ncid, NC_REQ_ALL, NULL, NULL);
	check_io(ret, ncmpi_wait_all);

	ret = dtf_transfer(file_path, ncid);
	check_error(!ret, dtf_transfer);

	ret = ncmpi_close(ncid);
	check_io(ret, ncmpi_close);

	if (cycle) dtf_readtime_end();
	
	file->var_read_buffers = arrays;

	// Check validity of read values
	for (i = 0; i < num_var; i++) {
		struct data_buf *read_buf = &file->var_read_buffers[i];
		
		ret = compare_buffer(pd, read_buf, cycle, SCALE_WEIGHT);
		check_error(!ret, compare_buffer);
	}

	return ret;
}

int write_anal(PD *pd, char *dir_path, int cycle)
{
	int ncid = -1;
	int i, ret = 0, first_run = 0;
	int num_var = LVARS_ANAL;
	char file_path[MAX_PATH_LEN] = { 0 };
	struct file_info *file = &pd->files[ANAL];
	struct data_buf *arrays = file->var_write_buffers;

	fmt_filename(cycle, pd->ens_id, 6, dir_path, ".anal.nc", file_path);
	prepare_file(file, pd->ens_comm, file_path, FILE_OPEN_W, &ncid);

	if (!arrays) {
		arrays = (struct data_buf *)malloc(sizeof(struct data_buf) * num_var);
		check_error(arrays, malloc);
		memset(arrays, 0, sizeof(struct data_buf) * num_var);
		file->nvar_write_buf = num_var;
		first_run = 1;
	}

	if (cycle) dtf_time_start();

	for (i = 0; i < num_var; i++) {
		struct var_pair *var = NULL;
		struct data_buf *array = &arrays[i];
		int varid = array->varid;
		int ndims = array->ndims;
		
		MPI_Offset *start, *count;
		float *data = array->data;

		if (first_run) {
			ret = ncmpi_inq_varid(ncid, anal_vars[i], &varid);
			check_io(ret, ncmpi_inq_varid);
		}
		
		var = &file->vars[varid];
		if (unlikely(strcmp(var->name, anal_vars[i]))) {
			int idx;
			fprintf(stderr, "[WARNING]: Unmatched var %s ID and idx\n", anal_vars[i]);
			idx = find_var(file, anal_vars[i]);
			check_error(idx >= NUM_AXIS_VARS, find_var);
			
			var = &file->vars[idx];
			array->varid = var->varid;
		}

		if (first_run) {
			int j;
			MPI_Offset *s_idx, *e_idx;
			MPI_Offset total_count = 1;

			ndims = var->ndims;
			start = (MPI_Offset *)malloc(sizeof(MPI_Offset) * ndims * 4);
			check_error(start, malloc);
			
			count = start + ndims;
			s_idx = count + ndims;
			e_idx = s_idx + ndims;

			for (j = 0; j < ndims; j++) {
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
				s_idx[j] = 0;
				e_idx[j] = count[j];
			}

			data = (float *)malloc(sizeof(float) * total_count);
			check_error(data, malloc);
			memset(data, -1, sizeof(float) * total_count);
			
			array->ndims = ndims;
			array->shape = start;
			array->idxes = s_idx;
			array->varid = varid;
			array->data = data;
			array->nelems = total_count;
		}
		else {
			start = array->shape;
			count = start + ndims;
		}

		ret = fill_buffer(array, (float)pd->world_rank, (float)cycle, LETKF_WEIGHT);
		check_error(!ret, fill_buffer);
		
		ret = ncmpi_iput_vara_float(ncid, var->varid, start, count, data, NULL);
		check_io(ret, ncmpi_iput_vara_float);
	}

	ret = ncmpi_wait_all(ncid, NC_REQ_ALL, NULL, NULL);
	check_io(ret, ncmpi_wait_all);
	
	ret = dtf_transfer(file_path, ncid);
	check_error(!ret, dtf_transfer);

	ret = ncmpi_close(ncid);
	check_io(ret, ncmpi_close);

	if (cycle) dtf_writetime_end();
	
	file->var_write_buffers = arrays;

	return ret;
}
