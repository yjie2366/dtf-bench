#include "util.h"
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <time.h>

extern MPI_Datatype subarray_type[];

static int random_float_generator(float *array, MPI_Offset num)
{
	int i;

	if (!array) {
		fprintf(stderr, "[ERROR]: array argument cannot be NULL\n");
		return -EINVAL;
	}
	
	srand48(time(NULL));
	for (i = 0; i < num; i++) {
		array[i] = (float)drand48() * (rand() % 1000);
	}

	return 0;
}

int get_subarray_type(struct var_pair *var)
{
	if (var->ndims == 2) return XY;

	if (strchr(var->dim_name[2], 'o')) return OCEAN;
	if (strchr(var->dim_name[2], 'l')) return LAND;
	if (strchr(var->dim_name[2], 'u')) return URBAN;
	if (strchr(var->dim_name[2], 'h')) return ZHXY2;
	 
	return ZXY2;
}

/* write axes and associated coord vars */
static int write_axis_vars(PD *pd, int file_idx, int ncid)
{
	int i, ret = 0;
	struct file_info *file = &pd->files[file_idx];
	int num_axis = NUM_AXIS_VARS;
	int num_coords = file->nassct_coords;
	int arr_count = 0;

	float **arrays = file->axes_buffer;
	if (!arrays) {
		arrays = (float **)malloc(sizeof(float *) * (num_axis + num_coords));
		check_error(arrays, malloc);
		memset(arrays, 0, sizeof(float *) * (num_axis + num_coords));

		file->naxes_buf = num_axis + num_coords;
	}

	/* Write to AXIS variables */
	for (i = 0; i < num_axis; i++) {
		struct var_pair *var = &file->vars[i];
		MPI_Offset dim_len = file->dims[var->dims[0]].length;
		MPI_Offset start = 0;
		MPI_Offset count = 0;
		float *array = arrays[i];
		int j;
		
		// y axes
		if (pd->proc_rank_x == 0) {
			if (strchr(var->dim_name[0], 'y')) {
				if (file_idx == HIST) {
					if (strchr(var->dim_name[0], 'h')) {
						start = SYH_hist(pd);
						count = CYH_hist(pd);
					}
					else {
						start = SY_hist(pd);
						count = JMAX(pd);
					}
				}
				else if (file_idx == ANAL) {
					start = JSGA(pd);
					count = JEB(pd) - JSB(pd) + 1;
				}

				if (!array) {
					array = (float *)malloc(sizeof(float) * count);
					check_error(array, malloc);
				}

				for (j = 0; j < count; j++) {
					array[j] = (float)(dim_len * (pd->ens_id - 1) + 
							count * pd->proc_rank_y + j);
				}
			}
		}

		// x axes
		if (pd->proc_rank_y == 0) {
			if (strchr(var->dim_name[0], 'x')) {
				if (file_idx == HIST) {
					if (strchr(var->dim_name[0], 'h')) {
						start = SXH_hist(pd);
						count = CXH_hist(pd);
					}
					else {
						start = SX_hist(pd);
						count = IMAX(pd);
					}
				} // ANAL
				else if (file_idx == ANAL) {
					start = ISGA(pd);
					count = IEB(pd) - ISB(pd) + 1;
				}

				if (!array) {
					array = (float *)malloc(sizeof(float) * count);
					check_error(array, malloc);
				}

				for (j = 0; j < count; j++) {
					array[j] = (float)(dim_len * (pd->ens_id - 1) + 
							count * pd->proc_rank_x + j);
				}
			}
		}		
		
		// z, Z, X, and Y axes
		if (pd->ens_rank == 0) {
			if ( strchr(var->dim_name[0], 'z') || strchr(var->dim_name[0], 'Z') ||
			strchr(var->dim_name[0], 'X') || strchr(var->dim_name[0], 'Y') ) {
				start = 0;
				count = dim_len;

				if (!array) {
					array = (float *)malloc(sizeof(float) * count);
					check_error(array, malloc);
				}

				for (j = 0; j < count; j++) {
					array[j] = (float)(dim_len * (pd->ens_id - 1) + j);
				}
			}
		}
		
		// put values to the variable
		if (count) {
			ret = ncmpi_iput_vara_float(ncid, var->varid, &start,
					&count, array, NULL);
			check_io(ret, ncmpi_iput_vara_float);
		
			arrays[i] = array;
			arr_count++;
		}
	}

	/* Write to AssociatedCoord variables */
	for (; i < (num_coords + num_axis); i++) {
		struct var_pair *var = &file->vars[i];
		int ndims = var->ndims;
		MPI_Offset *start = (MPI_Offset *)malloc(sizeof(MPI_Offset) * ndims * 2);
		check_error(start, malloc);

		MPI_Offset *count = start + ndims;
		float *array = arrays[i];
		MPI_Offset total_count = 1;
		int j;

		for (j = 0; j < ndims; j++) {
			MPI_Offset dim_len = file->dims[var->dims[j]].length;

			if (file_idx == ANAL) {
				switch (j) {
				case 0:
					start[j] = JSGA(pd);
					count[j] = JEB(pd) - JSB(pd) + 1;
					break;
				case 1:
					start[j] = ISGA(pd);
					count[j] = IEB(pd) - ISB(pd) + 1;
					break;
				case 2:
					start[j] = 0;
					count[j] = dim_len;
					break;
				default:
					fprintf(stderr, "[ERROR]: Invalid dim index\n");
					MPI_Abort(MPI_COMM_WORLD, EINVAL);
				}
			}
			// HIST
			else {
				if (strchr(var->dim_name[j], 'z')) {
					start[j] = 0;
					count[j] = dim_len;
				}
				else if (strchr(var->dim_name[j], 'y')) {
					if (strchr(var->dim_name[j], 'h')) {
						start[j] = SYH_hist(pd);
						count[j] = CYH_hist(pd);
					}
					else {
						start[j] = SY_hist(pd);
						count[j] = JMAX(pd);
					}
				}
				else if (strchr(var->dim_name[j], 'x')) {
					if (strchr(var->dim_name[j], 'h')) {
						start[j] = SXH_hist(pd);
						count[j] = CXH_hist(pd);
					}
					else {
						start[j] = SX_hist(pd);
						count[j] = IMAX(pd);
					}
				}
			}
			total_count *= count[j];
		}

		if (!array) {
			array = (float *)malloc(sizeof(float) * total_count);
			check_error(array, malloc);
		}

		for (j = 0; j < total_count; j++) {
			array[j] = (float)(total_count * (pd->ens_id - 1) + j);
		}
		
		ret = ncmpi_iput_vara_float(ncid, var->varid, start,
				count, array, NULL);
		check_io(ret, ncmpi_iput_vara_float);

		arrays[i] = array;
		arr_count++;
	
		if (start) free(start);	
	}

	ret = ncmpi_wait_all(ncid, NC_REQ_ALL, NULL, NULL);
	check_io(ret, ncmpi_wait_all);

	file->axes_buffer = arrays;

	return ret;
}

static int write_data_vars(PD *pd, int file_idx, int ncid, int cycle)
{
	int i, ret = 0;
	int offset = (file_idx == ANAL) ? ANAL_DATA_VARS_OFFSET : HIST_DATA_VARS_OFFSET;
	struct file_info *file = &pd->files[file_idx];
	
	float **var_arrays = file->var_write_buffers;
	if (!var_arrays) {
		var_arrays = (float **)malloc(sizeof(float *) * file->ndata_vars);
		check_error(var_arrays, malloc);
		memset(var_arrays, 0, sizeof(float *) * file->ndata_vars);

		file->nvar_write_buf = file->ndata_vars;
	}

	/* Write to data variables */
	for (i = offset; i < file->nvars; i++) {
		int j;
		int arr_idx = i - offset;
		struct var_pair *var = &file->vars[i];
		MPI_Offset total_count = 1;
		float *array = var_arrays[arr_idx];

		MPI_Offset *start = (MPI_Offset *)calloc(var->ndims * 2, sizeof(MPI_Offset));
		check_error(start, calloc);
		MPI_Offset *count = start + var->ndims;

		for (j = 0; j < var->ndims; j++) {
			if (!strcmp(var->dim_name[j], "time")) {
				start[j] = 0;
				count[j] = 1;
			}
			else if (strchr(var->dim_name[j], 'z')) { // contains z
				start[j] = 0;
				count[j] = file->dims[var->dims[j]].length;
			}
			else if (strchr(var->dim_name[j], 'y')) {
				if (file_idx == ANAL) {
					start[j] = JSGA(pd);
					count[j] = JEB(pd) - JSB(pd) + 1;
				}
				else {
					start[j] = SY_hist(pd);
					count[j] = JMAX(pd);
				}
			}
			else if (strchr(var->dim_name[j], 'x')) {
				if (file_idx == ANAL) {
					start[j] = ISGA(pd);
					count[j] = IEB(pd) - ISB(pd) + 1;
				}
				else {
					start[j] = SX_hist(pd);
					count[j] = IMAX(pd);
				}
			}
			total_count *= count[j];
		}

		if (!array) {
			array = (float *)malloc(sizeof(float) * total_count);
			check_error(array, malloc);
		}

		for (j = 0; j < total_count; j++) {
			array[j] = (float)(total_count * (pd->ens_id - 1) + 
					total_count * cycle + j);
		}

		ret = ncmpi_bput_vara_float(ncid, var->varid, start, count, array, NULL);
		check_io(ret, ncmpi_bput_vara_float);

		var_arrays[arr_idx] = array;

		free(start);
	}

	ret = ncmpi_wait_all(ncid, NC_REQ_ALL, NULL, NULL);
	check_io(ret, ncmpi_wait_all);
	
	file->var_write_buffers = var_arrays;
	return ret;
}

static int write_time_var(PD *pd, int ncid, int cycle)
{
	static double time = TIME_INT;
	double time_s, time_e;
	int ret = 0, varid_time, varid_bnds;
	MPI_Offset index[2] = { 0 };

	time *= cycle;
	time_s = time - TIME_INT;
	time_e = time;

	ret = ncmpi_inq_varid(ncid, "time", &varid_time);
	check_io(ret, ncmpi_inq_varid);

	ret = ncmpi_inq_varid(ncid, "time_bnds", &varid_bnds);
	check_io(ret, ncmpi_inq_varid);

	index[0] = 0;
	ret = ncmpi_put_var1_double_all(ncid, varid_time, index, &time);
	check_io(ret, ncmpi_put_var1_double_all);

	index[1] = 0;
	ret = ncmpi_put_var1_double_all(ncid, varid_bnds, index, &time_s);
	check_io(ret, ncmpi_put_var1_double_all);

	index[1] = 1;
	ret = ncmpi_put_var1_double_all(ncid, varid_bnds, index, &time_e);
	check_io(ret, ncmpi_put_var1_double_all);

	return ret;
}

int write_hist(PD *pd, char *dir_path, int cycle)
{
	int ret = 0;
	int ncid = -1;
	MPI_Offset buf_size = 0;
	char file_path[MAX_PATH_LEN] = { 0 };
	struct file_info *file = &pd->files[HIST];
	
	fmt_filename(cycle, pd->ens_id, 6, dir_path, ".hist.nc", file_path);
	
	ret = create_dirs(file_path);
	check_error(!ret, create_dirs);
	
	prepare_file(file, pd->ens_comm, file_path, FILE_CREATE, &ncid);
	
	buf_size = get_databuf_size(pd, HIST);
	check_error(buf_size>0, get_databuf_size);

	ret = ncmpi_buffer_attach(ncid, buf_size);
	check_io(ret, ncmpi_buffer_attach);

	ret = ncmpi_enddef(ncid);
	check_io(ret, ncmpi_enddef);
	
	if (cycle) dtf_time_start();

	write_axis_vars(pd, HIST, ncid);
	write_time_var(pd, ncid, cycle);
	write_data_vars(pd, HIST, ncid, cycle);

	ret = dtf_transfer(file_path, ncid);
	check_error(!ret, dtf_transfer);

	ret = ncmpi_buffer_detach(ncid);
	check_io(ret, ncmpi_buffer_detach);

	ret = ncmpi_close(ncid);
	check_io(ret, ncmpi_close);
	
	if (cycle) dtf_time_end();

	return ret;
}

int write_anal(PD *pd, char *dir_path, int cycle)
{
	int ret;
	int ncid = -1;
	MPI_Offset buf_size = 0;
	char file_path[MAX_PATH_LEN] = { 0 };
	struct file_info *file = &pd->files[ANAL];

	// Generate file path
	fmt_filename(cycle, pd->ens_id, 6, dir_path, ".anal.nc", file_path);

	ret = create_dirs(file_path);
	check_error(!ret, create_dirs);

	/* Initial cycle : if the file exists, delete it*/
	if (cycle == 0 && pd->ens_rank == 0) {
		if ((ret = open(file_path, O_CREAT | O_EXCL, S_IRUSR | S_IWUSR)) > 0) {
			ret = unlink(file_path);
			check_error(!ret, unlink);
			
			ret = close(ret);
			check_error(!ret, close);
		}
	}
	
	prepare_file(file, pd->ens_comm, file_path, FILE_CREATE, &ncid);

	buf_size = get_databuf_size(pd, ANAL);
	check_error(buf_size>0, get_databuf_size);

	ret = ncmpi_buffer_attach(ncid, buf_size);
	check_io(ret, ncmpi_buffer_attach);

	ret = ncmpi_enddef(ncid);
	check_io(ret, ncmpi_enddef);

	if (cycle) dtf_time_start();

	write_axis_vars(pd, ANAL, ncid);
	write_data_vars(pd, ANAL, ncid, cycle);

	ret = dtf_transfer(file_path, ncid);
	check_error(!ret, dtf_transfer);
	
	ret = ncmpi_buffer_detach(ncid);
	check_io(ret, ncmpi_buffer_detach);

	ret = ncmpi_close(ncid);
	check_io(ret, ncmpi_close);

	if (cycle) dtf_time_end();

	return 0;
}

int read_anal(PD *pd, char *dir_path, int cycle)
{
	static int first_time = 0;

	int ret, ncid = -1;
	int i;
	char file_path[MAX_PATH_LEN] = { 0 };
	struct file_info *file = &pd->files[ANAL];
	float **arrays = file->var_read_buffers;
	
	// Generate file path
	fmt_filename(cycle, pd->ens_id, 6, dir_path, ".anal.nc", file_path);

	if (!arrays) {
		arrays = malloc(sizeof(float *) * file->ndata_vars);
		check_error(arrays, malloc);
		memset(arrays, 0, sizeof(float *) * file->ndata_vars);
		
		first_time = 1;
		file->nvar_read_buf = file->ndata_vars;
	}
	else {
		first_time = 0;
	}

	prepare_file(file, pd->ens_comm, file_path, FILE_OPEN_R, &ncid);

	if (cycle) dtf_time_start();

	for (i = ANAL_DATA_VARS_OFFSET; i < file->nvars; i++) {
		int arr_idx = i - ANAL_DATA_VARS_OFFSET;
		struct var_pair *var = &file->vars[i];
		int ndims = var->ndims;
		float *array = arrays[arr_idx];
		MPI_Offset ntypes = 0;
		MPI_Datatype dtype = MPI_DATATYPE_NULL;

		if ((first_time && (var->rflag & VAR_READ_ONCE)) ||
				(var->rflag & VAR_READ_ALWAYS) || !cycle) {
			MPI_Offset *start = malloc(sizeof(MPI_Offset) *ndims * 2);
			check_error(start, malloc);
			MPI_Offset *count = start + ndims;
			MPI_Offset total_size = 0;

			int type = get_subarray_type(var);
	
			dtype = subarray_type[type];
			ntypes = (type == XY)? IMAX(pd)*JMAX(pd) : 1;
			
			start[0] = JS_inG(pd);
			start[1] = IS_inG(pd);
			start[2] = 0;
			count[0] = JMAX(pd);
			count[1] = IMAX(pd);
			switch (type) {
				case XY:
					break;
				case ZXY2:
				case ZHXY2:
					count[2] = KMAX;
					break;
				case URBAN:
					count[2] = UKMAX;
					break;
				case LAND:
					count[2] = LKMAX;
					break;
				case OCEAN:
					count[2] = OKMAX;
					break;
				default:
					fprintf(stderr, "[ERROR] Invalid datatype\n");
					MPI_Abort(MPI_COMM_WORLD, EINVAL);
			}

			total_size = IMAX(pd) * JMAX(pd);
			if (type != XY) total_size *= count[2];

			if (!array) {
				array = malloc(sizeof(float) * total_size);
				check_error(array, malloc);
			}
			memset(array, 0, sizeof(float) * total_size);

			ret = ncmpi_iget_vara_float(ncid, var->varid, start, count, array, NULL);
			check_io(ret, ncmpi_iget_vara_float);

			arrays[arr_idx] = array;
			free(start);
		}
	}
	file->var_read_buffers = arrays;

	ret = ncmpi_wait_all(ncid, NC_REQ_ALL, NULL, NULL);
	check_io(ret, ncmpi_wait_all);
	
	ret = dtf_transfer(file_path, ncid);
	check_error(!ret, dtf_transfer);

	ret = ncmpi_close(ncid);
	check_io(ret, ncmpi_close);

	if (cycle) dtf_time_end();

	return 0;
}
