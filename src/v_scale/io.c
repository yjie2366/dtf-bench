#include "util.h"
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <time.h>

extern MPI_Datatype subarray_type[];

/* write axes and associated coord vars */
static int write_axis_vars(PD *pd, int file_idx, int cycle, int ncid)
{
	int first_run = 0;
	int i, ret = 0;
	struct file_info *file = &pd->files[file_idx];
	int num_axis = NUM_AXIS_VARS;
	int num_coords = file->nassct_coords;
	struct data_buf *arrays = file->axes_buffer;
	
	if (!arrays) {
		file->naxes_buf = num_axis + num_coords;
		arrays = (struct data_buf *)malloc(sizeof(struct data_buf) * file->naxes_buf);
		check_error(arrays, malloc);
		memset(arrays, 0, sizeof(struct data_buf) * file->naxes_buf);
		first_run = 1;
	}

	/* Write to AXIS variables */
	for (i = 0; i < num_axis; i++) {
		struct var_pair *var = &file->vars[i];
		MPI_Offset dim_len = file->dims[var->dims[0]].length;
		struct data_buf *array = &arrays[i];
		int ndims = var->ndims;
		int varid = var->varid;
		float *data = array->data;
		MPI_Offset start = 0, count = 0;

		/* Axies are 1-D variables */
		if (first_run) {
			MPI_Offset *shape = NULL, *idxes = NULL;
		
			// y axes
			if (pd->proc_rank_x == 0 && (strchr(var->dim_name[0], 'y'))) {
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
			}

			// x axes
			if (pd->proc_rank_y == 0 && (strchr(var->dim_name[0], 'x'))) {
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
			}		
			
			// z, Z, X, and Y axes
			if (pd->ens_rank == 0) {
				if (    strchr(var->dim_name[0], 'z') ||
					strchr(var->dim_name[0], 'Z') ||
					strchr(var->dim_name[0], 'X') ||
					strchr(var->dim_name[0], 'Y')	) {
					start = 0;
					count = dim_len;
				}
			}

			if (count) {
				shape = (MPI_Offset *)malloc(sizeof(MPI_Offset)*ndims*4);
				check_error(shape, malloc);
				idxes = shape + ndims * 2;
			
				shape[0] = start;
				shape[1] = count;
				/* axes variables don't care halo
				 * because they are never read by LETKF
				 *  [start_idx, end_idx)
				 */
				idxes[0] = 0;		 // s_idx
				idxes[1] = count;	 // e_idx
				
				data = (float *)malloc(sizeof(float) * count);
				check_error(data, malloc);
				memset(data, -1, sizeof(float) * count);
			}

			array->data = data;
			array->shape = shape;
			array->idxes = idxes;
			array->ndims = ndims;
			array->varid = varid;
			array->nelems = count;
		}
		else {
			if (array->shape) {
				start = array->shape[0];
				count = array->shape[1];
			}
		}

		/* 
		 * MIND that some process does not have data for this variable
		*/
		if (data) {
			ret = fill_buffer(array, (float)pd->world_rank, (float)cycle, SCALE_WEIGHT);
			check_error(!ret, fill_buffer);

			ret = ncmpi_iput_vara_float(ncid, varid, &start,
					&count, data, NULL);
			check_io(ret, ncmpi_iput_vara_float);
		}
	}

	/* Write to AssociatedCoord variables */
	for (; i < (num_coords + num_axis); i++) {
		int j;
		struct var_pair *var = &file->vars[i];
		int ndims = var->ndims;
		int varid = var->varid;
		struct data_buf *array = &arrays[i];
		MPI_Offset *start, *count;
		float *data = array->data;

		/* initialize shape of variable in PnetCDF file */
		if (first_run) {
			MPI_Offset *s_idx, *e_idx;
			MPI_Offset total_count = 1;

			start = (MPI_Offset *)malloc(sizeof(MPI_Offset) * ndims * 4);
			check_error(start, malloc);
			memset(start, 0, sizeof(MPI_Offset) * ndims * 4);

			count = start + ndims;
			s_idx = count + ndims;
			e_idx = s_idx + ndims;

			for (j = 0; j < ndims; j++) {

				MPI_Offset dim_len = file->dims[var->dims[j]].length;

				if (file_idx == ANAL) {
					switch (j) {
					case 0:
						start[j] = JSGA(pd);
						count[j] = JEB(pd) - JSB(pd) + 1;
						s_idx[j] = (!pd->proc_rank_y) ? JHALO : 0;
						e_idx[j] = (pd->proc_rank_y == (pd->proc_num_y - 1))?
							count[j] - JHALO : count[j];
						break;
					case 1:
						start[j] = ISGA(pd);
						count[j] = IEB(pd) - ISB(pd) + 1;
						s_idx[j] = (!pd->proc_rank_x) ? IHALO : 0;
						e_idx[j] = (pd->proc_rank_x == (pd->proc_num_x - 1))?
							count[j] - IHALO : count[j];
						break;
					case 2:
						start[j] = 0;
						count[j] = dim_len;
						s_idx[j] = 0;
						e_idx[j] = count[j];
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
					s_idx[j] = 0;
					e_idx[j] = count[j];
				}

				total_count *= count[j];
			}

			data = (float *)malloc(sizeof(float) * total_count);
			check_error(data, malloc);
			memset(data, -1, sizeof(float) * total_count);
			
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

		ret = fill_buffer(array, (float)pd->world_rank, (float)cycle, SCALE_WEIGHT);
		check_error(!ret, fill_buffer);

		ret = ncmpi_iput_vara_float(ncid, var->varid, start,
				count, data, NULL);
		check_io(ret, ncmpi_iput_vara_float);
	}
	
	ret = ncmpi_wait_all(ncid, NC_REQ_ALL, NULL, NULL);
	check_io(ret, ncmpi_wait_all);

	file->axes_buffer = arrays;

	return ret;
}

static int write_data_vars(PD *pd, int file_idx, int ncid, int cycle)
{
	int i, first_run = 0, ret = 0;
	int offset = (file_idx == ANAL) ? ANAL_DATA_VARS_OFFSET : HIST_DATA_VARS_OFFSET;
	struct file_info *file = &pd->files[file_idx];
	struct data_buf *arrays = file->var_write_buffers;

	if (!arrays) {
		arrays = (struct data_buf *)malloc(sizeof(struct data_buf) * file->ndata_vars);
		check_error(arrays, malloc);
		memset(arrays, 0, sizeof(float *) * file->ndata_vars);

		file->nvar_write_buf = file->ndata_vars;
		first_run = 1;
	}

	/* Write to data variables */
	for (i = offset; i < file->nvars; i++) {
		int arr_idx = i - offset;
		struct var_pair *var = &file->vars[i];
		struct data_buf *array = &arrays[arr_idx];

		int j;
		int ndims = var->ndims;
		int varid = var->varid;
		MPI_Offset *start = array->shape;
		MPI_Offset *count = start + ndims;
		float *data = array->data;

		// Keep in mind that processes at the edge of ensemble 
		// may have bigger buffer because of the included HALO area
		if (first_run) {
			MPI_Offset *s_idx, *e_idx;
			MPI_Offset total_count = 1;

			start = (MPI_Offset *)malloc(sizeof(MPI_Offset) * ndims * 4);
			check_error(start, malloc);
			count = start + ndims;
			s_idx = count + ndims;
			e_idx = s_idx + ndims;

			for (j = 0; j < ndims; j++) {
				if (!strcmp(var->dim_name[j], "time")) {
					start[j] = 0;
					count[j] = 1;
					s_idx[j] = 0;
					e_idx[j] = 1;
				}
				else if (strchr(var->dim_name[j], 'z')) { // contains z
					start[j] = 0;
					count[j] = file->dims[var->dims[j]].length;
					s_idx[j] = 0;
					e_idx[j] = count[j];
				}
				else if (strchr(var->dim_name[j], 'y')) {
					if (file_idx == ANAL) {
						start[j] = JSGA(pd);
						count[j] = JEB(pd) - JSB(pd) + 1;
						s_idx[j] = (!pd->proc_rank_y) ? JHALO : 0;
						e_idx[j] = (pd->proc_rank_y == (pd->proc_num_y - 1))?
							count[j]-JHALO : count[j];
					}
					else {
						start[j] = SY_hist(pd);
						count[j] = JMAX(pd);
						s_idx[j] = 0;
						e_idx[j] = count[j];
					}
				}
				else if (strchr(var->dim_name[j], 'x')) {
					if (file_idx == ANAL) {
						start[j] = ISGA(pd);
						count[j] = IEB(pd) - ISB(pd) + 1;
						s_idx[j] = (!pd->proc_rank_x) ? IHALO : 0;
						e_idx[j] = (pd->proc_rank_x == (pd->proc_num_x - 1))?
							count[j] - IHALO : count[j];
					}
					else {
						start[j] = SX_hist(pd);
						count[j] = IMAX(pd);
						s_idx[j] = 0;
						e_idx[j] = count[j];
					}
				}

				total_count *= count[j];
			}

			data = (float *)malloc(sizeof(float) * total_count);
			check_error(data, malloc);
			memset(data, -1, sizeof(float) * total_count);

			array->shape = start;
			array->idxes = s_idx;
			array->data = data;
			array->ndims = ndims;
			array->varid = varid;
			array->nelems = total_count;
		}

		ret = fill_buffer(array, (float)pd->world_rank, (float)cycle, SCALE_WEIGHT);
		check_error(!ret, fill_buffer);

		if (!pd->world_rank) {
			int idx;
			for (idx = 0; idx < num_var; idx++) {
				int j;
				fprintf(stderr, "Variable: %s\n", hist_vars[idx]);
				for (j = 0; j < arrays[idx].nelems; j++) {
					fprintf(stderr, "%f ", arrays[idx].data[j]);
				}
				fprintf(stderr, "\n");
			}
		}

		ret = ncmpi_bput_vara_float(ncid, var->varid, start, count, data, NULL);
		check_io(ret, ncmpi_bput_vara_float);
	}

	ret = ncmpi_wait_all(ncid, NC_REQ_ALL, NULL, NULL);
	check_io(ret, ncmpi_wait_all);
	
	file->var_write_buffers = arrays;

	return 0;
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
	
	// Create data folder to store files in the top directory
	// if it does not exist 
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

	write_axis_vars(pd, HIST, cycle, ncid);
	write_time_var(pd, ncid, cycle);
	write_data_vars(pd, HIST, ncid, cycle);

	ret = dtf_transfer(file_path, ncid);
	check_error(!ret, dtf_transfer);

	ret = ncmpi_buffer_detach(ncid);
	check_io(ret, ncmpi_buffer_detach);

	ret = ncmpi_close(ncid);
	check_io(ret, ncmpi_close);
	
	if (cycle) dtf_writetime_end();

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

	// Create data folder to store files in the top directory
	// if it does not exist 
	ret = create_dirs(file_path);
	check_error(!ret, create_dirs);

	prepare_file(file, pd->ens_comm, file_path, FILE_CREATE, &ncid);

	buf_size = get_databuf_size(pd, ANAL);
	check_error(buf_size>0, get_databuf_size);

	ret = ncmpi_buffer_attach(ncid, buf_size);
	check_io(ret, ncmpi_buffer_attach);

	ret = ncmpi_enddef(ncid);
	check_io(ret, ncmpi_enddef);

	if (cycle) dtf_time_start();

	write_axis_vars(pd, ANAL, cycle, ncid);
	write_data_vars(pd, ANAL, ncid, cycle);

	ret = dtf_transfer(file_path, ncid);
	check_error(!ret, dtf_transfer);
	
	ret = ncmpi_buffer_detach(ncid);
	check_io(ret, ncmpi_buffer_detach);

	ret = ncmpi_close(ncid);
	check_io(ret, ncmpi_close);

	if (cycle) dtf_writetime_end();

	return 0;
}

#define VARS_ANAL_READ		11

static char *anal_vars[VARS_ANAL_READ] = {
	"DENS", "MOMZ", "MOMX", "MOMY", "RHOT",
	"QV", "QC", "QR", "QI", "QS", "QG"
};

int read_anal(PD *pd, char *dir_path, int cycle)
{
	int first_run = 0;
	int i, ret, ncid = -1;
	char file_path[MAX_PATH_LEN] = { 0 };
	struct file_info *file = &pd->files[ANAL];
	struct data_buf *arrays = file->var_read_buffers;

	if (!arrays) {
		arrays = malloc(sizeof(struct data_buf) * VARS_ANAL_READ);
		check_error(arrays, malloc);
		memset(arrays, 0, sizeof(struct data_buf) * VARS_ANAL_READ);
		
		file->nvar_read_buf = VARS_ANAL_READ;
		first_run = 1;
	}

	// Generate file path
	fmt_filename(cycle, pd->ens_id, 6, dir_path, ".anal.nc", file_path);
	prepare_file(file, pd->ens_comm, file_path, FILE_OPEN_R, &ncid);

	if (cycle) dtf_time_start();

	for (i = 0; i < VARS_ANAL_READ; i++) {
		struct var_pair *var = NULL;
		struct data_buf *array = &arrays[i];
		float *data = array->data;
		int ndims = array->ndims, varid = array->varid;

		MPI_Offset ntypes = 0;
		MPI_Datatype dtype = MPI_DATATYPE_NULL;
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
		
		/* 
		 * OCEAN, LAND and URBAN variables are not read
		 * after the initial cycle; so we ignore them
		 */
		if (first_run) {
			MPI_Offset *s_idx, *e_idx;
			MPI_Offset total_count = 1;

			ndims = var->ndims;

			start = malloc(sizeof(MPI_Offset) * ndims * 4);
			check_error(start, malloc);
			memset(start, 0, sizeof(MPI_Offset) * ndims * 4);

			count = start + ndims;
			s_idx = count + ndims;
			e_idx = s_idx + ndims;

			if (ndims == 2) {
				start[0] = JS_inG(pd) - JHALO;
				start[1] = IS_inG(pd) - IHALO;

				count[0] = JA(pd);
				count[1] = IA(pd);

				s_idx[0] = JHALO;
				s_idx[1] = IHALO;

				e_idx[0] = s_idx[0] + JMAX(pd);
				e_idx[1] = s_idx[1] - IMAX(pd);

				total_count = IA(pd) * JA(pd);
			}
			else {
				start[0] = JS_inG(pd);
				start[1] = IS_inG(pd);
				start[2] = 0;

				count[0] = JMAX(pd);
				count[1] = IMAX(pd);
				count[2] = KMAX;

				s_idx[0] = JHALO;
				s_idx[1] = IHALO;
				s_idx[2] = KHALO;

				e_idx[0] = s_idx[0] + count[0];
				e_idx[1] = s_idx[1] + count[1];
				e_idx[2] = s_idx[2] + count[2];

				total_count = IA(pd) * JA(pd) * KA;
			}

			data = malloc(sizeof(float) * total_count);
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

		if (ndims == 2) {
			dtype = subarray_type[XY];
			ntypes = IA(pd) * JA(pd);
		}
		else {
			dtype = (strcmp(var->dim_name[2], "zh")) ? 
				subarray_type[ZXY2] :
				subarray_type[ZHXY2];
			ntypes = 1;
		}

		ret = ncmpi_iget_vara(ncid, varid, start, count, data,
				ntypes, dtype, NULL);
		check_io(ret, ncmpi_iget_vara);
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
	for (i = 0; i < VARS_ANAL_READ; i++) {
		struct data_buf *read_buf = &file->var_read_buffers[i];
		
		ret = compare_buffer(pd, read_buf, cycle, LETKF_WEIGHT);
		check_error(!ret, compare_buffer);
	}
	
	return 0;
}
