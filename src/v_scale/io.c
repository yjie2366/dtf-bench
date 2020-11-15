#include "util.h"
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <time.h>

/* write axes and associated coord vars */
static int write_axis_vars(PD *pd, int file_idx, int cycle, int ncid)
{
	int i, ret = 0;
	struct file_info *file = &pd->files[file_idx];
	int num_axis = NUM_AXIS_VARS;
	int num_coords = file->nassct_coords;
	struct data_buf *arrays = file->var_write_buffers;
	
	/* Write to AssociatedCoord variables */
	for (i = 0; i < (num_coords + num_axis); i++) {
		struct data_buf *array = &arrays[i];
		int ndims = array->ndims;
		int varid = array->varid;
		MPI_Offset *start = NULL;
	       	MPI_Offset *count = NULL;
		float *data = array->data;

		if (data) {
			start = array->shape;
			count = start + ndims;

			cycle_file_start(pd);

			ret = ncmpi_iput_vara_float(ncid, varid, start, count, data, NULL);
			check_io(ret, ncmpi_iput_vara_float);

			cycle_file_wend(pd, cycle);
		}
	}
	
	cycle_file_start(pd);
	
	ret = ncmpi_wait_all(ncid, NC_REQ_ALL, NULL, NULL);
	check_io(ret, ncmpi_wait_all);
	
	cycle_file_wend(pd, cycle);

	return ret;
}

static int write_data_vars(PD *pd, int file_idx, int ncid, int cycle)
{
	int i, ret = 0;
	int offset = (file_idx == ANAL) ? ANAL_DATA_VARS_OFFSET : HIST_DATA_VARS_OFFSET;
	struct file_info *file = &pd->files[file_idx];
	struct data_buf *arrays = file->var_write_buffers;

	/* Write to data variables */
	for (i = offset; i < file->nvars; i++) {
		struct data_buf *array = &arrays[i];
		int ndims = array->ndims;
		int varid = array->varid;
		MPI_Offset *start = array->shape;
		MPI_Offset *count = start + ndims;
		float *data = array->data;

		// Keep in mind that processes at the edge of ensemble 
		// may have bigger buffer because of the included HALO area
		cycle_file_start(pd);

		ret = ncmpi_bput_vara_float(ncid, varid, start, count, data, NULL);
		check_io(ret, ncmpi_bput_vara_float);

		cycle_file_wend(pd, cycle);
	}
	
	cycle_file_start(pd);

	ret = ncmpi_wait_all(ncid, NC_REQ_ALL, NULL, NULL);
	check_io(ret, ncmpi_wait_all);
	
	cycle_file_wend(pd, cycle);

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

	cycle_file_start(pd);

	index[0] = 0;
	ret = ncmpi_put_var1_double_all(ncid, varid_time, index, &time);
	check_io(ret, ncmpi_put_var1_double_all);

	index[1] = 0;
	ret = ncmpi_put_var1_double_all(ncid, varid_bnds, index, &time_s);
	check_io(ret, ncmpi_put_var1_double_all);

	index[1] = 1;
	ret = ncmpi_put_var1_double_all(ncid, varid_bnds, index, &time_e);
	check_io(ret, ncmpi_put_var1_double_all);

	cycle_file_wend(pd, cycle);

	return ret;
}

int write_hist(PD *pd, int cycle)
{
	int ret = 0;
	int ncid = -1;
	struct file_info *file = &pd->files[HIST];
	char *file_path = file->file_names[cycle];
	MPI_Offset buf_size = file->databuf_sz;
	
	prepare_file(file, pd->ens_comm, file_path, FILE_CREATE, &ncid);
	
	ret = ncmpi_buffer_attach(ncid, buf_size);
	check_io(ret, ncmpi_buffer_attach);

	ret = ncmpi_enddef(ncid);
	check_io(ret, ncmpi_enddef);
	
	write_axis_vars(pd, HIST, cycle, ncid);
	write_time_var(pd, ncid, cycle);
	write_data_vars(pd, HIST, ncid, cycle);

	cycle_transfer_start(pd);

	ret = dtf_transfer(file_path, ncid);
	check_error(!ret, dtf_transfer);

	cycle_transfer_wend(pd, cycle);

	report_put_size(pd, HIST, ncid);

	ret = ncmpi_buffer_detach(ncid);
	check_io(ret, ncmpi_buffer_detach);

	ret = ncmpi_close(ncid);
	check_io(ret, ncmpi_close);

	return ret;
}

int write_anal(PD *pd, int cycle)
{
	int ret;
	int ncid = -1;
	struct file_info *file = &pd->files[ANAL];
	char *file_path = file->file_names[cycle];
	MPI_Offset buf_size = file->databuf_sz;

	prepare_file(file, pd->ens_comm, file_path, FILE_CREATE, &ncid);

	ret = ncmpi_buffer_attach(ncid, buf_size);
	check_io(ret, ncmpi_buffer_attach);

	ret = ncmpi_enddef(ncid);
	check_io(ret, ncmpi_enddef);

	write_axis_vars(pd, ANAL, cycle, ncid);
	write_data_vars(pd, ANAL, ncid, cycle);

	cycle_transfer_start(pd);

	ret = dtf_transfer(file_path, ncid);
	check_error(!ret, dtf_transfer);
	
	cycle_transfer_wend(pd, cycle);
	
	report_put_size(pd, ANAL, ncid);

	ret = ncmpi_buffer_detach(ncid);
	check_io(ret, ncmpi_buffer_detach);

	ret = ncmpi_close(ncid);
	check_io(ret, ncmpi_close);

	return 0;
}

extern struct io_vars anal_vars[];

int read_anal(PD *pd, int cycle)
{
	int i, ret, ncid = -1;
	struct file_info *file = &pd->files[ANAL];
	struct data_buf *arrays = file->var_read_buffers;
	char *file_path = file->file_names[cycle];

	prepare_file(file, pd->ens_comm, file_path, FILE_OPEN_R, &ncid);

	for (i = 0; i < NUM_IOVARS_ANAL; i++) {
		int idx = anal_vars[i].idx;
		struct var_pair *var = &file->vars[idx];
		struct data_buf *array = &arrays[idx];
		float *data = NULL;
		int ndims, varid;
		MPI_Offset ntypes = 1;
		MPI_Datatype dtype = MPI_DATATYPE_NULL;
		MPI_Offset *start, *count;

		/* 
		 * OCEAN, LAND and URBAN variables are not read
		 * after the initial cycle; so we ignore them
		 */
		if (!cycle) {
			int size[3] = { 0 };
			int sub_size[3] = { 0 };
			int sub_off[3] = { 0 };
			int j;

			MPI_Offset total_count = 1;

			ndims = var->ndims;
			varid = var->varid;

			start = (MPI_Offset *)malloc(sizeof(MPI_Offset) * ndims * 2);
			check_error(start, malloc);
			memset(start, 0, sizeof(MPI_Offset) * ndims * 2);
			count = start + ndims;

			for (j = 0; j < ndims; j++) {
				if (strchr(var->dim_name[j], 'z')) {
					start[j] = 0;
					count[j] = KMAX;
					size[j] = KA;
					sub_size[j] = KMAX;
					sub_off[j] = KHALO;

					if (strchr(var->dim_name[j], 'h'))
						sub_off[j] -= 1;
				}
				else if (strchr(var->dim_name[j], 'y')) {
					start[j] = JS_inG(pd);
					count[j] = JMAX(pd);
					size[j] = JA(pd);
					sub_size[j] = JMAX(pd);
					sub_off[j] = JHALO;
				}
				else if (strchr(var->dim_name[j], 'x')) {
					start[j] = IS_inG(pd);
					count[j] = IMAX(pd);
					size[j] = IA(pd);
					sub_size[j] = IMAX(pd);
					sub_off[j] = IHALO;
				}
			}

			ntypes = 1;

			MPI_Type_create_subarray(ndims, size, sub_size,
					sub_off, MPI_ORDER_C, MPI_FLOAT, &dtype);
			MPI_Type_commit(&dtype);

			total_count = IA(pd) * JA(pd) * KA;

			data = (float *)malloc(sizeof(float) * total_count);
			check_error(data, malloc);
			memset(data, 0, sizeof(float) * total_count);

			array->ndims = ndims;
			array->shape = start;
			array->varid = varid;
			array->data = data;
			array->nelems = total_count;
			array->dtype = dtype;
			array->ntypes = ntypes;
		}
		else {
			ndims = array->ndims;
			varid = array->varid;
			data = array->data;
			ntypes = array->ntypes;
			dtype = array->dtype;
			start = array->shape;
			count = start + ndims;
		}

		cycle_file_start(pd);

		ret = ncmpi_iget_vara(ncid, varid, start, count, data, ntypes, dtype, NULL);
		check_io(ret, ncmpi_iget_vara);

		cycle_file_rend(pd, cycle);
	}

	cycle_file_start(pd);

	ret = ncmpi_wait_all(ncid, NC_REQ_ALL, NULL, NULL);
	check_io(ret, ncmpi_wait_all);

	cycle_file_rend(pd, cycle);

	cycle_transfer_start(pd);

	ret = dtf_transfer(file_path, ncid);
	check_error(!ret, dtf_transfer);

	cycle_transfer_rend(pd, cycle);

	report_get_size(pd, ANAL, ncid);

	ret = ncmpi_close(ncid);
	check_io(ret, ncmpi_close);

	return 0;
}
