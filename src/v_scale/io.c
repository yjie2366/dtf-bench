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
	struct data_buf *arrays = file->var_write_buffers;
	int num_axis = NUM_AXIS_VARS;
	int num_coords = file->nassct_coords;
	
	/* Write to AXIS variables */
	for (i = 0; i < num_axis + num_coords; i++) {
		struct data_buf *array = &arrays[i];
		int ndims = array->ndims;
		float *data = array->data;
		MPI_Offset *start = array->shape;
	       	MPI_Offset *count = array->shape + ndims;

		int varid = array->varid;

		/* 
		 * MIND that some process does not have data
		 * for 1D axis variable
		*/
		if (data) {
			cycle_file_start(pd);

			ret = ncmpi_iput_vara_float(ncid, varid, start,
					count, data, NULL);
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
	int ret = 0; int ncid = -1;
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

//	report_put_size(pd, HIST, ncid);

	ret = ncmpi_buffer_detach(ncid);
	check_io(ret, ncmpi_buffer_detach);

	ret = ncmpi_close(ncid);
	check_io(ret, ncmpi_close);

	return ret;
}

int write_anal(PD *pd, int cycle)
{
	int ret; int ncid = -1;
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

//	report_put_size(pd, ANAL, ncid);

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

	fprintf(stderr, "READ file %s\n", file_path);
	prepare_file(file, pd->ens_comm, file_path, FILE_OPEN_R, &ncid);

	for (i = 0; i < NUM_IOVARS_ANAL; i++) {
		int idx = anal_vars[i].idx;
		struct data_buf *array = &arrays[idx];
		float *data = array->data;
		int ndims = array->ndims;
		int varid = array->varid;

		MPI_Offset *start = array->shape;
	       	MPI_Offset *count = start + ndims;
		MPI_Offset ntypes = array->ntypes;
		MPI_Datatype dtype = array->dtype;

		/* 
		 * OCEAN, LAND and URBAN variables are not read
		 * after the initial cycle; so we ignore them
		 */
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

//	report_get_size(pd, ANAL, ncid);

	ret = ncmpi_close(ncid);
	check_io(ret, ncmpi_close);

	return 0;
}
