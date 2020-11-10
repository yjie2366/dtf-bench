#include "util.h"

extern struct io_vars anal_vars[];
extern struct io_vars hist_vars[];

int read_hist(PD *pd, int cycle)
{
	int i, ncid;
	int ret = 0;
	struct file_info *file = &pd->files[HIST];
	struct data_buf *arrays = file->var_read_buffers;
	char *file_path = file->file_names[cycle];
	
	prepare_file(file, pd->ens_comm, file_path, FILE_OPEN_R, &ncid);

	for (i = 0; i < NUM_IOVARS_HIST; i++) {
		int idx = hist_vars[i].idx;
		struct data_buf *array = &arrays[idx];
		int varid = array->varid;
		int ndims = array->ndims;
		MPI_Offset *start = array->shape;
	       	MPI_Offset *count = start + ndims;;
		float *data = array->data;

		int real_id;
		ncmpi_inq_varid(ncid, hist_vars[i].name, &real_id);
		if (real_id != varid) {
			fprintf(stderr, "read ID: %d index: %d\n", real_id, varid);
			MPI_Abort(MPI_COMM_WORLD, EINVAL);
		}

		cycle_file_start(pd);

		ret = ncmpi_iget_vara_float(ncid, varid, start, count, data, NULL);
		check_io(ret, ncmpi_iget_vara_float);

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
	
//	report_get_size(pd, HIST, ncid);

	ret = ncmpi_close(ncid);
	check_io(ret, ncmpi_close);

	return ret;
}

int read_anal(PD *pd, int cycle)
{
	int i, ncid;
	int ret = 0;
	struct file_info *file = &pd->files[ANAL];
	struct data_buf *arrays = file->var_read_buffers;
	char *file_path = file->file_names[cycle];

	prepare_file(file, pd->ens_comm, file_path, FILE_OPEN_R, &ncid);

	for (i = 0; i < NUM_IOVARS_ANAL; i++) {
		int idx = anal_vars[i].idx;
		struct data_buf *array = &arrays[idx];
		int varid = array->varid;
		int ndims = array->ndims;
		float *data = array->data;

		MPI_Offset *start = array->shape;
	       	MPI_Offset *count = start + ndims;

		// DEBUG
		int real_id;

		ncmpi_inq_varid(ncid, anal_vars[i].name, &real_id);
		if (real_id != varid) {
			fprintf(stderr, "read ID: %d index: %d\n", real_id, varid);
			MPI_Abort(MPI_COMM_WORLD, EINVAL);
		}

		cycle_file_start(pd);

		ret = ncmpi_iget_vara_float(ncid, varid, start, count, data, NULL);
		check_io(ret, ncmpi_iget_vara_float);

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

	return ret;
}

int write_anal(PD *pd, int cycle)
{
	int ncid = -1;
	int i, ret = 0;
	struct file_info *file = &pd->files[ANAL];
	struct data_buf *arrays = file->var_write_buffers;
	char *file_path = file->file_names[cycle];

	fprintf(stderr, "WRITE file %s\n", file_path);
	prepare_file(file, pd->ens_comm, file_path, FILE_OPEN_W, &ncid);

	for (i = 0; i < NUM_IOVARS_ANAL; i++) {
		int idx = anal_vars[i].idx;
		struct data_buf *array = &arrays[idx];
		int varid = array->varid;
		int ndims = array->ndims;

		MPI_Offset *start = array->shape;
	       	MPI_Offset *count = start + ndims;

		float *data = array->data;
		cycle_file_start(pd);

		ret = ncmpi_iput_vara_float(ncid, varid, start, count, data, NULL);
		check_io(ret, ncmpi_iput_vara_float);

		cycle_file_wend(pd, cycle);
	}

	cycle_file_start(pd);

	ret = ncmpi_wait_all(ncid, NC_REQ_ALL, NULL, NULL);
	check_io(ret, ncmpi_wait_all);
	
	cycle_file_wend(pd, cycle);

	cycle_transfer_start(pd);

	ret = dtf_transfer(file_path, ncid);
	check_error(!ret, dtf_transfer);

	cycle_transfer_wend(pd, cycle);
	
//	report_put_size(pd, ANAL, ncid);

	ret = ncmpi_close(ncid);
	check_io(ret, ncmpi_close);

	return ret;
}
