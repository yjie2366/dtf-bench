#include "util.h"

extern struct io_vars anal_vars[];
extern struct io_vars hist_vars[];

int read_hist(PD *pd, int cycle)
{
	int i, ncid, ret = 0;
	struct file_info *file = &pd->files[HIST];
	struct data_buf *arrays = file->var_read_buffers;
	char *file_path = file->file_names[cycle];
	
	prepare_file(file, pd->ens_comm, file_path, FILE_OPEN_R, &ncid);

#ifdef TIMING
	double get_t = 0.0;
	double wait_t = 0.0;
	double trans_t = 0.0;
	double close_t = 0.0;
#endif
	for (i = 0; i < NUM_IOVARS_HIST; i++) {
		int idx = hist_vars[i].idx;
		struct data_buf *array = &arrays[idx];
		struct var_pair *var = &file->vars[idx];
		int varid, ndims;
		MPI_Offset *start, *count;
		float *data = NULL;

		if (!cycle) {
			int j;
			MPI_Offset total_count = 1;
			
			ndims = var->ndims;
			varid = var->varid;
			
			start = (MPI_Offset *)malloc(sizeof(MPI_Offset) * ndims * 2);
			check_error(start, malloc);
			count = start + ndims;

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

				total_count *= count[j];
			}
			
			data = (float *)malloc(sizeof(float) * total_count);
			check_error(data, malloc);
			memset(data, 0, total_count * sizeof(float));
		
			array->data = data;
			array->shape = start;
			array->ndims = ndims;
			array->varid = varid;
			array->nelems = total_count;
		}
		else {
			data = array->data;
			ndims = array->ndims;
			varid = array->varid;
			start = array->shape;
			count = start + ndims;
		}

		cycle_file_start(pd);
#ifdef TIMING
		double get_ts = MPI_Wtime();
#endif
		ret = ncmpi_iget_vara_float(ncid, varid, start, count, data, NULL);
		check_io(ret, ncmpi_iget_vara_float);
#ifdef TIMING
		double get_te = MPI_Wtime();
		get_t += get_te - get_ts;
#endif
		cycle_file_rend(pd, cycle);
	}

	cycle_file_start(pd);
#ifdef TIMING
	double wait_ts = MPI_Wtime();
#endif
	ret = ncmpi_wait_all(ncid, NC_REQ_ALL, NULL, NULL);
	check_io(ret, ncmpi_wait_all);
#ifdef TIMING
	double wait_te = MPI_Wtime();
	wait_t += wait_te - wait_ts;
#endif
	cycle_file_rend(pd, cycle);

	cycle_transfer_start(pd);
#ifdef TIMING
	double trans_ts = MPI_Wtime();
#endif
	ret = dtf_transfer(file_path, ncid);
	check_error(!ret, dtf_transfer);
#ifdef TIMING
	double trans_te = MPI_Wtime();
	trans_t += trans_te - trans_ts;
#endif	
	cycle_transfer_rend(pd, cycle);

#ifdef TIMING
	report_get_size(pd, HIST, ncid);
	double close_ts = MPI_Wtime();
#endif
	ret = ncmpi_close(ncid);
	check_io(ret, ncmpi_close);
#ifdef TIMING
	double close_te = MPI_Wtime();
	close_t = close_te - close_ts;

	fprintf(stderr, "CYCLE[%d] LETKF read HIST:\n", cycle);
	fprintf(stderr, "iget: %f wait_all: %f  transfer: %f  close: %f\n",
			get_t, wait_t, trans_t, close_t);
	fprintf(stderr, "----- CYCLE[%d] LETKF READ HIST END -----\n", cycle);
#endif
	return ret;
}

int read_anal(PD *pd, int cycle)
{
	int i, ncid, ret = 0;
	struct file_info *file = &pd->files[ANAL];
	struct data_buf *arrays = file->var_read_buffers;
	char *file_path = file->file_names[cycle];
	
	prepare_file(file, pd->ens_comm, file_path, FILE_OPEN_R, &ncid);
	
#ifdef TIMING
	double get_t = 0.0;
	double wait_t = 0.0;
	double trans_t = 0.0;
	double close_t = 0.0;
#endif
	for (i = 0; i < NUM_IOVARS_ANAL; i++) {
		int j;
		int idx = anal_vars[i].idx;
		struct data_buf *array = &arrays[idx];
		struct var_pair *var = &file->vars[idx];
		int varid;
		int ndims;
		float *data = NULL;
		MPI_Offset *start, *count;

		/* first cycle initialize read buffer */
		if (!cycle) {
			MPI_Offset total_count = 1;

			ndims = var->ndims;
			varid = var->varid;

			start = (MPI_Offset *)malloc(sizeof(MPI_Offset) * ndims * 2);
			check_error(start, malloc);
			count = start + ndims;

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
			}

			data = (float *)malloc(sizeof(float) * total_count);
			check_error(data, malloc);
			memset(data, 0, sizeof(float) * total_count);

			array->ndims = ndims;
			array->shape = start;
			array->varid = varid;
			array->data = data;
			array->nelems = total_count;
		}
		else {
			data = array->data;
			varid = array->varid;
			ndims = array->ndims;
			start = array->shape;
			count = start + ndims;
		}

		cycle_file_start(pd);
#ifdef TIMING
		double get_ts = MPI_Wtime();
#endif
		ret = ncmpi_iget_vara_float(ncid, varid, start, count, data, NULL);
		check_io(ret, ncmpi_iget_vara_float);
#ifdef TIMING
		double get_te = MPI_Wtime();
		get_t += get_te - get_ts;
#endif
		cycle_file_rend(pd, cycle);
	}

	cycle_file_start(pd);
#ifdef TIMING
	double wait_ts = MPI_Wtime();
#endif
	ret = ncmpi_wait_all(ncid, NC_REQ_ALL, NULL, NULL);
	check_io(ret, ncmpi_wait_all);
#ifdef TIMING
	double wait_te = MPI_Wtime();
	wait_t += wait_te - wait_ts;
#endif
	cycle_file_rend(pd, cycle);

	cycle_transfer_start(pd);
#ifdef TIMING
	double trans_ts = MPI_Wtime();
#endif
	ret = dtf_transfer(file_path, ncid);
	check_error(!ret, dtf_transfer);
#ifdef TIMING
	double trans_te = MPI_Wtime();
	trans_t += trans_te - trans_ts;
#endif	
	cycle_transfer_rend(pd, cycle);

#ifdef TIMING
	report_get_size(pd, ANAL, ncid);
	double close_ts = MPI_Wtime();
#endif
	ret = ncmpi_close(ncid);
	check_io(ret, ncmpi_close);
#ifdef TIMING
	double close_te = MPI_Wtime();
	close_t += close_te - close_ts;
	
	fprintf(stderr, "CYCLE: %d LETKF read ANAL:\n", cycle);
	fprintf(stderr, "iget: %f wait_all: %f  transfer: %f  close: %f\n",
			get_t, wait_t, trans_t, close_t);
	fprintf(stderr, "----- CYCLE[%d] LETKF READ ANAL END -----\n", cycle);
#endif
	return ret;
}

int write_anal(PD *pd, int cycle)
{
	int ncid = -1;
	int i, ret = 0;
	struct file_info *file = &pd->files[ANAL];
	struct data_buf *arrays = file->var_write_buffers;
	char *file_path = file->file_names[cycle];

	prepare_file(file, pd->ens_comm, file_path, FILE_OPEN_W, &ncid);

#ifdef TIMING
	double put_t = 0.0;
	double wait_t = 0.0;
	double trans_t = 0.0;
	double close_t = 0.0;
#endif
	for (i = 0; i < NUM_IOVARS_ANAL; i++) {
		int idx = anal_vars[i].idx;
		struct data_buf *array = &arrays[idx];
		int varid = array->varid;
		int ndims = array->ndims;
		float *data = array->data;
		MPI_Offset *start = array->shape;
		MPI_Offset *count = start + ndims;

		cycle_file_start(pd);
#ifdef TIMING
		double put_ts = MPI_Wtime();
#endif
		ret = ncmpi_iput_vara_float(ncid, varid, start, count, data, NULL);
		check_io(ret, ncmpi_iput_vara_float);
#ifdef TIMING
		double put_te = MPI_Wtime();
		put_t += put_te - put_ts;
#endif
		cycle_file_wend(pd, cycle);
	}

	cycle_file_start(pd);
#ifdef TIMING
	double wait_ts = MPI_Wtime();
#endif
	ret = ncmpi_wait_all(ncid, NC_REQ_ALL, NULL, NULL);
	check_io(ret, ncmpi_wait_all);
#ifdef TIMING
	double wait_te = MPI_Wtime();
	wait_t += wait_te - wait_ts;
#endif
	cycle_file_wend(pd, cycle);

	cycle_transfer_start(pd);

#ifdef TIMING
	double trans_ts = MPI_Wtime();
#endif
	ret = dtf_transfer(file_path, ncid);
	check_error(!ret, dtf_transfer);
#ifdef TIMING
	double trans_te = MPI_Wtime();
	trans_t += trans_te - trans_ts;
#endif
	cycle_transfer_wend(pd, cycle);

#ifdef TIMING	
	report_put_size(pd, ANAL, ncid);
	double close_ts = MPI_Wtime();
#endif
	ret = ncmpi_close(ncid);
	check_io(ret, ncmpi_close);
#ifdef TIMING
	double close_te = MPI_Wtime();
	close_t += close_te - close_ts;
	
	fprintf(stderr, "CYCLE[%d]: LETKF write ANAL:\n", cycle);
	fprintf(stderr, "iput: %f wait_all: %f  transfer: %f  close: %f\n",
			put_t, wait_t, trans_t, close_t);
	fprintf(stderr, "----- CYCLE[%d] LETKF WRITE ANAL END -----\n", cycle);
#endif
	return ret;
}
