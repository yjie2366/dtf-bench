#include "util.h"
#include <math.h>

extern char *comp_name;

struct dim_pair anal_dims[NUM_ANALDIMS] = {
	{ "z"   , KMAX, 	-1 },
	{ "zh"  , KMAX+1, 	-1 },
	{ "oz"  , OKMAX, 	-1 },
	{ "ozh" , OKMAX+1,	-1 },
	{ "lz"  , LKMAX, 	-1 },
	{ "lzh" , LKMAX+1, 	-1 },
	{ "uz"  , UKMAX, 	-1 },
	{ "uzh" , UKMAX+1,	-1 },
	{ "x"   , IHALO*2,	-1 },
	{ "xh"  , IHALO*2,	-1 },
	{ "y"   , JHALO*2,	-1 },
	{ "yh"  , JHALO*2, 	-1 },
	{ "CZ"  , KA, 		-1 },
	{ "FZ"  , KA+1,	 	-1 },
	{ "FDZ" , KA-1, 	-1 },
	{ "OCZ" , OKMAX, 	-1 },
	{ "OFZ" , OKMAX+1, 	-1 },
	{ "LCZ" , LKMAX, 	-1 },
	{ "LFZ" , LKMAX+1, 	-1 },
	{ "UCZ" , UKMAX, 	-1 },
	{ "UFZ" , UKMAX+1, 	-1 },
	{ "CX"  , IHALO*2, 	-1 },
	{ "CY"  , JHALO*2, 	-1 },
	{ "FX"  , IHALO*2+1, 	-1 },
	{ "FY"  , JHALO*2+1,	-1 },
	{ "FDX" , IHALO*2-1,	-1 },
	{ "FDY" , JHALO*2-1,	-1 },
	{ "CXG" , IHALO*2, 	-1 },
	{ "CYG" , JHALO*2, 	-1 },
	{ "FXG" , IHALO*2+1, 	-1 },
	{ "FYG" , JHALO*2+1, 	-1 },
	{ "FDXG", IHALO*2-1, 	-1 },
	{ "FDYG", JHALO*2-1, 	-1 }
};

struct dim_pair hist_dims[NUM_HISTDIMS] = {
	{ "z"   , KMAX, 	-1 },
	{ "zh"  , KMAX+1, 	-1 },
	{ "oz"  , OKMAX, 	-1 },
	{ "ozh" , OKMAX+1,	-1 },
	{ "lz"  , LKMAX, 	-1 },
	{ "lzh" , LKMAX+1, 	-1 },
	{ "uz"  , UKMAX, 	-1 },
	{ "uzh" , UKMAX+1,	-1 },
	{ "x"   , 0,		-1 },
	{ "xh"  , 1,	 	-1 },
	{ "y"	, 0,		-1 },
	{ "yh"  , 1,		-1 },
	{ "CZ"  , KA, 		-1 },
	{ "FZ"  , KA+1,	 	-1 },
	{ "FDZ" , KA-1, 	-1 },
	{ "OCZ" , OKMAX, 	-1 },
	{ "OFZ" , OKMAX+1, 	-1 },
	{ "LCZ" , LKMAX, 	-1 },
	{ "LFZ" , LKMAX+1, 	-1 },
	{ "UCZ" , UKMAX, 	-1 },
	{ "UFZ" , UKMAX+1, 	-1 },
	{ "CX"  , IHALO*2, 	-1 },
	{ "CY"  , JHALO*2, 	-1 },
	{ "FX"  , IHALO*2+1, 	-1 },
	{ "FY"  , JHALO*2+1,	-1 },
	{ "FDX" , IHALO*2-1,	-1 },
	{ "FDY" , JHALO*2-1,	-1 },
	{ "CXG" , IHALO*2, 	-1 },
	{ "CYG" , JHALO*2, 	-1 },
	{ "FXG" , IHALO*2+1, 	-1 },
	{ "FYG" , JHALO*2+1, 	-1 },
	{ "FDXG", IHALO*2-1, 	-1 },
	{ "FDYG", JHALO*2-1, 	-1 },
	{ "time", NC_UNLIMITED, -1 },
	{ "nv"	, 2,		-1 }
};

/* ANAL vars SCALE read and LETKF read/write */
struct io_vars anal_vars[NUM_IOVARS_ANAL] = {
	{ "DENS",	-1 },
	{ "MOMX",	-1 },
	{ "MOMY",	-1 }, 
	{ "MOMZ",	-1 },
	{ "RHOT",	-1 },
	{ "QV",		-1 }, 
	{ "QC",		-1 },
	{ "QR",		-1 },
	{ "QI", 	-1 },
	{ "QS", 	-1 },
	{ "QG", 	-1 }
};

/* HIST vars LETKF read */
struct io_vars hist_vars[NUM_IOVARS_HIST] = {
	{ "U",		-1 }, 
	{ "V", 		-1 },
	{ "W", 		-1 },
	{ "T", 		-1 },
	{ "PRES",	-1 },
	{ "QV", 	-1 },
	{ "QC",		-1 },
	{ "QR",		-1 },
	{ "QI",		-1 },
	{ "QS",		-1 },
	{ "QG",		-1 },
	{ "RH",		-1 },
	{ "height", 	-1 },
	{ "topo",	-1 },
	{ "SFC_PRES",	-1 },
	{ "PREC",	-1 },
	{ "U10",	-1 },
	{ "V10",	-1 },
	{ "T2",		-1 },
	{ "Q2",		-1 }
};

MPI_Offset get_dim_length(struct file_info *file, char *dim_name)
{
	int i, ndims = file->ndims;

	for (i = 0; i < ndims; i++) {
		struct dim_pair *dim = &file->dims[i];
		if (!strcmp(dim_name, dim->name)) return dim->length;
	}

	return -1;
}

static void init_iovars(PD *pd)
{
	int i;
	struct file_info *file_anal = &pd->files[ANAL];
	struct file_info *file_hist = &pd->files[HIST];

	for (i = 0; i < NUM_IOVARS_ANAL; i++) {
		int idx;
		
		idx = find_var(file_anal, anal_vars[i].name);
		check_error(idx >= 0, find_var);
		anal_vars[i].idx = idx;
	}

	for (i = 0; i < NUM_IOVARS_HIST; i++) {
		int idx;
		
		idx = find_var(file_hist, hist_vars[i].name);
		check_error(idx >= 0, find_var);
		hist_vars[i].idx = idx;
	}
}

/* Only initiate write buffers */
static void init_file_buffers(PD *pd)
{
	int fi;

	init_iovars(pd);

	for (fi = 0; fi < pd->nfiles; fi++) {
		int i;
		struct file_info *file = &pd->files[fi];
		struct data_buf *var_read_buffers, *var_write_buffers;
		int num_vars = file->nvars;

		var_read_buffers = (struct data_buf *)malloc(sizeof(struct data_buf) * num_vars);
		check_error(var_read_buffers, malloc);
		memset(var_read_buffers, 0, sizeof(struct data_buf) * num_vars);

		var_write_buffers = (struct data_buf *)malloc(sizeof(struct data_buf) * num_vars);
		check_error(var_write_buffers, malloc);
		memset(var_write_buffers, 0, sizeof(struct data_buf) * num_vars);

		/* HIST+ANAL Write buffers: All the variables*/
		if (strstr(comp_name, "scale")) {
			for (i = 0; i < file->nvars; i++) {
				struct var_pair *var = &file->vars[i];
				struct data_buf *wbuf = &var_write_buffers[i];
				int ndims = var->ndims;
				MPI_Offset *shape = NULL;
				float *data = NULL;
				MPI_Offset nelems = 1;
				int varid = i;

				if (i < NUM_AXIS_VARS) {
					MPI_Offset dim_len = get_dim_length(file, var->dim_name[0]);
					MPI_Offset start = 0, count = 0;
					int j;

					if (pd->proc_rank_x == 0 && (strchr(var->dim_name[0], 'y'))) {
						if (fi == HIST) {
							if (strchr(var->dim_name[0], 'h')) {
								start = SYH_hist(pd);
								count = CYH_hist(pd);
							}
							else {
								start = SY_hist(pd);
								count = JMAX(pd);
							}
						}
						else if (fi == ANAL) {
							start = JSGA(pd);
							count = JEB(pd) - JSB(pd) + 1;
						}
					}

					if (pd->proc_rank_y == 0 && (strchr(var->dim_name[0], 'x'))) {
						if (fi == HIST) {
							if (strchr(var->dim_name[0], 'h')) {
								start = SXH_hist(pd);
								count = CXH_hist(pd);
							}
							else {
								start = SX_hist(pd);
								count = IMAX(pd);
							}
						} // ANAL
						else if (fi == ANAL) {
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

					nelems = count;

					if (count) {
						shape = (MPI_Offset *)malloc(sizeof(MPI_Offset) * ndims * 2);
						check_error(shape, malloc);

						shape[0] = start;
						shape[1] = count;

						data = (float *)malloc(sizeof(float) * count);
						check_error(data, malloc);
						memset(data, 0, sizeof(float) * count);
						
						for (j = 0; j < count; j++)
							data[j] = i + 1;
					}
				}
				else {
					int j;
					MPI_Offset *start = NULL, *count = NULL;

					/* time* vars in HIST are ignored */
					if (strstr(var->name, "time")) {
						wbuf->data = NULL;
						wbuf->shape = NULL;
						continue;
					}

					shape = (MPI_Offset *)malloc(sizeof(MPI_Offset) * ndims * 2);
					check_error(shape, malloc);
					memset(shape, 0, sizeof(MPI_Offset) * ndims * 2);

					start = shape;
					count = start + ndims;

					for (j = 0; j < ndims; j++) {
						if (strstr(var->dim_name[j], "time")) {
							start[j] = 0;
							count[j] = 1;
						}
						else if (strchr(var->dim_name[j], 'z')) {
							start[j] = 0;
							count[j] = get_dim_length(file, var->dim_name[j]);
						}
						else if (strchr(var->dim_name[j], 'y')) {
							if (fi == ANAL) {
								start[j] = JSGA(pd);
								count[j] = JEB(pd) - JSB(pd) + 1;
							}
							else {
								if (strchr(var->dim_name[j], 'h')) {
									start[j] = SYH_hist(pd);
									count[j] = CYH_hist(pd);
								}
								else {
									start[j] = SY_hist(pd);
									count[j] = JMAX(pd);
								}
							}
						}
						else if (strchr(var->dim_name[j], 'x')) {
							if (fi == ANAL) {
								start[j] = ISGA(pd);
								count[j] = IEB(pd) - ISB(pd) + 1;
							}
							else {
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
						nelems *= count[j];
					}

					data = (float *)malloc(sizeof(float) * nelems);
					check_error(data, malloc);
					memset(data, 0, sizeof(float) * nelems);

					for (j = 0; j < nelems; j++)
						data[j] = i + 1;
				}

				wbuf->data = data;
				wbuf->shape = shape;
				wbuf->ndims = ndims;
				wbuf->varid = varid;
				wbuf->nelems = nelems;
			}
		}
		/* LETKF write buffers: only anal_vars */
		else if (strstr(comp_name, "letkf")) {
			for (i = 0; i < NUM_IOVARS_ANAL; i++) {
				int j;
				int idx = anal_vars[i].idx;
				struct var_pair *var = &file->vars[idx];
				struct data_buf *wbuf = &var_write_buffers[idx];
				int ndims = var->ndims;
				MPI_Offset *start, *count;
				float *data = NULL;
				MPI_Offset nelems = 1;
				int varid = idx;

				start = (MPI_Offset *)malloc(sizeof(MPI_Offset) * ndims * 2);
				check_error(start, malloc);
				count = start + ndims;

				for (j = 0; j < ndims; j++) {
					if (strchr(var->dim_name[j], 'z')) {
						start[j] = 0;
						count[j] = KMAX;
					}
					else if (strchr(var->dim_name[j], 'y')) {
						start[j] = JS_inG(pd);
						count[j] = JMAX(pd);
					}
					else if (strchr(var->dim_name[j], 'x')) {
						start[j] = IS_inG(pd);
						count[j] = IMAX(pd);
					}

					nelems *= count[j];
				}

				data = (float *)malloc(sizeof(float) * nelems);
				check_error(data, malloc);
				memset(data, 0, sizeof(float) * nelems);
				for (j = 0; j < nelems; j++)
					data[j] = idx + 1;

				wbuf->data = data;
				wbuf->shape = start;
				wbuf->ndims = ndims;
				wbuf->varid = varid;
				wbuf->nelems = nelems;
			}
		}
		else {
			fprintf(stderr, "[ERROR] Invalid component name\n");
			MPI_Abort(MPI_COMM_WORLD, EINVAL);
		}

		file->var_read_buffers = var_read_buffers;
		file->var_write_buffers = var_write_buffers;
	}
}

static void init_fileinfo(PD *pd)
{
	int ret, i;
	char data_path[MAX_PATH_LEN] = { 0 };
	char *env = NULL;
	int len_path;

	/* Get path where data files will be stored */
	if ((env = getenv("DATA_PATH")) != NULL) {
		memcpy(data_path, env, MAX_PATH_LEN);
	}
	else {
		memcpy(data_path, DATA_PATH, MAX_PATH_LEN);
	}

	/* Add trailing slash behind path */
	len_path = strlen(data_path);
	if(data_path[len_path-1] != '/' && len_path < (MAX_PATH_LEN-1)) {
		data_path[len_path++] = '/';
		data_path[len_path] = '\0';
	}

	/* Create data folder if it does not exist */
	ret = create_dirs(data_path);
	check_error(!ret, create_dirs);

	/* Init file info */
	pd->files = (struct file_info *)calloc(pd->nfiles, sizeof(struct file_info));
	check_error(pd->files, calloc);

	for (i = 0; i < pd->nfiles; i++) {
		int j;
		FILE *fp = NULL;
		struct file_info *file = &pd->files[i];
		char info_fname[MAX_PATH_LEN] = { 0 };
		struct dim_pair *src_dim = NULL;

		file->nvars = 0;

		if (i == ANAL) {
			strcpy(info_fname, INFO_PATH(anal));
			src_dim = anal_dims;
			file->ndims = NUM_ANALDIMS;
			file->nvars = NUM_ANALVARS;
			file->nassct_coords = NUM_ANAL_ASSOCIATECOORD_VARS;
			file->ndata_vars = NUM_ANAL_DATA_VARS;
		}
		else if (i == HIST) {
			strcpy(info_fname, INFO_PATH(hist));
			src_dim = hist_dims;
			file->ndims = NUM_HISTDIMS;
			file->nvars = NUM_HISTVARS;
			file->nassct_coords = NUM_HIST_ASSOCIATECOORD_VARS;
			file->ndata_vars = NUM_HIST_DATA_VARS;
		}
		else {
			fprintf(stderr, "[ERROR] Unknown file type\n");
			MPI_Abort(MPI_COMM_WORLD, EINVAL);
		}

		/* Initialize dimension information */
		file->dims = (struct dim_pair *)malloc(file->ndims * sizeof(struct dim_pair));
		check_error(file->dims, malloc);

		/* Copy anal_dim/hist_dim --> file->dims */
		memcpy(file->dims, src_dim, sizeof(struct dim_pair) * file->ndims);
		
		/* Initialize variable information */
		file->vars = (struct var_pair *)malloc(file->nvars * sizeof(struct var_pair));
		check_error(file->vars, malloc);
		memset(file->vars, 0, sizeof(file->nvars * sizeof(struct var_pair)));

		fp = fopen(info_fname, "r");
		check_error(fp, fopen);
	
		/* get variable name, type, ndims and mapped dims name */
		for (j = 0; j < file->nvars; j++) {
			struct var_pair *var = &file->vars[j];
			int iter;
			int num_var_dims = 0;
			nc_type var_type;
			char var_name[NC_MAX_NAME+1] = { 0 };
			
			ret = fscanf(fp, "%s %d %d", var_name, &var_type, &num_var_dims);
			check_error(ret == 3, fscanf);
			
			strcpy(var->name, var_name);
			var->type = var_type;
			var->ndims = num_var_dims;

			var->dim_name = (char (*)[NC_MAX_NAME+1])malloc(
					(NC_MAX_NAME + 1) * num_var_dims);
			check_error(var->dim_name, malloc);
			memset(var->dim_name, 0, (NC_MAX_NAME+1)*num_var_dims);
			
			for (iter = 0; iter < num_var_dims; iter++) {
				ret = fscanf(fp, "%s", var->dim_name[iter]);
				check_error(ret == 1, fscanf);
			}

			/* Will be filled by pnetcdf later */
			var->dims = (int *)calloc(num_var_dims, sizeof(int));
			check_error(var->dims, calloc);

		}

		ret = fclose(fp);
		check_error(ret != EOF, fclose);

		/* calculate data buffer size for SCALE pnetcdf bput_*() */
		file->databuf_sz = get_databuf_size(pd, i);
		check_error(file->databuf_sz > 0, get_databuf_size);

		/* Fill pnetcdf file name for each cycle */
		file->file_names = (char (*)[NC_MAX_NAME+1])malloc((NC_MAX_NAME+1)*pd->cycles);
		check_error(file->file_names, malloc);

		char *suffix = (i == ANAL) ? ".anal.nc" : ".hist.nc";
		for (j = 0; j < pd->cycles; j++) {
			fmt_filename(j, pd->ens_id, 6, data_path, suffix, file->file_names[j]);
		}
	}

	init_file_buffers(pd);
}

void init_pd(int argc, char **argv, PD **p_pd)
{
	int i, ret = 0;
	int opt, proc_per_ens, color;
	int px = 0, py = 0; // 2D-proc map
	PD *pd = NULL;

	if (!p_pd) {
		fprintf(stderr, "[ERROR] Invalid address of PD");
		MPI_Abort(MPI_COMM_WORLD, EINVAL);
	}

	pd = (struct proc_data *)malloc(sizeof(struct proc_data));
	check_error(pd, malloc);
	memset(pd, 0, sizeof(PD));

	for (opt = 1; opt < argc; opt++) {
		char tmp[MAX_PARAM_LEN] = { 0 };

		/* An option */
		if (argv[opt][0] == '-') {
			int j = 1;
			int tmp_len;
			char *tmp_p = tmp;
			char key[512] = { 0 };
			char value[512] = { 0 };

			/* Ignore all the leading dashes */
			while (argv[opt][j] == '-') j++;
			tmp_len = strlen(argv[opt]) - j;
			memcpy(tmp, argv[opt] + j, tmp_len);
			tmp[tmp_len] = '\0';
			
			while(*tmp_p != '=') {
				tmp_p++;
				if (*tmp_p == '\0') {
					fprintf(stderr, "[ERROR] Invalid option format: %s\n", tmp);
					MPI_Abort(MPI_COMM_WORLD, EINVAL);
				}
			}
			*tmp_p = '\0';
			tmp_p++;
			memcpy(key, tmp, strlen(tmp));
			memcpy(value, tmp_p, strlen(tmp_p));

			if (!(strcmp(TOLOWER(key), "imax"))) {
				pd->imax = atoi(value);
			}
			else if (!(strcmp(TOLOWER(key), "jmax"))) {
				pd->jmax = atoi(value);
			}
			else if (!(strcmp(TOLOWER(key), "member"))) {
				pd->num_ens = atoi(value);
			}
			else if (!(strcmp(TOLOWER(key), "cycles"))){
				pd->cycles = atoi(value);
			}
			else if (!(strcmp(TOLOWER(key), "px"))){
				px = atoi(value);
			}
			else if (!(strcmp(TOLOWER(key), "py"))){
				py = atoi(value);
			}
			else {
				fprintf(stderr, "[ERROR] Unknown options! %s\n", tmp);
				MPI_Abort(MPI_COMM_WORLD, EINVAL);
			}
		}
	}
	
	if (pd->imax == 0) pd->imax = 32;
	if (pd->jmax == 0) pd->jmax = 20;
	if (pd->num_ens == 0) pd->num_ens = 2;
	if (pd->cycles == 0) pd->cycles = 4;
	pd->kmax = KMAX;
	pd->nfiles = 2;

	MPI_Comm_size(MPI_COMM_WORLD, &pd->world_size);
	MPI_Comm_rank(MPI_COMM_WORLD, &pd->world_rank);

	if (pd->world_size % pd->num_ens) {
		fprintf(stderr, "[ERROR]: number of processes should be dividable by 'member'\n");
		MPI_Abort(MPI_COMM_WORLD, EINVAL);
	}

	proc_per_ens = pd->world_size / pd->num_ens;
	color = pd->world_rank / proc_per_ens;
	pd->ens_id = color + 1;

	ret = MPI_Comm_split(MPI_COMM_WORLD, color, pd->world_rank, &pd->ens_comm);
	check_mpi(ret, MPI_Comm_split);

	MPI_Comm_rank(pd->ens_comm, &pd->ens_rank);
	MPI_Comm_size(pd->ens_comm, &pd->ens_size);
	
	if (!px || !py) {
		px = upper(sqrt(pd->ens_size));
		while (pd->ens_size % px) px++;
		py = pd->ens_size / px;
	}
	else if ((px * py) != pd->ens_size) {
		fprintf(stderr, "[ERROR] px [%d] * py [%d] is not equal to"
			" #process %d in each ensemble\n", px, py, pd->ens_size);
		MPI_Abort(MPI_COMM_WORLD, EINVAL);
	}

	pd->proc_num_x = px;
	pd->proc_num_y = py;
	pd->proc_rank_x = pd->ens_rank % pd->proc_num_x;
	pd->proc_rank_y = pd->ens_rank / pd->proc_num_x;

	if (!pd->world_rank) {
		fprintf(stdout, "%s: Number of Processes: %d Grid Size (IMAX * JMAX): %ld * %ld\n"
				"\tNumber of Ensembles: %d Number of Cycles: %d\n"
				"\tProcess Coordinate (X*Y): %d * %d\n",
				comp_name, pd->world_size, pd->imax, pd->jmax, pd->num_ens,
				pd->cycles, pd->proc_num_x, pd->proc_num_y);
	}

	/* Add imax and jmax to length of each PnetCDF variables */
	for (i = 0; i < NUM_ANALDIMS; i++) {
		if ((strchr(anal_dims[i].name, 'X')) || (strchr(anal_dims[i].name, 'x'))) {
			anal_dims[i].length += (pd->imax * pd->proc_num_x);
		}
		else if ((strchr(anal_dims[i].name, 'Y')) || (strchr(anal_dims[i].name, 'y'))) {
			anal_dims[i].length += (pd->jmax * pd->proc_num_y);
		}

		if ((strchr(hist_dims[i].name, 'X')) || (strchr(hist_dims[i].name, 'x'))) {
			hist_dims[i].length += (pd->imax * pd->proc_num_x);
		}
		else if ((strchr(hist_dims[i].name, 'Y')) || (strchr(hist_dims[i].name, 'y'))) {
			hist_dims[i].length += (pd->jmax * pd->proc_num_y);
		}
	}

	init_fileinfo(pd);

	/*
	 Include the cycle 0, which needs to initialize data structs
	 allocate memory for time of RW, time of R and time of W
	*/
	struct timing *t = &pd->time;
	t->cycle_transfer_time = (double *) malloc(sizeof(double) * pd->cycles * 6);
	check_error(t->cycle_transfer_time, malloc);
	memset(t->cycle_transfer_time, 0, sizeof(double) * pd->cycles * 6);

	t->cycle_transfer_rtime = t->cycle_transfer_time + pd->cycles;
	t->cycle_transfer_wtime = t->cycle_transfer_rtime + pd->cycles;
	t->cycle_file_time = t->cycle_transfer_wtime + pd->cycles;
	t->cycle_file_rtime = t->cycle_file_time + pd->cycles;
	t->cycle_file_wtime = t->cycle_file_rtime + pd->cycles;
	t->trans_checkpoint = 0.0;
	t->file_checkpoint = 0.0;

	*p_pd = pd;
}

static int free_datatype(MPI_Datatype *type)
{
	int ret;
	if (*type != MPI_DATATYPE_NULL) {
		int combiner, num_a, num_d, num_i;
		ret = MPI_Type_get_envelope(*type, &num_i, &num_a, &num_d, &combiner);
		check_mpi(ret, MPI_Type_get_envelope);

		if (combiner != MPI_COMBINER_NAMED) {
			ret = MPI_Type_free(type);
			check_mpi(ret, MPI_Type_free);
		}

		*type = MPI_DATATYPE_NULL;
	}
	return 0;
}

int finalize_pd(PD *pd)
{
	int i;

	for (i = 0; i < pd->nfiles; i++) {
		int j;
		struct file_info *file = &pd->files[i];

		if (file->var_write_buffers) {
			for (j = 0; j < file->nvars; j++) {
				struct data_buf *buf = &file->var_write_buffers[j];
				if (buf->shape) free(buf->shape);
				if (buf->data) free(buf->data);
				if (buf->ntypes == 1) free_datatype(&buf->dtype);
			}
			free(file->var_write_buffers);
		}
		
		if (file->var_read_buffers) {
			for (j = 0; j < file->nvars; j++) {
				struct data_buf *buf = &file->var_read_buffers[j];
				if (buf->shape) free(buf->shape);
				if (buf->data) free(buf->data);
				if (buf->ntypes == 1) free_datatype(&buf->dtype);
			}
			free(file->var_read_buffers);
		}

		if (file->dims) free(file->dims);
		if (file->vars) {
			for (j = 0; j < file->nvars; j++) {
				struct var_pair *var = &file->vars[j];
				if (var->dim_name) free(var->dim_name);
				if (var->dims) free(var->dims);
			}
			free(file->vars);
		}
		if (file->file_names) free(file->file_names);
	}

	if (pd->time.cycle_transfer_time) free(pd->time.cycle_transfer_time);
	MPI_Comm_free(&pd->ens_comm);
	free(pd->files);
	free(pd);

	return 0;
}

static double stand_devi(double myval, double sum, int nranks)
{
	int err;
	double tmpsum=0;
	double mean = sum/(double)nranks;
	double tmp = (myval - mean)*(myval - mean);

	if(myval == 0) tmp = 0;

	err = MPI_Reduce(&tmp, &tmpsum, 1, MPI_DOUBLE, MPI_SUM,0, MPI_COMM_WORLD);
	check_mpi(err, MPI_Reduce);

	return sqrt(tmpsum/(double)nranks);
}

void output_stat(PD *pd, char *comp_name)
{
	struct timing *time = &pd->time;
	int i, ret; int num_cycles = pd->cycles;
	double *cycle_transfer_time = NULL, *cycle_transfer_rtime, *cycle_transfer_wtime;
	double *cycle_file_time = NULL, *cycle_file_rtime, *cycle_file_wtime;
	double *std_devi_t = NULL, *std_devi_rt, *std_devi_wt;
	double *std_devi_f = NULL, *std_devi_rf, *std_devi_wf;

	/* allocate buffers for cycle time data */
	cycle_transfer_time = (double *)malloc(sizeof(double) * num_cycles * 6);
	check_error(cycle_transfer_time, malloc);
	memset(cycle_transfer_time, 0, sizeof(double) * num_cycles * 6);

	cycle_transfer_rtime = cycle_transfer_time + num_cycles;
	cycle_transfer_wtime = cycle_transfer_rtime + num_cycles;
	cycle_file_time = cycle_transfer_wtime + num_cycles;
	cycle_file_rtime = cycle_file_time + num_cycles;
	cycle_file_wtime = cycle_file_rtime + num_cycles;

	std_devi_t = (double *)malloc(sizeof(double) * num_cycles * 6);
	check_error(std_devi_t, malloc);
	memset(std_devi_t, 0, sizeof(double) * num_cycles * 6);
	std_devi_rt = std_devi_t + num_cycles;
	std_devi_wt = std_devi_rt + num_cycles;
	std_devi_f = std_devi_wt + num_cycles;
	std_devi_rf = std_devi_f + num_cycles;
	std_devi_wf = std_devi_rf + num_cycles;

	/* get sum of cycle time of all the processes
	 * (for calculating std devi and avg)
	 */
	ret = MPI_Allreduce(time->cycle_transfer_time, cycle_transfer_time, num_cycles,
			MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	check_mpi(ret, MPI_Allreduce);

	ret = MPI_Allreduce(time->cycle_transfer_rtime, cycle_transfer_rtime, num_cycles,
			MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	check_mpi(ret, MPI_Allreduce);

	ret = MPI_Allreduce(time->cycle_transfer_wtime, cycle_transfer_wtime, num_cycles,
			MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	check_mpi(ret, MPI_Allreduce);

	ret = MPI_Allreduce(time->cycle_file_time, cycle_file_time, num_cycles,
			MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	check_mpi(ret, MPI_Allreduce);

	ret = MPI_Allreduce(time->cycle_file_rtime, cycle_file_rtime, num_cycles,
			MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	check_mpi(ret, MPI_Allreduce);

	ret = MPI_Allreduce(time->cycle_file_wtime, cycle_file_wtime, num_cycles,
			MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	check_mpi(ret, MPI_Allreduce);

	for (i = 0; i < num_cycles; i++) {
		std_devi_t[i] = stand_devi(time->cycle_transfer_time[i], cycle_transfer_time[i], pd->world_size);
		std_devi_rt[i] = stand_devi(time->cycle_transfer_rtime[i], cycle_transfer_rtime[i], pd->world_size);
		std_devi_wt[i] = stand_devi(time->cycle_transfer_wtime[i], cycle_transfer_wtime[i], pd->world_size);
		std_devi_f[i] = stand_devi(time->cycle_file_time[i], cycle_file_time[i], pd->world_size);
		std_devi_rf[i] = stand_devi(time->cycle_file_rtime[i], cycle_file_rtime[i], pd->world_size);
		std_devi_wf[i] = stand_devi(time->cycle_file_wtime[i], cycle_file_wtime[i], pd->world_size);
	}

	if (!pd->world_rank) {
		double t_ct = 0.0, t_rct = 0.0, t_wct = 0.0;
		double f_ct = 0.0, f_rct = 0.0, f_wct = 0.0;

		for (i = 0; i < num_cycles; i++) {

			cycle_transfer_time[i]  /= (double)pd->world_size;
			cycle_transfer_rtime[i] /= (double)pd->world_size;
			cycle_transfer_wtime[i] /= (double)pd->world_size;
			cycle_file_time[i]  /= (double)pd->world_size;
			cycle_file_rtime[i] /= (double)pd->world_size;
			cycle_file_wtime[i] /= (double)pd->world_size;
		
			fprintf(stderr, "TRANSFER [%s] Cycle[%d]: Rd time: %.4f(%.4f)"
					" Wr time: %.4f(%.4f) Total: %.4f(%.4f)\n"
					"FILE [%s] Cycle[%d]: Rd time: %.4f(%.4f)"
					" Wr time: %.4f(%.4f) Total: %.4f(%.4f)\n",
					comp_name, i, cycle_transfer_rtime[i], std_devi_rt[i],
					cycle_transfer_wtime[i], std_devi_wt[i],
					cycle_transfer_time[i], std_devi_t[i],
					comp_name, i, cycle_file_rtime[i], std_devi_rf[i],
					cycle_file_wtime[i], std_devi_wf[i],
					cycle_file_time[i], std_devi_f[i]);

			/* ignore the initial cycle */
			if (i) {
				t_ct += cycle_transfer_time[i];
				t_rct += cycle_transfer_rtime[i];
				t_wct += cycle_transfer_wtime[i];
				f_ct += cycle_file_time[i];
				f_rct += cycle_file_rtime[i];
				f_wct += cycle_file_wtime[i];
			}
		}
		fprintf(stderr, "TRANSFER [%s] AVG Rd time: %.4f Wr time: %.4f Total: %.4f\n"
				"FILE [%s] AVG Rd time: %.4f Wr time: %.4f Total: %.4f\n",
				comp_name, t_rct / (double)(pd->cycles - 1),
				t_wct / (double)(pd->cycles - 1),
				t_ct / (double)(pd->cycles - 1), 
				comp_name, f_rct / (double)(pd->cycles - 1),
				f_wct / (double)(pd->cycles - 1),
				f_ct / (double)(pd->cycles - 1));
	}

	if (cycle_transfer_time) free(cycle_transfer_time);
	if (std_devi_t) free(std_devi_t);
}
