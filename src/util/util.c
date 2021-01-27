#include <string.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <dirent.h>
#include "util.h"

extern struct dim_pair anal_dims[];
extern struct dim_pair hist_dims[];

static float pseudo_rand = 1.0;

static void reset_seed(float c, float a, int seed)
{
	pseudo_rand = (seed % 1000) * (++c) + a;
}

static float get_float(float c, float w)
{
	return (pseudo_rand = (float)((int)pseudo_rand % 1000) * (++c) + w);
}

int fill_buffer(struct data_buf *buf, float c, float a, float w)
{
	int varid = buf->varid;
	int ndims = buf->ndims;
	MPI_Offset d1, d2, d3, d4;
	MPI_Offset *count = buf->shape + ndims;
	MPI_Offset *s_idx = buf->idxes;
	MPI_Offset *e_idx = s_idx + ndims;
	float *data = buf->data;

	reset_seed(c, a, varid);

	if (ndims == 1) {
		for (d1 = s_idx[0]; d1 < e_idx[0]; d1++)
			data[d1] = get_float(c, w);
	}
	else if (ndims == 2) {
		for (d1 = s_idx[0]; d1 < e_idx[0]; d1++)
		for (d2 = s_idx[1]; d2 < e_idx[1]; d2++) {
			MPI_Offset idx = d1 * count[1] + d2;
			data[idx] = get_float(c, w);
		}
	}
	else if (ndims == 3) {
		for (d1 = s_idx[0]; d1 < e_idx[0]; d1++)
		for (d2 = s_idx[1]; d2 < e_idx[1]; d2++)
		for (d3 = s_idx[2]; d3 < e_idx[2]; d3++) {
			MPI_Offset idx = d1 * count[1] * count[2] 
					+ d2 * count[2] + d3;
			data[idx] = get_float(c, w);
		}
	}
	else if (ndims == 4) {
		for (d1 = s_idx[0]; d1 < e_idx[0]; d1++)
		for (d2 = s_idx[1]; d2 < e_idx[1]; d2++)
		for (d3 = s_idx[2]; d3 < e_idx[2]; d3++)
		for (d4 = s_idx[3]; d4 < e_idx[3]; d4++) {
			MPI_Offset idx = d1 * count[1] * count[2] * count[3]
					+ d2 * count[2] * count[3]
					+ d3 * count[3] + d4;
			data[idx] = get_float(c, w);
		}
	}
	else {
		fprintf(stderr, "[ERROR] Invalid number of dimension %d\n", ndims);
		errno = EINVAL;

		return errno;
	}
	
	return 0;	
}

int compare_buffer(PD *pd, struct data_buf *buf, int cycle, float weight)
{
	int j;
	int varid = buf->varid;
	int ndims = buf->ndims;
	MPI_Offset d1, d2, d3, d4, cnt = 0;

	MPI_Offset *count = buf->shape + ndims;
	MPI_Offset *s_idx = buf->idxes;
	MPI_Offset *e_idx = s_idx + ndims;
	float *data = buf->data;
	float *expect = NULL;
	int rank = pd->world_rank;

	MPI_Offset nelems = 1;

	if (!buf || !data) {
		fprintf(stderr, "[ERROR] Invalid buffer!\n");
		errno = EINVAL;

		return errno;
	}

	if ((ndims > 4) || (ndims < 2)) {
		fprintf(stderr, "[ERROR] Comparing %d-D array unsupported\n", ndims);
		errno = EINVAL;

		return errno;
	}
	
	for (j = 0; j < buf->ndims; j++)
		nelems *= (e_idx[j] - s_idx[j]);

	expect = (float *)malloc(sizeof(float) * nelems);
	check_error(expect, malloc);

	reset_seed(rank, cycle, varid);
	for (j = 0; j < nelems; j++) expect[j] = get_float(rank, weight);

	if (ndims == 2) {
		for (d1 = s_idx[0]; d1 < e_idx[0]; d1++)
		for (d2 = s_idx[1]; d2 < e_idx[1]; d2++) {
			MPI_Offset idx = d1 * count[1] + d2;
			float exp = expect[cnt++];

			if (exp != data[idx]) {
				fprintf(stderr, "[ERROR] Unmatched value of var ID "
						"%d [%ld, %ld]\n",
						varid, d1, d2);
				fprintf(stderr, "Expect: %f Acquired: %f\n",
					exp, data[idx]);
				errno = EINVAL;
				return errno;
			}
		}
	}
	else if (ndims == 3) {
		for (d1 = s_idx[0]; d1 < e_idx[0]; d1++)
		for (d2 = s_idx[1]; d2 < e_idx[1]; d2++)
		for (d3 = s_idx[2]; d3 < e_idx[2]; d3++) {
			MPI_Offset idx =  d1 * count[1] * count[2] 
					+ d2 * count[2] + d3;
			float exp = expect[cnt++];

			if (exp != data[idx]) {
				fprintf(stderr, "[ERROR] Unmatched value of var ID "
						"%d [%ld, %ld, %ld]\n",
						varid, d1, d2, d3);
				fprintf(stderr, "Expect: %f Acquired: %f\n",
					exp, data[idx]);
				errno = EINVAL;
				return errno;
			}
		}
	}
	else if (ndims == 4) {
		for (d1 = s_idx[0]; d1 < e_idx[0]; d1++)
		for (d2 = s_idx[1]; d2 < e_idx[1]; d2++)
		for (d3 = s_idx[2]; d3 < e_idx[2]; d3++)
		for (d4 = s_idx[3]; d4 < e_idx[3]; d4++) {
			MPI_Offset idx =  d1 * count[1] * count[2] * count[3]
					+ d2 * count[2] * count[3]
					+ d3 * count[3] + d4;
			float exp = expect[cnt++];

			if (exp != data[idx]) {
				fprintf(stderr, "[ERROR] Unmatched value of var ID "
						"%d [%ld, %ld, %ld]\n",
						varid, d1, d2, d3);
				fprintf(stderr, "Expect: %f Acquired: %f\n",
					exp, data[idx]);
				errno = EINVAL;
				return errno;
			}
		}
	}

	free(expect);

	return 0;
}

int create_dirs(char *path)
{
	int offset = 0;
	char tmp[MAX_PATH_LEN] = { 0 };
	char *pt = NULL;
	DIR *d = NULL;

	if (strlen(path) >= MAX_PATH_LEN) {
		fprintf(stderr, "[ERROR] %s Path %s is too long to process.\n", __func__, path);
		errno = EINVAL;

		return errno;
	}
	
	/* If directory exists, return immediately */
	d = opendir(path);
	if (d) {
		closedir(d);
		return 0;
	}

	/* Otherwise, mkdir -p */
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
	name[off] = '\0';
}

MPI_Offset get_databuf_size(PD *pd, int file_idx)
{
	struct file_info *file = NULL;
	MPI_Offset total_size = 0;
	int idx;

	if (!pd->files) {
		fprintf(stderr, "[ERROR]: files info are not filled yet\n");
		errno = EINVAL;

		return -errno;
	}
	
	file = &pd->files[file_idx];
	if (!file->dims || !file->vars) {
		fprintf(stderr, "[ERROR]: dims and vars are not filled yet\n");
		errno = EINVAL;

		return -errno;
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
			else {
				fprintf(stderr, "[ERROR] Unrecoginized dimension: %s\n", name);
				errno = EINVAL;

				return -errno;
			}
			total_dim_size *= dim_size;
		}
		total_size += total_dim_size;
	}

	return total_size * sizeof(float);
}

int prepare_file(struct file_info *file, MPI_Comm comm, char *file_path, int flag, int *ncid)
{
	int i, ret, fid;
	
	if (flag & FILE_CREATE) {
		ret = ncmpi_create(comm, file_path, NC_CLOBBER | NC_64BIT_OFFSET,
					MPI_INFO_NULL, &fid);
		check_io(ret, ncmpi_create);
	}
	else {
		int mode = (flag & FILE_OPEN_W) ? NC_WRITE : NC_NOWRITE;

		ret = ncmpi_open(comm, file_path, mode, MPI_INFO_NULL, &fid);
		check_io(ret, ncmpi_open);
	}
	
	for (i = 0; i < file->ndims; i++) {
		struct dim_pair *dim = &file->dims[i];

		if (flag & FILE_CREATE) {
			ret = ncmpi_def_dim(fid, dim->name, dim->length, &dim->dimid);
			check_io(ret, ncmpi_def_dim);
		}
		else {
			ret = ncmpi_inq_dimid(fid, dim->name, &dim->dimid);
			check_io(ret, ncmpi_inq_dimid);
		}
	}

	for (i = 0; i < file->nvars; i++) {
		int j;
		struct var_pair *var = &file->vars[i];

		for (j = 0; j < var->ndims; j++) {
			ret = ncmpi_inq_dimid(fid, var->dim_name[j], &var->dims[j]);
			check_io(ret, ncmpi_inq_dimid);
		}

		if (flag & FILE_CREATE) {
			ret = ncmpi_def_var(fid, var->name, var->type, var->ndims,
					var->dims, &var->varid);
			check_io(ret, ncmpi_def_var);
		}
		else {
			ret = ncmpi_inq_varid(fid, var->name, &var->varid);
			check_io(ret, ncmpi_inq_varid);
		}
	}

	*ncid = fid;

	return 0;
}

int find_var(struct file_info *file, char *var_name)
{
	int idx = -1, j;

	for (j = NUM_AXIS_VARS; j < file->nvars; j++) {
		char *name = file->vars[j].name;

		if (strcmp(name, var_name)) continue;
		idx = j; 
		break;
	}

	return idx;
}

void cycle_transfer_start(PD *pd)
{
	if (!pd) {
		fprintf(stderr, "[ERROR] PD pointer is NULL!\n");
		return;
	}
	pd->time.trans_checkpoint = MPI_Wtime();
}

void cycle_transfer_end(PD *pd, int cycle)
{
	double t = MPI_Wtime();
	double dur = t - pd->time.trans_checkpoint;

	if (!pd) {
		fprintf(stderr, "[ERROR] PD pointer is NULL!\n");
		return;
	}
	if (cycle < 0) {
		fprintf(stderr, "[ERROR] Illegal cycle ID\n");
		return;
	}
	
	pd->time.cycle_transfer_time[cycle] += dur;
	pd->time.trans_checkpoint = t;
}

void cycle_transfer_rend(PD *pd, int cycle)
{
	double t = MPI_Wtime();
	double dur = t - pd->time.trans_checkpoint;

	if (!pd) {
		fprintf(stderr, "[ERROR] %s: PD pointer is NULL!\n", __func__);
		return;
	}
	if (cycle < 0) {
		fprintf(stderr, "[ERROR] %s: Illegal cycle ID\n", __func__);
		return;
	}
	
	pd->time.cycle_transfer_rtime[cycle] += dur;
	pd->time.cycle_transfer_time[cycle] += dur;
	pd->time.trans_checkpoint = t;
}

void cycle_transfer_wend(PD *pd, int cycle)
{
	double t = MPI_Wtime();
	double dur = t - pd->time.trans_checkpoint;

	if (!pd) {
		fprintf(stderr, "[ERROR] %s: PD pointer is NULL!\n", __func__);
		return;
	}
	if (cycle < 0) {
		fprintf(stderr, "[ERROR] %s: Illegal cycle ID\n", __func__);
		return;
	}
	
	pd->time.cycle_transfer_wtime[cycle] += dur;
	pd->time.cycle_transfer_time[cycle] += dur;
	pd->time.trans_checkpoint = t;
}

void cycle_file_start(PD *pd)
{
	if (!pd) {
		fprintf(stderr, "[ERROR] PD pointer is NULL!\n");
		return;
	}
	pd->time.file_checkpoint = MPI_Wtime();
}

void cycle_file_end(PD *pd, int cycle)
{
	double t = MPI_Wtime();
	double dur = t - pd->time.file_checkpoint;

	if (!pd) {
		fprintf(stderr, "[ERROR] PD pointer is NULL!\n");
		return;
	}
	if (cycle < 0) {
		fprintf(stderr, "[ERROR] Illegal cycle ID\n");
		return;
	}
	
	pd->time.cycle_file_time[cycle] += dur;
	pd->time.file_checkpoint = t;
}

void cycle_file_rend(PD *pd, int cycle)
{
	double t = MPI_Wtime();
	double dur = t - pd->time.file_checkpoint;

	if (!pd) {
		fprintf(stderr, "[ERROR] %s: PD pointer is NULL!\n", __func__);
		return;
	}
	if (cycle < 0) {
		fprintf(stderr, "[ERROR] %s: Illegal cycle ID\n", __func__);
		return;
	}
	
	pd->time.cycle_file_rtime[cycle] += dur;
	pd->time.cycle_file_time[cycle] += dur;
	pd->time.file_checkpoint = t;
}

void cycle_file_wend(PD *pd, int cycle)
{
	double t = MPI_Wtime();
	double dur = t - pd->time.file_checkpoint;

	if (!pd) {
		fprintf(stderr, "[ERROR] %s: PD pointer is NULL!\n", __func__);
		return;
	}
	if (cycle < 0) {
		fprintf(stderr, "[ERROR] %s: Illegal cycle ID\n", __func__);
		return;
	}
	
	pd->time.cycle_file_wtime[cycle] += dur;
	pd->time.cycle_file_time[cycle] += dur;
	pd->time.file_checkpoint = t;
}

void report_put_size(PD *pd, int file_idx, int ncid)
{
	int ret;
	MPI_Offset io_size;
	char *msg = (file_idx == ANAL) ? "ANAL write amount" : "HIST write amount";

	ret = ncmpi_inq_put_size(ncid, &io_size);
	check_io(ret, ncmpi_inq_put_size);

	if (!pd->world_rank) {
		fprintf(stderr, "%s: %ld\n", msg, io_size);
	}
}

void report_get_size(PD *pd, int file_idx, int ncid)
{
	int ret;
	MPI_Offset io_size;
	char *msg = (file_idx == ANAL) ? "ANAL read amount" : "HIST read amount";

	ret = ncmpi_inq_get_size(ncid, &io_size);
	check_io(ret, ncmpi_inq_get_size);

	if (!pd->world_rank) {
		fprintf(stderr, "%s: %ld\n", msg, io_size);
	}
}
