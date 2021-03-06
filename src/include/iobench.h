#ifndef _IOBENCH_H_
#define _IOBENCH_H_

#define FILE_CREATE	0x04
#define FILE_OPEN_W	0x08
#define FILE_OPEN_R	0x10

enum { XY = 0, ZXY2, ZHXY2, OCEAN, LAND, URBAN };
enum file_type { ANAL = 0, HIST };

struct dim_pair {
	char name[64];
	MPI_Offset length;
	int dimid;	// WRITTEN BY NCMPI CALL
};

struct var_pair {
	int *dims; 		// WRITTEN BY NCMPI CALL
	int varid; 		// WRITTEN BY NCMPI CALL
	nc_type type;		// READ FROM FILE
	unsigned int ndims;	// READ FROM FILE
	char **dim_name; // READ FROM FILE
	char name[64];		// READ FROM FILE
};

struct io_vars {
	char *name;
	int idx;
};

struct data_buf {
	// which variable this buffer belongs to
	int varid;
	// For getting the layout of shape and idxes
	int ndims;
	// Including HALO areas
	int nelems;

	// For SCALE's read_anal()
	int ntypes;
	MPI_Datatype dtype;
	/* shape in PnetCDF file
	 *
	 * [0 ... ndims-1] = start
	 * [ndims ... ndims*2-1] = count
	 */
	MPI_Offset *shape;

	/* NOTE: this region is not used for perf branch
	 * region without halo in data buffer
	 *
	 * [0 ... ndims-1] = start_idxes
	 * [ndims ... ndims*2-1] = end_idxes
	 */
	MPI_Offset *idxes;
	/* where data is stored */
	float *data;
};

struct file_info {
	int ndims;   // Total number of PnetCDF dimensions
	int nvars;     // Total number of PnetCDF variables
	int nassct_coords;  // Number of associated coordinate variables
	int ndata_vars;     // Number of data variables
	//struct data_buf *axes_buffer;
	MPI_Offset databuf_sz; // for bput_*
	struct data_buf *var_read_buffers;
	struct data_buf *var_write_buffers;
	struct dim_pair *dims;
	struct var_pair *vars;
	char **file_names;  // Path to output file for each cycle
};

struct timing {
	double trans_checkpoint;
	double file_checkpoint;
	double *cycle_transfer_time;
	double *cycle_transfer_rtime;
	double *cycle_transfer_wtime;
	double *cycle_file_time;
	double *cycle_file_rtime;
	double *cycle_file_wtime;
};

typedef struct proc_data {
	struct file_info *files;
	MPI_Offset imax;     	/* Size of data for the process at coord x without HALO*/
	MPI_Offset jmax;	/* Size of data for the process at coord y without HALO*/
	MPI_Offset kmax;	/* Size of data for the process at coord z without HALO*/
	MPI_Offset lkmax;
	MPI_Offset ukmax;
	MPI_Offset okmax;
	int num_ens;  		/* Number of ensembles */
	int cycles;
	int nfiles;
	int world_rank;
	int world_size;
	int ens_rank;
	int ens_size;
	int ens_id;
	int proc_num_x;
	int proc_num_y;
	int proc_rank_x;
	int proc_rank_y;
	MPI_Comm ens_comm;
	struct timing time;
} PD;

int create_dirs(char *path);
void fmt_filename(int cycle, int id, int total_chrs, char *prefix, char *suffix, char **p_name);
int prepare_file(struct file_info *file, MPI_Comm comm, char *file_path, int flag, int *ncid);
MPI_Offset get_databuf_size(PD *pd, int file_idx);
int fill_buffer(struct data_buf *buf, float c, float a, float w);
int compare_buffer(PD *pd, struct data_buf *buf, int cycle, float weight);
int find_var(struct file_info *file, char *var_name);
int init_data_buf(struct data_buf **buf, int num);
void cycle_transfer_start(PD *pd);
void cycle_transfer_end(PD *pd, int cycle);
void cycle_transfer_rend(PD *pd, int cycle);
void cycle_transfer_wend(PD *pd, int cycle);

void cycle_file_start(PD *pd);
void cycle_file_end(PD *pd, int cycle);
void cycle_file_rend(PD *pd, int cycle);
void cycle_file_wend(PD *pd, int cycle);

void report_put_size(PD *pd, int file_idx, int ncid);
void report_get_size(PD *pd, int file_idx, int ncid);

#endif // _IOBENCH_H_
