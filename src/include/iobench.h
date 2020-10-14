#ifndef _IOBENCH_H_
#define _IOBENCH_H_

#include "constants.h"

#define VAR_READ_NEVER		0x00
#define VAR_READ_ONCE		0x01
#define VAR_READ_ALWAYS    	0x02

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
	int rflag;
	int wflag;
	nc_type type;		// READ FROM FILE
	unsigned int ndims;	// READ FROM FILE
	char (*dim_name)[NC_MAX_NAME+1]; // READ FROM FILE
	char name[64];		// READ FROM FILE
};

/* TODO: replace each buffer in file_info with this struct */
struct data_buffer {
	int num_buf;
	MPI_Offset *size_list;
	float **buffer;
};

struct file_info {
	int ndims;   // Total number of dimensions
	int nvars;     // Total number of variables
	int nassct_coords;  // Number of associated coordinate variables
	int ndata_vars;     // Number of data variables
	int naxes_buf;  // Number of axes buffer
	int nvar_read_buf;   // Number of data read buffer
	int nvar_write_buf;  // Number of data write buffer
	float **axes_buffer;
	float **var_read_buffers;
	float **var_write_buffers;
	struct dim_pair *dims;
	struct var_pair *vars;
};

typedef struct proc_data {
	struct file_info *files;
	MPI_Offset imax;     	/* Size of data for the process at coord x without HALO*/
	MPI_Offset jmax;	/* Size of data for the process at coord y without HALO*/	
	MPI_Offset kmax;	/* Size of data for the process at coord z without HALO*/
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
} PD;

#endif // _IOBENCH_H_
