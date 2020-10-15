#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include "util.h"

int main(int argc, char **argv)
{
	int cycle, ret;
	int len_path;
	char *env = NULL;
	char data_path[MAX_PATH_LEN] = { 0 };
	PD *pd = NULL;

	/* INITIALIZATION */
	MPI_Init(&argc, &argv);
	ret = dtf_init(DTF_INIT_FILE, "v_scale");
	check_error(!ret, dtf_init);

	pd = malloc(sizeof(struct proc_data));
	check_error(pd, malloc);

	init_pd(argc, argv, pd);

	if ((env = getenv("INIT_DATA_PATH")) != NULL) {
		memcpy(data_path, env, MAX_PATH_LEN);
	}
	else {
		memcpy(data_path, DATA_PATH, MAX_PATH_LEN);
	}

	len_path = strlen(data_path);
	if(data_path[len_path-1] != '/' && len_path < (MAX_PATH_LEN-1)) {
		data_path[len_path++] = '/';
		data_path[len_path] = '\0';
	}
	
	/* Prepare init data */
	write_anal(pd, data_path, 0);

	MPI_Barrier(MPI_COMM_WORLD);

	/* Start Cycling */
	for (cycle = 1; cycle < pd->cycles; cycle++) {
		write_hist(pd, data_path, cycle);
		write_anal(pd, data_path, cycle);
		read_anal(pd, data_path, cycle); 
	}

	finalize_pd(pd);
	free(pd);

	ret = dtf_finalize();
	check_error(!ret, dtf_finalize);
	
	MPI_Finalize();

	return 0;

}


