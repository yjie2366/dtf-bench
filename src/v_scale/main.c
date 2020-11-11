#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include "util.h"

int write_hist(PD *pd, int cycle);
int write_anal(PD *pd, int cycle);
int read_anal(PD *pd, int cycle);

char *comp_name = "v_scale";

int main(int argc, char **argv)
{
	int cycle = 0, ret;
	PD *pd = NULL;

	/* INITIALIZATION */
	MPI_Init(&argc, &argv);
	ret = dtf_init(DTF_INIT_FILE, comp_name);
	check_error(!ret, dtf_init);

	pd = (struct proc_data *)malloc(sizeof(struct proc_data));
	check_error(pd, malloc);

	init_pd(argc, argv, pd);

	/* Start Cycling */
	for (cycle = 0; cycle < pd->cycles; cycle++) {
//		write_hist(pd, cycle);
		write_anal(pd, cycle);
		read_anal(pd, cycle); 
	}

	output_stat(pd, comp_name);
	finalize_pd(pd);
	free(pd);

	ret = dtf_finalize();
	check_error(!ret, dtf_finalize);
	
	MPI_Finalize();

	return 0;
}


