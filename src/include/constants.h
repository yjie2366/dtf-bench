// NOTE: ZX dim type is not included
#ifndef _CONSTANTS_H_
#define _CONSTANTS_H_

#define SCALE_WEIGHT	((float)0.0)
#define LETKF_WEIGHT 	((float)5.0)

#define NUM_DATATYPE	6
#define NUM_ANALDIMS	33
#define NUM_HISTDIMS	35

#define NUM_ANALVARS	143
#define NUM_HISTVARS	89

#define NUM_AXIS_VARS	51
#define NUM_ANAL_ASSOCIATECOORD_VARS 10
#define NUM_HIST_ASSOCIATECOORD_VARS 18

#define ANAL_DATA_VARS_OFFSET	(NUM_AXIS_VARS+NUM_ANAL_ASSOCIATECOORD_VARS)
// exclude time and time_bnds
#define HIST_DATA_VARS_OFFSET	(NUM_AXIS_VARS+NUM_HIST_ASSOCIATECOORD_VARS + 2)
#define NUM_ANAL_DATA_VARS 		(NUM_ANALVARS - ANAL_DATA_VARS_OFFSET)
#define NUM_HIST_DATA_VARS 		(NUM_HISTVARS - HIST_DATA_VARS_OFFSET)

/* Index & variable size */
#define KHALO	2
#define IHALO	2
#define JHALO	2

#define KMAX 	98
#define OKMAX 	1
#define LKMAX 	5
#define UKMAX 	5

#define IMAX(pd) (pd->imax)
#define JMAX(pd) (pd->jmax)

#define IMAXG(pd) (IMAX(pd) * pd->proc_num_x)
#define JMAXG(pd) (JMAX(pd) * pd->proc_num_y)

// ensemble size + halo zone
#define IAG(pd) (IMAXG(pd) + IHALO * 2)
#define JAG(pd) (JMAXG(pd) + JHALO * 2)

#define KA	KMAX + KHALO * 2
#define IA(pd) (IMAX(pd) + IHALO * 2)
#define JA(pd) (JMAX(pd) + JHALO * 2)

#define KS	KHALO
#define KE	(KMAX + KHALO - 1)
#define IS	IHALO
#define IE(pd) (IMAX(pd) + IHALO - 1)
#define JS	JHALO
#define JE(pd) (JMAX(pd) + JHALO - 1)
// Ocean
#define OKS	0
#define OKE	(OKMAX - 1)
// Land
#define LKS	0
#define LKE	(LKMAX - 1)
// Urban
#define UKS	0
#define UKE	(UKMAX - 1)

// global index
#define IS_inG(pd) (IHALO + pd->proc_rank_x * IMAX(pd))
#define IE_inG(pd) (IS_inG(pd) + IMAX(pd) - 1)
#define JS_inG(pd) (JHALO + pd->proc_rank_y * JMAX(pd))
#define JE_inG(pd) (JS_inG(pd) + JMAX(pd) - 1)

/* Global start indices */
#define ISGA(pd)  ({\
	int rank = pd->proc_rank_x, idx = 0; \
	if (rank != 0) { idx = IS_inG(pd); }\
	idx; })

#define IEGA(pd)  ({\
	int rank = pd->proc_rank_x, idx = IAG(pd) - 1 ; \
	if (rank != (pd->proc_num_x-1)) { idx = IE_inG(pd); }\
	idx; })

#define JSGA(pd)  ({\
	int rank = pd->proc_rank_y, idx = 0;\
	if (rank != 0) { idx = JS_inG(pd); }\
	idx; })

#define JEGA(pd)  ({\
	int rank = pd->proc_rank_y, idx = JAG(pd) - 1;\
	if (rank != (pd->proc_num_y-1)) { idx = JE_inG(pd); }\
	idx; })

/* Local start and end indices */
#define ISB(pd)	  ({\
	int rank = pd->proc_rank_x, idx = IS;\
	if (rank == 0) idx = 0;\
	idx; })

#define IEB(pd)	  ({\
	int rank = pd->proc_rank_x, edge = pd->proc_num_x-1, idx = IE(pd);\
	if (rank == edge) idx = IA(pd)-1;\
	idx; })

#define JSB(pd)   ({\
	int rank = pd->proc_rank_y, idx = JS;\
	if (rank == 0) idx = 0;\
	idx; })

#define JEB(pd)	  ({\
	int rank = pd->proc_rank_y, edge = pd->proc_num_y-1, idx = JE(pd);\
	if (rank == edge) idx = JA(pd)-1;\
	idx; })

/* For history file */
#define SX_hist(pd) (pd->proc_rank_x * IMAX(pd))
#define SY_hist(pd) (pd->proc_rank_y * JMAX(pd))

#define SXH_hist(pd) ({\
		int rank = pd->proc_rank_x, idx = SX_hist(pd);\
		if (rank) { idx += 1; }\
		idx; })
#define SYH_hist(pd) ({\
		int rank = pd->proc_rank_y, idx = SY_hist(pd);\
		if (rank) { idx += 1; }\
		idx; })

#define CXH_hist(pd) ({\
		int rank = pd->proc_rank_x, count = IMAX(pd);\
		if (rank == 0) { count += 1; }\
		count;})

#define CYH_hist(pd) ({\
		int rank = pd->proc_rank_y, count = JMAX(pd);\
		if (rank == 0) { count += 1; }\
		count;})


#define MAX_PATH_LEN		4096
#define MAX_NAME_LEN 		1024
#define MAX_PARAM_LEN		1024

#define TIME_INT		30

#endif
