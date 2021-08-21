script_dir=$(cd $(dirname ${BASH_SOURCE[0]}) &> /dev/null && pwd)

log_dir=${script_dir}/../log
if [ ! -e ${log_dir} ]; then
	mkdir -p ${log_dir} || exit 1
fi

args=()
group=($(id -nG))
master=
target=
ppn=
nnodes=
MCK_PATH=
MCK_MEM=

# Switch for LLIO 
enable_llio=0
enable_mckernel=0

bin_scale=${script_dir}/../bin/scale
bin_letkf=${script_dir}/../bin/letkf
info_anal=${script_dir}/../info/anal.info
info_hist=${script_dir}/../info/hist.info

usage()
{
cat << EOF
$0 [OPTION_1] [ARG_1] ...
    If no options and arguments explicitly specified,
    a job with default settings will be submitted.    

    -h  Print help message

    -i	Atmos latitude (IMAX)
    -j	Atmos longitude (JMAX)
    -k  Atmos Height (KMAX)
    -l  Land Height (LKMAX)
    -u  Urban Height (UKMAX)
    -o  Ocean Height (OKMAX)

    -t  Transfer mode [file|transfer]
    -n	Process-per-node (PPN)
    -p	Number of processes per component
    -x	Number of processes on x-coord
    -y	Number of processes on y-coord
    -m	Number of ensembles
    -c	Number of I/O cycles
    -a	Number of matcher processes in DTF
    -r	Number of runs in a row

    -b  Enable LLIO for DTF file-IO mode
    -d  Enable McKernel execution
EOF
}

# Need to add another case for Fugaku
case `hostname` in
	*ofp*)
		target=ofp
		;;
	fn01sv*)
		target=fugaku
		;;
	*)
		echo "[ERROR] Unsupported Machine"
		exit 1
		;;
esac

while getopts "n:p:i:j:k:l:u:o:c:m:a:r:x:y:t:d:hb" OPT; do
	case ${OPT} in
	n)
		# Number of processes per node
		ppn=${OPTARG}
		;;
	p)
		# Number of processes per comp
		nprocs=${OPTARG}
		;;
	i)
		args+=("-imax ${OPTARG}")
		;;
	j)
		args+=("-jmax ${OPTARG}")
		;;
	k)
		args+=("-kmax ${OPTARG}")
		;;
	l)
		args+=("-lkmax ${OPTARG}")
		;;
	u)
		args+=("-ukmax ${OPTARG}")
		;;
	o)
		args+=("-okmax ${OPTARG}")
		;;
	c)
		args+=("-cycles ${OPTARG}")
		;;
	m)
		member=${OPTARG}
		;;
	a)
		master=${OPTARG}
		;;
	t)
		args+=("-mode ${OPTARG}")
		;;
	r)
		args+=("-run ${OPTARG}")
		;;
	x)
		px=${OPTARG}
		args+=("-px ${px}")
		;;
	y)
		py=${OPTARG}
		args+=("-py ${py}")
		;;
	h)
		usage
		exit
		;;
	b)
		enable_llio=1
		args+=("-llio")
		;;
	d)
		enable_mckernel=1
		if [ "${target}" = "ofp" ]; then
			MCK_PATH=${OPTARG}
			MCK_MEM="48G@0,all@1"
		fi
		args+=("-mck ${MCK_PATH}")
		;;
	?)
		echo "[ERROR] Invalid option"
		exit 1
	esac
done

if [ -z "${master}" ]; then master=2; fi
if [ -z "${member}" ]; then member=2; fi

if [ -z "${nprocs}" ]; then
	if [ -z "${px}" -a -z "${py}" ]; then
		nprocs=4
	else
		nprocs=$((px*py*member))
	fi
fi
# Add number of process per comp into argument list
args+=("-member ${member}")
args+=("-master ${master}")
args+=("-np ${nprocs}")
total_nprocs=$((nprocs*2))

if [ -z "${ppn}" ]; then
	ppn=1; nnodes=${total_nprocs}
else
	nnodes=$((total_nprocs/ppn))
fi
args+=(-ppn ${ppn})

if [ "${target}" = "ofp" ]; then
	#rsc_args="rscgrp=debug-cache"
	rsc_args="rscgrp=regular-flat"
elif [ "${target}" = "fugaku" ]; then
	if [ ${enable_mckernel} -eq 0 ]; then
		if [ ${nnodes} -gt 385 ]; then
			#rsc_args="rscgrp=eap-large"
			rsc_args="rscgrp=large"
		else
			#rsc_args="rscgrp=eap-small"
			rsc_args="rscgrp=small"
		fi
	else
		# mckernel resource group
		rsc_args="rscunit=rscunit_ft01,rscgrp=dvsys-mck1,jobenv=linux2"
	fi
else
	echo "[ERROR] Unsupported Machine"
	exit 1
fi

# Batch script variables
batch_script="${script_dir}/batch.${target}.${nprocs}.${member}"
elapse_time="03:30:00"
jobname="d-${nprocs}-${member}"
runlog_dir="${log_dir}/run-${nprocs}-${master}-${member}"
if [ ! -d "runlog_dir" ]; then
	mkdir -p ${runlog_dir}
fi
output="${runlog_dir}/%n.%j.out"
error="${runlog_dir}/%n.%j.err"
stat="${runlog_dir}/%n.%j.stat"

cat <<- EOF > ${batch_script}
#!/bin/bash
#
#PJM -N "${jobname}"
#PJM -L "node=${nnodes}"
#PJM -L "${rsc_args}"
#PJM -L "elapse=${elapse_time}"
#PJM -g ${group[-2]}
#PJM -S
#PJM --spath ${stat}
#PJM -o ${output}
#PJM -e ${error}
#PJM --mpi "proc=${total_nprocs}"
`if [ ${enable_llio} -eq 1 ]; then
echo -e "#PJM --llio sharedtmp-size=10Gi"
else echo "#"
fi`
`if [ ${enable_mckernel} -eq 1 -a ${target} = "ofp" ]; then
echo "#PJM -x MCK=\"${MCK_PATH}\""
echo "#PJM -x MCK_MEM=\"${MCK_MEM}\""
else
echo "#PJM --mpi \"max-proc-per-node=${ppn}\""
fi`

sh ${script_dir}/exec.sh ${args[@]}

EOF

ret="$(pjsub ${batch_script})"
retval=$?
if [ "$retval" -eq 0 ]; then
	jobid=`echo $ret | grep -i "submitted" | awk '{print $6}'`
	echo "JOB ${jobid} submitted."

	pjwait ${jobid}
else
	echo "[ERROR] Error happened in pjsub with return value: ${retval}"
	exit
fi
