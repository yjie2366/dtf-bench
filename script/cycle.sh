script_dir=$(cd $(dirname ${BASH_SOURCE[0]}) &> /dev/null && pwd)
log_dir=${script_dir}/../log
args=()
group=($(id -nG))
master=
target=
ppn=
nnodes=
mck=${mck:=0}

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

    -n	Process-per-node
    -p	Number of processes per component
    -x	Number of processes on x-coord
    -y	Number of processes on y-coord
    -m	Number of ensembles
    -c	Number of I/O cycles
    -a	Number of matcher processes in DTF
    -r	Number of runs in a row
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

if [ ! -e ${log_dir} ]; then
	mkdir -p ${log_dir} || exit 1
fi

while getopts "n:p:i:j:k:l:u:o:c:m:a:r:x:y:h" OPT; do
	case ${OPT} in
	n)
		# Number of processes per node
		ppn=${OPTARG}
		;;
	p)
		# Number of processes per comp
		nprocs=${OPTARG}
		args+=("-np ${nprocs}")
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
		args+=("-member ${OPTARG}")
		;;
	a)
		master=${OPTARG}
		args+=("-master ${OPTARG}")
		;;
	r)
		args+=("-run ${OPTARG}")
		;;
	x)
		args+=("-px ${OPTARG}")
		;;
	y)
		args+=("-py ${OPTARG}")
		;;
	h)
		usage
		exit
		;;
	?)
		echo "[ERROR] Invalid option"
		exit 1
	esac
done

if [ -z "${nprocs}" ]; then nprocs=4; fi
if [ -z "${ppn}" ]; then ppn=1; nnodes=${nprocs}; fi
if [ -z "${master}" ]; then master=2; fi

total_nprocs=$nprocs
nnodes=$((nprocs/ppn))

if [ "${target}" = "ofp" ]; then
	rsc_args="rscgrp=debug-flat"
	#rsc_args="rscgrp=regular-cache"
elif [ "${target}" = "fugaku" ]; then
	total_nprocs=$((nprocs*2))
	if [ ${mck} -eq 0 ]; then
		if [ ${nnodes} -gt 385 ]; then
			rsc_args="rscgrp=eap-large"
		else
			rsc_args="rscgrp=eap-small"
		fi
	else
		# mckernel resource group
		rsc_args="rscunit=rscunit_ft01,rscgrp=dvsys-mck1,jobenv=linux2"
	fi
else
	echo "[ERROR] Unsupported Machine"
	exit 1
fi

batch_script="${script_dir}/batch.${target}.${nprocs}"
elapse_time="00:30:00"
jobname="d-${nprocs}"

cat <<- EOF > ${batch_script}
#!/bin/bash
#
#PJM -N "${jobname}"
#PJM -L "node=${nnodes}"
#PJM -L "${rsc_args}"
#PJM -L "elapse=${elapse_time}"
#PJM -g ${group[-1]}
#PJM -S
#PJM --spath ${log_dir}/%n.%j.stat
#PJM -o ${log_dir}/%n.%j.out
#PJM -e ${log_dir}/%n.%j.err
#PJM --mpi "proc=${total_nprocs}"
#	PJM --mpi "max-proc-per-node=${ppn}"

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

runlog_dir="${log_dir}/run-${nprocs}-${master}"
if [ ! -d "runlog_dir" ]; then
	mkdir -p ${runlog_dir}
fi

find ${log_dir} -maxdepth 1 -name "${jobname}.${jobid}.*" -exec mv {} ${runlog_dir} \;

