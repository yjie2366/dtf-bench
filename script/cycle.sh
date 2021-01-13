script_dir=$(cd $(dirname ${BASH_SOURCE[0]}) &> /dev/null && pwd)
log_dir=${script_dir}/../log
args=()
group=($(id -nG))
master=
target=
ppn=2
nnodes=

# Need to add another case for Fugaku
case `hostname` in
	*ofp*)
		target=ofp
		;;
	*)
		target=unknown
		;;
esac

if [ ! -e ${log_dir} ]; then
	mkdir -p ${log_dir} || exit 1
fi

while getopts "n:p:i:j:c:m:a:u:x:y:" OPT; do
	case ${OPT} in
	n)
		nnodes=${OPTARG}
		args+=("-node ${nnodes}")
		;;
	p)
		nprocs=${OPTARG}
		args+=("-np ${nprocs}")
		;;
	i)
		args+=("-imax ${OPTARG}")
		;;
	j)
		args+=("-jmax ${OPTARG}")
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
	u)
		args+=("-run ${OPTARG}")
		;;
	x)
		args+=("-px ${OPTARG}")
		;;
	y)
		args+=("-py ${OPTARG}")
		;;
	?)
		echo "[ERROR] Invalid Option"
		exit 1
	esac
done

if [ -z "${nprocs}" ]; then nprocs=4; fi
if [ -z "${nnodes}" ]; then nnodes=${nprocs}; fi
if [ -z "${master}" ]; then master=2; fi

if [ "${target}" = "ofp" ]; then
#	rsc="debug-cache"
	rsc="regular-flat"
elif [ "${target}" = "fugaku" ]; then
	rsc="eap-small"
	if [ $((nprocs*2)) -gt 385 ]; then
		rsc="eap-large"
	fi
else
	echo "[ERROR] Unsupported Machine"
	exit 1
fi

batch_script="${script_dir}/batch.${target}.${nprocs}"
elapse_time="00:30:00"

cat <<- EOF > ${batch_script}
#!/bin/bash
#
#PJM -N "d-${nprocs}"
#PJM -L "node=$((nnodes*2))"
#PJM -L "rscgrp=${rsc}"
#PJM -L "elapse=${elapse_time}"
#PJM -g ${group[-1]}
#PJM -S
#PJM --spath ${log_dir}/%n.%j.stat
#PJM -o ${log_dir}/%n.%j.out
#PJM -e ${log_dir}/%n.%j.err
#PJM --mpi "proc=$((nprocs*2))"

sh ${script_dir}/exec.sh ${args[@]}

EOF

ret=$(pjsub ${batch_script})
jobid=`echo $ret | grep -i "submitted" | awk '{print $6}'`

echo "JOB ${jobid} submitted."

pjwait ${jobid}

runlog_dir="${log_dir}/run-${nprocs}-${master}"
if [ ! -d "runlog_dir" ]; then
	mkdir -p ${runlog_dir}
fi

find ${log_dir} -maxdepth 1 -name "d-${nprocs}.${jobid}.*" -exec mv {} ${runlog_dir} \;

