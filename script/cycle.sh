script_dir=$(cd $(dirname ${BASH_SOURCE[0]}) &> /dev/null && pwd)
log_dir=${script_dir}/../log
target="fugaku"
args=()
group=($(id -nG))
master=

if [ ! -e ${log_dir} ]; then
	mkdir -p ${log_dir} || exit 1
fi

while getopts "n:i:j:c:m:o:a:u:" OPT; do
	case ${OPT} in
	n)
		nprocs=${OPTARG}
		args+=("-np ${OPTARG}")
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
	o)
		target=${OPTARG}
		target=${target,,*}
		if [ "${target}" != "ofp" -a "${target}" != "fugaku" ]; then
			echo "[ERROR] Invalid Target Machine"
			exit 1
		fi
		;;
	?)
		echo "[ERROR] Invalid Option"
		exit 1
	esac
done

if [ -z "${nprocs}" ]; then
	nprocs=4
fi

if [ -z "${master}" ]; then 
	master=2
fi

if [ "${target}" = "ofp" ]; then
	rsc="regular-cache"
elif [ "${target}" = "fugaku" ]; then
	rsc="eap-small"
	if [ $((nprocs*2)) -gt 385 ]; then
		rsc="eap-large"
	fi	
else
	echo "[ERROR] Unsupported Machine"
	exit 1
fi

batch_script="${script_dir}/batch.${target}"

cat <<- EOF > ${batch_script}
#!/bin/bash
#
#PJM -N "d-${nprocs}"
#PJM -L "node=$((nprocs * 2))"
#PJM -L "rscgrp=${rsc}"
#PJM -L "elapse=00:60:00"
#PJM -g ${group[-1]}
#PJM -S
#PJM --spath ${log_dir}/%n.%j.stat
#PJM -o ${log_dir}/%n.%j.out
#PJM -e ${log_dir}/%n.%j.err
#	PJM --mpi "proc=$((nprocs * 2))"
#PJM --mpi "max-proc-per-node=1"

sh ${script_dir}/exec.sh ${args[@]}

EOF

ret=$(pjsub ${batch_script})
jobid=`echo $ret | grep -i "submitted" | awk '{print $6}'`
runlog_dir="${log_dir}/run-${nprocs}-${master}"

pjwait ${jobid}
find ${log_dir} -maxdepth 1 -name "d-${nprocs}.${jobid}.*" -exec mv {} ${runlog_dir} \;

