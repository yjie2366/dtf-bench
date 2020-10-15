script_dir=$(cd $(dirname ${BASH_SOURCE[0]}) &> /dev/null && pwd)
target="ofp"
args=()

while getopts "n:i:j:c:m:o:" OPT; do
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

if [ -z ${nprocs} ]; then
	nprocs=4
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
#PJM -N "dtf-bench"
#PJM -L "node=$((nprocs * 2))"
#PJM -L "rscgrp=${rsc}"
#PJM -L "elapse=00:60:00"
#PJM -g `id -nG | cut -d" " -f2`
#PJM -o ${script_dir}/../log/%n.%j.out
#PJM -e ${script_dir}/../log/%n.%j.err
#PJM --mpi "max-proc-per-node=1"

sh ${script_dir}/exec.sh ${args[@]}
EOF

pjsub ${batch_script}
