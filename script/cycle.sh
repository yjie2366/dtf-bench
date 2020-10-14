script_dir=$(cd $(dirname ${BASH_SOURCE[0]}) &> /dev/null && pwd)
target="ofp"
args=()

while getopts "n:i:j:c:m:o:" OPT; do
	case ${OPT} in
	n)
		nprocs=$(( ${OPTARG} * 2 ))
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
	echo "[ERROR] Number of processes should be specified"
	exit 1
fi

if [ "${target}" = "ofp" ]; then
	rsc="regular-cache"
elif [ "${target}" = "fugaku" ]; then
	rsc="eap-large"
else
	echo "[ERROR] Unsupported Machine"
	exit 1
fi

batch_script="${script_dir}/batch.${target}"

cat <<- EOF > ${batch_script}
#!/bin/bash
#
#PJM -N "dtf-bench"
#PJM -L "node=${nprocs}"
#PJM -L "rscgrp=${rsc}"

sh ${script_dir}/exec.sh ${args[@]}
EOF

pjsub ${batch_script}
