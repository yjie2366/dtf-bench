#!/bin/sh

filename=
num_proc=

while getopts "f:p:" OPT; do
	case ${OPT} in
	f)
		filename=${OPTARG}
		;;
	p)
		num_proc=${OPTARG}
		;;
	*)
		echo "[ERROR] Invalid option"
		exit
		;;
	esac
done

if [ -z "$filename" ]; then filename="err_verbose"; fi

if [ ! -e "$filename" ]; then
	echo "[ERROR] File $filename is non-existent"
	exit
fi

if [ -z "$num_proc" ]; then
	num_proc=4
	echo "[WARNING]: Number of processes for each component is not given."
fi

log_dir="`pwd`/log_${filename}"
if [ ! -d "$log_dir" ]; then
	mkdir -p $log_dir
fi

for comp in "scale" "letkf"; do
	for n in `seq 0 $((num_proc-1))`; do
		grep "v_${comp} ${n} " ${filename} > ${log_dir}/${comp}-${n}.log
	done
done
