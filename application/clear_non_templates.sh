#!/bin/bash
#for di in ls -d output_max/!(\(template\)); do rm -rf $di; done
# dir=output_max_withoutMP_SC_RT_kapp_mapping
dir=output_max_withoutMP
cd $dir
for di in $(ls -d *); do
	if [[ "$di" != "(template)" && "$di" != "run.sh" ]]; then
		rm -rf $di
	fi
done