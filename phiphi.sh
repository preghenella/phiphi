#! /usr/bin/env bash

strun=0
nruns=5000
nevents=100000
njobs=16

workdir=runAna
mkdir $workdir || exit 1

for run in $(seq $strun $(($strun + $nruns - 1))); do

    ### wait if there are too many jobs running
    while true; do
	bkgjobs=$(jobs | grep phiphi | wc -l | xargs)
	if [ $bkgjobs -lt $njobs ]; then
	    break
	fi
	echo "[---] sleep while waiting for a free job slot"
	sleep 60
    done
    
    runid=$(printf "%05d" $run)
    seed=$((123456789 + $run * 2))
    
    echo "[---] starting run: $runid (nevents=$nevents, seed=$seed)"
    ./phiphi -n $nevents -s $seed -o $workdir/phiphi.$runid.root &> /dev/null & 

done
echo "[---] waiting for all runs to finish"
wait
exit 0
