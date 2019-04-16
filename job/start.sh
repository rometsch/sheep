#!/usr/bin/env bash

# BINAC: Change dir
if [[ ${PBS_O_WORKDIR} ]]
then
	cd $PBS_O_WORKDIR
fi

source job/env

CMD="./fargo3d setups/$SETUP/$SETUP.par"

if (( $NCPU > 1 )); then
	mpirun -n $NCPU $CMD
else
	$CMD
fi
