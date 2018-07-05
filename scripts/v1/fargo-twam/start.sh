#!/usr/bin/env bash
# sheep version: 1

source env_sheep

# BINAC: Change dir
if [[ ${PBS_O_WORKDIR} ]]
then
	cd $PBS_O_WORKDIR
	# Load cuda when on a gpu node
	if [[ $(hostname -s) = gpu* ]]
	then
		module load devel/cuda/8.0
		unset CUDA_VISIBLE_DEVICES
	fi
fi

mpirun -n $NCPU fargo setups/$SETUP/in.par
