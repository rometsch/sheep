#!/usr/bin/env bash

source job/env

# apply all patches in the patches directory
for PATCH in $(find patches -mindepth 1 -maxdepth 1 -type d)
do
	$PATCH/apply.sh | tee $LOGDIR/patches.log
done

# make fargo3d with the right compile flags
FLAGS="UNITS=CGS RESCALE=0"

# build parallel or not
if (( $NCPU > 1 )); then
	FLAGS="$FLAGS PARALLEL=1"
else
	FLAGS="$FLAGS PARALLEL=0"
fi

# use GPU or not
if (( $NGPU > 0 )); then
	FLAGS="$FLAGS GPU=1"
else
	FLAGS="$FLAGS GPU=0"
fi

make SETUP=$SETUP $FLAGS | tee $LOGDIR/make.log

