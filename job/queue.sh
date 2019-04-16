#!/usr/bin/env bash

source job/env

ABSOLUTE_PATH="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/$(basename "${BASH_SOURCE[0]}")"
PROJECT_NAME="$(basename $ABSOLUTE_PATH)"
UUID="$(cat job/uuid.txt)"

echo $PROJECT_NAME
echo $ABSOLUTE_PATH

TOQUEUE=$1
shift 1
PASSEDARGS=$@

echo $PASSEDARGS

case $(hostname -s) in
	login0*)
		if (( $NGPU > 0 )); then
			RESOURCESTR="nodes=1:ppn=$NGPU:gpus=$NGPU:exclusive_process"
		else
			RESOURCESTR="nodes=1:ppn=$NCPU"
		fi
		# We are on binac
		qsub -N "${PROJECT_NAME}-ID-$UUID" \
			 -e $LOGDIR/err.log \
			 -o $LOGDIR/out.log \
			 -m bae \
			 -M "thomas.rometsch@uni-tuebingen.de" \
			 -l walltime=$HOURS:00:00 \
			 -l $RESOURCESTR \
			 -q $QUEUE \
			 $TOQUEUE $PASSEDARGS
		;;
  	cpt-*)
		if (( $NGPU > 0 )); then
			bsub -J $PROJECT_NAME -o $LOGDIR/out.log -e $LOGDIR/err.log -q gpu$DEVICE $PWD/$TOQUEUE $PASSEDARGS
		else
			$PWD/$TOQUEUE $PASSEDARGS 1>$LOGDIR/out.log 2>$LOGDIR/err.log &
		fi
  		;;
	*)
		$TOQUEUE $PASSEDARGS
		;;
esac
