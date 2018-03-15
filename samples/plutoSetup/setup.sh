#!/usr/bin/env bash
# This script runs the pluto setup.py script 
# to auto-generate the makefile and then
# runs make.

# Load modules on binac.
if [ $(hostname) == "login01" ]
then
	eval `/usr/share/Modules/bin/modulecmd-wrapper bash load compiler/gnu/5.2`
	eval `/usr/share/Modules/bin/modulecmd-wrapper bash load devel/cuda/8.0`
fi
python2 $PLUTO_DIR/setup.py --auto-update
echo "MAKING" > ../status.txt
make &> ../log/make.log
make clean >> ../log/make.log
