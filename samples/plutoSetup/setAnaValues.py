#!/usr/bin/env python3
#------------------------------------------------------------
# brief  :	Handle pluto analysis definitions
# author :  Thomas Rometsch
# date   :	17-11-16
#------------------------------------------------------------
# Parse the definitions and names of analysis values from
# easy to use definitions in this file to the definition
# files of the pluto code.
# The convention that analysis variables start with
# *AN_* is used here. Though this is applied automatically
# and variables should be defined without it.
# This convention is used to distinguish definition of variable
# name (which are regenerated) from other definitions in the
# Analysis section of the header file which are copied.
#------------------------------------------------------------
import os;
from collections import OrderedDict;

# Definition of file names.
pluto_def = "definitions.h";
#plutocuda_def = "definitionsCuda.cuh";

# Definition of scalar variables, all are supposed to be of type double
vars = [
	 { "name" : "M_DISK",	"op" : "OP_SUM",	"len" : 1 }
	,{ "name" : "KE_R", 	"op" : "OP_SUM",	"len" : 1 }
	,{ "name" : "KE_TH", 	"op" : "OP_SUM",	"len" : 1 }
	,{ "name" : "KE_PHI", 	"op" : "OP_SUM",	"len" : 1 }
#	,{ "name" : "N_CELLS_RHO_MIN", "op" : "OP_SUM",	"len" : 1 }
	,{ "name" : "RHO_MIN", 	"op" : "OP_MIN",	"len" : 1 }
	,{ "name" : "RHO_MAX", 	"op" : "OP_MAX",	"len" : 1 }
	,{ "name" : "J_DISK", 	"op" : "OP_SUM",	"len" : 3 }
	,{ "name" : "F",		"op" : "OP_SUM",	"len" : 3 }
];

N_ana_vals = 0;
for var in vars:
	N_ana_vals += var["len"];

# Load the definition file and find the block for user definitions.
lines = [];
with open(pluto_def, "r") as f:
	lines = f.read().splitlines();
Nbeg = -1;
Nend = -1;
for n,line in enumerate(lines):
	if "[usrBegAna]" in line:
		Nbeg = n;
	if "[usrEndAna]" in line:
		Nend = n;
		break;
if Nbeg >= Nend:
	raise AssertionError("End of user block (line {}) found before begin of user block (line {})".format(Nend,Nbeg));

# Generate a string with all variable names 
# which is to be printed in the header of the output file
# and a string with all operators.
header_var_line = "\"";
ops_var_line = "";
is_first_var = True;
for var in vars:
	if var["len"] == 1:	
		header_var_line += "{}\\t".format(var["name"]);
		ops_var_line += "{},".format(var["op"]);
	else:
		for n in range(var["len"]):
			header_var_line += "{}_{}\\t".format(var["name"],n);
			ops_var_line += "{},".format(var["op"]);
header_var_line = header_var_line.rstrip("\\t") + "\"";
ops_var_line = ops_var_line.rstrip(",");


# Flags for num of ana vals and header
is_numvals_set = False;
is_headerstr_set = False;
is_opsstr_set = False;

# Generate new definition file
newlines = lines[:Nbeg+1];
# Copy all other definitions and set number of values
for line in lines[Nbeg+1:Nend]:
	if "NUMBER_OF_ANALYSIS_VALUES" in line:
		newlines += [ "#define  NUMBER_OF_ANALYSIS_VALUES   {}".format(N_ana_vals) ];
		is_numvals_set = True;	
	elif "ANA_VAL_HEADER" in line:
		newlines += [ "#define  ANA_VAL_HEADER   {:s}".format(header_var_line)];
		is_headerstr_set = True;
	elif "ANA_VAL_OPERATORS" in line:
		newlines += [ "#define  ANA_VAL_OPERATORS   {:s}".format(ops_var_line)];
		is_opsstr_set = True;
	else:
		try:
			name = line.split()[1];
			if name[:3] != "AN_":
				newlines += [line];
		except IndexError:
			continue;

# Write numvals and or headerstr if not done yet.
if not is_numvals_set:
	newlines += [ "#define  NUMBER_OF_ANALYSIS_VALUES   {}".format(N_ana_vals) ];
if not is_headerstr_set:
	newlines += [ "#define  AUTO_ANA_VAL_HEADER   {:s}".format(header_var_line)];
if not is_opsstr_set:
	newlines += [ "#define  ANA_VAL_OPERATORS   {:s}".format(ops_var_line)];

# Generate the names and indices of the variables.
IDX_var = 0;
for var in vars:
	newlines += [ "#define  AN_{:20s}  {:d}".format(var["name"], IDX_var) ];
	if var["len"] > 1:
		len_name = "{}_LEN".format(var["name"]);
		newlines += [ "#define  AN_{:20s}  {:d}".format(len_name, var["len"] ) ];
	IDX_var += var["len"];

# Append the rest of the original def file
newlines += lines[Nend:];

with open(pluto_def, "w") as of:
	for line in newlines:
		of.write(line+"\n")

