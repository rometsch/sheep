#!/usr/bin/env python3
# Replace a parameter value in a text based config file.
# The only assumption is that the parameter name and the value
# are separated by whitespaces and followed by whitespaces or newlines.
import re
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("file", help="config file to change parameter in")
parser.add_argument("param", help="parameter to be changed")
parser.add_argument("value", help="new value")
parser.add_argument("-o", help="output file path")
args = parser.parse_args()

infile_path = args.file
outfile_path = args.o if args.o is not None else args.file
param_name = args.param
new_param_value = args.value

ptrn = re.compile(r"({}[ \t]+)([\S]+(?=\s))".format(param_name))

with open(infile_path, 'r') as infile:
    outstr, Nrepl = re.subn(ptrn, r"\g<1>{}".format(new_param_value), infile.read())

if Nrepl > 1:
    raise AssertionError("Replaced {} instead of 1 occurances of parameter '{}'".format(Nrepl, param_name))

if Nrepl == 0:
    raise AssertionError("Could not find parameter '{}' in '{}'".format(param_name, infile_path))

with open(outfile_path, 'w') as outfile:
    outfile.write(outstr)
