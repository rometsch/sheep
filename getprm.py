#!/usr/bin/env python3
# Print a parameter value in a text based config file.
# The only assumption is that the parameter name and the value
# are separated by whitespaces and followed by whitespaces or newlines.
import re
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("file", help="config file to change parameter in")
parser.add_argument("param", help="parameter to be changed")
args = parser.parse_args()

infile_path = args.file
param_name = args.param

ptrn = re.compile(r"\A([ \t]*{}[ \t]+)([\S]+(?=\s))".format(param_name))

with open(infile_path, 'r') as infile:
    for line in infile:
        m = re.search(ptrn, line)
        if m:
            print(m.group(2))
