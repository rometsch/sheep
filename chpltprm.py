#!/usr/bin/env python3
# Replace a parameter value in a Fargo3d planet config file.
import re
import argparse

# Parameters allowed in the planet config.
prm_list = ["name",
            "distance",
            "mass",
            "accretion",
            "feelsdisk",
            "feelsothers"]
syntax = { key : n for n,key in enumerate(prm_list) } 

parser = argparse.ArgumentParser()
parser.add_argument("file", help="config file to change parameter in")
parser.add_argument("planetnumber", type=int, help="which planet to edit, n = 1,2,3,...")
parser.add_argument("param", choices=prm_list, help="parameter to be changed")
parser.add_argument("value", help="new value")
parser.add_argument("-o", help="output file path")
args = parser.parse_args()

infile_path = args.file
outfile_path = args.o if args.o is not None else args.file
param_name = args.param
new_param_value = args.value

# parse the config and separate header and planet config
header_lines = []
planet_lines = []
with open(infile_path, 'r') as infile:
    for line in infile:
        if len(line.strip()) == 0 or line.strip()[0] == "#":
            header_lines.append(line)
        else:
            planet_lines.append(line)

# Parse the planet lines
Nplanets = len(planet_lines)
planets = []
for line in planet_lines:
    planets.append( line.strip().split() )

    
if args.planetnumber > Nplanets or args.planetnumber < 1:
    raise ValueError("Planet {} selected but only {} in planet config.".format(args.planetnumber, Nplanets))

# Change the parameter
planets[args.planetnumber-1][syntax[args.param]] = args.value

# Reconstruct the parameter file
with open(outfile_path, 'w') as outfile:
    for line in header_lines:
        outfile.write(line)
    for planet in planets:
        line = "\t".join(planet) + "\n"
        outfile.write(line)
