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

def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("file", help="config file to change parameter in")
    parser.add_argument("planetnumber", type=int, help="which planet to edit, n = 1,2,3,...")
    parser.add_argument("param", choices=prm_list, help="parameter to be changed")
    parser.add_argument("value", help="new value")
    parser.add_argument("-o", help="output file path")
    args = parser.parse_args()


    change_planet_param(args.file, args.planetnumber,
                        args.param, args.value,
                        outfile_path = args.o)

def change_planet_param(infile_path, planet_number, param, value,
                        outfile_path=None):

    if param not in syntax:
        raise ValueError("param {} not supported, only {} are supported".format(param, ", ".join(prm_list)))
    
    outfile_path = outfile_path if outfile_path is not None else infile_path
                        
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
        ent = line.strip().split()
        if len(ent) != len(prm_list):
            raise ValueError("File syntax ({} columns) does not fit expected syntax ({} columns).".format(len(ent), len(prm_list))) 
        planets.append( ent )
        
    if planet_number > Nplanets or planet_number < 1:
        raise ValueError("Planet {} selected but only {} in planet config.".format(planet_number, Nplanets))

    
    # Change the parameter
    planets[planet_number-1][syntax[param]] = value

    # Reconstruct the parameter file
    with open(outfile_path, 'w') as outfile:
        for line in header_lines:
            outfile.write(line)
        for planet in planets:
            line = "\t".join(planet) + "\n"
            outfile.write(line)

if __name__=="__main__":
    main()
