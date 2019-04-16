#!/usr/bin/env bash

# Obtain the directory of the scripts
# relative to the workspace
SD="$(basename $(dirname $(realpath $0)))"

# Preprate the simulation
$SD/prep.sh

# Queue a job to run the simulation
$SD/queue.sh $SD/start.sh
