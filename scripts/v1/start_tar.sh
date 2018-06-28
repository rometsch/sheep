#!/usr/bin/env bash

filename="$1"
name="${filename%.*}"

mkdir $name

tar -xzvf $1 -C $name

cd $name
./sheep_setup
./sheep_queue sheep_start
cd ..

mv $filename $name
