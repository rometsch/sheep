#!/usr/bin/env bash

filename="$1"
name="${filename%.*}"

mkdir $name

tar -xzvf $1 -C $name

cd $name
./make_sheep
./queue_sheep start_sheep
cd ..

mv $filename $name
