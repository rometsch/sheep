#!/usr/bin/env bash
# sheep v1.2

filename="$1"
name="${filename%.*}"

rm -rf $name
mkdir $name

tar -xzvf $1 -C $name

cd $name
./run
cd ..

mv $filename $name
