#!/bin/bash

nodes=24
name="$nodes"n
nodefile=nodes.txt

./run.generic.sh 1 12 $name samplerun 2 g $nodes $nodefile 
