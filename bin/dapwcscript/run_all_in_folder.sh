#!/bin/bash
folder=$6
roundCount=$7
tpp=$1
ppn=$2
mem=$3
numNodes=$4
nodesFile=$5

parentFolder=pwc_round_$roundCount
mkdir $parentFolder
subFolder=$parentFolder"_"
for confName in $folder*.properties; do
	echo $confName
	clusterNum=$(echo $confName| cut -d'_' -f 3 | cut -d '.' -f 1)
	echo $clusterNum
	mkdir $parentFolder/$subFolder"cluster_"$clusterNum
	./run.dapwc_init.sh $tpp $ppn $parentFolder"/"$subFolder"cluster" $clusterNum $mem g $numNodes $nodesFile $confName
done
