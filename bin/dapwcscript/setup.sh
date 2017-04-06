#!/bin/bash
#used to setp the initial scripts
distFile=$PWD/$1
numPoints=$2
numClusters=$3
confDAMDS=$PWD/conf/damds/config.properties
confDAPWC=$PWD/conf/dapwc/config.properties
edited=.edited

. init.properties
awk -F"=" -v newval="=$distFile"  '/DistanceMatrixFile/{$2=newval;print;next}1' $confDAMDS | awk -F"=" -v newval="=$numPoints"  '/numPoints/{$2=newval;print;next}1' > $confDAMDS$edited

awk -F"=" -v newval="=$distFile"  '/DistanceMatrixFile/{$2=newval;print;next}1' $confDAPWC | awk -F"=" -v newval="=$numPoints"  '/NumberDataPoints/{$2=newval;print;next}1' | awk -F"=" -v newval="=$numClusters"  '/MaxNcent/{$2=newval;print;next}1' > $confDAPWC$edited


echo "Done initial setp running DAPWC"
./init.sh $numNodes $nodesFile $tpp $ppn $mem $confDAPWC$edited $confDAMDS$edited $distFile
