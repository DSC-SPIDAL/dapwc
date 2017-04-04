#!/bin/bash
#edit the input properties file to the cluster extractor
pointsFile=$1
clusterDirPatt=$2
clusterFilePatt=$3
outDir=$4
numClusters=$5
pattern=$6

cp=$HOME/.m2/repository/edu/indiana/soic/spidal/dapwc/1.0-ompi1.8.1/dapwc-1.0-ompi1.8.1.jar
java -cp $cp edu.indiana.soic.spidal.tools.CollateClustersSimple $pointsFile $clusterDirPatt $clusterFilePatt $outDir $numClusters $pattern
