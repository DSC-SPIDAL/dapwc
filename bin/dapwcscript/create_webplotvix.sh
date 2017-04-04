#!/bin/bash
clusterFile=$1
pointsFile=$2
newFileName=$3
cp=$HOME/.m2/repository/edu/indiana/soic/spidal/dapwc/1.0-ompi1.8.1/dapwc-1.0-ompi1.8.1.jar
java -cp $cp edu.indiana.soic.spidal.tools.CreatePlotSimple $pointsFile $clusterFile $newFileName
#paste $2 $1 $1 | cut -f1-4,7,9 | sed 's/\t/ /g' > $3
