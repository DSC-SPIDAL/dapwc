#!/bin/bash
#edit the input properties file to the cluster extractor
distFile=$1
clusterFile=$2
clusters=$3
newclusters=$4
outDir=$5
#sed -i "s/\(distFile *= *\).*/\1$3/" $1 > $2
#awk -F"=" -v newval="=$3"  '/distFile/{$2=newval;print;next}1' $1 | awk -F"=" -v newval1="=$4" '/clusterFile/{$2=newval1;print;next}1' | awk -F"=" -v newval="=$5" '/numPoints/{$2=newval;print;next}1' | awk -F"=" -v newval="=$6" '/clusters/{$2=newval;print;next}1' | awk -F"=" -v newval="=$7" '/newclusters/{$2=newval;print;next}1'| awk -F"=" -v newval="=$8" '/joins/{$2=newval;print;next}1' | awk -F"=" -v newval="=$outdir" '/outDir/{$2=newval;print;next}1' > $2

cp=/N/u/pswickra/.m2/repository/edu/indiana/soic/spidal/dapwc/1.0-ompi1.8.1/dapwc-1.0-ompi1.8.1.jar
java -cp $cp edu.indiana.soic.spidal.tools.ClusterExtractorSimple $clusterFile $distFile $outDir $clusters $newclusters false
