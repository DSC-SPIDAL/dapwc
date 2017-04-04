#!/bin/bash

cp=$HOME/.m2/repository/com/google/guava/guava/15.0/guava-15.0.jar:$HOME/.m2/repository/commons-cli/commons-cli/1.2/commons-cli-1.2.jar:$HOME/.m2/repository/habanero-java-lib/habanero-java-lib/0.1.1/habanero-java-lib-0.1.1.jar:$HOME/.m2/repository/ompi/ompijavabinding/1.10.1/ompijavabinding-1.10.1.jar:$HOME/.m2/repository/edu/indiana/soic/spidal/dapwc/1.0-ompi1.8.1/dapwc-1.0-ompi1.8.1.jar:$HOME/.m2/repository/edu/indiana/soic/spidal/common/1.0/common-1.0.jar

wd=`pwd`
x='x'

cps=6
spn=4
cpn=$(($cps*$spn))

tpp=$1
ppn=$2
nodes=$7

memmultype=g
xmx=1
if [ "$5" ];then
    xmx=$5
    memmultype=$6
fi

bindToUnit=core
if [ $tpp -gt 1 ];then
    bindToUnit=none
fi

if [ $ppn -gt $cpn ];then
    bindToUnit=hwthread
fi

pat=$3_$4

TimingFile=$pat/timing.txt
SummaryFile=$pat/summary.txt
ClusterFile=$pat/cluster.txt

PWC_OPS="-DTimingFile=$TimingFile -DSummaryFile=$SummaryFile -DClusterFile=$ClusterFile"


opts="-XX:+UseG1GC -Xms256m -Xmx"$xmx"$memmultype"

#kill.java.sh $7
echo "Running $pat on `date`" >> $pat/status.txt

$BUILD/bin/mpirun --hostfile $8 --mca btl ^tcp --report-bindings --map-by ppr:$ppn:node:SPAN  --bind-to $bindToUnit --rank-by core  -np $(($nodes*$ppn)) java $opts  $PWC_OPS -cp $cp edu.indiana.soic.spidal.dapwc.Program -c $9 -n $nodes -t $tpp 2>&1 | tee $pat/out.txt
echo "Finished $pat on `date`" >> $pat/status.txt

