#!/bin/bash

#SBATCH -A skamburu
#SBATCH -N 4
#SBATCH --tasks-per-node=12
#SBATCH --time=00:30:00 

cp=/N/u/skamburu/.m2/repository/com/google/guava/guava/15.0/guava-15.0.jar:/N/u/skamburu/.m2/repository/commons-cli/commons-cli/1.2/commons-cli-1.2.jar:/N/u/skamburu/.m2/repository/habanero-java-lib/habanero-java-lib/0.1.1/habanero-java-lib-0.1.1.jar:/N/u/skamburu/.m2/repository/ompi/ompijavabinding/1.8.1/ompijavabinding-1.8.1.jar:/N/u/skamburu/.m2/repository/edu/indiana/soic/spidal/dapwc/1.0-ompi1.8.1/dapwc-1.0-ompi1.8.1.jar

x='x'
#opts="-XX:+UseConcMarkSweepGC -XX:ParallelCMSThreads=4 -Xms2G -Xmx2G"
opts="-XX:+UseG1GC -Xms512m -Xmx512m"

tpn=1
echo $SLURM_JOB_NUM_NODES 

/N/u/skamburu/projects/software/openmpi-1.8.1/build/bin/mpirun --report-bindings --mca btl ^tcp java $opts -cp $cp edu.indiana.soic.spidal.dapwc.Program -c config.properties -n $SLURM_JOB_NUM_NODES -t $tpn 
echo "Finished $pat on `date`" >> status.txt

