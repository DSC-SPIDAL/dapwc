package edu.indiana.soic.spidal.dapwc;

import edu.indiana.soic.spidal.common.BinaryReader2D;
import edu.indiana.soic.spidal.mpi.MpiOps;
import mpi.MPI;
import mpi.MPIException;

import static edu.rice.hj.Module0.finalizeHabanero;
import static edu.rice.hj.Module0.initializeHabanero;

public class PWCParallelism {
    public static void SetupParallelism(String[] args) throws MPIException {
        initializeHabanero();
        MPI.Init(args);

        PWCUtility.worldProcsComm = MPI.COMM_WORLD;
        PWCUtility.worldProcsRank = PWCUtility.worldProcsComm.getRank();
        PWCUtility.worldProcsCount = PWCUtility.worldProcsComm.getSize();

        PWCUtility.worldProcsPerNode = PWCUtility.worldProcsCount / PWCUtility
				.nodeCount;

        if ((PWCUtility.worldProcsPerNode * PWCUtility.nodeCount) != PWCUtility
				.worldProcsCount) {
            PWCUtility.printAndThrowRuntimeException("Inconsistent MPI counts" +
					" Nodes " + PWCUtility.nodeCount + " Size " + PWCUtility
					.worldProcsCount);
        }

        PWCUtility.mpiOps = new MpiOps(PWCUtility.worldProcsComm);
        PWCUtility.ParallelPattern =
                "---------------------------------------------------------\nMachine:" +
                        MPI.getProcessorName() + " " + PWCUtility.threadCount +
                        "x" + PWCUtility.worldProcsPerNode + "x" +
                        PWCUtility.nodeCount;
        if (PWCUtility.worldProcsRank == 0) {
            // TODO - distance type - short
            PWCUtility.SALSAPrint(0, " Distance Data Type: " + (Program
					.config.getDataTypeSize() == 2 ? "short" : Program.config
					.getDataTypeSize()));
            PWCUtility.SALSAPrint(0, PWCUtility.ParallelPattern);
        }

    } // End SetupParallelism

    public static void TearDownParallelism() throws MPIException {
        // Finalize threads
        finalizeHabanero();

        // End MPI
        MPI.Finalize();
    } // End TearDownParallelism

    public static void SetParallelDecomposition() {
        //	First divide points among processes
        Range[] processRanges = RangePartitioner.Partition(PWCUtility
				.PointCount_Global, PWCUtility.worldProcsCount);
        Range processRange = processRanges[PWCUtility.worldProcsRank]; // The
		// answer for this process

        PWCUtility.PointStart_Process = processRange.getStartIndex();
        PWCUtility.PointCount_Process = processRange.getLength();
        PWCUtility.PointCount_Largest = Integer.MIN_VALUE;

        for (Range r : processRanges) {
            PWCUtility.PointCount_Largest = Math.max(r.getLength(),
					PWCUtility.PointCount_Largest);
        }

        // We need points per process for all processes as used by load
		// balancing algorithm and then to reorder distance matrix consistently
        PWCUtility.PointsperProcess = new int[PWCUtility.worldProcsCount];
        PWCUtility.PointsperThreadperProcess = new int[PWCUtility
				.worldProcsCount][];

        for (int i = 0; i < PWCUtility.worldProcsCount; i++) {
            PWCUtility.PointsperProcess[i] = processRanges[i].getLength();
            Range[] threadprocessRanges = RangePartitioner.Partition
					(processRanges[i], PWCUtility.threadCount);
            PWCUtility.PointsperThreadperProcess[i] = new int[PWCUtility
					.threadCount];
            for (int j = 0; j < PWCUtility.threadCount; j++) {
                PWCUtility.PointsperThreadperProcess[i][j] =
						threadprocessRanges[j].getLength();
            }
        }

        //	Now divide points among threads for this process
        Range[] threadRanges = RangePartitioner.Partition
				(processRanges[PWCUtility.worldProcsRank], PWCUtility
						.threadCount);
        PWCUtility.PointsperThread = new int[PWCUtility.threadCount];
        PWCUtility.StartPointperThread = new int[PWCUtility.threadCount];

        for (int j = 0; j < PWCUtility.threadCount; j++) {
            PWCUtility.PointsperThread[j] = threadRanges[j].getLength();
            PWCUtility.StartPointperThread[j] = threadRanges[j].getStartIndex();
        }
    }

    // read data from file to memory
    // set starting position and end number for data points assigned to each
	// thread.
    // These are in startindex and lenindex respectively
    // Return data stored as
    // Distance [ClusterCenter,ClusterIndex] with all
	// ClusterIndex<=ClusterCenter stored in order for each ClusterCenter.
    // Diagonal values are stored (as zero)
    // Total space m(m+1)/2 if lower triangular store (checkerboard =1)
    public static void ReadDataFromFile(String fname) {
        //	First divide points among processes
        Range[] processRanges = RangePartitioner.Partition(PWCUtility
				.PointCount_Global, PWCUtility.worldProcsCount);
        Range processRange = processRanges[PWCUtility.worldProcsRank];

        PWCUtility.PointCount_Process = processRange.getLength();
        PWCUtility.PointStart_Process = processRange.getStartIndex();
        PWCUtility.PointCount_Largest = Integer.MIN_VALUE;

        for (Range r : processRanges) {
            PWCUtility.PointCount_Largest = Math.max(r.getLength(),
					PWCUtility.PointCount_Largest);
        }

        //	Now divide points among threads
        Range[] threadRanges = RangePartitioner.Partition(processRange,
				PWCUtility.threadCount);
        PWCUtility.StartPointperThread = new int[PWCUtility.threadCount];
        PWCUtility.PointsperThread = new int[PWCUtility.threadCount];

        for (int i = 0; i < PWCUtility.threadCount; i++) {
            PWCUtility.StartPointperThread[i] = threadRanges[i].getStartIndex();
            PWCUtility.PointsperThread[i] = threadRanges[i].getLength();
        }

        // Note - read binary distance data from file
        final edu.indiana.soic.spidal.common.Range range =
                new edu.indiana.soic.spidal.common.Range(
                        PWCUtility.PointStart_Process,
                        (PWCUtility.PointStart_Process + PWCUtility
								.PointCount_Process
                                - 1));
        PWCUtility.PointDistances = BinaryReader2D.readRowRange(fname, range,
				PWCUtility.PointCount_Global, PWCUtility.endianness, true,
				null);

    } //  End routine controlling reading of data

    public static int OwnerforThisPoint(int GlobalPointIndex) { // Return
		// process number for GlobalPointIndex

        int startpoint = 0;
        for (int mpiloop = 0; mpiloop < PWCUtility.worldProcsCount; mpiloop++) {
            int endpoint = startpoint + PWCUtility.PointsperProcess[mpiloop];
            if (GlobalPointIndex < endpoint) {
                return mpiloop;
            }
            startpoint = endpoint;
        }
        PWCUtility.printAndThrowRuntimeException(
                " Illegal Point in No Process " + (new Integer(GlobalPointIndex)).toString());
        return -1;

    } // End InThisProcess(int GlobalPointIndex)

}