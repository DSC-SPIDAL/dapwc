package edu.indiana.soic.spidal.dapwc;

import edu.indiana.soic.spidal.mpi.MpiOps;
import mpi.MPI;
import mpi.MPIException;

public class PWCParallelism
{
	public static void SetupParallelism(String[] args, String parentAffinityString) throws MPIException
	{
		//  Set up MPI
        MPI.Init(args);
		PWCUtility.MPI_communicator = MPI.COMM_WORLD; //initializing MPI world communicator

		PWCUtility.MPI_Rank = PWCUtility.MPI_communicator.getRank(); // Rank of this process
		PWCUtility.MPI_Size = PWCUtility.MPI_communicator.getSize(); // Number of MPI Processes
        PWCUtility.mpiOps = new MpiOps(PWCUtility.MPI_communicator);

		// Set up MPI
		PWCUtility.MPIperNodeCount = PWCUtility.MPI_Size / PWCUtility.NodeCount;

		if ((PWCUtility.MPIperNodeCount * PWCUtility.NodeCount) != PWCUtility.MPI_Size)
		{
            PWCUtility.printAndThrowRuntimeException("Inconsistent MPI counts Nodes " + PWCUtility.NodeCount + " Size " + PWCUtility.MPI_Size);
		}

        if (PWCUtility.bindThreads){
            int availableComputeUnits = getAvailableComputeUnits(parentAffinityString);
            if (availableComputeUnits < PWCUtility.ThreadCount){
                PWCUtility.printAndThrowRuntimeException("The process is bind to compute units [" + getHumanReadableString(parentAffinityString) + "], which is less than the requested number of threads (" + PWCUtility.ThreadCount + ")");
            }
            PWCUtility.bindThreadToCore = getBindToCoreIds(parentAffinityString, PWCUtility.ThreadCount);

        }

		PWCUtility.ParallelPattern =
				"---------------------------------------------------------\nMachine:" + MPI.getProcessorName() + " " +
						PWCUtility.ThreadCount + "x" + PWCUtility.MPIperNodeCount + "x" + PWCUtility.NodeCount;
		if (PWCUtility.MPI_Rank == 0)
		{
            // TODO - distance type - short
			PWCUtility.SALSAPrint(0, " Distance Data Type: " + (Program.config.getDataTypeSize() == 2 ? "short" : Program.config.getDataTypeSize()));
			PWCUtility.SALSAPrint(0, PWCUtility.ParallelPattern);
		}

	} // End SetupParallelism

    private static int[] getBindToCoreIds(String affinityString, int threadCount) {
        int [] bindings = new int[threadCount];
        int idx = 0;
        for (int i = affinityString.length() - 1; i >= 0; --i) {
            if (affinityString.charAt(i) == '1') {
                bindings[idx] = ((affinityString.length() - 1) - i);
                ++idx;
            }
        }
        return bindings;
    }

    public static String getHumanReadableString(String affinityMask) {
        String humanReadable = "";
        for (int i = 0; i < affinityMask.length(); i++) {
            if (affinityMask.charAt(i) == '1'){
                humanReadable += ((affinityMask.length() - 1) - i) + " ";
            }

        }
        return humanReadable;
    }

    private static int getAvailableComputeUnits(String parentAffinityString) {
        int count = 0;
        for (char c : parentAffinityString.toCharArray()){
            if (c == '1') ++count;
        }
        return count;
    }

    public static void TearDownParallelism() throws MPIException {
        // End MPI
        MPI.Finalize();
	} // End TearDownParallelism

	public static void SetParallelDecomposition()
	{
		//	First divide points among processes
		Range[] processRanges = RangePartitioner.Partition(PWCUtility.PointCount_Global, PWCUtility.MPI_Size);
		Range processRange = processRanges[PWCUtility.MPI_Rank]; // The answer for this process

		PWCUtility.PointStart_Process = processRange.getStartIndex();
		PWCUtility.PointCount_Process = processRange.getLength();
		PWCUtility.PointCount_Largest = Integer.MIN_VALUE;

		for (Range r : processRanges)
		{
			PWCUtility.PointCount_Largest = Math.max(r.getLength(), PWCUtility.PointCount_Largest);
		}

		// We need points per process for all processes as used by load balancing algorithm and then to reorder distance matrix consistently
		PWCUtility.PointsperProcess = new int[PWCUtility.MPI_Size];
		PWCUtility.PointsperThreadperProcess = new int[PWCUtility.MPI_Size][];

		for (int i = 0; i < PWCUtility.MPI_Size; i++)
		{
			PWCUtility.PointsperProcess[i] = processRanges[i].getLength();
			Range[] threadprocessRanges = RangePartitioner.Partition(processRanges[i], PWCUtility.ThreadCount);
			PWCUtility.PointsperThreadperProcess[i] = new int[PWCUtility.ThreadCount];
			for (int j = 0; j < PWCUtility.ThreadCount; j++)
			{
				PWCUtility.PointsperThreadperProcess[i][j] = threadprocessRanges[j].getLength();
			}
		}

		//	Now divide points among threads for this process
		Range[] threadRanges = RangePartitioner.Partition(processRanges[PWCUtility.MPI_Rank], PWCUtility.ThreadCount);
		PWCUtility.PointsperThread = new int[PWCUtility.ThreadCount];
		PWCUtility.StartPointperThread = new int[PWCUtility.ThreadCount];

		for (int j = 0; j < PWCUtility.ThreadCount; j++)
		{
			PWCUtility.PointsperThread[j] = threadRanges[j].getLength();
			PWCUtility.StartPointperThread[j] = threadRanges[j].getStartIndex();
		}
	}

	// read data from file to memory
	// set starting position and end number for data points assigned to each thread. 
	// These are in startindex and lenindex respectively
	// Return data stored as
	// Distance [ClusterCenter,ClusterIndex] with all ClusterIndex<=ClusterCenter stored in order for each ClusterCenter. 
	// Diagonal values are stored (as zero)
	// Total space m(m+1)/2 if lower triangular store (checkerboard =1)
	public static void ReadDataFromFile(String fname)
	{
		//	First divide points among processes
		Range[] processRanges = RangePartitioner.Partition(PWCUtility.PointCount_Global, PWCUtility.MPI_Size);
		Range processRange = processRanges[PWCUtility.MPI_Rank];

		PWCUtility.PointCount_Process = processRange.getLength();
		PWCUtility.PointStart_Process = processRange.getStartIndex();
		PWCUtility.PointCount_Largest = Integer.MIN_VALUE;

		for (Range r : processRanges)
		{
			PWCUtility.PointCount_Largest = Math.max(r.getLength(), PWCUtility.PointCount_Largest);
		}

		//	Now divide points among threads
		Range[] threadRanges = RangePartitioner.Partition(processRange, PWCUtility.ThreadCount);
		PWCUtility.StartPointperThread = new int[PWCUtility.ThreadCount];
		PWCUtility.PointsperThread = new int[PWCUtility.ThreadCount];

		for (int i = 0; i < PWCUtility.ThreadCount; i++)
		{
			PWCUtility.StartPointperThread[i] = threadRanges[i].getStartIndex();
			PWCUtility.PointsperThread[i] = threadRanges[i].getLength();
		}

		// Note - read binary distance data from file
        PWCUtility.PointDistances = Matrix.readRowRange(fname, processRange, PWCUtility.PointCount_Global,PWCUtility.dataTypeSize, PWCUtility.endianness, PWCUtility.isMemoryMapped);

	} //  End routine controlling reading of data

	public static int OwnerforThisPoint(int GlobalPointIndex)
	{ // Return process number for GlobalPointIndex

		int startpoint = 0;
		for (int mpiloop = 0; mpiloop < PWCUtility.MPI_Size; mpiloop++)
		{
			int endpoint = startpoint + PWCUtility.PointsperProcess[mpiloop];
			if (GlobalPointIndex < endpoint)
			{
				return mpiloop;
			}
			startpoint = endpoint;
		}
		PWCUtility.printAndThrowRuntimeException(
                " Illegal Point in No Process " + (new Integer(GlobalPointIndex)).toString());
		return -1;

	} // End InThisProcess(int GlobalPointIndex)

}