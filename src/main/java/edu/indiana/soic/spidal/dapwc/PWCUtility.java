package edu.indiana.soic.spidal.dapwc;

import com.google.common.base.Stopwatch;
import com.google.common.base.Strings;
import mpi.Intracomm;
import mpi.MPIException;
import edu.indiana.soic.spidal.mpi.MpiOps;

import java.io.IOException;
import java.io.PrintWriter;
import java.nio.ByteOrder;
import java.nio.charset.Charset;
import java.nio.file.*;
import java.util.ArrayList;
import java.util.concurrent.TimeUnit;

public class PWCUtility
{
    public static final double INV_SHORT_MAX = 1.0/Short.MAX_VALUE;
	public static int PointCount_Global = 0; // Total number of points summed over all threads and processes
	public static int PointCount_Process = 0; // Total number of points summed over all threads in this process
	public static int PointCount_Largest = 0; // Largest number of points in all processes
	public static int PointStart_Process = 0; //    First data point in this process

	// Parallel Parameters
	public static int MPI_Rank = 0; // Rank of process
	public static int MPI_Size = 1; // Number of MPI Processes
	public static Intracomm MPI_communicator = null; //MPI communicator

	public static String ParallelPattern = " "; // Pattern of parallel execution (e.g. 8x1x8 indicates [threads/process][processes/node][nodes])
	public static String PatternLabel = " "; // Title line for print

	//  Within a job data points will be divided into MPI_Size parts -- each part is assigned to a separate MPI Process
	public static int[] PointsperProcess = null; //how many data points each process will take care
	public static int[][] PointsperThreadperProcess = null; // Number of data points in each process-thread

	//	Within a process, data points will be divided into ThreadCount segments( each thread a segment), the array keep the size of each segment
	public static int[] PointsperThread = null; //how many data points each thread will take care
	public static int[] StartPointperThread = null; //the starting point that a thread will take care

	public static int ThreadCount = 1; // maximum number of parallel threads in a process
	public static int NodeCount = 1; // maximum number of separate nodes in run
	public static int MPIperNodeCount = 1; // Number of MPI processes per node

	public static int MPIIOStrategy = 0; // MPI I/O Strategy

	// TODO - conversion fix - make an appropriate data structure to hold distance values
	/*public static Matrix<TDistance> PointDistances;*/

	// Timing Parameters
    public static Stopwatch mainTimer;
	public static Stopwatch PreciseTimer; //    Hold Precise Timing
	public static int NumberofSubTimings = 0; // Number of subtimings
	public static Stopwatch[] SubTimers; // Timing Objects
    public static long mainDuration = 0; // duration measured by mainTimer in milliseconds
	public static double HPDuration = 0.0; // Time with Precision
	public static double[] SubDurations; // Hold partial timing
	public static String[] SubTimingNames; //  Labels of partial timing
	public static boolean[] SubTimingEnable;
	public static int MPIREDUCETiming = -1;
	public static int MPIREDUCETiming1 = -1;
	public static int MPISENDRECEIVETiming = -1;
	public static int MPISENDRECEIVEEigenTiming = -1;
	public static int MPIBROADCASTTiming = -1;
	public static int ThreadTiming = -1;

	//  These are general parameters for C# codes
	public static ArrayList<String> CosmicOutput = new ArrayList<>(1000); // Monitoring Output
	public static boolean ConsoleDebugOutput = true; // If true send Monitoring output to console
	public static int DebugPrintOption = 2; // Control Printing (= 0 None, ==1 Summary, = 2 Full)

	//  Center Finding Parameters
	public static int NumberofCenters = 8; // Number of centers to be found with each method
	public static double[] BucketFractions = new double[] {0.15, 0.4, 0.75}; // Fractions to be used in Bucket method -- centers are those with most neighbors in a radius determined by BucketFraction
	public static int NumberofBuckets = BucketFractions.length; // Number of Buckets
	public static int addMDS = 1; // Specify MDS versions of center finding; = 0 ignore Currently all nonzero values treated in same way
	public static String addMDSfile = ""; // File with MDS Information
	public static String Labelfile = ""; // File with Label and Length Information
	public static String ClusterNumberfile = ""; // File with Cluster Numbers; can be same as addMDSFile
	public static int CenterPointsPerCenterTypeInOutput = 3; // number of center points to include for each center type in the output plot
	public static String CenterPlotFile = ""; // output plot file with centers
	public static double MinimumDistanceCut = -0.1; // Insist Sequence distances bigger than this; negative values remove test
	public static int LinkCountinCenterFinding = 0; // Minimum Number of links
	public static int LengthCut1 = -1; // Insists first sequence has length bigger thanm this
	public static int LengthCut2 = -1; // Insist comparison sequence has length bigger than this
    public static MpiOps mpiOps;
    public static int dataTypeSize;
    public static ByteOrder endianness;
    public static boolean isMemoryMapped;
    public static int repetitions;
    public static short[] PointDistances; // 1st dimension is rowIdx, 2nd is colIdx
    public static boolean timingCompleted = false;

    public static void printException(Exception e) {
        System.out.println("SALSA Error " + e.getMessage());
    }

    public static void printAndThrowRuntimeException(RuntimeException e) {
        System.out.println("SALSA Error " + e.getMessage());
        throw e;
    }

    public static void printAndThrowRuntimeException(String message) {
        System.out.println("SALSA Error " + message);
        throw new RuntimeException(message);
    } // end printAndThrowRuntimeException

    // PrintOption = 0 Essential Printout
    // PrintOption = 1 Summary Printout
    // PrintOption = 2 Only if full print out requested
    public static void SALSAPrint(int PrintOption, String StufftoPrint) {
        if (MPI_Rank != 0) {
            return;
        }
        if (DebugPrintOption < PrintOption) {
            return;
        }
        CosmicOutput.add(StufftoPrint);

        if (ConsoleDebugOutput) {
            System.out.println(StufftoPrint);
        }
    } // End SALSAPrint

    public static String PrintFixedInteger(int data, int digits) {
        String returned = (new Integer(data)).toString();
        for (int charloop = returned.length(); charloop < digits; charloop++) {
            returned = " " + returned;
        }
        return returned;

    } // End PrintFixedInteger

    public static String LeftPrintFixedInteger(int data, int digits) {
        String returned = (new Integer(data)).toString();
        for (int charloop = returned.length(); charloop < digits; charloop++) {
            returned = returned + " ";
        }
        return returned;

    } // End LeftPrintFixedInteger

    public static void InitializeTiming(int InputNumberofTimers) {
        NumberofSubTimings = InputNumberofTimers;
        SubTimers = new Stopwatch[NumberofSubTimings]; // Timing Objects
        SubDurations = new double[NumberofSubTimings]; // Hold partial timing
        SubTimingEnable = new boolean[NumberofSubTimings];
        SubTimingNames = new String[NumberofSubTimings];

        for (int itimer = 0; itimer < NumberofSubTimings; itimer++) {
            SubTimers[itimer] = Stopwatch.createUnstarted();
            SubDurations[itimer] = 0.0;
            SubTimingEnable[itimer] = true;
        }
        PreciseTimer = Stopwatch.createStarted();
        mainTimer = Stopwatch.createStarted();

    } // End InitializeTiming

    public static void SetUpSubTimer(int TimingIndex, String TimingLabel) {
        if (TimingIndex >= NumberofSubTimings) {
            printAndThrowRuntimeException("Error in Timing Index " + (new Integer(TimingIndex)).toString() + " Max " +
                                                  (new Integer(NumberofSubTimings)).toString());
            return;
        }
        SubTimingNames[TimingIndex] = TimingLabel;
        SubTimingEnable[TimingIndex] = true;
        SubDurations[TimingIndex] = 0.0;
    } // End SetUpSubTimer

    public static void SetUpMPISubTimers(int StartTimingIndex, String MPISetLabel) {
        if ((StartTimingIndex + 5) >= NumberofSubTimings) {
            printAndThrowRuntimeException(
                    "Error in  MPI Timing Index " + (new Integer(StartTimingIndex + 2)).toString() + " Max " +
                            (new Integer(NumberofSubTimings)).toString());
            return;
        }
        SetUpSubTimer(StartTimingIndex, MPISetLabel + "MPI Reduce");
        SetUpSubTimer(StartTimingIndex + 1, MPISetLabel + "MPI SendRec");
        SetUpSubTimer(StartTimingIndex + 2, MPISetLabel + "MPI SendRec Eig");
        SetUpSubTimer(StartTimingIndex + 3, MPISetLabel + "MPI Bcast");
        SetUpSubTimer(StartTimingIndex + 4, MPISetLabel + "MPI Global Reductions");
        SetUpSubTimer(StartTimingIndex + 5, MPISetLabel + "Thread Global Reductions");

        if ((!MPISetLabel.equals("")) && (!MPISetLabel.equals("Lib "))) {
            return;
        }
        MPIREDUCETiming = StartTimingIndex;
        MPISENDRECEIVETiming = StartTimingIndex + 1;
        MPISENDRECEIVEEigenTiming = StartTimingIndex + 2;
        MPIBROADCASTTiming = StartTimingIndex + 3;
        MPIREDUCETiming1 = StartTimingIndex + 4;
        ThreadTiming = StartTimingIndex + 5;

    } // End SetUpMPISubTimers

    public static void InterimTiming() {
        PreciseTimer.stop();
        HPDuration += PreciseTimer.elapsed(TimeUnit.MICROSECONDS);
        PreciseTimer.reset();
        PreciseTimer.start();

    } // end InterimTiming

    public static void EndTiming() {
        timingCompleted = true;
        mainTimer.stop();
        PreciseTimer.stop();
        HPDuration += PreciseTimer.elapsed(TimeUnit.MICROSECONDS);
        mainDuration = mainTimer.elapsed(TimeUnit.MILLISECONDS);
        mainTimer.reset();
        PreciseTimer.reset();

    } // end EndTiming

    public static void StartSubTimer(int TimingIndex) {
        if (TimingIndex < 0) {
            return;
        }

        if (SubTimingEnable[TimingIndex]) {
            SubTimers[TimingIndex].start();
        }

    } // End StartSubTimer

    public static void StopSubTimer(int TimingIndex) {
        if (TimingIndex < 0) {
            return;
        }

        if (SubTimingEnable[TimingIndex]) {
            SubTimers[TimingIndex].stop();
            SubDurations[TimingIndex] += SubTimers[TimingIndex].elapsed(TimeUnit.MICROSECONDS);
            SubTimers[TimingIndex].reset();
        }

    } // End StopSubTimer


    public static int synchronizeMPIVariable(int sync) throws MPIException {
        if (PWCUtility.MPI_Size > 1) {
            PWCUtility.StartSubTimer(PWCUtility.MPIBROADCASTTiming);
            // Note - MPI Call - Broadcast - int
            sync = PWCUtility.mpiOps.broadcast(sync, 0);
            PWCUtility.StopSubTimer(PWCUtility.MPIBROADCASTTiming);
        }
        return sync;
    }

    public static boolean synchronizeMPIVariable(boolean sync) throws MPIException {
        if (PWCUtility.MPI_Size > 1) {
            PWCUtility.StartSubTimer(PWCUtility.MPIBROADCASTTiming);
            // Note - MPI Call - Broadcast - boolean
            sync = PWCUtility.mpiOps.broadcast(sync, 0);
            PWCUtility.StopSubTimer(PWCUtility.MPIBROADCASTTiming);
        }
        return sync;
    }

    public static void writeClusterResults(String file, ArrayList<String> lines) {
        if (Strings.isNullOrEmpty(file)) {
            PWCUtility.printAndThrowRuntimeException(new IllegalArgumentException(Constants.ERR_EMPTY_FILE_NAME));
        }
        Path filePath = Paths.get(file);
        OpenOption mode = StandardOpenOption.CREATE;
        try (PrintWriter writer = new PrintWriter(Files.newBufferedWriter(filePath, Charset.defaultCharset(), mode),
                                                  true)) {
            writer.println();
            for (String line : lines) {
                writer.println(line);
            }
            writer.close();
        } catch (IOException e) {
            System.err.format("Failed writing cluster results due to I/O exception: %s%n", e);
        }
    } // End writeClusterResults

    public static String formatElapsedMillis(long elapsed){
        String format = "%dd:%02dH:%02dM:%02dS:%03dmS";
        short millis = (short)(elapsed % (1000.0));
        elapsed = (elapsed - millis) / 1000; // remaining elapsed in seconds
        byte seconds = (byte)(elapsed % 60.0);
        elapsed = (elapsed - seconds) / 60; // remaining elapsed in minutes
        byte minutes =  (byte)(elapsed % 60.0);
        elapsed = (elapsed - minutes) / 60; // remaining elapsed in hours
        byte hours = (byte)(elapsed % 24.0);
        long days = (elapsed - hours) / 24; // remaining elapsed in days
        return String.format(format, days, hours, minutes,  seconds, millis);
    }

} // End PWCUtility