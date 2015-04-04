package edu.indiana.soic.spidal.dapwc;

import com.google.common.base.Optional;
import com.google.common.base.Strings;
import com.google.common.io.Files;
import edu.indiana.soic.spidal.configuration.ConfigurationMgr;
import edu.indiana.soic.spidal.configuration.sections.PairwiseClusteringSection;
import edu.indiana.soic.spidal.general.IntArray;
import mpi.MPI;
import mpi.MPIException;
import net.openhft.affinity.AffinitySupport;
import org.apache.commons.cli.*;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.math.BigDecimal;
import java.math.RoundingMode;
import java.nio.ByteOrder;
import java.nio.charset.Charset;
import java.nio.file.OpenOption;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.Arrays;
import java.util.concurrent.TimeUnit;
import java.util.regex.Pattern;

import static edu.rice.hj.Module0.launchHabaneroApp;
import static edu.rice.hj.Module1.forallChunked;

/**
 * The {@code Program} class
 *
 * ConvergenceTest decides if to stop EM Computation based on change of epsilon
 * calculated in Dist.calculateEpsi
 *
 * shouldSplit invokes vectorclass.getMinimumEigenvalue to see if to split clusters
 * and which one to split
 *
 * split generates an extra cluster and sets initial Malpha(k) -- not epsilonalpha(k)
 *
 * initializeDist.Temperature calculates the maximum temperature from average distance
 *
 * ReadDataFromFile reads a symmetric distance matrix into memory
 *
 * getArrayindexFromMatrix is a general utility to find location in memory of a particular distance
 *
 * Main reads arguments, sets initial conditions, calls ReadDataFromFile and then getDist
 * Finally it outputs results
 *
 * Below t runs over threads, alpha runs over points, k runs over clusters
 * Dist.oldepsi[alpha][k] 		Old value of epsilon from previous iteration
 * Dist.epsi[alpha][k]		Value of epsilon calculated on current iteration
 * Dist.epsidiff[k] 	Change of epsilon on this iteration used for convergence
 * Dist.Malpha_k_[alpha][k]	Value of Malpha[k] at start of iteration
 * Dist.partialsum_C_k_[t][k] 	Contribution to C[k] from thread t
 * Dist.partialsum_A_k_[t][k] 	Contribution to A[k] from thread t
 * Dist.A_k_[k]			Value of A[k]
 * Dist.Balpha_k_[alpha][k] 	Value of Balpha[k]
 * Dist.C_k_[k] 			Value of C[k]
 *
 * vectorclass.initialvector[alpha]		Initial vector in data point space for Power method (same vector used for each cluster)
 * vectorclass.Ax[k][alpha]			Full iteration vector (in data point space) for cluster k
 * vectorclass.oldAx[k][alpha]			Previous full iteration vector (in data point space) for cluster k
 * vectorclass.eigens[k]				Current eigenvalue per cluster k for powermethod
 * vectorclass.maxeigens[k]			Maximum eigenvalue per cluster k for powermethod
 * vectorclass.eigenconverged[k]			Indicator of status of eigenvalue determination per cluster k
 *
 * Dist.Ncent			The current number of clusters initialized to 2
 * Program.ClustertoSplit			The cluster number that will be picked up to be split
 * Program.T			The current temperature
 * Program.iter			Iteration number counting EM convergence steps
 * Program.oldepsiset = 0              =0 says epsi has no previous value
 *
 * Program.lenindex[t]		Number of data points each thread t takes care of
 * Program.startindex[t]		Index of the first datapoint that thread t takes care
 * DistanceMatrix[Program.getArrayindexFromMatrix[ClusterCenter,ClusterIndex]] is distance between ClusterCenter and ClusterIndex (DistanceMatrix passed as argument)
 *
 * Program.MAX_VALUE		maximum number of parallel threads in a process (read as an argument)
 * PWCUtility.PointCount_Global		Number of Points (read as an argument) summed over threads and processes
 * Program.dataFile		Input data file name (read as an argument)
 * Program.maxNcent		Maximum number of cluster centers (read as an argument)
 * PWCUtility.pattern			String Pattern of parallel execution (e.g. 8x1x8 indicates [threads/process][processes/node][nodes])
 * 				(read as an argument)
 * Program.splitorexpandit		Execution parameter, scale factor of input data seta (read as an argument)
 *
 * Program.outputFile		Output (Excel) file name (preset)
 *
 * Program.cachelinesize = 32	Fixed definition of Cache line size in units of double
 * Program.max_change = 0.001	Convergence test condition for change in epsilon
 *
 * Program.Tmin			Minimal temperature (preset)
 * Program.Tinitial		Initial (largest) temperature calculated in initializeDist.Temperature
 * Program.CoolingFactor = 0.95	Cooling Factor in Annealing (preset)
 * Program.eigenvaluechange = .001		Convergence test in power method for eigenvalues
 * Program.eigenvectorchange = .001	Convergence test in power method for eigenvector
 * Program.poweriterations = 50		Maximum iteration count in power method
 *
 * MPI version
 * We arrange points in any order and assign them to processes/threads by simple one dimension decomposition
 * Points assigned to each process p run from PointIndex_Global1[p] to PointIndex_Global2[p] with PointTotal[p] points in all
 * PointTotal[p] = PointIndex_Global2[p] - PointIndex_Global1[p] + 1
 *
 * A) First Calculate Malpha(k) in parallel over alpha (unless we have this preset either for initial iteration or for split cluster)
 * B) Then Calculate C(k) as sums over alpha of Malpha(k). This is straightforwardly parallel over alpha with final reduction
 *
 * C) Then we need to calculate Balpha(k) as sums over beta (all Data points). This is dominant computation -- parallel over beta
 * It needs Malpha(k) Mbeta(k) Distance between alpha and beta, it calculates Balpha(k)and Bbeta(k)
 *
 *
 * D) Then we calculate A(k) by summing weighted Balpha(k) in parallel over alpha
 * E) This allows epsilon-alpha(k) and oldepsilon-alpha(k) to be calculated
 *
 * F) Now we check for convergence. This is a simple parallel reduction operation
 * If not converged, repeat from A)
 * G) If converged, check if possible to split (we have Tmin and maximum cluster counts)
 * If possible to split, calculate eigenvalues by power method
 * J) If should split, reset cluster numbers and Malpha(k) and go to A)
 * If split not needed, then reduce temperaturte and go to A)
 *
 * Calculation of critical eigenvalues is done in two almost identical steps considering Lmax maximum of C and then maximum of Lmax-C
 * These are done sepately for each cluster
 *
 * Power Method algorithm
 * This needs to iterate a computation Matrix * powervector where
 * powervector-alpha(k) is calculated in parallel over alpha similarly to calculation of Balpha(k)
 * H) It needs A(k), Malpha(k) Mbeta(k) Distance between alpha and beta, Balpha(k), Bbeta(k)
 * I) This is followed by a simple convergence test for power vector
 *
 * We could chose approach below which minimizes memory use but probably simpler approach that stores full matrix and not triangular part of it
 * allows better use of cache
 *
 * Storage of distances reflects order of processing terms that depend on two data points ClusterCenter,ClusterIndex. ClusterCenter,ClusterIndex components are calculated in home of ClusterCenter
 * For each data point ClusterCenter, there are indices index-nn defined into CheckeredDistance[index-nn]
 * CD-FirstPoint1[ClusterCenter] is first data point label ClusterIndex of type 1 associated with ClusterCenter and it is stored in CD-index1[ClusterCenter]'th Position in CheckeredDistance
 * CD-FirstPoint2[ClusterCenter] is first data point label ClusterIndex of type 2 associated with ClusterCenter and it is stored in CD-index2[ClusterCenter]'th Position in CheckeredDistance
 * CD-LastPoint2[ClusterCenter] is last data point label ClusterIndex of type 2 associated with ClusterCenter
 * CD-NumPoints[ClusterCenter] is number (type 1 and 2) of points ClusterIndex associated with ClusterCenter in CheckeredDistance
 * values of ClusterIndex start at CD-FirstPoint1[ClusterCenter] and increment by 2 as long as ClusterIndex < CD-FirstPoint2[ClusterCenter]
 * Then they start again at CD-FirstPoint2[ClusterCenter] and increment by 2 as long as ClusterIndex <= CD-LastPoint2[ClusterCenter]
 *
 * Block Computation is done for steps C) and H). This does parallelism minimizing memory access to Distance Matrix and b
 * Block is arranged to enhance communication and cache use
 * Other steps are straightforwardly parallel
 * There are two blocks -- one associated with Home index alpha and one with Away index beta
 * In C) and H), each Block has Malpha(k) and Balpha(k)
 */
public class Program
{
    private static Options programOptions = new Options();
    static {
        programOptions.addOption(String.valueOf(Constants.CMD_OPTION_SHORT_C),Constants.CMD_OPTION_LONG_C, true,
                                 Constants.CMD_OPTION_DESCRIPTION_C);
        programOptions.addOption(String.valueOf(Constants.CMD_OPTION_SHORT_N),Constants.CMD_OPTION_LONG_N,true,
                                 Constants.CMD_OPTION_DESCRIPTION_N);
        programOptions.addOption(String.valueOf(Constants.CMD_OPTION_SHORT_T),Constants.CMD_OPTION_LONG_T,true,
                                 Constants.CMD_OPTION_DESCRIPTION_T);
    }
	public static int ClusterCountOutput = 1; // Control Label Output = 0 not at all, = 1 at each count, = 2 also find centers
	public static int[] ClusterAssignments = null; // This gives for all points their cluster assignments (set in OutputClusterLabels)

	public static int cachelinesize = 32;

	public static String ControlFileName = ""; // Control File Name
	public static int ProcessingOption = 0; // Processing Option
	public static int PrintInterval = 5; // Output at this interval

	public static int FirstClusterNumber = 0; // Number (0 or 1) of first cluster
	public static boolean ContinuousClustering = false; // If true use the Ken Rose Continuous Clustering
	public static int maxNcent = 0; // maximum number of cluster center
	public static int InitialNcent = 1; // Initial value of Ncent
	public static int minimumclustercount = 25; // absorb clusters with fewer points than this

	public static double InitialCoolingFactor = 0.95; // InitialCooling Factor in Annealing
	public static double FineCoolingFactor = 0.995; // Refined Cooling Factor in Annealing
	public static double ConvergingCoolingFactor = 0.98; // Refined Cooling Factor in final iterations
	public static int Waititerations = 10; // Wait this number of Temperature iterations before splitting

	public static int Iterationatend = 1000; // Finish up EM Loop with at most this number of iterations
	public static int ConvergenceLoopLimit = 50; // Limit on EM Convergence for each step at a given temperature
	public static int CountPMExtraIterations = 0; // Count in Pairwise converging P-M for Continuous Clustering
	public static int EMIterationStepCountMax = -1; // Count Maximum Iterations per loop

	public static double Epsi_max_change = 0.001; //converge test condition for change in epsi
	public static double FreezingLimit = 0.002; // In finish stop when all freezing measures are < FreezingLimit

	public static boolean ConvergeIntermediateClusters = true; // If true converge intermediate cluster cases
	public static double ToosmalltoSplit = 50.0; // Size of Cluster measured by C_k_ that should not be split

	public static int PowerIterationLimit = 200; // Maximum iteration count in power method
	public static double MinEigtest = -0.01; // Factor of Pass 0 eigenvalue used to test negative full eigenvalue (can be negative).  Test easier to satisfy as MinEigtest gets bigger
	public static int EMlimit_Duplication = 200; // Limit on EM loops in Duplication mode
	public static int Eigenvalue_Methodology = 3; // Specify matrix whose second derivative is to be calculated
	// =0 Do not look at second derivative for splitting
	// =1 original simple method using unaltered second derivative matrix
	// =2  estimate of Duplicated second derivative matrix valid for continuousclustering true
	// =3 EM estimate of Duplicated second derivative matrix
	// =4 Sophisticated estimate of Duplicated second derivative matrix (used in test only as 3 faster and same)
	public static double eigenvaluechange = 0.001; // convergence test in power method for eigenvalues
	public static double eigenvectorchange = 0.001; // convergence test in power method for eigenvector
	public static int MaxNumberSplitClusters = 3; // System will split upto this number simultaneously
	public static int CountTotalPowerIterations = 0; // Count Power Iterations
	public static int CountTotalPowerErrors = 0; // Count Power Iterations

	public static double ExpectedChange = 0.0; // Expected Change in Cluster Count in Perturbation at Split
	public static int TestExpectedChange = -1; // Cluster Index for Expected Change in Cluster Count in Perturbation at Split; if = -1 Nothing Set
	public static int ExpectedChangeMethod = -1; // Method for Expected Change in Cluster Count in Perturbation at Split; 0 Change in C; 1 Special Ncent=1 case; =2 size of split vector
	public static double PreviousC = 0.0; // Original Cluster Count in Perturbation at Split
	public static double[] deltaCoverC = new double[5]; //  Array to average splitting changes
	public static int[] CountdeltaCoverC = new int[5]; //  Count to average splitting changes
	public static int Countmethodzero = 0; // Count Method 0 Perturbations
	public static int Countmethodtwo = 0; // Count Method 2 Perturbations

	public static boolean PerformEigenTest = false; // Set special Eigenvalue Testing
	public static double Epsi_max_change_Duplication = 0.001; //converge test condition for change in epsi for Duplication mode

	public static int PerturbationVehicle = 0; // = 0 Perturb Epsilon ; if 1 Perturb Malpha_k_
	public static double SplitPerturbationFactor = 1.0; // Normalization of Perturbation for Split

	public static int JiggleOption = 0; //    Control Perturbation of System every JiggleOption iterations
	public static int JigglePerturbation = 1; //    Control How Perturbation implemented
	//  JigglePerturbation = 1 do randomly
	//  JigglePerturbation = 2 do using lowest second derivative eigenvalue
	public static double JigglePerturbationFactor = 10.0; // Normalization of Perturbation for Jiggle

    public static int dataTypeSize = 2; // 2 for short
    public boolean isBigEndian = false; // true for Java style binary data and false for C# style binary data

    //Config Settings
    public static PairwiseClusteringSection config;

    /**
     * Pairwise (non-vector) Parallel clustering based on Deterministic Annealing algorithm
     * @param args command line arguments to the program, which should include
     *             -c "path to config file" -t "number of threads" -n "number of nodes"
     *             The options may also be given as longer names
     *             --configFile, --threadCount, and --nodeCount respectively
     */
	public static void main(String[] args) throws MPIException {
        Optional<CommandLine> parserResult = parseCommandLineArguments(args, programOptions);
        if (!parserResult.isPresent()){
            System.out.println(Constants.ERR_PROGRAM_ARGUMENTS_PARSING_FAILED);
            new HelpFormatter().printHelp(Constants.PROGRAM_NAME, programOptions);
            return;
        }

        CommandLine cmd = parserResult.get();
        if (!(cmd.hasOption(Constants.CMD_OPTION_LONG_C) &&
                cmd.hasOption(Constants.CMD_OPTION_LONG_N) &&
                cmd.hasOption(Constants.CMD_OPTION_LONG_T))){
            System.out.println(Constants.ERR_INVALID_PROGRAM_ARGUMENTS);
            new HelpFormatter().printHelp(Constants.PROGRAM_NAME, programOptions);
            return;
        }

		//  Read Metadata using this as source of other metadata
		ReadControlFile(cmd);

		//  Set up Parallelism
        //  Set up MPI and threads parallelism
        try {
            PWCParallelism.SetupParallelism(args, Long.toBinaryString(AffinitySupport.getAffinity()));
        } catch (MPIException e) {
            PWCUtility.printException(e);
            return; // End program on error
        }

		// Set up Decomposition of USED points
		PWCParallelism.SetParallelDecomposition();

        PWCUtility.PatternLabel = String.format("==== mpi-cluster(%1$s) ==== Threads:%2$s Clusters:%3$s PointCount_Global:%4$s ==== %5$s (%6$s) ==== ", PWCUtility.ParallelPattern, (new Integer(PWCUtility.ThreadCount)).toString(), (new Integer(Program.maxNcent)).toString(), PWCUtility.PointCount_Global, config.DistanceMatrixFile, new java.util.Date());

		PWCUtility.SALSAPrint(0, "\n" + PWCUtility.PatternLabel);

		// Initial Processing Complete
		PWCUtility.mpiOps.barrier(); // Make certain all processes have processed original data before writing updated

		//  read data into memory
		PWCParallelism.ReadDataFromFile(config.DistanceMatrixFile);

		// Start real work
		Program.ClusterAssignments = new int[PWCUtility.PointCount_Global]; // Set whenever clusters output

        // TODO - fix PlotTools and uncomment this if block
/*
		if (Program.ProcessingOption == 101)
		{ // Find centers of Previously determined Clusters

			int NumberofClusters = 0;
			edu.indiana.soic.spidal.Boxspidal.general.Box<Integer> boxNumberofClusters = new edu.indiana.soic.spidal.Boxspidal.general.Box<Integer>(NumberofClusters);
			Program.ReadClusterNumbers(PWCUtility.ClusterNumberfile, boxNumberofClusters, Program.ClusterAssignments, 0, PWCUtility.PointCount_Global);
			NumberofClusters = boxNumberofClusters.content;

			PWCUtility.SALSAPrint(0, PWCUtility.PatternLabel);
			PWCUtility.SALSAPrint(0, "Labels File: " + config.ClusterFile);
			String file = "CenterFile-M" + (new Integer(Program.maxNcent)).toString() + "-C" + NumberofClusters + ".txt";
			String directory1 = (new File(config.ClusterFile)).getParent();
			String CenterFileName = Paths.get(directory1, file).toString();
			String place = (new File(config.DistanceMatrixFile)).getParent();
			String MDSFileName = Paths.get(place, PWCUtility.addMDSfile).toString();
			String FullLabelFileName = PWCUtility.Labelfile;

			PWCUtility.LengthCut1 = 250; // Not used as yet
			PWCUtility.LengthCut2 = 250; // Cut on Second Sequence Length
			PWCUtility.MinimumDistanceCut = 0.001;
			PWCUtility.LinkCountinCenterFinding = 30;

			PWCUtility.SALSAPrint(0, "Find Centers for " + NumberofClusters + " Cluster Numbers read from " + PWCUtility.ClusterNumberfile + " Labels from " + FullLabelFileName + " MDS Read from " + MDSFileName + "\nLength Cuts " + (new Integer(PWCUtility.LengthCut1)).toString() + " " + (new Integer(PWCUtility.LengthCut2)).toString() + " Minimum Sequence Distance " + String.format("%1$.3f", PWCUtility.MinimumDistanceCut) + " Minimum Link Count " + (new Integer(PWCUtility.LinkCountinCenterFinding)).toString());

			//  Create Centers
			FindCenters.FindGroupCenters(Program.ClusterAssignments, NumberofClusters, 0, CenterFileName, MDSFileName, FullLabelFileName);
			// Produce pviz file
			PWCUtility.MPI_communicator.Barrier();
			if (PWCUtility.MPI_Rank == 0)
			{
				String plotDescription = PWCUtility.PointCount_Global + "_points_into_" + NumberofClusters + "_clusters_with_centers";
				String clusterNumberFile = config.ClusterNumberFile;
				PlotTools.CreatePlotWithCenters(CenterFileName, MDSFileName, clusterNumberFile, PWCUtility.CenterPointsPerCenterTypeInOutput, PWCUtility.CenterPlotFile, plotDescription);
			}
			PWCUtility.MPI_communicator.Barrier();

			if ((PWCUtility.MPIIOStrategy > 0) || (PWCUtility.MPI_Rank == 0))
			{
				PWCUtility.writeClusterResults(config.SummaryFile, PWCUtility.CosmicOutput);
			}
			PWCParallelism.TearDownParallelism(); //  Finalize MPI
			return;
		}
*/

		for (int CountDCoverC = 0; CountDCoverC < 5; CountDCoverC++)
		{
			Program.deltaCoverC[CountDCoverC] = 0.0;
			Program.CountdeltaCoverC[CountDCoverC] = 0;
		}

		// Actually Determine Clusters
		Program.Eigenvalue_Methodology = 3;
		if (Program.ContinuousClustering)
		{
			Program.Eigenvalue_Methodology = 2;
		}
		CountPMExtraIterations = 0;

		//  Set up Timing
		PWCUtility.InitializeTiming(7);
		PWCUtility.SetUpMPISubTimers(1, "");
		PWCUtility.SetUpSubTimer(0, "Splitting");

		//  Do Clustering
		Dist runitall = new Dist();
		runitall.getDist();

		//  End Timing
		PWCUtility.EndTiming();

		PWCUtility.SALSAPrint(0, "\n----------------------------- Final Output\n" + PWCUtility.PatternLabel);
		PWCUtility.SALSAPrint(0, "Labels File: " + config.ClusterFile);
		PWCUtility.SALSAPrint(0, "Timing Output: " + config.TimingFile);
		PWCUtility.SALSAPrint(0, "Initial Number of Centers: " + (new Integer(Program.InitialNcent)).toString());
		PWCUtility.SALSAPrint(0, "Maximum Number of Centers: " + (new Integer(Program.maxNcent)).toString());
		PWCUtility.SALSAPrint(0, "Actual Number of Centers: " + (new Integer(Dist.RunningPWC.Ncent)).toString());
		PWCUtility.SALSAPrint(0, "Minimum Number Points in Final Clusters " + (new Integer(Program.minimumclustercount)).toString());
		PWCUtility.SALSAPrint(0, "Converge Intermediate Clusters " + Program.ConvergeIntermediateClusters);
		PWCUtility.SALSAPrint(0, "Initial Cooling Factor in Annealing: " + String.format("%1$.4f", Program.InitialCoolingFactor));
		PWCUtility.SALSAPrint(0, "Refined Cooling Factor in Annealing: " + String.format("%1$.4f", Program.FineCoolingFactor));
		PWCUtility.SALSAPrint(0, "Converging(Final) Cooling Factor in Annealing: " + String.format("%1$.4f", Program.ConvergingCoolingFactor));
		PWCUtility.SALSAPrint(0, "Continuous Clustering " + Program.ContinuousClustering);
		PWCUtility.SALSAPrint(0, "Eigenvalue Methodology: " + (new Integer(Program.Eigenvalue_Methodology)).toString());
		PWCUtility.SALSAPrint(0, "Do not split Clusters smaller than this: " + (new Double(Program.ToosmalltoSplit)).toString());
		PWCUtility.SALSAPrint(0, "Pass 1 Eigenvalue Fractional Test: " + (new Double(Program.MinEigtest)).toString());
		PWCUtility.SALSAPrint(0, "Wait stages between splits: " + (new Integer(Program.Waititerations)).toString());
		PWCUtility.SALSAPrint(0, "Maximum Number of Simultaneous Cluster Splits " + (new Integer(Program.MaxNumberSplitClusters)).toString());
		PWCUtility.SALSAPrint(0, "Jiggle and Split Perturbation Vehicle: " + (new Integer(Program.PerturbationVehicle)).toString());
		PWCUtility.SALSAPrint(0, "Split and Split Perturbation Factor: " + (new Double(Program.SplitPerturbationFactor)).toString());
		PWCUtility.SALSAPrint(0, "Jiggle Option:" + (new Integer(Program.JiggleOption)).toString());
		PWCUtility.SALSAPrint(0, (new Integer(Program.JigglePerturbation)).toString() + " Jiggle Perturbation Method");
		PWCUtility.SALSAPrint(0, (new Double(Program.JigglePerturbationFactor)).toString() + " Jiggle Perturbation Factor");

		// Calculate Cluster Statistics
		int[] counts = new int[Dist.RunningPWC.Ncent];

		if (ClusterCountOutput > 0)
		{
			// Note Program.ClusterAssignments 
			Program.OutputClusterLabels(counts);
			if (ClusterCountOutput == 2)
			{
				String file = "CenterFile-M" + (new Integer(Program.maxNcent)).toString() + "-C" + (new Integer(Dist.RunningPWC.Ncent)).toString() + ".txt";
				String directory1 = (new File(config.ClusterFile)).getParent();
				String CenterFileName = Paths.get(directory1, file).toString();
				String place = (new File(config.DistanceMatrixFile)).getParent();
				String MDSFileName = Paths.get(place, PWCUtility.addMDSfile).toString();
				String FullLabelFileName = PWCUtility.Labelfile;
				FindCenters.FindGroupCenters(Program.ClusterAssignments, Dist.RunningPWC.Ncent, 0, CenterFileName, MDSFileName, FullLabelFileName);
			}

		}
		else
		{
			for (int ClusterIndex = 0; ClusterIndex < Dist.RunningPWC.Ncent; ClusterIndex++)
			{
				counts[ClusterIndex] = 0;
			}
		}

		PWCUtility.SALSAPrint(0, "\n******************\nT " + String.format("%1$.4f", Dist.RunningPWC.Temperature) + " Cluster " + (new Integer(Dist.RunningPWC.Ncent)).toString() + " Iter " + (new Integer(Dist.EMIterationCount)).toString() + " Extra Iter " + (new Integer(Dist.Extra_EMIterationCount)).toString() + " Max Iterations per Step " + (new Integer(Program.EMIterationStepCountMax)).toString());
		PWCUtility.SALSAPrint(0, "Extra P-M Iterations " + (new Integer(Program.CountPMExtraIterations)).toString() + " Total Number of Power Iterations " + (new Integer(Program.CountTotalPowerIterations)).toString() + " with errors " + (new Integer(Program.CountTotalPowerErrors)).toString());
		PWCUtility.SALSAPrint(0, "Iterations at End " + (new Integer(Program.Iterationatend)).toString() + " Iteration Limit to converge a fixed EM iteration " + Program.ConvergenceLoopLimit + " Count in Pairwise converging P-M for Continuous Clustering " + (new Integer(Program.CountPMExtraIterations)).toString() + " Maximum Iterations per loop " + (new Integer(Program.EMIterationStepCountMax)).toString());
		PWCUtility.SALSAPrint(0, "Freezing Convergence Limit " + String.format("%1$.5f", Program.FreezingLimit) + " Convergence test condition for change in epsilon " + String.format("%1$.5f", Program.Epsi_max_change));

		String nextline = " Center Averages (Counts) [Width] Width here and earlier is sum of d(i,j) weighted with M and divided by squared of occupation count C";
		for (int ClusterIndex = 0; ClusterIndex < Dist.RunningPWC.Ncent; ClusterIndex++)
		{
			double meandist = -2.0 * Dist.RunningPWC.A_k_[ClusterIndex];
			nextline += String.format("%1$.2f", Dist.RunningPWC.C_k_[ClusterIndex]) + " (" + (new Integer(counts[ClusterIndex])).toString() + ") [" + String.format("%1$.5f", meandist) + "] ";
		}
		PWCUtility.SALSAPrint(0, nextline);

		// Calculate Distance Correlation Matrix
		runitall.CorrelationCalculation();
		PWCUtility.SALSAPrint(0, "\nDistance Correlation Function");
		double successmeasure1 = 0.0;
		double AverageFreeze = 0.0;
		double MaxFreeze = 0.0;
		for (int ClusterIndex1 = 0; ClusterIndex1 < Dist.RunningPWC.Ncent; ClusterIndex1++)
		{
			AverageFreeze += Dist.RunningPWC.FreezingMeasure_k_[ClusterIndex1];
			MaxFreeze = Math.max(MaxFreeze, Dist.RunningPWC.FreezingMeasure_k_[ClusterIndex1]);

			String anotherline = String.valueOf(ClusterIndex1 + Program.FirstClusterNumber);
			for (int ClusterIndex2 = 0; ClusterIndex2 < Dist.RunningPWC.Ncent; ClusterIndex2++)
			{
				double tmp = Dist.DistanceCorrelation[ClusterIndex1][ClusterIndex2] / Math.sqrt(Dist.DistanceCorrelation[ClusterIndex1][ClusterIndex1] * Dist.DistanceCorrelation[ClusterIndex2][ClusterIndex2]);
				anotherline += "\t" + String.format("%1$.2f", tmp);
				if (ClusterIndex1 != ClusterIndex2)
				{
					successmeasure1 += tmp;
				}

			}
			anotherline += "   Count " + (new Integer(counts[ClusterIndex1])).toString();
			PWCUtility.SALSAPrint(0, anotherline);
		}

		successmeasure1 = successmeasure1 / ((Dist.RunningPWC.Ncent - 1) * Dist.RunningPWC.Ncent);
		PWCUtility.SALSAPrint(0, "\nOff Diagonal Average " + String.format("%1$.2f", successmeasure1));

		AverageFreeze = AverageFreeze / Dist.RunningPWC.Ncent;
		PWCUtility.SALSAPrint(0, "\nAverage Freezing Coefficient " + String.format("%1$.6f", AverageFreeze) + " Maximum Freezing Coefficient " + String.format("%1$.6f", MaxFreeze));

		//  Output Statistics on Perturbations at Cluster Splitting
		PWCUtility.SALSAPrint(0, "\nCluster Splitting Perturbation Statistics    Method 0: " + (new Integer(Program.Countmethodzero)).toString() + " Method 2: " + (new Integer(Program.Countmethodtwo)).toString());
		if (Program.CountdeltaCoverC[0] > 0)
		{
			PWCUtility.SALSAPrint(0, "Actual Delta C over C Method 1 (Ncent =1)    " + String.format("%1$.6f", Program.deltaCoverC[0] / Program.CountdeltaCoverC[0]) + " Count " + (new Integer(Program.CountdeltaCoverC[0])).toString());
		}
		if (Program.CountdeltaCoverC[2] > 0)
		{
			PWCUtility.SALSAPrint(0, "Well Predicted Delta C/C Method 1 (Ncent =1) " + String.format("%1$.6f", Program.deltaCoverC[2] / Program.CountdeltaCoverC[2]) + " Count " + (new Integer(Program.CountdeltaCoverC[2])).toString());
		}
		if (Program.CountdeltaCoverC[1] > 0)
		{
			PWCUtility.SALSAPrint(0, "Actual Delta C over C Method 0 or 2          " + String.format("%1$.6f", Program.deltaCoverC[1] / Program.CountdeltaCoverC[1]) + " Count " + (new Integer(Program.CountdeltaCoverC[1])).toString());
		}
		if (Program.CountdeltaCoverC[3] > 0)
		{
			PWCUtility.SALSAPrint(0, "Well Predicted Delta C over C Method 0 or 2  " + String.format("%1$.6f", Program.deltaCoverC[3] / Program.CountdeltaCoverC[3]) + " Count " + (new Integer(Program.CountdeltaCoverC[3])).toString());
		}
		if (Program.CountdeltaCoverC[4] > 0)
		{
			PWCUtility.SALSAPrint(0, "Naive Predicted Delta C over C Method 0 or 2 " + String.format("%1$.6f", Program.deltaCoverC[4] / Program.CountdeltaCoverC[4]) + " Count " + (new Integer(Program.CountdeltaCoverC[4])).toString());
		}

		/* Compute the duration between the initial and the end time ignoring print out. */
		PWCUtility.SALSAPrint(0, "\nTotal Time excluding I/O  " + PWCUtility.formatElapsedMillis(PWCUtility.mainDuration) + " " + String.format("%1$.0f", PWCUtility.HPDuration * .001));
		nextline = "Partial Times ";
		for (int itimer = 0; itimer < PWCUtility.NumberofSubTimings; itimer++)
		{
			if (PWCUtility.SubTimingEnable[itimer])
			{
				double tmp = PWCUtility.SubDurations[itimer] / PWCUtility.HPDuration;
				nextline += PWCUtility.SubTimingNames[itimer] + " " + String.format("%1$.0f", PWCUtility.SubDurations[itimer] * .001) + " " + String.format("%1$.4f", tmp) + " ";
			}
		}
		PWCUtility.SALSAPrint(0, nextline);

		if ((PWCUtility.MPIIOStrategy > 0) || (PWCUtility.MPI_Rank == 0))
		{
			PWCUtility.writeClusterResults(config.SummaryFile, PWCUtility.CosmicOutput);
            WriteTimingFile(config.TimingFile, PWCUtility.mainDuration, PWCUtility.HPDuration, PWCUtility.ThreadCount,
                            PWCUtility.MPIperNodeCount, PWCUtility.NodeCount, PWCUtility.PointCount_Process,
                            PWCUtility.PointCount_Global, Program.maxNcent, PWCUtility.SubDurations[0] * 0.001,
                            PWCUtility.SubDurations[0] / PWCUtility.HPDuration, PWCUtility.SubDurations[1] * 0.001,
                            PWCUtility.SubDurations[1] / PWCUtility.HPDuration, PWCUtility.SubDurations[2] * 0.001,
                            PWCUtility.SubDurations[2] / PWCUtility.HPDuration, PWCUtility.SubDurations[3] * 0.001,
                            PWCUtility.SubDurations[3] / PWCUtility.HPDuration, PWCUtility.SubDurations[4] * 0.001,
                            PWCUtility.SubDurations[4] / PWCUtility.HPDuration, config.DistanceMatrixFile,
                            new java.util.Date(), MPI.getProcessorName());

		}
		PWCParallelism.TearDownParallelism(); //  Finalize MPI
	}
	// End Main

    /**
     * Parse command line arguments
     * @param args Command line arguments
     * @param opts Command line options
     * @return An <code>Optional&lt;CommandLine&gt;</code> object
     */
    private static Optional<CommandLine> parseCommandLineArguments(String [] args, Options opts){

        CommandLineParser optParser = new GnuParser();

        try {
            return Optional.fromNullable(optParser.parse(opts, args));
        } catch (ParseException e) {
            System.out.println(e);
        }
        return Optional.fromNullable(null);
    }


	public static void OutputClusterLabels(int[] OccupationCounts) throws MPIException {
		// Generate Cluster Labels (cluster number) from final epsi
		final int[] labels = new int[PWCUtility.PointCount_Process];
		int[][] partialsum_OccupationCounts = new int[PWCUtility.ThreadCount][];

		for (int ThreadNo = 0; ThreadNo < PWCUtility.ThreadCount; ThreadNo++)
		{
			partialsum_OccupationCounts[ThreadNo] = new int[Dist.RunningPWC.Ncent + cachelinesize];
		}

		//  Parallel Section setting cluster labels
        // Note - parallel for
        launchHabaneroApp(() -> {
            forallChunked(0, PWCUtility.ThreadCount - 1, (threadIndex) ->
                    {
                        if (PWCUtility.bindThreads){
                            AffinitySupport.setAffinity(1L<<PWCUtility.bindThreadToCore[threadIndex]);
                        }
                        Arrays.fill(partialsum_OccupationCounts[threadIndex], 0, Dist.RunningPWC.Ncent, 0);
                        int indexlen = PWCUtility.PointsperThread[threadIndex];
                        int beginpoint = PWCUtility.StartPointperThread[threadIndex] - PWCUtility.PointStart_Process;
                        for (int index = beginpoint; index < indexlen + beginpoint; index++) {
                            double distmin = 9999999999999.0;
                            int knear = 0;
                            for (int ClusterIndex = 0; ClusterIndex < Dist.RunningPWC.Ncent; ClusterIndex++) {
                                if (Dist.RunningPWC.Epsilonalpha_k_[index][ClusterIndex] < distmin) {
                                    distmin = Dist.RunningPWC.Epsilonalpha_k_[index][ClusterIndex];
                                    knear = ClusterIndex;
                                }
                            }
                            labels[index] = knear + Program.FirstClusterNumber;
                            partialsum_OccupationCounts[threadIndex][knear]++;
                        }
                    }
            );
        });


        int[] LocalOccupationCounts = new int[Dist.RunningPWC.Ncent];
		Arrays.fill(LocalOccupationCounts, 0, Dist.RunningPWC.Ncent, 0);
		for (int ThreadNo = 0; ThreadNo < PWCUtility.ThreadCount; ThreadNo++)
		{
			for (int ClusterIndex = 0; ClusterIndex < Dist.RunningPWC.Ncent; ClusterIndex++)
			{
				LocalOccupationCounts[ClusterIndex] += partialsum_OccupationCounts[ThreadNo][ClusterIndex];
			}
		}


		if (PWCUtility.MPI_Size > 1)
		{
			PWCUtility.StartSubTimer(PWCUtility.MPIREDUCETiming);
            // Note - MPI Call - Allreduce - int - sum
            PWCUtility.mpiOps.allReduce(LocalOccupationCounts, MPI.SUM);
			PWCUtility.StopSubTimer(PWCUtility.MPIREDUCETiming);
		}

        System.arraycopy(LocalOccupationCounts, 0, OccupationCounts, 0, Dist.RunningPWC.Ncent);
        if (PWCUtility.PreciseTimer.isRunning()){ // Check is necessary as this gets called after EndTiming() as well
            PWCUtility.PreciseTimer.stop(); // temporarily stop timer
            PWCUtility.HPDuration += PWCUtility.PreciseTimer.elapsed(TimeUnit.MICROSECONDS);
        }

		String directory = (new File(config.ClusterFile)).getParent();
		String file = Files.getNameWithoutExtension(config.ClusterFile) + "-M" + Program.maxNcent + "-C" + Dist.RunningPWC.Ncent + Files.getFileExtension(config.ClusterFile);

		String ClusternumberFileName = Paths.get(directory, file).toString();

		if (PWCUtility.MPIIOStrategy > 0)
		{
			WriteClusterFile(ClusternumberFileName, labels, PWCUtility.PointCount_Process, PWCUtility.PointStart_Process, false);
		}


		int MPItag = 100;
		if (PWCUtility.MPI_Rank == 0)
		{
			WriteClusterFile(ClusternumberFileName, labels, PWCUtility.PointCount_Process, PWCUtility.PointStart_Process, false);

            // Note - parallel for
            launchHabaneroApp(() -> {
                forallChunked(0, PWCUtility.ThreadCount - 1, (threadIndex) ->
                        {
                            if (PWCUtility.bindThreads){
                                AffinitySupport.setAffinity(1L<<PWCUtility.bindThreadToCore[threadIndex]);
                            }
                            int indexlen = PWCUtility.PointsperThread[threadIndex];
                            int beginpoint = PWCUtility.StartPointperThread[threadIndex] - PWCUtility.PointStart_Process;
                            for (int index = beginpoint; index < indexlen + beginpoint; index++) {
                                Program.ClusterAssignments[index + PWCUtility.PointStart_Process] = labels[index] - Program.FirstClusterNumber;
                            }
                        }
                );
            });

            for (int MPISource = 1; MPISource < PWCUtility.MPI_Size; MPISource++)
			{
                // Note - MPI Call - Receive - MPIPacket
				MPIPacket fromsource = PWCUtility.mpiOps.receive(MPISource, MPItag, MPIPacket.Type.Integer);
				if (PWCUtility.MPIIOStrategy == 0)
				{
					WriteClusterFile(ClusternumberFileName, fromsource::getMArrayIntAt, fromsource.getNumberOfPoints(), fromsource.getFirstPoint(), true);
				}

				for (int index = 0; index < fromsource.getNumberOfPoints(); index++)
				{
					Program.ClusterAssignments[index + fromsource.getFirstPoint()] = fromsource.getMArrayIntAt(index) - Program.FirstClusterNumber;
				}
			}
		} // End Root Process that receives cluster assignments from afar
		else
		{
			MPIPacket tosend = MPIPacket.newIntegerPacket(PWCUtility.PointCount_Process);
			tosend.setFirstPoint(PWCUtility.PointStart_Process);
			tosend.setNumberOfPoints(PWCUtility.PointCount_Process);
            for (int index = 0; index < PWCUtility.PointCount_Process; index++) {
                tosend.setMArrayIntAt(index,labels[index]);
            }
            // Note - MPI Call - Send - MPIPacket
			PWCUtility.mpiOps.send(tosend, 0, MPItag);
		}
        // Note - MPI call - Broadcast - int[]
		PWCUtility.mpiOps.broadcast(Program.ClusterAssignments, 0);
		PWCUtility.mpiOps.barrier();
        if (!PWCUtility.timingCompleted){ // If PWCUtility.EndTiming() has been called then no need to reset the precise timer
            /* Restart Timer - requires reset and start */
            PWCUtility.PreciseTimer.reset();
            PWCUtility.PreciseTimer.start();
        }
	}

	public static void ReadControlFile(CommandLine cmd)
	{
        config = ConfigurationMgr.LoadConfiguration(cmd.getOptionValue(Constants.CMD_OPTION_LONG_C)).pairwiseClusteringSection;
        Program.ControlFileName = cmd.getOptionValue(Constants.CMD_OPTION_LONG_C);
        PWCUtility.NodeCount = Integer.parseInt(cmd.getOptionValue(Constants.CMD_OPTION_LONG_N));
        PWCUtility.ThreadCount = Integer.parseInt(cmd.getOptionValue(Constants.CMD_OPTION_LONG_T));

        // Override section's node and thread counts with command line values if different
        if (config.getNodeCount() != PWCUtility.NodeCount) {
            config.setNodeCount(PWCUtility.NodeCount);
        }
        if (config.getThreadCount() != PWCUtility.ThreadCount) {
            config.setThreadCount(PWCUtility.ThreadCount);
        }

	    PWCUtility.PointCount_Global = config.NumberDataPoints;
	    Program.ProcessingOption = config.ProcessingOption;
	    Program.maxNcent = config.MaxNcent;
	    PWCUtility.ThreadCount = config.ThreadCount;     // Number of Threads
	    PWCUtility.NodeCount = config.NodeCount;       // Number of Nodes
	    PWCUtility.MPIIOStrategy = config.MPIIOStrategy;   // Controls strategy of file handling with MPI =0 is ONE FILE
	    Program.ToosmalltoSplit = config.TooSmallToSplit;
	    Program.Waititerations = config.WaitIterations;
	    Program.MinEigtest = config.MinEigTest;
	    Program.Epsi_max_change = config.EpsiMaxChange;     //converge test condition for change in epsi
	    Program.InitialCoolingFactor = config.InitialCoolingFactor;  // InitialCooling Factor in Annealing
	    Program.FineCoolingFactor = config.FineCoolingFactor;    // Refined Cooling Factor in Annealing
	    Program.eigenvaluechange = config.EigenValueChange;    // convergence test in power method for eigenvalues
	    Program.eigenvectorchange = config.EigenVectorChange;   // convergence test in power method for eigenvector
	    Program.Iterationatend = config.IterationAtEnd;       // Finish up EM Loop with this number of iterations
	    Program.ConvergenceLoopLimit = config.ConvergenceLoopLimit; // Limit on EM Convergence
	    Program.FreezingLimit = config.FreezingLimit;       // In finish stop when all freezing measures are < FreezingLimit
	    Program.PowerIterationLimit = config.PowerIterationLimit;   // Maximum iteration count in power method
	    Program.ConvergeIntermediateClusters = config.ConvergeIntermediateClusters;
	    PWCUtility.DebugPrintOption = config.DebugPrintOption;
	    PWCUtility.ConsoleDebugOutput = config.ConsoleDebugOutput;
	    Program.ContinuousClustering = config.ContinuousClustering;
	    PWCUtility.Labelfile = config.LabelFile;
	    PWCUtility.addMDSfile = config.AddMdsFile;
	    PWCUtility.ClusterNumberfile = config.ClusterNumberFile;
	    PWCUtility.addMDS = config.AddMds;
	    PWCUtility.BucketFractions = config.BucketFractions;
	    PWCUtility.NumberofBuckets = PWCUtility.BucketFractions != null ? PWCUtility.BucketFractions.length : 0;
	    PWCUtility.NumberofCenters = config.NumberOfCenters;
	    PWCUtility.CenterPointsPerCenterTypeInOutput = config.CenterPointsPerCenterTypeInOuput;
	    PWCUtility.CenterPlotFile = config.CenterPlotFile;
        PWCUtility.dataTypeSize = config.getDataTypeSize();
        PWCUtility.endianness = config.isBigEndian() ? ByteOrder.BIG_ENDIAN : ByteOrder.LITTLE_ENDIAN;
        PWCUtility.isMemoryMapped = config.isMemoryMapped();
        PWCUtility.bindThreads = config.isBindThreads();
    }
    public static void WriteClusterFile(String file, int[] labels, int dataPoints, int startPosition, boolean append){
        WriteClusterFile(file, i -> labels[i], dataPoints, startPosition, append);
    }
    public static void WriteClusterFile(String file, Integer[] labels, int dataPoints, int startPosition, boolean append){
        WriteClusterFile(file, i -> labels[i], dataPoints, startPosition, append);
    }
    public static void WriteClusterFile(String file, IntArray labels, int dataPoints, int startPosition, boolean append)
	{
        if (Strings.isNullOrEmpty(file)) {
            PWCUtility.printAndThrowRuntimeException(new IllegalArgumentException(Constants.ERR_EMPTY_FILE_NAME));
        }
        Path filePath = Paths.get(file);
        OpenOption mode = append ? StandardOpenOption.APPEND : StandardOpenOption.CREATE;
        try (PrintWriter writer = new PrintWriter(
                java.nio.file.Files.newBufferedWriter(filePath, Charset.defaultCharset(), mode),
                true)) {
            for (int i = 0; i < dataPoints; i++)
            {
                writer.println((i + startPosition) + "\t" + labels.get(i));
            }
            writer.close();
        } catch (IOException e) {
            System.err.format("Failed writing cluster results due to I/O exception: %s%n", e);
        }
	}

    public static void WriteTimingFile(String fileName, long duration, double HPDuration, int ThreadCount, int MPIperNodeCount, int NodeCount, int PointCount_Process, int PointCount_Global, int maxNcent, double PTsplitting, double Ratio1, double MPIReduce, double Ratio2, double MPISRBasic, double Ratio3, double MPISREigen, double Ratio4, double MPIBdcast, double Ratio5, String DataFileName, java.util.Date CurrentTime, String ProcessorName)
    {
        String header =
                "Duration(ms)\tHPDuration(us)" +
                        "\tThread#\tMPIperNode\tNode#\tPattern\tParallelism\tPoint#/Process\tPoint#/Global\tmaxNcent" +
                        "\tPTsplitting\tRatio1\tMPIReduce\tRatio2\tMPISRBasic\tRatio3\tMPISREigen\tRatio4\tMPIBdcast" +
                        "\tRatio5\tDataFileName\tCurrentTime\tProcessorName";

        String format = "%1$s\t%2$s\t%3$s\t%4$s\t%5$s\t%6$s\t%7$s\t%8$s\t%9$s\t%10$s\t%11$s\t%12$s\t%13$s\t%14$s\t%15$s\t%16$s\t%17$s\t%18$s\t%19$s\t%20$s\t%21$s\t%22$s\t%23$s";
        String line =
                String.format(format,
                              Math.round(duration),
                              Math.round(HPDuration),
                              ThreadCount,
                              MPIperNodeCount,
                              NodeCount,
                              String.format("%1$sx%2$sx%3$s", ThreadCount, MPIperNodeCount, NodeCount),
                              (ThreadCount * MPIperNodeCount * NodeCount),
                              PointCount_Process,
                              PointCount_Global,
                              maxNcent,
                              Math.round(PTsplitting),
                              new BigDecimal(Ratio1).setScale(4, RoundingMode.UP).doubleValue(),
                              Math.round(MPIReduce),
                              new BigDecimal(Ratio2).setScale(4, RoundingMode.UP).doubleValue(),
                              Math.round(MPISRBasic),
                              new BigDecimal(Ratio3).setScale(4, RoundingMode.UP).doubleValue(),
                              Math.round(MPISREigen),
                              new BigDecimal(Ratio4).setScale(4, RoundingMode.UP).doubleValue(),
                              Math.round(MPIBdcast),
                              new BigDecimal(Ratio5).setScale(4, RoundingMode.UP).doubleValue(),
                              DataFileName, CurrentTime,
                              ProcessorName); // Name of input data file -  MPIBdcast vs HPDuration -  MPI broadcast -  MPISREigen vs HPDuration -  MPI SR Eigen -  MPISRBasic vs HPDuration -  MPI SR Basic -  MPIReduce vs HPDuration -  MPI reduce -  PTsplitting vs HPDuration -  Partial times splitting -  cluster# -  Global points -  Local points -  Pattern -  Node# aquired -  Process# per node -  Thread# per process -  High performance timer -  Total time

        Path filePath = Paths.get(fileName);
        try {
            if (!java.nio.file.Files.exists(filePath))
            {
                java.nio.file.Files.write(filePath, header.getBytes(), StandardOpenOption.CREATE);
            }
            java.nio.file.Files.write(filePath, (System.lineSeparator() + line).getBytes(), StandardOpenOption.APPEND);
        } catch (IOException e) {
            System.err.format("Failed writing timing file due to I/O exception: %s%n", e);
        }
    }

	// Read Cluster Numbers
	public static int  ReadClusterNumbers(String fname, int[] ClusterNumbers, int BeginPoint, int PointstoRead)
	{
        if (Strings.isNullOrEmpty(fname)) {
            PWCUtility.printAndThrowRuntimeException(new IllegalArgumentException(Constants.ERR_EMPTY_FILE_NAME));
        }

		boolean success = false;
		String line = "NotSet";
		int count = 0;
		int NumberofClusters = 0;

        try (BufferedReader br = java.nio.file.Files.newBufferedReader(Paths.get(fname), Charset.defaultCharset()))
		{
            Pattern pattern = Pattern.compile("[\t ]");
            while ((line = br.readLine()) != null)
            {
                if (Strings.isNullOrEmpty(line)) {
                    continue; // continue on empty lines - "while" will break on null anyway;
                }
                String[] splits = pattern.split(line.trim());
                int NumberPosition = 1;
                if (splits.length > 2)
                {
                    NumberPosition = 4;
                }
                int index = Integer.parseInt(splits[0]);
                if (index < BeginPoint)
                {
                    continue;
                }
                if (index >= BeginPoint + PointstoRead)
                {
                    break;
                }
                int readnumber = Integer.parseInt(splits[NumberPosition]);
                ClusterNumbers[index - BeginPoint] = readnumber - Program.FirstClusterNumber;
                NumberofClusters = Math.max(NumberofClusters, readnumber);
                if (readnumber < Program.FirstClusterNumber)
                {
                    PWCUtility.printAndThrowRuntimeException(
                            "Illegal Cluster Number on Cluster Number / MDS file Index " +
                                    (new Integer(index)).toString() + " Cluster Number Read " + readnumber +
                                    " First Cluster # " + (new Integer(Program.FirstClusterNumber)).toString());
                }
                count++;
            }
            if (count != PointstoRead)
            {
                PWCUtility.printAndThrowRuntimeException(
                        "Illegal count on Cluster Number file " + (new Integer(count)).toString() + " " +
                                (new Integer(PointstoRead)).toString());
            }
            success = true;
		}
		catch (Exception e)
		{
            PWCUtility.printAndThrowRuntimeException("Failed reading Cluster Number data: count " + count + " Line" +
                                       " " + line + "\n" + e);
		}
		if (!success)
		{
			PWCUtility.printAndThrowRuntimeException(
                    "Cluster Number File read error: count " + (new Integer(count)).toString() + " Line " + line +
                            "\n" + fname);
		}

		if (Program.FirstClusterNumber == 0)
		{
			++NumberofClusters;
		}
        return NumberofClusters;

	} // End Read Cluster Numbers
}