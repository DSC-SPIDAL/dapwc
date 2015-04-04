package edu.indiana.soic.spidal.dapwc;

//	Class Dist **********************************************************
// ------------------------------------------------------------------------------
//	getDist controls complete pairwise computation
//	It sets up problem and calls parallel threads PairwiseThread. 

//	PairwiseThread does major pairwise computation looping over all EM iterations

//	In first iteration it is assumed that Malpha(k) is set  and epsilonalpha(k) is calculated
//	Other iterations use value of epsilonalpha(k) to first calculate a new value of Malpha(k)
//	One calls calculatEpsi to find A(k) Balpha(k) from C(k) Malpha(k) and distance matrix
//	Then Program.convergenceTest in Program decides if current EM iteration converged
//	Program.ShouldWeSplit is called to see if clusters should be split
//	Program.DoTheSplit is called to split chosen cluster
//	Dist.Temperature is reduced and EM Iteration continues
//	
//	calculateEpsi finds A(k) Balpha(k) from C(k) Malpha(k) and distance matrix
//	These give the new value of epsilonalpha(k)

// ------------------------------------------------------------------------------
// public void getDist()
// Implements logic of EM and Temperature loop and decision on splitting
// 
// ------------------------------------------------------------------------------
// public static void PrintIteration(string looptype)
// Print summary of EM iteration with a given temperature value
// 
// ------------------------------------------------------------------------------
// public void PairwiseThread()
// Control calculation of key functions starting by calculating M and then calling
// calculateEpsi(Malpha_k_, A_k_, Balpha_k_, C_k_, Epsilonalpha_k_, Dist.Ncent)
// 
// ------------------------------------------------------------------------------
// public void calculateEpsi(double[][] localMalpha_k_, double[] localA_k_, double[][] localBalpha_k_, double[] localC_k_, double[][] localepsi, int localNcent)
//	Perform multiple parallel steps calculating A_k_, Balpha_k_, epsi and differences
//	Value of M is assumed
// Complicated MPI parallelism broadcasting M in parts with a ring
// 
// ------------------------------------------------------------------------------
// public static int convergenceTest(double ChangeLimit)
// The change of sum(points) Delta(epsi)/# Points should be less than argument ChangeLimit for each cluster
// epsiDiff calculated in CalculateEpsi while Epsi are updated
//  Return 0 if not converged; 1 if converged and can continue; 2 if converged but no refinement
// 
// Test if iteration has converged for given temperature
// 
// ------------------------------------------------------------------------------
// public void SaveCurrentTask()
// {   // Save current task that should be split but we need to converge first
// public void RestorePreviousTask()
// {   // Restore previous task that should be split but we needed to converge current cluster configuration first
// 
// Currently before splitting we converge current cluster configuration to zero Temperature before saving
// These routines allow one to save configuration before split and return to it
// 
// ------------------------------------------------------------------------------
// JiggleClusters() // Jiggle Values of Epsilon -- leaving other Cluster Parameters
// 
// Either do totally randomly or along most negative eigenvector Based on Program.JigglePerturbation
// Eigenvalue analysis called from JiggleClusters()
// 
// Set new values ofd epsilon(alpha, k).
// Calculate change in C and M for output only
// Set flag that M has to be calculated
// 
// ------------------------------------------------------------------------------
// public bool shouldweSplit()
// Decide if to split based on negative eigenvalue of second derivative matrix
// MinimumEigenvalue is MINIMUM eigenvalue
// ClustertoSplit is cluster number of this
// D is distance matrix
//  Skip first time this happens as initialize to two clusters
// 
// Control by Program.Eigenvalue_Methodology
// =0 Do not look at second derivative for splitting
// =1 original simple method using unaltered second derivative matrix (use for Continuous Clustering)
// =2 quick and dirty estimate of Duplicated second derivative matrix
// =3 EM estimate of Duplicated second derivative matrix (DEFAULT). Each cluster has weight increased to 2 and system iterated to convergence
// =4 Sophisticated estimate of Duplicated second derivative matrix (used in test only as 3 faster and same)
// 
// Program.PerformEigenTest is defaulted to False.
// If true perform a test of approach
// 
// shouldweSplit() calls eigenvalue routine
// 
// ------------------------------------------------------------------------------
// public void dothesplit()
//  Do a split on the identified cluster found in shouldweSplit()
//  Malpha_k_ is value of M indexed by (point,cluster)
//  New clusters stored in last cluster and old position
//  Cluster count is incremented by one
//  ClustertoSplit is cluster number to split
//	Embarassingly Parallel over Processes -- no MPI calls needed
// 
// public static int PerturbationVehicle = 0; // = 0 Perturb Epsilon ; if 1 Perturb Malpha_k_
// public static double SplitPerturbationFactor = 1.0; // Normalization of Perturbation for Split
// 
// M is halved and either M or epsilon perturbed
// epsilon uses eigenvector for perturbation direction (assumes currently that calculated)
// 
// -------------------------------------------------------------------------------
// public static void initializeTemperature()
//	Find initial Temperature that can be used to start with one cluster
// Just need an upper bound as whips through this part ofd program
// 
// ------------------------------------------------------------------------------
// Calculate Distance Correlation Function
// public void CorrelationCalculation()
// Use at end of run for final clusters

// Calculate Pairwise Algorithm

import edu.rice.hj.api.SuspendableException;
import net.openhft.affinity.AffinitySupport;

import static edu.rice.hj.Module0.launchHabaneroApp;
import static edu.rice.hj.Module1.forallChunked;

public class ClusteringSolution
{
	public double[][] Old_Epsilonalpha_k_; // Previous value of Epsilon
	public double[][] Epsilonalpha_k_; // Epsilon
	public double[][] Best_Epsilonalpha_k_; // best Epsilon in loop over Clusters Duplicated
	public double[][] Master_Epsilonalpha_k_; // Master Epsilon in loop over Clusters Duplicated
	public double[][] Malpha_k_; // Probability that point in Cluster k
	public double[][] Previous_Malpha_k_; // Previous value of Malpha_k_ used to test converge of p(k) Malpha(k) loop
	public double[][] Balpha_k_; // B(alpha,k) = sum_N d(ClusterCenter, a) M_i(k) / C(k)

	public double[] A_k_; // A(k) = -0.5 sum_N Balpha(k) M_a(k) / C(k)
	public double[] C_k_; // Summation of C(k) over parallel threads
	public double[] Weight_k_; // weight -- typically 1 -- of each cluster
	public double[] FreezingMeasure_k_; // Freezing Measure for Cluster
	public double[] P_k_; // P(k) used in Continuous Clustering

	public double[][] Saved_Ax; // New power vector
	public double[][] Saved_oldAx; // Old power vector

	public double PairwiseHammy = 0.0; // Value of Hamiltonian
	public double OldHammy = 0.0; //  Previous value of Hamiltonian
	public int Ncent = 1; //  the current number of clusters
	public int ClustertoSplit = -1; //  the cluster label that will be split
	public double MinimumEigenvalue = 0.0; // Minimum full eigenvalue associated with ClustertoSplit
	public double ActualCoolingFactor = Program.InitialCoolingFactor; // Actual Cooling Factor
	public double Temperature = Dist.Tinitial; // The current temperature
	public int IterationSetAt = 0; // Iteration set at
	public boolean Axset = false; // If True Saved_Ax and Saved_oldAx set
	public boolean SolutionSet = false; // If True Solution Set

	public static int NumberofPointsinProcess = -1; // Number of Points in Process
	public static int MaximumNumberClusters = 0; // Maximum Number of Centers
	public static int cachelinesize = 0; // Cacheline increment

    public static void SetParameters(int NumberofPointsinProcessINPUT, int MaximumNumberClustersINPUT,
                                     int cachelinesizeINPUT) {
        if (NumberofPointsinProcess > 0) {
            return;
        }
        NumberofPointsinProcess = NumberofPointsinProcessINPUT;
        MaximumNumberClusters = MaximumNumberClustersINPUT;
        cachelinesize = cachelinesizeINPUT;
    }

    public ClusteringSolution() {
        if (NumberofPointsinProcess < 0) {
            PWCUtility.printAndThrowRuntimeException("NumberofPointsinProcess Unset");
        }
        Epsilonalpha_k_ = new double[NumberofPointsinProcess][];
        Old_Epsilonalpha_k_ = new double[NumberofPointsinProcess][];
        Best_Epsilonalpha_k_ = new double[NumberofPointsinProcess][];
        Master_Epsilonalpha_k_ = new double[NumberofPointsinProcess][];
        Malpha_k_ = new double[NumberofPointsinProcess][];
        Previous_Malpha_k_ = new double[NumberofPointsinProcess][];
        Balpha_k_ = new double[NumberofPointsinProcess][];

        Saved_Ax = new double[NumberofPointsinProcess][];
        Saved_oldAx = new double[NumberofPointsinProcess][];

        // Note - parallel for

        launchHabaneroApp(() -> {
            forallChunked(0, PWCUtility.ThreadCount - 1, (threadIndex) -> {
                if (PWCUtility.bindThreads){
                    AffinitySupport.setAffinity(1L<<PWCUtility.bindThreadToCore[threadIndex]);
                }
                int indexlen = PWCUtility.PointsperThread[threadIndex];
                int beginpoint = PWCUtility.StartPointperThread[threadIndex] - PWCUtility.PointStart_Process;
                for (int ProcessPointIndex = beginpoint; ProcessPointIndex < indexlen + beginpoint; ProcessPointIndex++) {
                    Epsilonalpha_k_[ProcessPointIndex] = new double[MaximumNumberClusters];
                    Old_Epsilonalpha_k_[ProcessPointIndex] = new double[MaximumNumberClusters];
                    Best_Epsilonalpha_k_[ProcessPointIndex] = new double[MaximumNumberClusters];
                    Master_Epsilonalpha_k_[ProcessPointIndex] = new double[MaximumNumberClusters];
                    Malpha_k_[ProcessPointIndex] = new double[MaximumNumberClusters];
                    Previous_Malpha_k_[ProcessPointIndex] = new double[MaximumNumberClusters];
                    Balpha_k_[ProcessPointIndex] = new double[MaximumNumberClusters];
                    Saved_Ax[ProcessPointIndex] = new double[MaximumNumberClusters + cachelinesize];
                    Saved_oldAx[ProcessPointIndex] = new double[MaximumNumberClusters + cachelinesize];
                    for (int ClusterIndex = 0; ClusterIndex < Program.maxNcent; ClusterIndex++) {
                        Epsilonalpha_k_[ProcessPointIndex][ClusterIndex] = 0.0;
                        Old_Epsilonalpha_k_[ProcessPointIndex][ClusterIndex] = 0.0;
                    }
                }
            });
        });

        C_k_ = new double[MaximumNumberClusters + cachelinesize]; // Final value of C(k)
        A_k_ = new double[MaximumNumberClusters + cachelinesize];
        P_k_ = new double[MaximumNumberClusters + cachelinesize];

        FreezingMeasure_k_ = new double[MaximumNumberClusters + cachelinesize];
        Weight_k_ = new double[MaximumNumberClusters + cachelinesize];

        Axset = false;
    }

    public static void CopySolution(ClusteringSolution From, ClusteringSolution To) {
        To.PairwiseHammy = From.PairwiseHammy; // Value of Hamiltonian
        To.OldHammy = From.OldHammy; //  Previous value of Hamiltonian
        To.Ncent = From.Ncent; //the current number of clusters
        To.ClustertoSplit = From.ClustertoSplit; //the cluster label that will be split
        To.MinimumEigenvalue = From.MinimumEigenvalue; // Minimum full eigenvalue associated with ClustertoSplit
        To.ActualCoolingFactor = From.ActualCoolingFactor; // Actual Cooling Factor
        To.Temperature = From.Temperature; //the current temperature
        To.Axset = From.Axset;
        To.SolutionSet = From.SolutionSet;
        To.IterationSetAt = From.IterationSetAt;

        int NumberClusters = From.Ncent;
        for (int ClusterIndex = 0; ClusterIndex < NumberClusters; ClusterIndex++) {
            To.C_k_[ClusterIndex] = From.C_k_[ClusterIndex];
            To.A_k_[ClusterIndex] = From.A_k_[ClusterIndex];
            To.P_k_[ClusterIndex] = From.P_k_[ClusterIndex];

            To.FreezingMeasure_k_[ClusterIndex] = From.FreezingMeasure_k_[ClusterIndex];
            To.Weight_k_[ClusterIndex] = From.Weight_k_[ClusterIndex];
        }

        // Note - parallel for
        launchHabaneroApp(() -> {
            forallChunked(0, PWCUtility.ThreadCount - 1, (threadIndex) -> {
                if (PWCUtility.bindThreads){
                    AffinitySupport.setAffinity(1L<<PWCUtility.bindThreadToCore[threadIndex]);
                }
                int indexlen = PWCUtility.PointsperThread[threadIndex];
                int beginpoint = PWCUtility.StartPointperThread[threadIndex] - PWCUtility.PointStart_Process;
                for (int ProcessPointIndex = beginpoint; ProcessPointIndex < indexlen + beginpoint; ProcessPointIndex++) {
                    for (int ClusterIndex = 0; ClusterIndex < NumberClusters; ClusterIndex++) {
                        To.Epsilonalpha_k_[ProcessPointIndex][ClusterIndex] =
                                From.Epsilonalpha_k_[ProcessPointIndex][ClusterIndex];
                        To.Old_Epsilonalpha_k_[ProcessPointIndex][ClusterIndex] =
                                From.Old_Epsilonalpha_k_[ProcessPointIndex][ClusterIndex];
                        To.Best_Epsilonalpha_k_[ProcessPointIndex][ClusterIndex] =
                                From.Best_Epsilonalpha_k_[ProcessPointIndex][ClusterIndex];
                        To.Master_Epsilonalpha_k_[ProcessPointIndex][ClusterIndex] =
                                From.Master_Epsilonalpha_k_[ProcessPointIndex][ClusterIndex];
                        To.Malpha_k_[ProcessPointIndex][ClusterIndex] = From.Malpha_k_[ProcessPointIndex][ClusterIndex];
                        To.Previous_Malpha_k_[ProcessPointIndex][ClusterIndex] =
                                From.Previous_Malpha_k_[ProcessPointIndex][ClusterIndex];
                        To.Balpha_k_[ProcessPointIndex][ClusterIndex] = From.Balpha_k_[ProcessPointIndex][ClusterIndex];
                        if (From.Axset) {
                            To.Saved_Ax[ProcessPointIndex][ClusterIndex] = From.Saved_Ax[ProcessPointIndex][ClusterIndex];
                            To.Saved_oldAx[ProcessPointIndex][ClusterIndex] =
                                    From.Saved_oldAx[ProcessPointIndex][ClusterIndex];
                        }
                    }
                }
            });
        });
    } // End CopySolution

    public static void RemoveCluster(ClusteringSolution Changing, int RemovedIndex) {
        --Changing.Ncent;
        if (RemovedIndex == Changing.Ncent) {
            return;
        }

        for (int ClusterIndex = RemovedIndex; ClusterIndex < Changing.Ncent; ClusterIndex++) {
            Changing.C_k_[ClusterIndex] = Changing.C_k_[ClusterIndex + 1];
            Changing.A_k_[ClusterIndex] = Changing.A_k_[ClusterIndex + 1];
            Changing.P_k_[ClusterIndex] = Changing.P_k_[ClusterIndex + 1];

            Changing.FreezingMeasure_k_[ClusterIndex] = Changing.FreezingMeasure_k_[ClusterIndex + 1];
            Changing.Weight_k_[ClusterIndex] = Changing.Weight_k_[ClusterIndex + 1];
        }

        // Note - parallel for
        launchHabaneroApp(() -> {
            forallChunked(0, PWCUtility.ThreadCount - 1, (threadIndex) -> {
                if (PWCUtility.bindThreads){
                    AffinitySupport.setAffinity(1L<<PWCUtility.bindThreadToCore[threadIndex]);
                }
                int indexlen = PWCUtility.PointsperThread[threadIndex];
                int beginpoint = PWCUtility.StartPointperThread[threadIndex] - PWCUtility.PointStart_Process;
                for (int ProcessPointIndex = beginpoint; ProcessPointIndex < indexlen + beginpoint; ProcessPointIndex++) {
                    for (int ClusterIndex = 0; ClusterIndex < Changing.Ncent; ClusterIndex++) {
                        Changing.Epsilonalpha_k_[ProcessPointIndex][ClusterIndex] =
                                Changing.Epsilonalpha_k_[ProcessPointIndex][ClusterIndex] + 1;
                        Changing.Old_Epsilonalpha_k_[ProcessPointIndex][ClusterIndex] =
                                Changing.Old_Epsilonalpha_k_[ProcessPointIndex][ClusterIndex + 1];
                        Changing.Best_Epsilonalpha_k_[ProcessPointIndex][ClusterIndex] =
                                Changing.Best_Epsilonalpha_k_[ProcessPointIndex][ClusterIndex + 1];
                        Changing.Master_Epsilonalpha_k_[ProcessPointIndex][ClusterIndex] =
                                Changing.Master_Epsilonalpha_k_[ProcessPointIndex][ClusterIndex + 1];
                        Changing.Malpha_k_[ProcessPointIndex][ClusterIndex] =
                                Changing.Malpha_k_[ProcessPointIndex][ClusterIndex + 1];
                        Changing.Previous_Malpha_k_[ProcessPointIndex][ClusterIndex] =
                                Changing.Previous_Malpha_k_[ProcessPointIndex][ClusterIndex + 1];
                        Changing.Balpha_k_[ProcessPointIndex][ClusterIndex] =
                                Changing.Balpha_k_[ProcessPointIndex][ClusterIndex + 1];
                        if (Changing.Axset) {
                            Changing.Saved_Ax[ProcessPointIndex][ClusterIndex] =
                                    Changing.Saved_Ax[ProcessPointIndex][ClusterIndex + 1];
                            Changing.Saved_oldAx[ProcessPointIndex][ClusterIndex] =
                                    Changing.Saved_oldAx[ProcessPointIndex][ClusterIndex + 1];
                        }
                    }
                }
            });
        });
    } // End RemoveCluster

    public static void SetAxinSolution(ClusteringSolution Solution) {
        Solution.Axset = true;
        int NumberClusters = Solution.Ncent;

        // Note - parallel for
        launchHabaneroApp(() -> {
            forallChunked(0, PWCUtility.ThreadCount - 1, (threadIndex) -> {
                if (PWCUtility.bindThreads){
                    AffinitySupport.setAffinity(1L<<PWCUtility.bindThreadToCore[threadIndex]);
                }
                int indexlen = PWCUtility.PointsperThread[threadIndex];
                int beginpoint = PWCUtility.StartPointperThread[threadIndex] - PWCUtility.PointStart_Process;
                for (int ProcessPointIndex = beginpoint; ProcessPointIndex < indexlen + beginpoint; ProcessPointIndex++) {
                    for (int ClusterIndex = 0; ClusterIndex < NumberClusters; ClusterIndex++) {
                        Solution.Saved_Ax[ProcessPointIndex][ClusterIndex] =
                                vectorclass.Ax[ProcessPointIndex][ClusterIndex];
                        Solution.Saved_oldAx[ProcessPointIndex][ClusterIndex] =
                                vectorclass.oldAx[ProcessPointIndex][ClusterIndex];
                    }
                }
            });
        });
    }

    public static void RestoreAxfromSolution(ClusteringSolution Solution) {
        if (!Solution.Axset) {
            PWCUtility.printAndThrowRuntimeException("Axset false but restore requested");
        }

        int NumberClusters = Solution.Ncent;

        // Note - parallel for
        try {
            forallChunked(0, PWCUtility.ThreadCount - 1, (threadIndex) -> {
                if (PWCUtility.bindThreads){
                    AffinitySupport.setAffinity(1L<<PWCUtility.bindThreadToCore[threadIndex]);
                }
                int indexlen = PWCUtility.PointsperThread[threadIndex];
                int beginpoint = PWCUtility.StartPointperThread[threadIndex] - PWCUtility.PointStart_Process;
                for (int ProcessPointIndex = beginpoint; ProcessPointIndex < indexlen + beginpoint; ProcessPointIndex++) {
                    for (int ClusterIndex = 0; ClusterIndex < NumberClusters; ClusterIndex++) {
                        vectorclass.Ax[ProcessPointIndex][ClusterIndex] =
                                Solution.Saved_Ax[ProcessPointIndex][ClusterIndex];
                        vectorclass.oldAx[ProcessPointIndex][ClusterIndex] =
                                Solution.Saved_oldAx[ProcessPointIndex][ClusterIndex];
                    }
                }
            });
        } catch (SuspendableException e) {
            PWCUtility.printAndThrowRuntimeException(e.getMessage());
        }
    }

} // End ClusteringSolution