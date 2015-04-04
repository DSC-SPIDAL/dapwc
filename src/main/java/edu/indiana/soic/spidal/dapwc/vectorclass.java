package edu.indiana.soic.spidal.dapwc;

import mpi.MPIException;
import net.openhft.affinity.AffinitySupport;

import java.util.Arrays;

import static edu.rice.hj.Module0.launchHabaneroApp;
import static edu.rice.hj.Module1.forallChunked;

//	Class vectorclass **************************************************
//	Linear Algebra including detailed computation of second derivative matrix

//	Parameters is present in all major classes and consists of get and set routines 
//	for a few overall parameters used to communicate between threads. 
//	The choice of parameters is NOT optimal

//	vectorDot calculates simple sequential dot product of two vectors given vectors and length
//	Currently not used as need parallel implementation

//	getMinimumEigenvalue finds the minimum eigenvalue of Second Derivative Matrix C using Power Method
//	The logic is similar to that calculating first and secobd derivatives 
//	It calls PairwiseThread_SecDrv twice to calculate first the maximum eigenvalue Lmax of C
//	Then the maximum eigenvalue of (Lmax Identity Matrix -C)
//	Parameters of call are not very sensible and most could be dropped

//	PairwiseThread_SecDrv does major computations in power method.
//	It calculates Matrix.current iteration vector
//	For two matrices needed by getMinimumEigenvalue
//	It is controlled by Parameters passed by CCR

// Linear Algebra

public class vectorclass
{

	public static double[] initialvector = null; // Initial vector for Power method
	public static double[][] Ax = null; // New power vector
	public static double[][] oldAx = null; // Previous Power Vector
	public static double[] Eigenvalues_Current = null; // Current eigenvalue per cluster for powermethod
	public static double[] Eigenvalues_Pass0 = null; // Maximum eigenvalue per cluster for powermethod
	public static int[] eigenconverged = new int[Program.maxNcent + Program.cachelinesize]; // Indicator of status of eigenvalue determination per cluster

	public static boolean MandBset; // Set true if external Malpha_k_ and Balpha_k_ have been saved
	public static MPISecPacket[] MandBRepository;
	public static int MaxlengthMandB;

	// Find the minimum eigenvalue of second derivative matrix -- called from shouldSplit
	public final void getEigenvalues(int Methodology, double[][] localMalpha_k_, double[] localA_k_, double[][] localBalpha_k_, double[] localC_k_, int localNcent) throws MPIException {
		int cachelinesize = Program.cachelinesize;

		if (Ax == null)
		{ // Initialize arrays on first call

			Ax = new double[PWCUtility.PointCount_Process][];
			oldAx = new double[PWCUtility.PointCount_Process][];
			for (int PointIndex = 0; PointIndex < PWCUtility.PointCount_Process; PointIndex++)
			{
				Ax[PointIndex] = new double[Program.maxNcent + cachelinesize];
				oldAx[PointIndex] = new double[Program.maxNcent + cachelinesize];
			}

			Eigenvalues_Current = new double[Program.maxNcent + cachelinesize];
			Eigenvalues_Pass0 = new double[Program.maxNcent + cachelinesize];

			initialvector = new double[PWCUtility.PointCount_Process];
			double fudge = 1.0 / PWCUtility.PointCount_Global;

			GlobalReductions.FindDoubleSum Find_initnorm = new GlobalReductions.FindDoubleSum(PWCUtility.ThreadCount);

            // Note - parallel for
            launchHabaneroApp(() -> {
                forallChunked(0, PWCUtility.ThreadCount - 1, (threadIndex) ->
                        {
                            if (PWCUtility.bindThreads){
                                AffinitySupport.setAffinity(1L << PWCUtility.bindThreadToCore[threadIndex]);
                            }
                            int indexlen = PWCUtility.PointsperThread[threadIndex];
                            int beginpoint = PWCUtility.StartPointperThread[threadIndex] - PWCUtility.PointStart_Process;
                            for (int ProcessPointIndex = beginpoint; ProcessPointIndex < indexlen + beginpoint; ProcessPointIndex++) {
                                initialvector[ProcessPointIndex] = 1.0 + 0.5 * (ProcessPointIndex + PWCUtility.PointStart_Process) * fudge;
                                Find_initnorm.addapoint(threadIndex,
                                        initialvector[ProcessPointIndex] * initialvector[ProcessPointIndex]);
                            }
                        }
                );
            });

            Find_initnorm.sumoverthreadsandmpi();
			final double initnorm_global = 1.0 / Math.sqrt(Find_initnorm.Total);

            // Note - parallel for
            launchHabaneroApp(() -> {
                forallChunked(0, PWCUtility.ThreadCount - 1, (threadIndex) ->
                        {
                            if (PWCUtility.bindThreads){
                                AffinitySupport.setAffinity(1L<<PWCUtility.bindThreadToCore[threadIndex]);
                            }
                            //	Start Delegate Code normalizing initialvector
                            int indexlen = PWCUtility.PointsperThread[threadIndex];
                            int beginpoint = PWCUtility.StartPointperThread[threadIndex] - PWCUtility.PointStart_Process;
                            for (int ProcessPointIndex = beginpoint; ProcessPointIndex < indexlen + beginpoint; ProcessPointIndex++) {
                                initialvector[ProcessPointIndex] *= initnorm_global;
                            }
                        }
                );
            });

        } // End Initialize arrays on first call

		MaxlengthMandB = localNcent * PWCUtility.PointCount_Largest;
		MandBRepository = new MPISecPacket[PWCUtility.MPI_Size];
		for (int localrank = 0; localrank < PWCUtility.MPI_Size; localrank++)
		{
			MandBRepository[localrank] = new MPISecPacket(MaxlengthMandB);
		}
		MandBset = false;

		Dist.PlaceforEigsforMultipleCluster = -1;
		// Set in Eigenvalues_Pass0 the maximum eigenvalue per center of second derivative matrix   on first pass      
		if (Methodology > 2)
		{
			for (int ClusterIndex = 0; ClusterIndex < localNcent; ClusterIndex++)
			{
				if (Dist.ClusterSelected[ClusterIndex] == 1)
				{
					Dist.PlaceforEigsforMultipleCluster = ClusterIndex;
					break;
				}
			}
		}

		//	Invoke PairwiseThread_SecDrv with Flag = 0
		PairwiseThread_SecDrv(Methodology, 0, localMalpha_k_, localA_k_, localBalpha_k_, localC_k_, localNcent);

		// Set in Eigenvalues_Pass0 the maximum eigenvalue per center of second derivative matrix   on first pass      
		int[] savestatus = new int[Dist.RunningPWC.Ncent];
		for (int ClusterIndex = 0; ClusterIndex < localNcent; ClusterIndex++)
		{
			if ((Methodology > 2))
			{
				if (ClusterIndex != Dist.PlaceforEigsforMultipleCluster)
				{
					continue;
				}
				if ((eigenconverged[ClusterIndex] <= 0) || (Eigenvalues_Current[ClusterIndex] <= 0.0))
				{
					Eigenvalues_Pass0[ClusterIndex] = Eigenvalues_Current[ClusterIndex];
					eigenconverged[ClusterIndex] = -1;
					Eigenvalues_Current[ClusterIndex] = -1.0 + Eigenvalues_Pass0[ClusterIndex];
					++Program.CountTotalPowerErrors;
					MandBset = false;
					return;
				}
			}
			Eigenvalues_Pass0[ClusterIndex] = Eigenvalues_Current[ClusterIndex];
			savestatus[ClusterIndex] = eigenconverged[ClusterIndex];
			if ((eigenconverged[ClusterIndex] <= 0) || (Eigenvalues_Current[ClusterIndex] <= 0.0))
			{
				++Program.CountTotalPowerErrors;
			}
		}

		//	Invoke PairwiseThread_SecDrv with Flag = 1
		PairwiseThread_SecDrv(Methodology, 1, localMalpha_k_, localA_k_, localBalpha_k_, localC_k_, localNcent);
		for (int ClusterIndex = 0; ClusterIndex < localNcent; ClusterIndex++)
		{
			if ((Methodology > 2) && (ClusterIndex != Dist.PlaceforEigsforMultipleCluster))
			{
				continue;
			}
			if (savestatus[ClusterIndex] > 0 && eigenconverged[ClusterIndex] > 0)
			{
				eigenconverged[ClusterIndex] += savestatus[ClusterIndex];
			}
			if (savestatus[ClusterIndex] > 0 && eigenconverged[ClusterIndex] < 0)
			{
				++Program.CountTotalPowerErrors;
			}
		}
		MandBset = false;

	} // End getEigenvalues

	// Compute the crrelation matrix compared to PairwiseThread. Rest is same as PairwiseThread
	// PassIndicator = 0 Find MAXIMUM eigenvalue  LMAX[kk] in EACH Cluster sector kk
	// PassIndicator = 1 Find MINIMUM eigenvalue in each cluster sector by finding Maximum eigenvalue of (LMAX[kk]-Second Derivative Matrix)

	// Methodology = 1 -- Calculate Power Vector for all Clusters (one per Cluster) from Original Second Derivative Matrix
	// Methodology = 2 -- Calculate Power Vector for all Clusters (one per Cluster) from  modification of second derivative matrix to reflect replication 
	//                    for cases like ContinuousClustering=true that cluster can be divided into two with M halved for each case
	// Methodology = 3 -- Calculate Power Vector for ONE Clusters "Clustertoprocess" (Specified in ClusterSelected) from proper modification of second derivative matrix to reflect replication
	//                          Assuming that Clustertoprocess has its M value doubled in previous EM Iteration
	//                  2 and 3 are same except 3 assumes explicit facor of 2 weight in Z in calculation; Methodology=2 assumes no factor of 2
	// Methodology = 4 -- Calculate ONE Power Vector for set of Clusters (Could be just one but specified in ClusterSelected and usually at least two) from proper modification of second derivative matrix to reflect replication
	//                          Assuming NO doubled M values
	// Set Dist.ClusterSelected for all cases ( all zeros if Methodology 1,2)
	// Set ClustertoSplitin Methodology 3 CONSISTENT with Dist.ClusterSelected[ClustertoSplit]=1 being only ONE value
	// Set Dist.PlaceforEigsforMultipleCluster in Methodology 4 to indicate where results stored

	private void PairwiseThread_SecDrv(int Methodology, int PassIndicator, double[][] localMalpha_k_, double[] localA_k_, double[][] localBalpha_k_, double[] localC_k_, int localNcent) throws MPIException {
		double T = Dist.RunningPWC.Temperature;
        int[] UsethisCluster = new int[localNcent];
		double[] AxPattern = new double[localNcent];
		int REMOVEZEROEigs = 0;
		if (Methodology == 4)
		{
			REMOVEZEROEigs = 1;
		}

		for (int ClusterIndex = 0; ClusterIndex < localNcent; ClusterIndex++)
		{
			AxPattern[ClusterIndex] = 1.0;
			UsethisCluster[ClusterIndex] = Dist.ClusterSelected[ClusterIndex];
			if ((Methodology == 3) && (UsethisCluster[ClusterIndex] == 1))
			{
				eigenconverged[ClusterIndex] = -1;
			}
			if (Methodology < 3)
			{
				eigenconverged[ClusterIndex] = -1;
				if ((PassIndicator == 1) && (vectorclass.Eigenvalues_Pass0[ClusterIndex] <= 0.0))
				{
					UsethisCluster[ClusterIndex] = 0;
				}
			}
			if ((Methodology == 4) && (UsethisCluster[ClusterIndex] == 0))
			{
				REMOVEZEROEigs = 0;
			}
		}
		if (Methodology == 4)
		{
			eigenconverged[Dist.PlaceforEigsforMultipleCluster] = -1;
			int usedclusters = 0;
			for (int ClusterIndex = 0; ClusterIndex < localNcent; ClusterIndex++)
			{
				if (UsethisCluster[ClusterIndex] > 0)
				{
					++usedclusters;
				}
			}
			double vectormult = 1.0;
			double oddmultiplier = 1.0;
			if (usedclusters > 1)
			{
				int halfused = usedclusters / 2;
				if (halfused * 2 == usedclusters)
				{
					vectormult = 1.0 / Math.sqrt((double)usedclusters);
				}
				else
				{
					oddmultiplier = 1.0 + 1.0 / halfused;
					vectormult = 1.0 / Math.sqrt(2.0 * halfused + 3.0 + 1.0 / halfused);
				}
			}
			for (int ClusterIndex = 0; ClusterIndex < localNcent; ClusterIndex++)
			{
				if (Dist.ClusterSelected[ClusterIndex] == 0)
				{
					AxPattern[ClusterIndex] = 0.0;
				}
				else
				{
					AxPattern[ClusterIndex] = vectormult;
					vectormult = -vectormult;
					if (vectormult > 0.0)
					{
						AxPattern[ClusterIndex] *= oddmultiplier;
					}
				}
			}
		} // End Methodology = 4

        // Note - parallel for
        launchHabaneroApp(() -> {
            forallChunked(0, PWCUtility.ThreadCount - 1, (threadIndex) ->
                    {
                        if (PWCUtility.bindThreads){
                            AffinitySupport.setAffinity(1L<<PWCUtility.bindThreadToCore[threadIndex]);
                        }
                        //	Start Code initializing power vectors Ax oldAx
                        int indexlen = PWCUtility.PointsperThread[threadIndex];
                        int beginpoint = PWCUtility.StartPointperThread[threadIndex] - PWCUtility.PointStart_Process;
                        for (int ProcessPointIndex = beginpoint; ProcessPointIndex < indexlen + beginpoint; ProcessPointIndex++) {
                            for (int ClusterIndex = 0; ClusterIndex < localNcent; ClusterIndex++) {
                                if (Dist.ClusterSelected[ClusterIndex] == 0) {
                                    continue;
                                }
                                Ax[ProcessPointIndex][ClusterIndex] = initialvector[ProcessPointIndex] * AxPattern[ClusterIndex];
                                oldAx[ProcessPointIndex][ClusterIndex] = Ax[ProcessPointIndex][ClusterIndex];
                            }
                        }
                    }
            );
        });


        double Mfudge = 1.0;
		if (Methodology == 2)
		{
			Mfudge = 0.5;
		}
		final boolean Subtract_twiceDuplicatedCenter = (Methodology != 1) && (Methodology != 4);

		MPISecPacket fromafarMandB = new MPISecPacket(MaxlengthMandB);
		MPISecPacket toafarMandB = new MPISecPacket(MaxlengthMandB);
		MPISecPacket myownMandB = new MPISecPacket(MaxlengthMandB);
		double[] fromafarAxarray = null;
		double[] toafarAxarray = null;
		double[] myownAxarray = null;

		int oldMaxLength2 = Integer.MIN_VALUE;

		for (int NumPowerIterations = 0; NumPowerIterations < Program.PowerIterationLimit; NumPowerIterations++)
		{ // Loop over Power Iterations
			//	Check if any centers are still not converged
			int NumberofAVectorsUsed = 0;
			for (int ClusterIndex = 0; ClusterIndex < localNcent; ClusterIndex++)
			{
				if (UsethisCluster[ClusterIndex] == 1)
				{
					++NumberofAVectorsUsed;
				}
			}
			if (NumberofAVectorsUsed == 0)
			{
				break;
			}
			if ((Methodology == 4) && (vectorclass.eigenconverged[Dist.PlaceforEigsforMultipleCluster] > 0))
			{
				break;
			}

			// Increment Global Count of Iterations
			++Program.CountTotalPowerIterations;

			// Loop over MPI calls
			int sendtag = 0;
			int receivetag = 0;

			if (!MandBset)
			{
				fromafarMandB.Clear();
			}

			toafarMandB.Clear();
			myownMandB.Clear();

			int Maxlength2 = NumberofAVectorsUsed * PWCUtility.PointCount_Largest;

			if (oldMaxLength2 != Maxlength2)
			{
				fromafarAxarray = new double[Maxlength2];
				toafarAxarray = new double[Maxlength2];
				myownAxarray = new double[Maxlength2];
				oldMaxLength2 = Maxlength2;
			}
			else
			{
                Arrays.fill(fromafarAxarray, 0, Maxlength2, 0.0);
				Arrays.fill(toafarAxarray, 0, Maxlength2, 0.0);
				Arrays.fill(myownAxarray, 0, Maxlength2, 0.0);
			}

			myownMandB.setFirstPoint(PWCUtility.PointStart_Process);
			myownMandB.setNumberOfPoints(PWCUtility.PointCount_Process);
			if (PWCUtility.PointCount_Process < PWCUtility.PointCount_Largest)
			{
				int Acount = 0;
				for (int countC = 0; countC < localNcent; countC++)
				{
					int settozero = PWCUtility.PointCount_Process * localNcent + countC;
					myownMandB.setMArrayDoubleAt(settozero, 0.0);
					myownMandB.setBArrayDoubleAt(settozero, 0.0);
					if (UsethisCluster[countC] == 0)
					{
						continue;
					}
					myownAxarray[PWCUtility.PointCount_Process * NumberofAVectorsUsed + Acount] = 0.0;
					++Acount;
				}
			}
			int fromprocess = PWCUtility.MPI_Rank - 1;
			if (fromprocess < 0)
			{
				fromprocess = PWCUtility.MPI_Size - 1;
			}
			int toprocess = PWCUtility.MPI_Rank + 1;
			if (toprocess > PWCUtility.MPI_Size - 1)
			{
				toprocess = 0;
			}

			//	First communicationloop is local; then we have MPI_Size transfers of data in  a ring through processes             
			for (int MPICommunicationSteps = 0; MPICommunicationSteps < PWCUtility.MPI_Size; MPICommunicationSteps++)
			{
				if (MPICommunicationSteps == 1)
				{
					toafarMandB.setFirstPoint(PWCUtility.PointStart_Process);
					toafarMandB.setNumberOfPoints(PWCUtility.PointCount_Process);
					if (PWCUtility.PointCount_Process < PWCUtility.PointCount_Largest)
					{
						int Acount = 0;
						for (int ClusterIndex = 0; ClusterIndex < localNcent; ClusterIndex++)
						{
							int bigindex = PWCUtility.PointCount_Process * localNcent + ClusterIndex;
							toafarMandB.setMArrayDoubleAt(bigindex, myownMandB.getMArrayDoubleAt(bigindex));
							toafarMandB.setBArrayDoubleAt(bigindex,myownMandB.getBArrayDoubleAt(bigindex));
							if (UsethisCluster[ClusterIndex] == 0)
							{
								continue;
							}
							toafarAxarray[PWCUtility.PointCount_Process * NumberofAVectorsUsed + Acount] = myownAxarray[PWCUtility.PointCount_Process * NumberofAVectorsUsed + Acount];
							++Acount;
						}
					}
				}
				if (MPICommunicationSteps > 1)
				{
					toafarMandB.setFirstPoint(fromafarMandB.getFirstPoint());
					toafarMandB.setNumberOfPoints(fromafarMandB.getNumberOfPoints());
					if (PWCUtility.PointCount_Process < PWCUtility.PointCount_Largest)
					{
						int Acount = 0;
						for (int ClusterIndex = 0; ClusterIndex < localNcent; ClusterIndex++)
						{
							int bigindex = PWCUtility.PointCount_Process * localNcent + ClusterIndex;
							toafarMandB.setMArrayDoubleAt(bigindex,fromafarMandB.getMArrayDoubleAt(bigindex));
							toafarMandB.setBArrayDoubleAt(bigindex,fromafarMandB.getBArrayDoubleAt(bigindex));
							if (UsethisCluster[ClusterIndex] == 0)
							{
								continue;
							}
							toafarAxarray[PWCUtility.PointCount_Process * NumberofAVectorsUsed + Acount] = fromafarAxarray[PWCUtility.PointCount_Process * NumberofAVectorsUsed + Acount];
							++Acount;
						}
					}
				}
				if (MPICommunicationSteps > 0)
				{
					PWCUtility.StartSubTimer(PWCUtility.MPISENDRECEIVEEigenTiming);
					if (!MandBset)
					{
                        fromafarMandB = PWCUtility.mpiOps.sendReceive(toafarMandB, toprocess, sendtag, fromprocess, receivetag);
                        MPISecPacket.memberCopy(fromafarMandB, MandBRepository[MPICommunicationSteps]);
					}
					else
					{
						fromafarMandB = MandBRepository[MPICommunicationSteps];
					}
                    fromafarAxarray = PWCUtility.mpiOps.sendReceive(toafarAxarray, toprocess, sendtag, fromprocess, receivetag);
					PWCUtility.StopSubTimer(PWCUtility.MPISENDRECEIVEEigenTiming);
				}

				// Communication finished -- now update A vector
                final int MPICommunicationStepsLoopVar = MPICommunicationSteps;
                final int NumberofAVectorsUsedLoopVar = NumberofAVectorsUsed;
                final double [] toafarAxarrayLoopVar = toafarAxarray;
                final double [] fromafarAxarrayLoopVar = fromafarAxarray;
                final double [] myownAxarrayLoopVar = myownAxarray;
                final double MfudgeLoopVar = Mfudge;
                MPISecPacket fromafarMandBLoopVar = fromafarMandB;
                // Note - parallel for
                launchHabaneroApp(() -> {
                    forallChunked(0, PWCUtility.ThreadCount - 1, (threadIndex) ->
                                    //	Start Code calculating power vectors Ax oldAx
                                    // DiagonalTerm[ClusterIndex] = tmp / T;
                                    //  FirstTerm = -MalphaMu * MbetaLambda / T;
                                    // double MatrixElement = FirstTerm + ((MalphaMu * MbetaLambda * Mfudge * Mfudge) / (T * T)) * tmp;
                            {
                                if (PWCUtility.bindThreads){
                                    AffinitySupport.setAffinity(1L<<PWCUtility.bindThreadToCore[threadIndex]);
                                }
                                double[] DiagonalTerm = new double[localNcent];
                                int indexlen = PWCUtility.PointsperThread[threadIndex];
                                int beginpoint = PWCUtility.StartPointperThread[threadIndex] - PWCUtility.PointStart_Process;
                                for (int ProcessPointIndex = beginpoint; ProcessPointIndex < indexlen + beginpoint; ProcessPointIndex++) {
                                    int betatotal, betastart;
                                    if (MPICommunicationStepsLoopVar == 0) {
                                        betatotal = PWCUtility.PointCount_Process;
                                        betastart = PWCUtility.PointStart_Process;
                                        int Acount = 0;
                                        for (int ClusterIndex = 0; ClusterIndex < localNcent; ClusterIndex++) {
                                            int bigindex = ProcessPointIndex * localNcent + ClusterIndex;
                                            myownMandB.setMArrayDoubleAt(bigindex,
                                                    localMalpha_k_[ProcessPointIndex][ClusterIndex]);
                                            myownMandB.setBArrayDoubleAt(bigindex,
                                                    localBalpha_k_[ProcessPointIndex][ClusterIndex]);
                                            toafarMandB.setMArrayDoubleAt(bigindex,
                                                    localMalpha_k_[ProcessPointIndex][ClusterIndex]);
                                            toafarMandB.setBArrayDoubleAt(bigindex,
                                                    localBalpha_k_[ProcessPointIndex][ClusterIndex]);
                                            if (UsethisCluster[ClusterIndex] == 0) {
                                                continue;
                                            }
                                            int Aindex = ProcessPointIndex * NumberofAVectorsUsedLoopVar + Acount;
                                            ++Acount;
                                            toafarAxarrayLoopVar[Aindex] = oldAx[ProcessPointIndex][ClusterIndex];
                                            myownAxarrayLoopVar[Aindex] = oldAx[ProcessPointIndex][ClusterIndex];
                                            vectorclass.Ax[ProcessPointIndex][ClusterIndex] = 0.0;
                                        }
                                    } else {
                                        betatotal = fromafarMandBLoopVar.getNumberOfPoints();
                                        betastart = fromafarMandBLoopVar.getFirstPoint();
                                        if (MPICommunicationStepsLoopVar != (PWCUtility.MPI_Size - 1)) {
                                            int Acount = 0;
                                            for (int ClusterIndex = 0; ClusterIndex < localNcent; ClusterIndex++) {
                                                int bigindex = ProcessPointIndex * localNcent + ClusterIndex;
                                                toafarMandB.setMArrayDoubleAt(bigindex,
                                                        fromafarMandBLoopVar.getMArrayDoubleAt(bigindex));
                                                toafarMandB.setBArrayDoubleAt(bigindex,
                                                        fromafarMandBLoopVar.getBArrayDoubleAt(bigindex));
                                                if (UsethisCluster[ClusterIndex] == 0) {
                                                    continue;
                                                }
                                                int Aindex = ProcessPointIndex * NumberofAVectorsUsedLoopVar + Acount;
                                                ++Acount;
                                                toafarAxarrayLoopVar[Aindex] = fromafarAxarrayLoopVar[Aindex];
                                            }
                                        }
                                    }
                                    for (int ClusterIndex = 0; ClusterIndex < localNcent; ClusterIndex++) {
                                        if (UsethisCluster[ClusterIndex] == 0) {
                                            continue;
                                        }
                                        double tmp = localMalpha_k_[ProcessPointIndex][ClusterIndex] * MfudgeLoopVar;
                                        if (!Subtract_twiceDuplicatedCenter) {
                                            tmp -= tmp * tmp;
                                        }
                                        DiagonalTerm[ClusterIndex] = tmp;
                                    }
                                    for (int betalocal = 0; betalocal < betatotal; betalocal++) {
                                        int betafull = betastart + betalocal;
                                        double dijforthiscase = PWCUtility.PointDistances.getDistance(
                                                ProcessPointIndex + PWCUtility.PointStart_Process, betafull);
                                        boolean PointIndicesEqual = (betafull == (ProcessPointIndex + PWCUtility.PointStart_Process));
                                        for (int CenterVectorMu = 0; CenterVectorMu < localNcent; CenterVectorMu++) {
                                            double MalphaMu = localMalpha_k_[ProcessPointIndex][CenterVectorMu];
                                            if (UsethisCluster[CenterVectorMu] == 0) {
                                                continue;
                                            }
                                            for (int CenterVectorLambda = 0; CenterVectorLambda < localNcent; CenterVectorLambda++) {
                                                int ALambdacount = 0;
                                                if (UsethisCluster[CenterVectorLambda] == 0) {
                                                    continue;
                                                }
                                                int ALambdaindex = betalocal * NumberofAVectorsUsedLoopVar + ALambdacount;
                                                ALambdacount++;
                                                if ((Methodology < 4) && (CenterVectorLambda != CenterVectorMu)) {
                                                    continue;
                                                }
                                                double MbetaLambda, AxbetaLambda;
                                                if (MPICommunicationStepsLoopVar == 0) {
                                                    MbetaLambda = localMalpha_k_[betalocal][CenterVectorLambda];
                                                    AxbetaLambda = oldAx[betalocal][CenterVectorLambda];
                                                } else {
                                                    MbetaLambda = fromafarMandBLoopVar.getMArrayDoubleAt(
                                                            betalocal * localNcent + CenterVectorLambda);
                                                    AxbetaLambda = fromafarAxarrayLoopVar[ALambdaindex];
                                                }
                                                double FirstTerm = 0;
                                                if (PointIndicesEqual && (CenterVectorMu == CenterVectorLambda)) {
                                                    FirstTerm = DiagonalTerm[CenterVectorMu];
                                                }
                                                if ((Methodology == 4) && PointIndicesEqual && (CenterVectorMu != CenterVectorLambda)) {
                                                    FirstTerm = -MalphaMu * MbetaLambda;
                                                }
                                                double tmp = 0.0;
                                                for (int CentersSummed = 0; CentersSummed < localNcent; CentersSummed++) {
                                                    double BbetaCenterSummed, MbetaCenterSummed;
                                                    if (MPICommunicationStepsLoopVar == 0) {
                                                        MbetaCenterSummed = localMalpha_k_[betalocal][CentersSummed];
                                                        BbetaCenterSummed = localBalpha_k_[betalocal][CentersSummed];
                                                    } else {
                                                        MbetaCenterSummed = fromafarMandBLoopVar.getMArrayDoubleAt(
                                                                betalocal * localNcent + CentersSummed);
                                                        BbetaCenterSummed = fromafarMandBLoopVar.getBArrayDoubleAt(
                                                                betalocal * localNcent + CentersSummed);
                                                    }
                                                    double MMmultiplier;
                                                    double Mfudge_CenterCalculated = 1.0;
                                                    if (Methodology <= 3) {
                                                        if (CentersSummed == CenterVectorMu) {
                                                            Mfudge_CenterCalculated = MfudgeLoopVar;
                                                            if (Subtract_twiceDuplicatedCenter) {
                                                                MMmultiplier = 1.0;
                                                            } else {
                                                                MMmultiplier = (localMalpha_k_[ProcessPointIndex][CentersSummed] * Mfudge_CenterCalculated - 1.0) * (MbetaCenterSummed * Mfudge_CenterCalculated - 1.0);
                                                            }
                                                        } else {
                                                            if (Subtract_twiceDuplicatedCenter) {
                                                                MMmultiplier = 0.0;
                                                                continue;
                                                            } else {
                                                                MMmultiplier = localMalpha_k_[ProcessPointIndex][CentersSummed] * MbetaCenterSummed;
                                                            }
                                                        }
                                                    } else {
                                                        double DeltaSummedMu = 0.0;
                                                        double DeltaSummedLambda = 0.0;
                                                        if (CentersSummed == CenterVectorMu) {
                                                            DeltaSummedMu = 1.0;
                                                        }
                                                        if (CentersSummed == CenterVectorLambda) {
                                                            DeltaSummedLambda = 1.0;
                                                        }
                                                        MMmultiplier = (localMalpha_k_[ProcessPointIndex][CentersSummed] - DeltaSummedMu) * (MbetaCenterSummed - DeltaSummedLambda);
                                                    }
                                                    tmp = tmp + (-2.0 * localA_k_[CentersSummed] - BbetaCenterSummed - localBalpha_k_[ProcessPointIndex][CentersSummed] + dijforthiscase) * MMmultiplier / (localC_k_[CentersSummed] * Mfudge_CenterCalculated);
                                                }
                                                double MatrixElement = FirstTerm + ((MalphaMu * MbetaLambda * MfudgeLoopVar * MfudgeLoopVar) / T) * tmp;
                                                if (PassIndicator == 1) {
                                                    if (PointIndicesEqual && (CenterVectorMu == CenterVectorLambda)) {
                                                        int eigenposition = CenterVectorMu;
                                                        if (Methodology == 4) {
                                                            eigenposition = Dist.PlaceforEigsforMultipleCluster;
                                                        }
                                                        MatrixElement = Eigenvalues_Pass0[eigenposition] - MatrixElement;
                                                    } else {
                                                        MatrixElement = -MatrixElement;
                                                    }
                                                }
                                                vectorclass.Ax[ProcessPointIndex][CenterVectorMu] += MatrixElement * AxbetaLambda;
                                            }
                                        }
                                    }
                                }
                            }
                    );
                });
            } // End communicationloop

			MandBset = true;
			GlobalReductions.FindVectorDoubleSum Find_sum_t0 = new GlobalReductions.FindVectorDoubleSum(PWCUtility.ThreadCount, localNcent);
			GlobalReductions.FindVectorDoubleSum Find_sum_t1 = new GlobalReductions.FindVectorDoubleSum(PWCUtility.ThreadCount, localNcent);
			GlobalReductions.FindVectorDoubleSum Find_sum_t2 = new GlobalReductions.FindVectorDoubleSum(PWCUtility.ThreadCount, localNcent);

			//	sum over threads to get Power Eigenvalues
            // Note - parallel for
            final int REMOVEZEROEigsLoopVar = REMOVEZEROEigs;
            launchHabaneroApp(() -> {
                forallChunked(0, PWCUtility.ThreadCount - 1, (threadIndex) ->
                        {
                            if (PWCUtility.bindThreads){
                                AffinitySupport.setAffinity(1L<<PWCUtility.bindThreadToCore[threadIndex]);
                            }
                            double[] Accum_t0 = new double[localNcent];
                            double[] Accum_t1 = new double[localNcent];
                            double[] Accum_t2 = new double[localNcent];
                            for (int ClusterIndex = 0; ClusterIndex < localNcent; ClusterIndex++) {
                                Accum_t0[ClusterIndex] = 0.0;
                                Accum_t1[ClusterIndex] = 0.0;
                                Accum_t2[ClusterIndex] = 0.0;
                            }
                            int indexlen = PWCUtility.PointsperThread[threadIndex];
                            int beginpoint = PWCUtility.StartPointperThread[threadIndex] - PWCUtility.PointStart_Process;
                            for (int ProcessPointIndex = beginpoint; ProcessPointIndex < indexlen + beginpoint; ProcessPointIndex++) {
                                double SubtractAverage = 0.0;
                                if (REMOVEZEROEigsLoopVar == 1) {
                                    for (int ClusterIndex = 0; ClusterIndex < localNcent; ClusterIndex++) {
                                        SubtractAverage += Ax[ProcessPointIndex][ClusterIndex];
                                    }
                                    SubtractAverage = SubtractAverage / localNcent;
                                }
                                for (int ClusterIndex = 0; ClusterIndex < localNcent; ClusterIndex++) {
                                    if (UsethisCluster[ClusterIndex] == 0) {
                                        continue;
                                    }
                                    double tmp = Ax[ProcessPointIndex][ClusterIndex] - SubtractAverage;
                                    Ax[ProcessPointIndex][ClusterIndex] = tmp;
                                    double oldtmp = oldAx[ProcessPointIndex][ClusterIndex];
                                    Accum_t0[ClusterIndex] += tmp * oldtmp;
                                    Accum_t1[ClusterIndex] += oldtmp * oldtmp;
                                    Accum_t2[ClusterIndex] += tmp * tmp;
                                }
                            }
                            Find_sum_t0.addAPoint(threadIndex, Accum_t0);
                            Find_sum_t1.addAPoint(threadIndex, Accum_t1);
                            Find_sum_t2.addAPoint(threadIndex, Accum_t2);
                        }
                );
            });

            Find_sum_t0.sumOverThreadsAndMPI();
			Find_sum_t1.sumOverThreadsAndMPI();
			Find_sum_t2.sumOverThreadsAndMPI();

			double[] AxNormfactor = new double[localNcent];
			int[] ActiveAxvalue = new int[localNcent];
			for (int ClusterIndex = 0; ClusterIndex < localNcent; ClusterIndex++)
			{ // ClusterIndex loop checking eigenvalue and status determination for each cluster

				ActiveAxvalue[ClusterIndex] = UsethisCluster[ClusterIndex];
				if ((Methodology < 4) && (UsethisCluster[ClusterIndex] == 0))
				{
					continue;
				}
				if ((Methodology == 4) && (ClusterIndex != Dist.PlaceforEigsforMultipleCluster))
				{
					continue;
				}
				double[] globalsum_t012 = new double[3];
				for (int eiglist = 0; eiglist < 3; eiglist++)
				{
					globalsum_t012[eiglist] = 0.0;
				}

				for (int ClusterSummed = 0; ClusterSummed < localNcent; ClusterSummed++)
				{ // Cluster Summation for Methodology 4 (CenterSummed == ClusterIndex if Methodology < 4)
					if ((Methodology < 4) && (ClusterSummed != ClusterIndex))
					{
						continue;
					}

					if (Methodology < 4)
					{
						globalsum_t012[0] = Find_sum_t0.TotalVectorSum[ClusterSummed];
						globalsum_t012[1] = Find_sum_t1.TotalVectorSum[ClusterSummed];
						globalsum_t012[2] = Find_sum_t2.TotalVectorSum[ClusterSummed];
					} //  End Cluster Summation for Methodology < 4
					else
					{
						double[] globalsum_t012_temp = new double[3];
						globalsum_t012_temp[0] = Find_sum_t0.TotalVectorSum[ClusterSummed];
						globalsum_t012_temp[1] = Find_sum_t1.TotalVectorSum[ClusterSummed];
						globalsum_t012_temp[2] = Find_sum_t2.TotalVectorSum[ClusterSummed];
						for (int eiglist = 0; eiglist < 3; eiglist++)
						{
							globalsum_t012[eiglist] += globalsum_t012_temp[eiglist];
						}
					} //  End Cluster Summation for Methodology 4
				}

				double eigenvalue = globalsum_t012[0] / globalsum_t012[1];
				AxNormfactor[ClusterIndex] = 1.0 / Math.sqrt(globalsum_t012[2]);

				//	Check if converged
				//	Do this in one process ONLY 
				if ((NumPowerIterations > 10) && (eigenvalue > 0.0))
				{ // Arbitrary choice for Number of Power Iterations Cut
					int somethingtodo = 0;
					if (PWCUtility.MPI_Rank == 0)
					{ // Decisions can only be made in one process
						if (Math.abs(eigenvalue - Eigenvalues_Current[ClusterIndex]) > eigenvalue * Program.eigenvaluechange)
						{
							++somethingtodo;
						}
						double delta = globalsum_t012[2] - 2.0 * globalsum_t012[0] * eigenvalue + globalsum_t012[1] * eigenvalue * eigenvalue; // (Ax- Eigenvalue*Axold)**2
						if (Math.abs(delta) > eigenvalue * eigenvalue * Program.eigenvectorchange)
						{
							++somethingtodo;
						}
					} // End Test on Convergence
					somethingtodo = PWCUtility.synchronizeMPIVariable(somethingtodo);

					if (somethingtodo == 0)
					{
						vectorclass.eigenconverged[ClusterIndex] = 1 + NumPowerIterations;
						UsethisCluster[ClusterIndex] = 0;
					}
				}
				Eigenvalues_Current[ClusterIndex] = eigenvalue;

			} // End ClusterIndex loop checking eigenvalue status determination for each cluster

            // Note - parallel for
            launchHabaneroApp(() -> {
                forallChunked(0, PWCUtility.ThreadCount - 1, (threadIndex) ->
                        {
                            if (PWCUtility.bindThreads){
                                AffinitySupport.setAffinity(1L<<PWCUtility.bindThreadToCore[threadIndex]);
                            }
                            int indexlen = PWCUtility.PointsperThread[threadIndex];
                            int beginpoint = PWCUtility.StartPointperThread[threadIndex] - PWCUtility.PointStart_Process;
                            for (int PointIndex = beginpoint; PointIndex < indexlen + beginpoint; PointIndex++) {
                                for (int ClusterIndex = 0; ClusterIndex < localNcent; ClusterIndex++) {
                                    double Normfactor;
                                    if (Methodology == 4) {
                                        Normfactor = AxNormfactor[Dist.PlaceforEigsforMultipleCluster];
                                    } else {
                                        Normfactor = AxNormfactor[ClusterIndex];
                                    }
                                    if (ActiveAxvalue[ClusterIndex] == 0) {
                                        continue;
                                    }
                                    vectorclass.oldAx[PointIndex][ClusterIndex] = vectorclass.Ax[PointIndex][ClusterIndex] * Normfactor;
                                }
                            }
                        }
                );
            });
        } // End NumPowerIterations loop over Power Iteration Method

	} // End PairwiseThread_SecDrv

} // End vectorclass
 // End namespace cluster
