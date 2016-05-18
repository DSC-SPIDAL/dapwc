package edu.indiana.soic.spidal.dapwc;

import edu.rice.hj.api.SuspendableException;
import mpi.MPIException;

import java.util.Arrays;

import static edu.rice.hj.Module1.forallChunked;


public class Dist
{
	public static double[] Diff_Epsilon_k_ = null;

	public static double[][] DistanceCorrelation; // Distance Correlation
	public static int[] ClusterSelected; // NONZERO values are used to select Eigenvalues to Process

	public static ClusteringSolution RunningPWC; // The main Solution
	public static ClusteringSolution SavedPWC; // Saved Solution
	public static ClusteringSolution BestPWC; // Current Best Solution

	public static int PlaceforEigsforMultipleCluster = -1; // Where eigenvalues for complex situations stored
	public static int[] ListofClusterstoSplit; // List of Clusters to Split
	public static double[] EigsofClusterstoSplit; // List of Eigenvalues of Clusters to Split
	public static int Numberthatcanbesplit = 0; // Length of entries in ListofClusterstoSplit which is at most Program.MaxNumberSplitClusters

	public static double Tmin; //minimal temperature
	public static double Tinitial; //Initial temperature
	public static int ActualMaxNcent; // Maximum reduced if small clusters
	public static int CountValidityFailures = 0; // Count failures in validity checks of final solution
	public static boolean OnLastleg = false; // If True we are on final convergence step with full number of clusters

	public static int EMIterationCount = 0; // iteration number in EM method
	public static int EMIterationStepCount = 0;
	public static int Extra_EMIterationCount = 0; // Extra iterations in EM method
	public static int ActualWaititerations = 0;
	public static int countAfterFixingClusterCount = 0; // Count loops after end on number of centers
	public static int PMLoopUsed = 0; // Number of PM Loops
	public static int IterationNumberPrinted = -1; // Iteration Number Printed

	public static int oldepsiset = 0; // oldepsiset =1 if oldepsi set
	public static int needtocalculateMalpha_k_ = -1; // =0 Already Set; = 1 Calculate from EM; = -1 Initialize

	public static boolean HitConvergenceLoopLimit = false;
	public static boolean justconverging = false; // If true just converging -- No splits
	public static boolean taskswaiting = false; // If true just converging -- No splits
	public static int SplitFailures = 0; // Counts splitting failures

	public static boolean initialized = false;
	public static boolean HammyNotSet = true; // If True Hamiltonian Set


	/*
	 * The basic function in the program. It implements all steps of the pairwiseDA algorithm.
	 *  
	 */
	public final void getDist() throws MPIException {

		//allocate memory on first and indeed only call
		if (!initialized)
		{
			initialized = true;

			ClusteringSolution.SetParameters(PWCUtility.PointCount_Process, Program.maxNcent, Program.cachelinesize);
			RunningPWC = new ClusteringSolution();
			SavedPWC = new ClusteringSolution();
			BestPWC = new ClusteringSolution();

			Program.MaxNumberSplitClusters = Math.max(Program.MaxNumberSplitClusters, 1);
			Dist.ListofClusterstoSplit = new int[Program.MaxNumberSplitClusters];
			Dist.EigsofClusterstoSplit = new double[Program.MaxNumberSplitClusters];

			int cachelinesize = Program.cachelinesize;

			ClusterSelected = new int[Program.maxNcent + cachelinesize];


		} //end Initialization

		//  Do EM calculation
		Dist.ActualMaxNcent = Program.maxNcent;
		Dist.OnLastleg = false;

		// Set Initial Dist.Temperature
		initializeTemperature();
		Dist.RunningPWC.Temperature = Dist.Tinitial;
		Dist.Tmin = Dist.RunningPWC.Temperature / 100000.0;
		Dist.RunningPWC.ActualCoolingFactor = Program.InitialCoolingFactor;

		Dist.RunningPWC.Ncent = Program.InitialNcent;

		Dist.ActualWaititerations = Program.Waititerations; // Wait this number of Temperature iterations before splitting
		// This is changed upto 8 times input value if no cluster to split

		EMIterationCount = 0; // Initialize iteration count -- there is no limit
		EMIterationStepCount = 0; // This is number of counts in current step converging epsilon at given temperature
		Dist.countAfterFixingClusterCount = 0; // Counts iterations after maximum cluster count reached (this decreases temperature)
		int HammyViolations = 0; // Counts number of (illegal) consequitive INCREASES in Hamiltonian
		int CountBetweenJiggles = 0; // Count iterations in Jiggle with maximum JiggleOption (Only used if JiggleOption > 0)
		int CountBetweenSplits = 0; // Count iterations between splits. Limit i ActualWaitIterations

		//  Variables to control pending task spun off when converging current case
		taskswaiting = false; // If true we saved a Solution in SavedPWC and evolved that solution to define output at a particular cluster count
		boolean currenttaskfinished = false;
		boolean TakePreviousSplitDecision = false; // If true, we decided on a split but spun off a convergence task

		Dist.needtocalculateMalpha_k_ = -1; // =0 Malpha_k_ Already Set; = 1 Calculate Malpha_k_ from EM (usual); = -1 Initialize Malpha_k_

		//  Variables tracking convergence
		HitConvergenceLoopLimit = false; // If true, EMIterationStepCount has hit limit
		SplitFailures = 0;
		justconverging = false; // If true we are in a special loop freezing current case
		Dist.RunningPWC.OldHammy = Tinitial * PWCUtility.PointCount_Global * PWCUtility.PointCount_Global; // Set Old value of Hamiltonian

		//	Loop over EM calculations
		while (true)
		{

			for (int ClusterIndex = 0; ClusterIndex < Dist.RunningPWC.Ncent; ClusterIndex++)
			{
				Dist.RunningPWC.Weight_k_[ClusterIndex] = 1.0;
			}

			Dist.EMIterationCount++; // Increment EM Loop Count
			Dist.EMIterationStepCount++; // Iterate within a fixed temperature converging epilson
			if (Program.EMIterationStepCountMax >= 0)
			{
				Program.EMIterationStepCountMax = Math.max(Program.EMIterationStepCountMax, Dist.EMIterationStepCount);
			}
			else
			{
				Program.EMIterationStepCountMax = Dist.EMIterationStepCount;
			}

			PairwiseThread();
			Dist.RunningPWC.SolutionSet = true;
			Dist.RunningPWC.IterationSetAt = Dist.EMIterationCount;

			//	Now see if we are done -- take results from Rank 0
			int convergence = 0;
			convergence = Dist.convergenceTest(Program.Epsi_max_change);
			convergence = PWCUtility.synchronizeMPIVariable(convergence);

			//  If not converged at this temperature and cluster number just proceed with while(true) iteration
			if (convergence == 0)
			{
				continue;
			}

// Case we are converged = 1 in middle or = 2 at end of refinement
			Dist.EMIterationStepCount = 0;

// Calculate measures of Clustering and save features for later printing
			double[] Save_ToPrint_C_k_ = new double[Dist.RunningPWC.Ncent];
			double[] Save_ToPrint_FreezingMeasure_k_ = new double[Dist.RunningPWC.Ncent];
			double[] Save_ToPrint_A_k_ = new double[Dist.RunningPWC.Ncent];
			double[] Save_ToPrint_P_k_ = new double[Dist.RunningPWC.Ncent];

			Dist.RunningPWC.PairwiseHammy = 0.0;
			for (int ClusterIndex = 0; ClusterIndex < Dist.RunningPWC.Ncent; ClusterIndex++)
			{
				Dist.RunningPWC.PairwiseHammy += -2.0 * Dist.RunningPWC.A_k_[ClusterIndex] * Dist.RunningPWC.C_k_[ClusterIndex];
				Save_ToPrint_C_k_[ClusterIndex] = Dist.RunningPWC.C_k_[ClusterIndex];
				Save_ToPrint_FreezingMeasure_k_[ClusterIndex] = Dist.RunningPWC.FreezingMeasure_k_[ClusterIndex];
				Save_ToPrint_A_k_[ClusterIndex] = Dist.RunningPWC.A_k_[ClusterIndex];
				Save_ToPrint_P_k_[ClusterIndex] = Dist.RunningPWC.P_k_[ClusterIndex];
			}

			//  Set Best Solution and progress flag decreasing
			boolean decreasing;
			if (HammyNotSet)
			{
				HammyNotSet = false;
				decreasing = true;
			}
			else
			{
				decreasing = Dist.RunningPWC.PairwiseHammy <= Dist.RunningPWC.OldHammy;
			}
			decreasing = PWCUtility.synchronizeMPIVariable(decreasing);
			if (decreasing)
			{
				ClusteringSolution.CopySolution(Dist.RunningPWC, Dist.BestPWC);
			}

			double ChangeinHammy = Dist.RunningPWC.PairwiseHammy - Dist.RunningPWC.OldHammy;
			Dist.RunningPWC.OldHammy = Dist.RunningPWC.PairwiseHammy;

			//  Case when at end of stage (e.g. given number of clusters or of whole job). This still needs to be iterated to low Temperatures
			if (convergence == 2 || justconverging)
			{
				if (Dist.countAfterFixingClusterCount < Program.Iterationatend) //    do Iterationatend iterations after reaching the maximum cluster#
				{
					if (Dist.countAfterFixingClusterCount == 0)
					{ // First step of "justconverging stage"
						HammyViolations = 0;
						Dist.RunningPWC.ActualCoolingFactor = Program.ConvergingCoolingFactor;
					}
					int toobigfreezing = 0;
					for (int ClusterIndex = 0; ClusterIndex < Dist.RunningPWC.Ncent; ClusterIndex++)
					{
						if (Dist.RunningPWC.FreezingMeasure_k_[ClusterIndex] > Program.FreezingLimit)
						{
							++toobigfreezing;
						}
					}
					toobigfreezing = PWCUtility.synchronizeMPIVariable(toobigfreezing);
					Dist.countAfterFixingClusterCount++;

					String looptype = "Final Clean Up";
					if (taskswaiting)
					{
						looptype = "Intermediate Cluster Convergence Steps";
					}
					int printtype = -1;
					if ((Dist.countAfterFixingClusterCount == (Program.Iterationatend - 1)) || (toobigfreezing == 0) || (!decreasing))
					{
						printtype = Dist.RunningPWC.IterationSetAt;
					}
					PrintIteration(looptype, printtype);

					if (decreasing)
					{
						HammyViolations = 0;
					}
					else
					{
						++HammyViolations;
						if (HammyViolations < 5)
						{
							decreasing = true;
						}
					}
					decreasing = PWCUtility.synchronizeMPIVariable(decreasing);
					if (PWCUtility.MPI_Rank == 0)
					{
						currenttaskfinished = !((toobigfreezing > 0) && decreasing);
					}
					currenttaskfinished = PWCUtility.synchronizeMPIVariable(currenttaskfinished);

					// Iteration at end still going Reduce Temperature and proceed
					if (!currenttaskfinished)
					{
						Dist.RunningPWC.Temperature = Dist.RunningPWC.ActualCoolingFactor * Dist.RunningPWC.Temperature;
						++CountBetweenJiggles;
						if (CountBetweenJiggles > Program.JiggleOption)
						{
							CountBetweenJiggles = 0;
						}
						if ((CountBetweenJiggles == 0) && (Program.JiggleOption > 0))
						{ // Jiggle Parameters
							JiggleClusters();
						}
						continue; // Continue while over EM Loop
					}

					// We have finished this task
					if (!decreasing)
					{
						PWCUtility.SALSAPrint(1, "Stop as Hamiltonian Increasing with Change " + String.format("%1$.4E", ChangeinHammy));
					}
					else
					{
						PWCUtility.SALSAPrint(1, "Stop as Freezing Measures smaller than " + (new Double(Program.FreezingLimit)).toString());
					}

				} // End Convergence=2  or justconverging doing Final Program.Iterationatend iterations counted by countAfterFixingClusterCount

				else
				{ // Case when Program.Iterationatend iterations counted by countAfterFixingClusterCount exceeded
					PWCUtility.SALSAPrint(1, " Stop as Final Iteration Count larger than " + (new Integer(Program.Iterationatend)).toString());
					currenttaskfinished = true;
				} // End Convergence=2 or justconverging case where extra iteration count completed


				// This completes processing for convergence=2 and justconverging=true case
				if (currenttaskfinished)
				{
					if (Dist.RunningPWC.Temperature < Dist.Tmin)
					{
						PWCUtility.SALSAPrint(1, "Tmin Reached " + String.format("%1$.4E", Dist.Tmin));
					}
					if (Dist.HitConvergenceLoopLimit)
					{
						PWCUtility.SALSAPrint(1, "EM Convergence Loop Reached for fixed Temperature -- Terminate Run");
					}
					if (Dist.SplitFailures > 1)
					{
						PWCUtility.SALSAPrint(1, "Consequitive Split Failures");
					}
					if (!taskswaiting)
					{
						if (!Dist.CheckValidSolution())
						{ // Restart with small clusters removed
							OnLastleg = true;
							Dist.needtocalculateMalpha_k_ = 1;
							Dist.countAfterFixingClusterCount = 0;
							CountBetweenJiggles = 0;
							CountBetweenSplits = 0;
							HammyNotSet = true;
							justconverging = true;
							taskswaiting = false;
							currenttaskfinished = false;
							TakePreviousSplitDecision = false;
							continue;
						}
						break;
					}

					//  Resume Previous Task
					PWCUtility.SALSAPrint(1, "\nResuming Previous Task -- Implement Cluster Splitting");

					//  Output Converged Results
					if (Program.ClusterCountOutput > 0)
					{
						int[] counts = new int[Dist.RunningPWC.Ncent];
						Program.OutputClusterLabels(counts);
						PWCUtility.SALSAPrint(1, "Clusters Output " + (new Integer(Dist.RunningPWC.Ncent)).toString());
					}

					//  Restore Previous Task
					RestorePreviousTask();
					Dist.countAfterFixingClusterCount = 0; // Counts iterations after maximum cluster count reached
					Dist.needtocalculateMalpha_k_ = 1; // =0 Already Set; = 1 Calculate from EM (usual); = -1 Initialize
					CountBetweenJiggles = 0;
					CountBetweenSplits = 0;
					justconverging = false;
					taskswaiting = false;
					Dist.HitConvergenceLoopLimit = false;
					currenttaskfinished = false;
					TakePreviousSplitDecision = true;
					convergence = 1;
				} // end case when task finished

			} // End justconverging = true or convergence =2 case

			//  This section results either in EM Loop Continue (for ongoing refinement for a nonsplittable case) or EM Loop Break (to end) OR
			//  Restart previous task which has convergence=1 by definitio

			//  Convergence = 1 Case
			//	Converged so test for split -  -- take results from Rank 0 but rest do arithmetic
			//  For restarted task that decision has already been made
			boolean ResultofSplittingTest = false;
			if (!TakePreviousSplitDecision)
			{ // Need to decide if split to occur
				++CountBetweenSplits;
				if (CountBetweenSplits > Dist.ActualWaititerations)
				{
					CountBetweenSplits = 0;
				}
				if ((CountBetweenSplits > 0) && (Dist.ActualWaititerations > 0) && (Dist.RunningPWC.Ncent > 1))
				{
					PrintIteration("Ongoing Annealing", -1);
					Dist.RunningPWC.Temperature = Dist.RunningPWC.ActualCoolingFactor * Dist.RunningPWC.Temperature;
					++CountBetweenJiggles;
					if (CountBetweenJiggles > Program.JiggleOption)
					{
						CountBetweenJiggles = 0;
					}
					if ((CountBetweenJiggles == 0) && (Program.JiggleOption > 0))
					{ // Jiggle Parameters
						JiggleClusters();
					}
					continue; // EM Loop for compulsory iterations between splits
				}

				//  Decide if to split
				PWCUtility.StartSubTimer(0);
				ResultofSplittingTest = shouldweSplit();
                ResultofSplittingTest = PWCUtility.synchronizeMPIVariable(ResultofSplittingTest);
                Dist.RunningPWC.ClustertoSplit = PWCUtility.synchronizeMPIVariable(Dist.RunningPWC.ClustertoSplit);
				PWCUtility.StopSubTimer(0);
				PWCUtility.InterimTiming();

				//  Diagnostic Output for splitting
				if (ResultofSplittingTest || (Dist.EMIterationCount % Program.PrintInterval == 0)) // *********************
				{
					diagnosticsplitprint(ResultofSplittingTest, Save_ToPrint_A_k_, Save_ToPrint_C_k_, Save_ToPrint_FreezingMeasure_k_, Save_ToPrint_P_k_);
				}

			} // End !TakePreviousSplitDecision -- case when splitting needs to be determined -- rather than taken from old task

			else
			{ // Case where splitting determined earlier
				ResultofSplittingTest = true;
			}

			//	If split indicated perform this                          
			if (ResultofSplittingTest)
			{
				if (Dist.RunningPWC.Ncent > 1)
				{
					if (!TakePreviousSplitDecision)
					{
						if (Program.ConvergeIntermediateClusters)
						{ // Start iteration to converge current cluster output
							SaveCurrentTask();
							justconverging = true;
							taskswaiting = true;
							TakePreviousSplitDecision = true;
							Dist.ActualWaititerations = Program.Waititerations;
							Dist.RunningPWC.ActualCoolingFactor = Program.ConvergingCoolingFactor;
							Dist.countAfterFixingClusterCount = 0;
							CountBetweenSplits = 0;
							CountBetweenJiggles = 0;
							continue;
						}
						if ((Program.ClusterCountOutput > 0) && Program.ConvergeIntermediateClusters)
						{
							int[] counts = new int[Dist.RunningPWC.Ncent];
							Program.OutputClusterLabels(counts);
						}
					}
				}

				//  Need to perform already determined split
				TakePreviousSplitDecision = false;
				if (Dist.Numberthatcanbesplit > 1)
				{
					String nextline = "";
					for (int splitlistloop = 0; splitlistloop < Dist.Numberthatcanbesplit; splitlistloop++)
					{
						Dist.RunningPWC.ClustertoSplit = Dist.ListofClusterstoSplit[splitlistloop];
						nextline += (new Integer(Dist.RunningPWC.ClustertoSplit)).toString() + " " + String.format("%1$.4E", Dist.EigsofClusterstoSplit[splitlistloop]) + " * ";
						dothesplit();
					}
					PWCUtility.SALSAPrint(1, " Multiple Clusters Split Center Eigenvalue " + nextline);
					Dist.Numberthatcanbesplit = 0;
				}
				else
				{
					dothesplit();
				}
				Dist.SplitFailures = 0;
				Dist.ActualWaititerations = Program.Waititerations;
			} // End ResultofSplittingTest == true (either changed number of clusters or spun off a convergence task)

			//  Final portion of loop changing Temperature if needed
			if (!ResultofSplittingTest)
			{ //    Reduce T and continue iterating if converged and no split
				Dist.RunningPWC.Temperature = Dist.RunningPWC.ActualCoolingFactor * Dist.RunningPWC.Temperature;
				++CountBetweenJiggles;
				if (CountBetweenJiggles > Program.JiggleOption)
				{
					CountBetweenJiggles = 0;
				}
				if ((CountBetweenJiggles == 0) && (Program.JiggleOption > 0))
				{ // Jiggle Parameters
					JiggleClusters();
				}
			}
			else
			{ // Don't reduce T as clusters split
				Dist.RunningPWC.ActualCoolingFactor = Program.FineCoolingFactor;
				CountBetweenSplits = 0;
				CountBetweenJiggles = 0;
			}

		} // End while EM Loop

		//  Check if solution best
		if (!Dist.BestPWC.SolutionSet)
		{
			return;
		}
		if (Dist.BestPWC.Ncent != Dist.RunningPWC.Ncent)
		{
			PWCUtility.SALSAPrint(1, " Best Solution Not Used as Number of Centers " + Dist.BestPWC.Ncent + " Different from Running Solution with " + Dist.RunningPWC.Ncent);
			return;
		}
		boolean changesolution = Dist.BestPWC.PairwiseHammy < Dist.RunningPWC.PairwiseHammy;
		changesolution = PWCUtility.synchronizeMPIVariable(changesolution);
		if (changesolution)
		{
			PWCUtility.SALSAPrint(1, " Solution at Iteration " + (new Integer(Dist.BestPWC.IterationSetAt)).toString() + " Chisq " + String.format("%1$.4E", Dist.BestPWC.PairwiseHammy) + " Taken rather than Iteration " + (new Integer(Dist.RunningPWC.IterationSetAt)).toString() + " Chisq " + String.format("%1$.4E", Dist.RunningPWC.PairwiseHammy));
			ClusteringSolution.CopySolution(Dist.BestPWC, Dist.RunningPWC);
		}
		Dist.PrintIteration(" Final Solution ", Dist.RunningPWC.IterationSetAt);

	} // End of getDist

	public static void PrintIteration(String looptype, int linecheck)
	{
		if (linecheck < 0)
		{
			if (Dist.EMIterationCount % Program.PrintInterval != 0)
			{
				return;
			}
		}
		else if (linecheck == Dist.IterationNumberPrinted)
		{
			return;
		}
		Dist.IterationNumberPrinted = Dist.RunningPWC.IterationSetAt;

		PWCUtility.InterimTiming();
		PWCUtility.SALSAPrint(1, "B) " + (new Integer(Dist.RunningPWC.Ncent)).toString() + " Iter " + (new Integer(Dist.EMIterationCount)).toString() + " T " + String.format("%1$.4E", Dist.RunningPWC.Temperature) + " PWHammy " + String.format("%1$.4E", Dist.RunningPWC.PairwiseHammy) + " PMLoop " + (new Integer(Dist.PMLoopUsed)).toString() + " " + looptype + " Time " + String.format("%1$.0f", PWCUtility.HPDuration));

		String nextline = "";
		for (int ClusterIndex = 0; ClusterIndex < Dist.RunningPWC.Ncent; ClusterIndex++)
		{
			double tmp = -2.0 * Dist.RunningPWC.A_k_[ClusterIndex];
			if (ClusterIndex != 0)
			{
				nextline += "* ";
			}
			nextline += (new Integer(ClusterIndex)).toString() + " C " + String.format("%1$.1f", Dist.RunningPWC.C_k_[ClusterIndex]) + " (Frz " + String.format("%1$.6f", Dist.RunningPWC.FreezingMeasure_k_[ClusterIndex]) + ") " + "[Wdth " + String.format("%1$.4f", tmp) + "] ";
			if (Program.ContinuousClustering)
			{
				nextline += "P " + String.format("%1$.4f", Dist.RunningPWC.P_k_[ClusterIndex]) + " ";
			}
		}
		PWCUtility.SALSAPrint(1, nextline);

	} // End Print Iteration from getDist

	public final void diagnosticsplitprint(boolean ResultofSplittingTest, double[] Local_A_k_, double[] Local_C_k_, double[] Local_FreezingMeasure_k_, double[] Local_P_k_)
	{
		String nextline1 = "A) " + (new Integer(Dist.RunningPWC.Ncent)).toString() + " Iter " + (new Integer(Dist.EMIterationCount)).toString() + " T " + String.format("%1$.4E", Dist.RunningPWC.Temperature) + " PWHammy " + String.format("%1$.4E", Dist.RunningPWC.PairwiseHammy) + " PMLoop " + (new Integer(Dist.PMLoopUsed)).toString() + " Time " + String.format("%1$.0f", PWCUtility.HPDuration);
		if (Dist.RunningPWC.ClustertoSplit >= 0)
		{
			if (ResultofSplittingTest)
			{
				nextline1 += " C# To Split " + (new Integer(Dist.RunningPWC.ClustertoSplit)).toString();
			}
			else
			{
				nextline1 += " No Split";
			}
			nextline1 += " " + ResultofSplittingTest + " Status " + vectorclass.eigenconverged[Dist.RunningPWC.ClustertoSplit] + " Pass 1 " + String.format("%1$.4E", vectorclass.Eigenvalues_Current[Dist.RunningPWC.ClustertoSplit]) + " Pass 0 " + String.format("%1$.4E", vectorclass.Eigenvalues_Pass0[Dist.RunningPWC.ClustertoSplit]);
		}
		PWCUtility.SALSAPrint(1, nextline1);
		nextline1 = "   ";
		double a_k_sum = 0.0;
		for (int ClusterIndex = 0; ClusterIndex < Dist.RunningPWC.Ncent; ClusterIndex++)
		{
			a_k_sum += Local_A_k_[ClusterIndex];
		}
		for (int ClusterIndex = 0; ClusterIndex < Dist.RunningPWC.Ncent; ClusterIndex++)
		{
			double tmp = Local_A_k_[ClusterIndex] / a_k_sum;
			String endofline = "]";
			if (Program.ContinuousClustering)
			{
				endofline = "] [P " + String.format("%1$.4f", Local_P_k_[ClusterIndex]) + "]";
			}
			if (ClusterIndex != 0)
			{
				nextline1 += " * ";
			}
			nextline1 += (new Integer(ClusterIndex)).toString() + " " + String.format("%1$.1f", Local_C_k_[ClusterIndex]) + " (Frz " + String.format("%1$.6f", Local_FreezingMeasure_k_[ClusterIndex]) + ") [Hfract" + String.format("%1$.4f", tmp) + endofline;
		}
		PWCUtility.SALSAPrint(1, nextline1);
		PWCUtility.SALSAPrint(1, " ");
	} // End diagnosticsplitprint

	public final void PairwiseThread() throws MPIException {
		Diff_Epsilon_k_ = new double[Dist.RunningPWC.Ncent];

		GlobalReductions.FindVectorDoubleSum Find_C_k_ = new GlobalReductions.FindVectorDoubleSum(PWCUtility.ThreadCount, Dist.RunningPWC.Ncent);
		GlobalReductions.FindVectorDoubleSum Find_FreezingMeasure_k_ = new GlobalReductions.FindVectorDoubleSum(PWCUtility.ThreadCount, Dist.RunningPWC.Ncent);
		GlobalReductions.FindDoubleSum Find_ChangeinM = new GlobalReductions.FindDoubleSum(PWCUtility.ThreadCount);
		GlobalReductions.FindDoubleSum Find_NormalizechangeinM = new GlobalReductions.FindDoubleSum(PWCUtility.ThreadCount);

		int pmloopmax = 1;
		if (Program.ContinuousClustering && (Dist.needtocalculateMalpha_k_ == 1) && (Dist.RunningPWC.Ncent > 1))
		{
			pmloopmax = 40;
		}
		double ChangeinM = 0.0;
		double NormalizechangeinM = 0.0;
		Dist.PMLoopUsed = 0;

		for (int pmloop = 0; pmloop < pmloopmax; pmloop++)
		{
			if (pmloop > 0)
			{
				ChangeinM = 0.0;
				Find_C_k_.zero();
				Find_FreezingMeasure_k_.zero();
				Find_ChangeinM.zero();
			}

            // Note - parallel for
            int pmloopmaxLoopVar = pmloopmax;
            int pmloopLoopVar = pmloop;
            try {
                forallChunked(0, PWCUtility.ThreadCount - 1, (threadIndex) ->
                        {
                            //	Start Code setting Malpha_k_ and partialsum_C_k_
                            int localNcent = Dist.RunningPWC.Ncent;
                            double[] Accum_C_k_ = new double[localNcent];
                            double[] Accum_FreezingMeasure_k_ = new double[localNcent];
                            for (int ClusterIndex = 0; ClusterIndex < localNcent; ClusterIndex++) {
                                Accum_C_k_[ClusterIndex] = 0.0;
                                Accum_FreezingMeasure_k_[ClusterIndex] = 0.0;
                            }
                            int indexlen = PWCUtility.PointsperThread[threadIndex];
                            int beginpoint = PWCUtility.StartPointperThread[threadIndex] - PWCUtility.PointStart_Process;
                            for (int ProcessPointIndex = beginpoint; ProcessPointIndex < indexlen + beginpoint; ProcessPointIndex++) {
                                //	calculate Malpha_k_ if needed (controlled by Dist.needtocalculateMalpha_k_)
                                //                  if ((Dist.needtocalculateMalpha_k_ == -1) && (Dist.Ncent > 2)) Dist.needtocalculateMalpha_k_ = 1;
                                if (Dist.needtocalculateMalpha_k_ == -1) {
                                    for (int ClusterIndex = 0; ClusterIndex < localNcent; ClusterIndex++) {
                                        //                          double fudge = 1.0 + Program.PerturbationFactor;
                                        //                         if (PointIndex % 2 == 0) fudge = 2.0 - fudge;
                                        Dist.RunningPWC.Malpha_k_[ProcessPointIndex][ClusterIndex] = 1.0 / Dist.RunningPWC.Ncent;
                                        //                          Malpha_k_[PointIndex][1] = 1.0 - 0.5 * fudge;
                                    }
                                }
                                if (Dist.needtocalculateMalpha_k_ == 1) {
                                    double Minepsi = 0.0;
                                    for (int ClusterIndex = 0; ClusterIndex < localNcent; ClusterIndex++) {
                                        double tmpepsi = Dist.RunningPWC.Epsilonalpha_k_[ProcessPointIndex][ClusterIndex];
                                        if (ClusterIndex == 0) {
                                            Minepsi = tmpepsi;
                                        } else {
                                            Minepsi = Math.min(Minepsi, tmpepsi);
                                        }
                                    }
                                    double tmp = 0.0;
                                    for (int ClusterIndex = 0; ClusterIndex < localNcent; ClusterIndex++) {
                                        double tmpepsi = Dist.RunningPWC.Epsilonalpha_k_[ProcessPointIndex][ClusterIndex];
                                        Dist.RunningPWC.Malpha_k_[ProcessPointIndex][ClusterIndex] = Math.exp(
                                                -(tmpepsi - Minepsi) / Dist.RunningPWC.Temperature);
                                        if (Program.ContinuousClustering) {
                                            Dist.RunningPWC.Malpha_k_[ProcessPointIndex][ClusterIndex] *= Dist.RunningPWC.P_k_[ClusterIndex];
                                        }
                                        tmp += Dist.RunningPWC.Weight_k_[ClusterIndex] * Dist.RunningPWC.Malpha_k_[ProcessPointIndex][ClusterIndex];
                                    }
                                    for (int ClusterIndex = 0; ClusterIndex < localNcent; ClusterIndex++) {
                                        if (tmp != 0) {
                                            Dist.RunningPWC.Malpha_k_[ProcessPointIndex][ClusterIndex] = Dist.RunningPWC.Malpha_k_[ProcessPointIndex][ClusterIndex] / tmp;
                                        } else {
                                            Dist.RunningPWC.Malpha_k_[ProcessPointIndex][ClusterIndex] = 1.0 / localNcent;
                                        }
                                        if (pmloopmaxLoopVar > 1) {
                                            if (pmloopLoopVar == 0) {
                                                Find_NormalizechangeinM.addapoint(threadIndex,
                                                        Dist.RunningPWC.Malpha_k_[ProcessPointIndex][ClusterIndex]);
                                            }
                                            if (pmloopLoopVar > 0) {
                                                Find_ChangeinM.addapoint(threadIndex, Math.abs(
                                                        Dist.RunningPWC.Previous_Malpha_k_[ProcessPointIndex][ClusterIndex] - Dist.RunningPWC.Malpha_k_[ProcessPointIndex][ClusterIndex]));
                                            }
                                            Dist.RunningPWC.Previous_Malpha_k_[ProcessPointIndex][ClusterIndex] = Dist.RunningPWC.Malpha_k_[ProcessPointIndex][ClusterIndex];
                                        }
                                    }
                                }
                                for (int ClusterIndex = 0; ClusterIndex < localNcent; ClusterIndex++) {
                                    double tmp = Dist.RunningPWC.Malpha_k_[ProcessPointIndex][ClusterIndex];
                                    Accum_C_k_[ClusterIndex] += tmp;
                                    Accum_FreezingMeasure_k_[ClusterIndex] += tmp * (1.0 - tmp);
                                }
                            }
                            Find_C_k_.addAPoint(threadIndex, Accum_C_k_);
                            Find_FreezingMeasure_k_.addAPoint(threadIndex, Accum_FreezingMeasure_k_);
                        }
                );
            } catch (SuspendableException e) {
                PWCUtility.printAndThrowRuntimeException(e.getMessage());
            }

            Find_C_k_.sumOverThreadsAndMPI();
			Find_FreezingMeasure_k_.sumOverThreadsAndMPI();
			if (pmloop > 0)
			{
				Find_ChangeinM.sumoverthreadsandmpi();
				ChangeinM = Find_ChangeinM.Total;
			}
			if ((pmloop == 0) && (pmloopmax > 1))
			{
				Find_NormalizechangeinM.sumoverthreadsandmpi();
				NormalizechangeinM = Find_NormalizechangeinM.Total;
			}

			//  Form C_k_ and P_k_ from sum over threads
			for (int ClusterCount = 0; ClusterCount < Dist.RunningPWC.Ncent; ClusterCount++)
			{
				Dist.RunningPWC.C_k_[ClusterCount] = Find_C_k_.TotalVectorSum[ClusterCount];
				Dist.RunningPWC.FreezingMeasure_k_[ClusterCount] = Find_FreezingMeasure_k_.TotalVectorSum[ClusterCount];
			}

			for (int ClusterCount = 0; ClusterCount < Dist.RunningPWC.Ncent; ClusterCount++)
			{
				double tmp = Dist.RunningPWC.C_k_[ClusterCount];
				if (tmp > 0.0)
				{
					Dist.RunningPWC.FreezingMeasure_k_[ClusterCount] = Dist.RunningPWC.FreezingMeasure_k_[ClusterCount] / tmp;
				}
				if (Program.ContinuousClustering)
				{
					Dist.RunningPWC.P_k_[ClusterCount] = Dist.RunningPWC.C_k_[ClusterCount] / PWCUtility.PointCount_Global;
				}

			}
			//  See Shift in C after split
			if (Program.TestExpectedChange >= 0)
			{
				double C1 = Dist.RunningPWC.C_k_[Program.TestExpectedChange];
				double C2 = Dist.RunningPWC.C_k_[Dist.RunningPWC.Ncent - 1];
				PWCUtility.SALSAPrint(0, "Cluster " + (new Integer(Program.TestExpectedChange)).toString() + " Method " + (new Integer(Program.ExpectedChangeMethod)).toString() + " Expected Change " + String.format("%1$.4E", Program.ExpectedChange) + " Previous C " + String.format("%1$.4E", Program.PreviousC) + " New " + String.format("%1$.6E", C1) + " " + String.format("%1$.6E", C2) + " Temp " + String.format("%1$.4e", Dist.RunningPWC.Temperature) + " " + (new Integer(Dist.EMIterationCount)).toString());
				Program.TestExpectedChange = -1;
				if (Program.ExpectedChangeMethod == 1)
				{
					Program.deltaCoverC[0] += Math.abs(C1 - C2) / Program.PreviousC;
					Program.CountdeltaCoverC[0]++;
				}
				else
				{
					Program.deltaCoverC[1] += Math.abs(C1 - C2) / Program.PreviousC;
					Program.CountdeltaCoverC[1]++;
				}


			}
			//  End iteration over M and p given epsilon
			if (pmloopmax > 1)
			{
				Dist.PMLoopUsed = pmloop;
				boolean pmloopend = (pmloop > 0) && (ChangeinM < 0.001 * NormalizechangeinM);
				pmloopend = PWCUtility.synchronizeMPIVariable(pmloopend);
				if (pmloopend)
				{
					break;
				}
				if (pmloop == pmloopmax - 1)
				{
					PWCUtility.SALSAPrint(1, " pmloop limit reached " + " Iter " + (new Integer(Dist.EMIterationCount)).toString() + " T " + String.format("%1$.4E", Dist.RunningPWC.Temperature) + " PWHammy " + String.format("%1$.4E", Dist.RunningPWC.PairwiseHammy) + " " + String.format("%1$.4E", ChangeinM) + " " + String.format("%1$.4E", NormalizechangeinM));
				}
			}
		} // End pmloop
		Program.CountPMExtraIterations += Dist.PMLoopUsed;

		calculateEpsi(Dist.RunningPWC.Malpha_k_, Dist.RunningPWC.A_k_, Dist.RunningPWC.Balpha_k_, Dist.RunningPWC.C_k_, Dist.RunningPWC.Epsilonalpha_k_, Dist.RunningPWC.Ncent);
		Dist.needtocalculateMalpha_k_ = 1; // as epsi are now set

    } // End PairwiseThread

	public int oldMaxLength = Integer.MIN_VALUE;
	public static MPIPacket fromafar = null;
	public static MPIPacket toafar = null;
	public static MPIPacket myown = null;


	//	Perform multiple parallel steps calculating A_k_, Balpha_k_, epsi and differences
	public final void calculateEpsi(double[][] localMalpha_k_, double[] localA_k_, double[][] localBalpha_k_, double[] localC_k_, double[][] localepsi, int localNcent) throws MPIException {

		int sendtag = 0;
		int receivetag = 0;

		int Maxlength = localNcent * PWCUtility.PointCount_Largest;

		if (oldMaxLength != Maxlength)
		{
			fromafar = MPIPacket.newDoublePacket(Maxlength);
			toafar = MPIPacket.newDoublePacket(Maxlength);
			myown = MPIPacket.newDoublePacket(Maxlength);
			oldMaxLength = Maxlength;
		}
		else
		{
			fromafar.Clear();
			toafar.Clear();
			myown.Clear();
		}

		myown.setFirstPoint(PWCUtility.PointStart_Process);
		myown.setNumberOfPoints(PWCUtility.PointCount_Process);
		if (PWCUtility.PointCount_Process < PWCUtility.PointCount_Largest)
		{
			for (int countC = 0; countC < localNcent; countC++)
			{
				myown.setMArrayDoubleAt(PWCUtility.PointCount_Process * localNcent + countC,0.0);
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
				toafar.setFirstPoint(PWCUtility.PointStart_Process);
				toafar.setNumberOfPoints(PWCUtility.PointCount_Process);
				if (PWCUtility.PointCount_Process < PWCUtility.PointCount_Largest)
				{
					for (int countC = 0; countC < localNcent; countC++)
					{
						int bigindex = PWCUtility.PointCount_Process * localNcent + countC;
						toafar.setMArrayDoubleAt(bigindex, myown.getMArrayDoubleAt(bigindex));
					}
				}
			}
			if (MPICommunicationSteps > 1)
			{
				toafar.setFirstPoint(fromafar.getFirstPoint());
				toafar.setNumberOfPoints(fromafar.getNumberOfPoints());
				if (PWCUtility.PointCount_Process < PWCUtility.PointCount_Largest)
				{
					for (int countC = 0; countC < localNcent; countC++)
					{
						int bigindex = PWCUtility.PointCount_Process * localNcent + countC;
						toafar.setMArrayDoubleAt(bigindex, fromafar.getMArrayDoubleAt(bigindex));
					}
				}
			}
			if (MPICommunicationSteps > 0)
			{
				PWCUtility.StartSubTimer(PWCUtility.MPISENDRECEIVETiming);
                // Note - MPI Call - SendRecv - MPIPacket
                fromafar = PWCUtility.mpiOps.sendReceive(toafar,toprocess,sendtag,fromprocess,receivetag, MPIPacket.Type.Double);
				PWCUtility.StopSubTimer(PWCUtility.MPISENDRECEIVETiming);
			}

            // Note - parallel for
            final int MPICommunicationStepsLoopVar = MPICommunicationSteps;
            try {
                forallChunked(0, PWCUtility.ThreadCount-1, (threadIndex) ->
                {
                    //  Loop over Home (point) indices in thread using non-home values from MPI or locally for communicationloop == 0
                    int betastart, betatotal;
                    int indexlen = PWCUtility.PointsperThread[threadIndex];
                    int beginpoint = PWCUtility.StartPointperThread[threadIndex] - PWCUtility.PointStart_Process;
                    for (int ProcessPointIndex = beginpoint; ProcessPointIndex < indexlen + beginpoint; ProcessPointIndex++)
                    {
                        if (MPICommunicationStepsLoopVar == 0)
                        {
                            Arrays.fill(localBalpha_k_[ProcessPointIndex], 0, localNcent, 0.0);
                            betatotal = PWCUtility.PointCount_Process;
                            betastart = PWCUtility.PointStart_Process;
                            for (int ClusterIndex = 0; ClusterIndex < localNcent; ClusterIndex++)
                            {
                                double tmp = localMalpha_k_[ProcessPointIndex][ClusterIndex];
                                int bigindex = ProcessPointIndex * localNcent + ClusterIndex;
                                myown.setMArrayDoubleAt(bigindex, tmp);
                                toafar.setMArrayDoubleAt(bigindex, tmp);
                            }
                        }
                        else
                        {
                            betatotal = fromafar.getNumberOfPoints();
                            betastart = fromafar.getFirstPoint();
                            if (MPICommunicationStepsLoopVar != (PWCUtility.MPI_Size - 1))
                            {
                                for (int ClusterIndex = 0; ClusterIndex < localNcent; ClusterIndex++)
                                {
                                    int bigindex = ProcessPointIndex * localNcent + ClusterIndex;
                                    toafar.setMArrayDoubleAt(bigindex,fromafar.getMArrayDoubleAt(bigindex));
                                }
                            }
                        }
                        for (int betalocal = 0; betalocal < betatotal; betalocal++)
                        {
                            double tmp;
                            int betafull = betastart + betalocal;
                            double dijforthiscase = PWCUtility.PointDistances[ProcessPointIndex*PWCUtility.PointCount_Global+betafull] * PWCUtility.INV_SHORT_MAX;
                            for (int ClusterIndex = 0; ClusterIndex < localNcent; ClusterIndex++)
                            {
                                if (MPICommunicationStepsLoopVar == 0)
                                {
                                    tmp = localMalpha_k_[betalocal][ClusterIndex];
                                }
                                else
                                {
                                    tmp = fromafar.getMArrayDoubleAt(betalocal * localNcent + ClusterIndex);
                                }
                                localBalpha_k_[ProcessPointIndex][ClusterIndex] += dijforthiscase * tmp / localC_k_[ClusterIndex];
                            }
                        }
                    }
                }
               );
            } catch (SuspendableException e) {
                PWCUtility.printAndThrowRuntimeException(e.getMessage());
            }
        } // End loop over communicationloop

		//  Now calculate quantities involving global sums

		//	Calculate full A(k)
		GlobalReductions.FindVectorDoubleSum Find_A_k_ = new GlobalReductions.FindVectorDoubleSum(PWCUtility.ThreadCount, localNcent);

        // Note - parallel for
        try {
            forallChunked(0, PWCUtility.ThreadCount-1, (threadIndex) ->
            {
                double[] LocalContribution_A_k_ = new double[localNcent];
                int indexlen = PWCUtility.PointsperThread[threadIndex];
                int beginpoint = PWCUtility.StartPointperThread[threadIndex] - PWCUtility.PointStart_Process;
                for (int ProcessPointIndex = beginpoint; ProcessPointIndex < indexlen + beginpoint; ProcessPointIndex++)
                {
                    for (int ClusterIndex = 0; ClusterIndex < localNcent; ClusterIndex++)
                    {
                        LocalContribution_A_k_[ClusterIndex] = localBalpha_k_[ProcessPointIndex][ClusterIndex] * localMalpha_k_[ProcessPointIndex][ClusterIndex];
                    }
                    Find_A_k_.addAPoint(threadIndex, LocalContribution_A_k_);
                }
            }
           );
        } catch (SuspendableException e) {
            PWCUtility.printAndThrowRuntimeException(e.getMessage());
        }

        Find_A_k_.sumOverThreadsAndMPI();

		for (int ClusterIndex = 0; ClusterIndex < localNcent; ClusterIndex++)
		{
			localA_k_[ClusterIndex] = -0.5 * Find_A_k_.TotalVectorSum[ClusterIndex] / localC_k_[ClusterIndex];
		}

		// Calculate new values of epsi and do partial sums of differences   
		GlobalReductions.FindVectorDoubleSum Find_EpsiDiff = new GlobalReductions.FindVectorDoubleSum(PWCUtility.ThreadCount, localNcent);

        // Note - parallel for
        try {
            forallChunked(0, PWCUtility.ThreadCount-1, (threadIndex) ->
            {
                double[] Local_EpsiDiff = new double[localNcent];
                int indexlen = PWCUtility.PointsperThread[threadIndex];
                int beginpoint = PWCUtility.StartPointperThread[threadIndex] - PWCUtility.PointStart_Process;
                for (int ProcessPointIndex = beginpoint; ProcessPointIndex < indexlen + beginpoint; ProcessPointIndex++)
                {
                    for (int ClusterIndex = 0; ClusterIndex < localNcent; ClusterIndex++)
                    {
                        double tmp = localBalpha_k_[ProcessPointIndex][ClusterIndex] + localA_k_[ClusterIndex];
                        localepsi[ProcessPointIndex][ClusterIndex] = tmp;
                        if (Dist.oldepsiset > 0)
                        {
                            Local_EpsiDiff[ClusterIndex] = Math.abs(Dist.RunningPWC.Old_Epsilonalpha_k_[ProcessPointIndex][ClusterIndex] - tmp);
                        }
                            Dist.RunningPWC.Old_Epsilonalpha_k_[ProcessPointIndex][ClusterIndex] = tmp;
                    }
                    if (Dist.oldepsiset > 0)
                    {
                        Find_EpsiDiff.addAPoint(threadIndex, Local_EpsiDiff);
                    }
                }
            }
           );
        } catch (SuspendableException e) {
            PWCUtility.printAndThrowRuntimeException(e.getMessage());
        }

        if (Dist.oldepsiset > 0)
		{ //    Calculate epsidiff which is sum for each center over data points
			Find_EpsiDiff.sumOverThreadsAndMPI();

            System.arraycopy(Find_EpsiDiff.TotalVectorSum, 0, Diff_Epsilon_k_, 0, localNcent);

		} // End computation of epsidiff
		++Dist.oldepsiset;

	} // End calculateEpsi


	// The change of sum(points) Delta(epsi)/# Points should be less than argument ChangeLimit for each cluster
	// epsiDiff calculated in CalculateEpsi while Epsi are updated
	//  Return 0 if not converged; 1 if converged and can continue; 2 if converged but no refinement
	public static int convergenceTest(double ChangeLimit) throws MPIException {

		Dist.HitConvergenceLoopLimit = false;
		if (Dist.oldepsiset <= 1)
		{
			return 0; // Skip first two cases -- WHY?
		}
		double epsidiffmax = 0;
		for (int ClusterCount = 0; ClusterCount < Dist.RunningPWC.Ncent; ClusterCount++)
		{
			if (Dist.Diff_Epsilon_k_[ClusterCount] > epsidiffmax)
			{
				epsidiffmax = Dist.Diff_Epsilon_k_[ClusterCount];
			}
		}
		epsidiffmax = epsidiffmax / PWCUtility.PointCount_Global;
		boolean epslimit = (epsidiffmax <= ChangeLimit);
		epslimit = PWCUtility.synchronizeMPIVariable(epslimit);

		if (epslimit)
		{
			Dist.EMIterationStepCount = 0;
			if (Dist.OnLastleg || Dist.justconverging)
			{
				return 2;
			}
			boolean templimit = (Dist.RunningPWC.Temperature < Dist.Tmin);
			templimit = PWCUtility.synchronizeMPIVariable(templimit);
			if ((Dist.RunningPWC.Ncent >= Dist.ActualMaxNcent) || templimit || (Dist.SplitFailures > 1))
			{
				Dist.OnLastleg = true;
				return 2;
			}
			else
			{
				return 1;
			}

		}
		if (Dist.EMIterationStepCount > Program.ConvergenceLoopLimit)
		{
			PWCUtility.SALSAPrint(0, "EM Convergence Loop Reached for fixed Temperature at Iteration " + (new Integer(Dist.EMIterationCount)).toString());
			Dist.EMIterationStepCount = 0;
			Dist.HitConvergenceLoopLimit = true;
			return 2;
		}
		else
		{
			return 0;
		}

	} // End convergenceTest

	public final void SaveCurrentTask()
	{ // Save current task that should be split but we need to converge first

		ClusteringSolution.SetAxinSolution(Dist.RunningPWC);
		ClusteringSolution.CopySolution(Dist.RunningPWC, Dist.SavedPWC);

    } // End saving current task that should be split but we need to converge first

	public final void RestorePreviousTask()
	{ // Restore previous task that should be split but we needed to converge current cluster configuration first

		ClusteringSolution.CopySolution(Dist.SavedPWC, Dist.RunningPWC);
		ClusteringSolution.RestoreAxfromSolution(Dist.RunningPWC);
		Dist.RunningPWC.Axset = false;

    } // End Restore previous task that should be split but we needed to converge current cluster configuration first

	public final void JiggleClusters() throws MPIException { // Jiggle Values of Epsilon -- leaving other Cluster Parameters


		if (Dist.RunningPWC.Ncent == 1)
		{
			return;
		}

		// Initialize
		java.util.Random Randobject = new java.util.Random();
		vectorclass vcjiggle = new vectorclass();

		if (Program.JigglePerturbation == 2)
		{
			for (int ClusterIndex = 0; ClusterIndex < Dist.RunningPWC.Ncent; ClusterIndex++)
			{
				Dist.RunningPWC.Weight_k_[ClusterIndex] = 1.0;
				Dist.ClusterSelected[ClusterIndex] = 1;
			}
			Dist.RunningPWC.ClustertoSplit = 1;
			Dist.PlaceforEigsforMultipleCluster = Dist.RunningPWC.ClustertoSplit;

			// Note this calculates Full eigenvector of entire matrix -- NOT as in stability analysis, the eigenvector in each cluster sector
			vcjiggle.getEigenvalues(4, Dist.RunningPWC.Malpha_k_, Dist.RunningPWC.A_k_, Dist.RunningPWC.Balpha_k_, Dist.RunningPWC.C_k_, Dist.RunningPWC.Ncent);
			double TrueMinimumEigenvalue1 = vectorclass.Eigenvalues_Pass0[Dist.PlaceforEigsforMultipleCluster] - vectorclass.Eigenvalues_Current[Dist.PlaceforEigsforMultipleCluster];

			PWCUtility.SALSAPrint(1, " Jiggle Cluster " + (new Integer(Dist.PlaceforEigsforMultipleCluster)).toString() + " Status " + (new Integer(vectorclass.eigenconverged[Dist.PlaceforEigsforMultipleCluster])).toString() + " Pass 1 " + String.format("%1$.4E", vectorclass.Eigenvalues_Current[Dist.PlaceforEigsforMultipleCluster]) + " Pass 0 " + String.format("%1$.4E", vectorclass.Eigenvalues_Pass0[Dist.PlaceforEigsforMultipleCluster]));

			if ((vectorclass.Eigenvalues_Pass0[Dist.PlaceforEigsforMultipleCluster] <= 0.0) || (vectorclass.eigenconverged[Dist.PlaceforEigsforMultipleCluster] <= 0))
			{
				return;
			}
		} // End Initialization of JiggleOption = 2

		// NewC_k_ and AverageMalpha_k_Change calculated for output only
		GlobalReductions.FindVectorDoubleSum Find_NewC_k_ = new GlobalReductions.FindVectorDoubleSum(PWCUtility.ThreadCount, Dist.RunningPWC.Ncent);
		GlobalReductions.FindDoubleSum Find_AverageMalpha_k_Change = new GlobalReductions.FindDoubleSum(PWCUtility.ThreadCount);

        // Note - parallel for
        try {
            forallChunked(0, PWCUtility.ThreadCount-1, (threadIndex) ->
            {
                double[] NewMalpha_k_ = new double[Dist.RunningPWC.Ncent];
                double[] partialsum_NewC_k_ = new double[Dist.RunningPWC.Ncent];
                int indexlen = PWCUtility.PointsperThread[threadIndex];
                int beginpoint = PWCUtility.StartPointperThread[threadIndex] - PWCUtility.PointStart_Process;
                for (int ProcessPointIndex = beginpoint; ProcessPointIndex < indexlen + beginpoint; ProcessPointIndex++)
                {
                    double AverageMalpha_k_Change = 0.0;
                    double tmp = 0.0;
                    for (int ClusterIndex = 0; ClusterIndex < Dist.RunningPWC.Ncent; ClusterIndex++)
                    {
                        double perturb;
                        if (Program.JigglePerturbation == 1)
                        {
                            perturb = Randobject.nextDouble();
                        }
                        else
                        {
                            perturb = vectorclass.oldAx[ProcessPointIndex][ClusterIndex];
                        }
                        Dist.RunningPWC.Epsilonalpha_k_[ProcessPointIndex][ClusterIndex] += perturb * Program.JigglePerturbationFactor * Dist.RunningPWC.Temperature;
                        NewMalpha_k_[ClusterIndex] = Math.exp(-Dist.RunningPWC.Epsilonalpha_k_[ProcessPointIndex][ClusterIndex] / Dist.RunningPWC.Temperature);
                        tmp += NewMalpha_k_[ClusterIndex];
                    }
                    for (int ClusterIndex = 0; ClusterIndex < Dist.RunningPWC.Ncent; ClusterIndex++)
                    {
                        double addtoC = NewMalpha_k_[ClusterIndex] / tmp;
                        partialsum_NewC_k_[ClusterIndex] += addtoC;
                        AverageMalpha_k_Change += Math.abs(addtoC - Dist.RunningPWC.Malpha_k_[ProcessPointIndex][ClusterIndex]);
                    }
                    Find_NewC_k_.addAPoint(threadIndex, partialsum_NewC_k_);
                    Find_AverageMalpha_k_Change.addapoint(threadIndex, AverageMalpha_k_Change);
                }
            }
           );
        } catch (SuspendableException e) {
            PWCUtility.printAndThrowRuntimeException(e.getMessage());
        }

        Find_AverageMalpha_k_Change.sumoverthreadsandmpi();
		double FullAverageMalpha_k_Change = Find_AverageMalpha_k_Change.Total / (Dist.RunningPWC.Ncent * PWCUtility.PointCount_Global);
		Dist.needtocalculateMalpha_k_ = 1;

		double[] NewC_k_ = new double[Dist.RunningPWC.Ncent];
		Find_NewC_k_.sumOverThreadsAndMPI();
        System.arraycopy(Find_NewC_k_.TotalVectorSum, 0, NewC_k_, 0, Dist.RunningPWC.Ncent);

		String nextline = " Jiggle M change " + String.format("%1$.5f", FullAverageMalpha_k_Change) + "  Old(Jiggled) Sizes ";
		for (int ClusterIndex = 0; ClusterIndex < Dist.RunningPWC.Ncent; ClusterIndex++)
		{
			nextline += String.format("%1$.1f", Dist.RunningPWC.C_k_[ClusterIndex]) + " (" + String.format("%1$.6f", NewC_k_[ClusterIndex]) + ") ";
		}
		PWCUtility.SALSAPrint(1, nextline);

    } // End JiggleCluster()

	// Decide if to split based on negative eigenvalue of second derivative matrix
	// MinimumEigenvalue is MINIMUM eigenvalue
	// ClustertoSplit is cluster number of this
	// D is distance matrix
	//  Skip first time this happens as initialize to two clusters
	public final boolean shouldweSplit() throws MPIException {
		//	Calculate Cluster with minimum eigenvalue -- find cluster number and eigenvalue (which could be positive)
		Dist.Numberthatcanbesplit = 0;
		int LimitonSplits = Math.min(Program.MaxNumberSplitClusters, Dist.ActualMaxNcent - Dist.RunningPWC.Ncent);

		Dist.RunningPWC.ClustertoSplit = -1;
		Dist.RunningPWC.MinimumEigenvalue = 99999999999.0;
		int ActualMethodology = Program.Eigenvalue_Methodology;
		if (Dist.RunningPWC.Ncent == 1)
		{
			ActualMethodology = 2;
		}
		if (ActualMethodology == 0)
		{
			return false;
		}

		vectorclass vc = new vectorclass();
		for (int ClusterIndex = 0; ClusterIndex < Dist.RunningPWC.Ncent; ClusterIndex++)
		{
			vectorclass.eigenconverged[ClusterIndex] = 0;
			Dist.ClusterSelected[ClusterIndex] = 1;
		}

		// Continuous Clustering is simple case of ActualMethodology  = 1
		if (ActualMethodology <= 2)
		{
			vc.getEigenvalues(ActualMethodology, Dist.RunningPWC.Malpha_k_, Dist.RunningPWC.A_k_, Dist.RunningPWC.Balpha_k_, Dist.RunningPWC.C_k_, Dist.RunningPWC.Ncent);
			for (int ClusterToRefine = 0; ClusterToRefine < Dist.RunningPWC.Ncent; ClusterToRefine++)
			{
				if (Dist.RunningPWC.C_k_[ClusterToRefine] <= Program.ToosmalltoSplit)
				{
					continue;
				}
				double TrueMinimumEigenvalue = vectorclass.Eigenvalues_Pass0[ClusterToRefine] - vectorclass.Eigenvalues_Current[ClusterToRefine];
				PWCUtility.SALSAPrint(1, "       Cluster " + (new Integer(ClusterToRefine)).toString() + " Status " + (new Integer(vectorclass.eigenconverged[ClusterToRefine])).toString() + " Pass 1 " + String.format("%1$.4E", vectorclass.Eigenvalues_Current[ClusterToRefine]) + " Pass 0 " + String.format("%1$.4E", vectorclass.Eigenvalues_Pass0[ClusterToRefine]) + " Diff " + String.format("%1$.2E", TrueMinimumEigenvalue) + " Size " + String.format("%1$.1f", Dist.RunningPWC.C_k_[ClusterToRefine]));

				if ((vectorclass.Eigenvalues_Pass0[ClusterToRefine] > 0.0) || (vectorclass.eigenconverged[ClusterToRefine] > 0))
				{
					if (TrueMinimumEigenvalue < Dist.RunningPWC.MinimumEigenvalue)
					{
						Dist.RunningPWC.MinimumEigenvalue = TrueMinimumEigenvalue;
						Dist.RunningPWC.ClustertoSplit = ClusterToRefine;
					}
					if (LimitonSplits > 1)
					{
						double eigtest1 = Program.MinEigtest * vectorclass.Eigenvalues_Pass0[ClusterToRefine];
						if (TrueMinimumEigenvalue < eigtest1)
						{ // Candidate for Split List
							if (Dist.Numberthatcanbesplit == 0) // Initialize Split List
							{
								Dist.Numberthatcanbesplit = 1;
								Dist.EigsofClusterstoSplit[0] = TrueMinimumEigenvalue;
								Dist.ListofClusterstoSplit[0] = ClusterToRefine;
							}
							else // Add to Split List
							{
								int position = Dist.Numberthatcanbesplit;
								for (int positionloop = 0; positionloop < Dist.Numberthatcanbesplit; positionloop++)
								{
									if (TrueMinimumEigenvalue >= Dist.EigsofClusterstoSplit[positionloop])
									{
										continue;
									}
									position = positionloop;
									break;
								}
								if (position >= LimitonSplits)
								{
									continue;
								}
								for (int positionloop = Dist.Numberthatcanbesplit - 1; positionloop >= position; positionloop--)
								{
									if (positionloop == (LimitonSplits - 1))
									{
										continue;
									}
									Dist.EigsofClusterstoSplit[positionloop + 1] = Dist.EigsofClusterstoSplit[positionloop];
									Dist.ListofClusterstoSplit[positionloop + 1] = Dist.ListofClusterstoSplit[positionloop];
								}
								Dist.Numberthatcanbesplit = Math.min(Dist.Numberthatcanbesplit + 1, LimitonSplits);
								Dist.EigsofClusterstoSplit[position] = TrueMinimumEigenvalue;
								Dist.ListofClusterstoSplit[position] = ClusterToRefine;
							}
						}
					}
				}
			}

			if (Program.PerformEigenTest && (Dist.RunningPWC.Ncent == 1))
			{
				// Eigenvalue Test for Metholodology 2 -- only sensible if Ncent = 1 as otherwise Method 2 inaccurate so test will fail
				Dist.RunningPWC.ClustertoSplit = 0;

                try {
                    forallChunked(0, PWCUtility.ThreadCount-1, (threadIndex) ->
                    {
                        int indexlen = PWCUtility.PointsperThread[threadIndex];
                        int beginpoint = PWCUtility.StartPointperThread[threadIndex] - PWCUtility.PointStart_Process;
                        for (int index = beginpoint; index < indexlen + beginpoint; index++)
                        {
                            Dist.RunningPWC.Epsilonalpha_k_[index][Dist.RunningPWC.Ncent] = Dist.RunningPWC.Epsilonalpha_k_[index][Dist.RunningPWC.ClustertoSplit];
                            Dist.RunningPWC.Old_Epsilonalpha_k_[index][Dist.RunningPWC.Ncent] = Dist.RunningPWC.Epsilonalpha_k_[index][Dist.RunningPWC.ClustertoSplit];
                        }
                    }
                   );
                } catch (SuspendableException e) {
                    PWCUtility.printAndThrowRuntimeException(e.getMessage());
                }


                Dist.RunningPWC.Ncent++;
				for (int ClusterIndex = 0; ClusterIndex < Dist.RunningPWC.Ncent; ClusterIndex++)
				{
					Dist.RunningPWC.Weight_k_[ClusterIndex] = 1.0;
					Dist.ClusterSelected[ClusterIndex] = 0;
				}
				Dist.ClusterSelected[0] = 1;
				Dist.ClusterSelected[1] = 1;
				Dist.PlaceforEigsforMultipleCluster = Dist.RunningPWC.ClustertoSplit;
				PairwiseThread();
				vc.getEigenvalues(4, Dist.RunningPWC.Malpha_k_, Dist.RunningPWC.A_k_, Dist.RunningPWC.Balpha_k_, Dist.RunningPWC.C_k_, Dist.RunningPWC.Ncent);
				double TrueMinimumEigenvalue1 = vectorclass.Eigenvalues_Pass0[Dist.PlaceforEigsforMultipleCluster] - vectorclass.Eigenvalues_Current[Dist.PlaceforEigsforMultipleCluster];

				String nextline = "  C)Test Cluster " + Dist.RunningPWC.ClustertoSplit + " T " + String.format("%1$.4E", Dist.RunningPWC.Temperature) + " PMLoop " + (new Integer(Dist.PMLoopUsed)).toString() + " PWHammy " + String.format("%1$.4E", Dist.RunningPWC.PairwiseHammy) + " Status " + (new Integer(vectorclass.eigenconverged[Dist.PlaceforEigsforMultipleCluster])).toString() + " Status " + (new Integer(vectorclass.eigenconverged[Dist.PlaceforEigsforMultipleCluster])).toString() + " Pass 1 " + String.format("%1$.4E", vectorclass.Eigenvalues_Current[Dist.PlaceforEigsforMultipleCluster]) + " Pass 0 " + String.format("%1$.4E", vectorclass.Eigenvalues_Pass0[Dist.PlaceforEigsforMultipleCluster]) + " C ";
				for (int ClusterIndex = 0; ClusterIndex < Dist.RunningPWC.Ncent; ClusterIndex++)
				{
					nextline += String.format("%1$.1f", Dist.RunningPWC.C_k_[ClusterIndex]) + " ";
				}
				PWCUtility.SALSAPrint(1, nextline);

				Dist.RunningPWC.Ncent--;
			} // End Test of Eigenvalues for Methodology 2

		} // End ActualMethodology <= 2

		else
		{
			// Implement exact Duplication separately for each cluster ActualMethodology = 3 or 4
            // Parallel reSetting of Epsilon

            // Note - parallel for
            try {
                forallChunked(0, PWCUtility.ThreadCount-1, (threadIndex) ->
                {
                    int indexlen = PWCUtility.PointsperThread[threadIndex];
                    int beginpoint = PWCUtility.StartPointperThread[threadIndex] - PWCUtility.PointStart_Process;
                    for (int index = beginpoint; index < indexlen + beginpoint; index++)
                    {
                        for (int ClusterIndex = 0; ClusterIndex < Dist.RunningPWC.Ncent; ClusterIndex++)
                        {
                            Dist.RunningPWC.Master_Epsilonalpha_k_[index][ClusterIndex] = Dist.RunningPWC.Epsilonalpha_k_[index][ClusterIndex];
                            Dist.RunningPWC.Best_Epsilonalpha_k_[index][ClusterIndex] = Dist.RunningPWC.Epsilonalpha_k_[index][ClusterIndex];
                            Dist.RunningPWC.Old_Epsilonalpha_k_[index][ClusterIndex] = Dist.RunningPWC.Epsilonalpha_k_[index][ClusterIndex];
                        }
                    }
                }
               );
            } catch (SuspendableException e) {
                PWCUtility.printAndThrowRuntimeException(e.getMessage());
            }

            double[] Save_C_k_ = new double[Dist.RunningPWC.Ncent];
            System.arraycopy(Dist.RunningPWC.C_k_, 0, Save_C_k_, 0, Dist.RunningPWC.Ncent);

			for (int ClusterToRefine = 0; ClusterToRefine < Dist.RunningPWC.Ncent; ClusterToRefine++)
			{
				if (Save_C_k_[ClusterToRefine] <= Program.ToosmalltoSplit)
				{
					continue;
				}
				//  Do EM calculation
				Dist.countAfterFixingClusterCount = 0; // Counts iterations after maximum cluster count reached
				Dist.needtocalculateMalpha_k_ = 1; // =0 Already Set; = 1 Calculate from EM (usual); = -1 Initialize

				//	Loop over EM calculations
				for (int EMLoop = 0; EMLoop < Program.EMlimit_Duplication; EMLoop++)
				{
					for (int ClusterIndex = 0; ClusterIndex < Dist.RunningPWC.Ncent; ClusterIndex++)
					{
						Dist.RunningPWC.Weight_k_[ClusterIndex] = 1.0;
					}
					Dist.RunningPWC.Weight_k_[ClusterToRefine] = 2.0;
					PairwiseThread();
					Extra_EMIterationCount++;

					//	Now see if we are done -- take results from Rank 0
					int convergence = 0;
					convergence = Dist.convergenceTest(Program.Epsi_max_change_Duplication);
					convergence = PWCUtility.synchronizeMPIVariable(convergence);

					if ((convergence > 0) || (EMLoop == Program.EMlimit_Duplication - 1))
					{
						String nextline = "     D)Special Iter " + (new Integer(EMLoop)).toString() + " T " + String.format("%1$.4E", Dist.RunningPWC.Temperature) + " PMLoop " + (new Integer(Dist.PMLoopUsed)).toString() + " PWHammy " + String.format("%1$.4E", Dist.RunningPWC.PairwiseHammy) + " Cluster " + (new Integer(ClusterToRefine)).toString() + " Sizes ";
						for (int ClusterIndex = 0; ClusterIndex < Dist.RunningPWC.Ncent; ClusterIndex++)
						{
							double tmp = -2.0 * Dist.RunningPWC.A_k_[ClusterIndex];
							nextline += String.format("%1$.1f", Dist.RunningPWC.C_k_[ClusterIndex]) + " (" + String.format("%1$.6f", Dist.RunningPWC.FreezingMeasure_k_[ClusterIndex]) + ") [" + String.format("%1$.1f", tmp) + "] ";
						}
						PWCUtility.SALSAPrint(1, nextline);
					}
					if (convergence > 0)
					{
						break;
					}

				} // End Loop over EMLoop

				for (int ClusterIndex = 0; ClusterIndex < Dist.RunningPWC.Ncent; ClusterIndex++)
				{
					vectorclass.eigenconverged[ClusterIndex] = 0;
					Dist.ClusterSelected[ClusterIndex] = 0;
				}
				Dist.ClusterSelected[ClusterToRefine] = 1;
				vc.getEigenvalues(3, Dist.RunningPWC.Malpha_k_, Dist.RunningPWC.A_k_, Dist.RunningPWC.Balpha_k_, Dist.RunningPWC.C_k_, Dist.RunningPWC.Ncent);
				double TrueMinimumEigenvalue = vectorclass.Eigenvalues_Pass0[ClusterToRefine] - vectorclass.Eigenvalues_Current[ClusterToRefine];

				PWCUtility.SALSAPrint(1, "       Cluster " + (new Integer(ClusterToRefine)).toString() + " Status " + (new Integer(vectorclass.eigenconverged[ClusterToRefine])).toString() + " Pass 1 " + String.format("%1$.4E", vectorclass.Eigenvalues_Current[ClusterToRefine]) + " Pass 0 " + String.format("%1$.4E", vectorclass.Eigenvalues_Pass0[ClusterToRefine]));

				if ((vectorclass.Eigenvalues_Pass0[ClusterToRefine] > 0.0) && (vectorclass.eigenconverged[ClusterToRefine] > 0))
				{
					if (TrueMinimumEigenvalue < Dist.RunningPWC.MinimumEigenvalue)
					{
						Dist.RunningPWC.MinimumEigenvalue = TrueMinimumEigenvalue;
						Dist.RunningPWC.ClustertoSplit = ClusterToRefine;
					}
				}

                // Parallel Saving of Epsilon
                // Note - parallel for
                int ClusterToRefineLoopVar = ClusterToRefine;
                try {
                    forallChunked(0, PWCUtility.ThreadCount-1, (threadIndex) ->
                    {
                        int indexlen = PWCUtility.PointsperThread[threadIndex];
                        int beginpoint = PWCUtility.StartPointperThread[threadIndex] - PWCUtility.PointStart_Process;
                        for (int ProcessPointIndex = beginpoint; ProcessPointIndex < indexlen + beginpoint; ProcessPointIndex++)
                        {
                            for (int ClusterIndex = 0; ClusterIndex < Dist.RunningPWC.Ncent; ClusterIndex++)
                            {
                                if (Dist.RunningPWC.ClustertoSplit == ClusterToRefineLoopVar)
                                {
                                    Dist.RunningPWC.Best_Epsilonalpha_k_[ProcessPointIndex][ClusterIndex] = Dist.RunningPWC.Epsilonalpha_k_[ProcessPointIndex][ClusterIndex];
                                }
                                    if (ClusterToRefineLoopVar < Dist.RunningPWC.Ncent - 1)
                                    {
                                        Dist.RunningPWC.Epsilonalpha_k_[ProcessPointIndex][ClusterIndex] = Dist.RunningPWC.Master_Epsilonalpha_k_[ProcessPointIndex][ClusterIndex];
                                        Dist.RunningPWC.Old_Epsilonalpha_k_[ProcessPointIndex][ClusterIndex] = Dist.RunningPWC.Master_Epsilonalpha_k_[ProcessPointIndex][ClusterIndex];
                                    }
                                    else
                                    {
                                        Dist.RunningPWC.Epsilonalpha_k_[ProcessPointIndex][ClusterIndex] = Dist.RunningPWC.Best_Epsilonalpha_k_[ProcessPointIndex][ClusterIndex];
                                        Dist.RunningPWC.Old_Epsilonalpha_k_[ProcessPointIndex][ClusterIndex] = Dist.RunningPWC.Best_Epsilonalpha_k_[ProcessPointIndex][ClusterIndex];
                                    }
                            }
                        }
                    }
                   );
                } catch (SuspendableException e) {
                    PWCUtility.printAndThrowRuntimeException(e.getMessage());
                }
            } // End case Methodology 3 or 4

			for (int ClusterIndex = 0; ClusterIndex < Dist.RunningPWC.Ncent; ClusterIndex++)
			{
				Dist.RunningPWC.Weight_k_[ClusterIndex] = 1.0;
			}

			if (Program.PerformEigenTest)
			{
				// Eigenvalue Test for Methodology 3
                // Note - parallel for
                try {
                    forallChunked(0, PWCUtility.ThreadCount-1, (threadIndex) ->
                    {
                        int indexlen = PWCUtility.PointsperThread[threadIndex];
                        int beginpoint = PWCUtility.StartPointperThread[threadIndex] - PWCUtility.PointStart_Process;
                        for (int ProcessPointIndex = beginpoint; ProcessPointIndex < indexlen + beginpoint; ProcessPointIndex++)
                        {
                            Dist.RunningPWC.Epsilonalpha_k_[ProcessPointIndex][Dist.RunningPWC.Ncent] = Dist.RunningPWC.Epsilonalpha_k_[ProcessPointIndex][Dist.RunningPWC.ClustertoSplit];
                            Dist.RunningPWC.Old_Epsilonalpha_k_[ProcessPointIndex][Dist.RunningPWC.Ncent] = Dist.RunningPWC.Epsilonalpha_k_[ProcessPointIndex][Dist.RunningPWC.ClustertoSplit];
                        }
                    }
                   );
                } catch (SuspendableException e) {
                    PWCUtility.printAndThrowRuntimeException(e.getMessage());
                }

                Dist.RunningPWC.Ncent++;
				for (int ClusterIndex = 0; ClusterIndex < Dist.RunningPWC.Ncent; ClusterIndex++)
				{
					Dist.RunningPWC.Weight_k_[ClusterIndex] = 1.0;
					Dist.ClusterSelected[ClusterIndex] = 0;
				}
				Dist.ClusterSelected[Dist.RunningPWC.ClustertoSplit] = 1;
				Dist.ClusterSelected[Dist.RunningPWC.Ncent - 1] = 1;

				Dist.PlaceforEigsforMultipleCluster = Dist.RunningPWC.ClustertoSplit;
				PairwiseThread();
				vc.getEigenvalues(4, Dist.RunningPWC.Malpha_k_, Dist.RunningPWC.A_k_, Dist.RunningPWC.Balpha_k_, Dist.RunningPWC.C_k_, Dist.RunningPWC.Ncent);
				double TrueMinimumEigenvalue1 = vectorclass.Eigenvalues_Pass0[Dist.PlaceforEigsforMultipleCluster] - vectorclass.Eigenvalues_Current[Dist.PlaceforEigsforMultipleCluster];

				String nextline = " C)Test Cluster " + Dist.RunningPWC.ClustertoSplit + " T " + String.format("%1$.4E", Dist.RunningPWC.Temperature) + " PMLoop " + (new Integer(Dist.PMLoopUsed)).toString() + " PWHammy " + String.format("%1$.4E", Dist.RunningPWC.PairwiseHammy) + (new Integer(Dist.PlaceforEigsforMultipleCluster)).toString() + " Status " + (new Integer(vectorclass.eigenconverged[Dist.PlaceforEigsforMultipleCluster])).toString() + " Pass 1 " + String.format("%1$.4E", vectorclass.Eigenvalues_Current[Dist.PlaceforEigsforMultipleCluster]) + " Pass 0 " + String.format("%1$.4E", vectorclass.Eigenvalues_Pass0[Dist.PlaceforEigsforMultipleCluster]) + " C ";
				for (int ClusterIndex = 0; ClusterIndex < Dist.RunningPWC.Ncent; ClusterIndex++)
				{
					nextline += String.format("%1$.1f", Dist.RunningPWC.C_k_[ClusterIndex]) + " ";
				}
				PWCUtility.SALSAPrint(1, nextline);

				Dist.RunningPWC.Ncent--;
			} // End Eigenvalue Test for Metholodology 3

		} // End case Methodology 3 or 4

		if (Dist.RunningPWC.ClustertoSplit < 0)
		{
			Dist.SplitFailures++;
			return false;
		}
		double eigtest = Program.MinEigtest * vectorclass.Eigenvalues_Pass0[Dist.RunningPWC.ClustertoSplit];
        return Dist.RunningPWC.MinimumEigenvalue <= eigtest;

    } // End shouldweSplit

	//  Do a split on the identified cluster
	//  Malpha_k_ is value of M indexed by (point,cluster)
	//  New clusters stored in last cluster and old position
	//  Cluster count is incremented by one
	//  ClustertoSplit is cluster number to split
	//	Embarassingly Parallel over Processes -- no MPI calls needed except in normalization
	public final void dothesplit() throws MPIException {
		double PerturbationNormFactor = 1.0;

		//  Calculate Normalization
		if (Program.PerturbationVehicle == 0)
		{
			GlobalReductions.FindDoubleSum SumoverShifts = new GlobalReductions.FindDoubleSum(PWCUtility.ThreadCount);
			GlobalReductions.FindDoubleSum ShiftNorm = new GlobalReductions.FindDoubleSum(PWCUtility.ThreadCount);
			GlobalReductions.FindDoubleSum EpsNorm = new GlobalReductions.FindDoubleSum(PWCUtility.ThreadCount);

            // Note - parallel for
            try {
                forallChunked(0, PWCUtility.ThreadCount-1, (threadIndex) ->
                {
                    int indexlen = PWCUtility.PointsperThread[threadIndex];
                    int beginpoint = PWCUtility.StartPointperThread[threadIndex] - PWCUtility.PointStart_Process;
                    for (int ProcessPointIndex = beginpoint; ProcessPointIndex < indexlen + beginpoint; ProcessPointIndex++)
                    {
                        double tmp = Dist.RunningPWC.Malpha_k_[ProcessPointIndex][Dist.RunningPWC.ClustertoSplit] * 0.5;
                        tmp = tmp * vectorclass.oldAx[ProcessPointIndex][Dist.RunningPWC.ClustertoSplit] / Dist.RunningPWC.Temperature;
                        if (Dist.RunningPWC.Ncent == 1)
                        {
                            tmp = Math.abs(tmp);
                        }
                            double tmp1 = vectorclass.oldAx[ProcessPointIndex][Dist.RunningPWC.ClustertoSplit] * Program.SplitPerturbationFactor;
                            double tmp2 = Dist.RunningPWC.Epsilonalpha_k_[ProcessPointIndex][Dist.RunningPWC.ClustertoSplit];
                            SumoverShifts.addapoint(threadIndex, tmp);
                            ShiftNorm.addapoint(threadIndex, tmp1 * tmp1);
                            EpsNorm.addapoint(threadIndex, tmp2 * tmp2);
                    }
                }
               );
            } catch (SuspendableException e) {
                PWCUtility.printAndThrowRuntimeException(e.getMessage());
            }
            SumoverShifts.sumoverthreadsandmpi();
			ShiftNorm.sumoverthreadsandmpi();
			EpsNorm.sumoverthreadsandmpi();

			Program.TestExpectedChange = Dist.RunningPWC.ClustertoSplit;
			Program.ExpectedChangeMethod = 0;
			if (Dist.RunningPWC.Ncent == 1)
			{
				Program.ExpectedChangeMethod = 1;
			}

			double CShiftNorm = SumoverShifts.Total;
			CShiftNorm = Math.abs(0.05 * Dist.RunningPWC.C_k_[Dist.RunningPWC.ClustertoSplit] / CShiftNorm);
			PerturbationNormFactor = 0.1 * Math.sqrt(EpsNorm.Total / ShiftNorm.Total);
			PWCUtility.SALSAPrint(0, "Perturb Cluster " + (new Integer(Program.TestExpectedChange)).toString() + " Sum Shift " + String.format("%1$.4E", SumoverShifts.Total) + " Shift Norm " + String.format("%1$.4E", ShiftNorm.Total) + " Eps Norm " + String.format("%1$.4E", EpsNorm.Total) + " CShift " + String.format("%1$.4E", CShiftNorm) + " Eps Shift " + String.format("%1$.4E", PerturbationNormFactor) + " Iter " + (new Integer(Dist.EMIterationCount)).toString());

			boolean testchange = CShiftNorm < PerturbationNormFactor;
			testchange = PWCUtility.synchronizeMPIVariable(testchange);
			Program.ExpectedChange = 0.05 * Dist.RunningPWC.C_k_[Dist.RunningPWC.ClustertoSplit] * Program.SplitPerturbationFactor;
			if (testchange)
			{
				PerturbationNormFactor = CShiftNorm;
			}
			else
			{
				if (Program.ExpectedChangeMethod == 0)
				{
					Program.ExpectedChangeMethod = 2;
				}
				Program.ExpectedChange *= PerturbationNormFactor / CShiftNorm;
			}
			if (Program.ExpectedChangeMethod == 0)
			{
				++Program.Countmethodzero;
			}
			if (Program.ExpectedChangeMethod == 2)
			{
				++Program.Countmethodtwo;
			}

			Program.PreviousC = Dist.RunningPWC.C_k_[Dist.RunningPWC.ClustertoSplit] * 0.5;
			if (Program.ExpectedChangeMethod != 1)
			{
				Program.deltaCoverC[4] += 2.0 * Math.abs(Program.ExpectedChange) / Program.PreviousC;
				Program.CountdeltaCoverC[4]++;
			}

			if (Program.ContinuousClustering)
			{
				GlobalReductions.FindDoubleSum NewC1 = new GlobalReductions.FindDoubleSum(PWCUtility.ThreadCount);
				GlobalReductions.FindDoubleSum NewC2 = new GlobalReductions.FindDoubleSum(PWCUtility.ThreadCount);

                // Note - parallel for
                final double PerturbationNormFactorLoopVar = PerturbationNormFactor;
                try {
                    forallChunked(0, PWCUtility.ThreadCount-1, (threadIndex) ->
                    {
                        int indexlen = PWCUtility.PointsperThread[threadIndex];
                        int beginpoint = PWCUtility.StartPointperThread[threadIndex] - PWCUtility.PointStart_Process;
                        for (int ProcessPointIndex = beginpoint; ProcessPointIndex < indexlen + beginpoint; ProcessPointIndex++)
                        {
                            double tmp = Dist.RunningPWC.Malpha_k_[ProcessPointIndex][Dist.RunningPWC.ClustertoSplit] * 0.5;
                            double perturb = vectorclass.oldAx[ProcessPointIndex][Dist.RunningPWC.ClustertoSplit] * Program.SplitPerturbationFactor * PerturbationNormFactorLoopVar;
                            double tmp2 = tmp * Math.exp(perturb / Dist.RunningPWC.Temperature);
                            double tmp1 = tmp * Math.exp(-perturb / Dist.RunningPWC.Temperature);
                            double NewBottom = 1.0 - 2.0 * tmp + tmp1 + tmp2;
                            NewC1.addapoint(threadIndex, tmp1 / NewBottom);
                            NewC2.addapoint(threadIndex, tmp2 / NewBottom);
                        }
                    }
                   );
                } catch (SuspendableException e) {
                    PWCUtility.printAndThrowRuntimeException(e.getMessage());
                }
                NewC1.sumoverthreadsandmpi();
				NewC2.sumoverthreadsandmpi();
				PWCUtility.SALSAPrint(0, "Real Expectation " + String.format("%1$.6E", NewC1.Total) + " " + String.format("%1$.6E", NewC2.Total) + " Temp " + String.format("%1$.4e", Dist.RunningPWC.Temperature));
				if (Program.ExpectedChangeMethod == 1)
				{
					Program.deltaCoverC[2] += Math.abs(NewC1.Total - NewC2.Total) / Program.PreviousC;
					Program.CountdeltaCoverC[2]++;
				}
				else
				{
					Program.deltaCoverC[3] += Math.abs(NewC1.Total - NewC2.Total) / Program.PreviousC;
					Program.CountdeltaCoverC[3]++;
				}
			}
		}

		Dist.RunningPWC.Ncent++;

		//  Parallel Section Splitting Cluster with delegate for action reading thread # from StartindexPort
        // Note - parallel for
        double PerturbationNormFactorLoopVar = PerturbationNormFactor;
        try {
            forallChunked(0, PWCUtility.ThreadCount-1, (threadIndex) ->
            {
                int indexlen = PWCUtility.PointsperThread[threadIndex];
                int beginpoint = PWCUtility.StartPointperThread[threadIndex] - PWCUtility.PointStart_Process;
                for (int ProcessPointIndex = beginpoint; ProcessPointIndex < indexlen + beginpoint; ProcessPointIndex++)
                {
                    if (Program.PerturbationVehicle == 1)
                    {
                        double newvalueofM = Dist.RunningPWC.Malpha_k_[ProcessPointIndex][Dist.RunningPWC.ClustertoSplit] * 0.5;
                        double fudge = 1.0 + Program.SplitPerturbationFactor;
                        if (ProcessPointIndex % 2 == 0)
                        {
                            fudge = 2.0 - fudge;
                        }
                            Dist.RunningPWC.Malpha_k_[ProcessPointIndex][Dist.RunningPWC.ClustertoSplit] = newvalueofM * fudge;
                            Dist.RunningPWC.Malpha_k_[ProcessPointIndex][Dist.RunningPWC.Ncent - 1] = newvalueofM * (2.0 - fudge);
                    }
                    else
                    {
                        double perturb = vectorclass.oldAx[ProcessPointIndex][Dist.RunningPWC.ClustertoSplit] * Program.SplitPerturbationFactor * PerturbationNormFactorLoopVar;
                        double original = Dist.RunningPWC.Epsilonalpha_k_[ProcessPointIndex][Dist.RunningPWC.ClustertoSplit];
                        Dist.RunningPWC.Epsilonalpha_k_[ProcessPointIndex][Dist.RunningPWC.ClustertoSplit] = original + perturb;
                        Dist.RunningPWC.Epsilonalpha_k_[ProcessPointIndex][Dist.RunningPWC.Ncent - 1] = original - perturb;
                    }
                }
            }
           );
        } catch (SuspendableException e) {
            PWCUtility.printAndThrowRuntimeException(e.getMessage());
        }

        if (Program.PerturbationVehicle == 1)
		{
			Dist.needtocalculateMalpha_k_ = 0;
		}
		else
		{
			// If continuous Clustering set P_k_ as Malpha_k_ formula involves P_k_
			// No need to perturb P_k_ -- they are just halved from last time around
			Dist.needtocalculateMalpha_k_ = 1;
			if (Program.ContinuousClustering)
			{
				Dist.RunningPWC.P_k_[Dist.RunningPWC.ClustertoSplit] = 0.5 * Dist.RunningPWC.P_k_[Dist.RunningPWC.ClustertoSplit];
				Dist.RunningPWC.P_k_[Dist.RunningPWC.Ncent - 1] = Dist.RunningPWC.P_k_[Dist.RunningPWC.ClustertoSplit];
			}
		}

		oldepsiset = 0;

	} // End dothesplit

	//	Find initial Temperature
	public static void initializeTemperature() throws MPIException {

		GlobalReductions.FindDoubleSum Find_avg1 = new GlobalReductions.FindDoubleSum(PWCUtility.ThreadCount);
		GlobalReductions.FindDoubleSum Find_avg2 = new GlobalReductions.FindDoubleSum(PWCUtility.ThreadCount);
		GlobalReductions.FindDoubleSum Find_avg3 = new GlobalReductions.FindDoubleSum(PWCUtility.ThreadCount);

        // Note - parallel for
        try {
            forallChunked(0, PWCUtility.ThreadCount-1, (threadIndex) ->
            {
                double DistanceSum = 0.0;
                double NumberSum = 0.0;
                double STDSum = 0.0;
                int indexlen = PWCUtility.PointsperThread[threadIndex];
                int beginpoint = PWCUtility.StartPointperThread[threadIndex] - PWCUtility.PointStart_Process;
                for (int ProcessPointIndex = beginpoint; ProcessPointIndex < indexlen + beginpoint; ProcessPointIndex++)
                {
                    int GlobalIndex = ProcessPointIndex + PWCUtility.PointStart_Process;
                    for (int PointIndex1 = 0; PointIndex1 < PWCUtility.PointCount_Global; PointIndex1++)
                    {
                        if (PointIndex1 == GlobalIndex)
                        {
                            continue;
                        }
						double placevalue = PWCUtility.PointDistances[ProcessPointIndex*PWCUtility.PointCount_Global+PointIndex1]* PWCUtility.INV_SHORT_MAX;
						DistanceSum += placevalue;
						STDSum += placevalue * placevalue;
						NumberSum += 1.0;
					}
                }
                Find_avg1.addapoint(threadIndex, DistanceSum);
                Find_avg2.addapoint(threadIndex, NumberSum);
                Find_avg3.addapoint(threadIndex, STDSum);
            }
           );
        } catch (SuspendableException e) {
            PWCUtility.printAndThrowRuntimeException(e.getMessage());
        }

        Find_avg1.sumoverthreadsandmpi();
		Find_avg2.sumoverthreadsandmpi();
		Find_avg3.sumoverthreadsandmpi();

		// Calculate global averages
		double[] avg_global = new double[3];

		avg_global[0] = Find_avg1.Total;
		avg_global[1] = Find_avg2.Total;
		avg_global[2] = Find_avg3.Total;

		//  Estimate of Initial Dist.Temperature is average distance
		//  Fudge factor of 2 over estimated critical temperature for first split
		double DistceMean = avg_global[0] / avg_global[1];
		double DistceSTD = Math.sqrt((avg_global[2] / avg_global[1]) - DistceMean * DistceMean);
		double EstimatedDimension = 2.0 * DistceMean * DistceMean / (DistceSTD * DistceSTD);
		double IndividualSigma = Math.sqrt(DistceMean / EstimatedDimension);
		Dist.Tinitial = 2.0 * DistceMean;

		PWCUtility.SALSAPrint(1, "Initial Temperature " + String.format("%1$.4f", Dist.Tinitial) + " # Distces " + (new Double(avg_global[1])).toString() + " Mean " + String.format("%1$.4f", DistceMean) + " STD " + String.format("%1$.4f", DistceSTD) + " Estimated Dimension " + String.format("%1$.3f", EstimatedDimension) + " IndividualSigma " + String.format("%1$.4f", IndividualSigma));
	} // End Initialize  Dist.InitializeTemperature


	// Calculate Distance Correlation Function
	public final void CorrelationCalculation() throws MPIException {
		int localNcent = Dist.RunningPWC.Ncent;
		Dist.DistanceCorrelation = new double[localNcent][localNcent];

// Summation of one row of Correlation Matrix
		GlobalReductions.FindVectorDoubleSum Find_CorrelationRow = new GlobalReductions.FindVectorDoubleSum(PWCUtility.ThreadCount, localNcent);

		// Loop over rows over Correlation Matrix
		for (int CorrelationrowIndex = 0; CorrelationrowIndex < Dist.RunningPWC.Ncent; CorrelationrowIndex++)
		{
			if (CorrelationrowIndex != 0)
			{
				Find_CorrelationRow.zero();
			}

            // Note - parallel for
            int CorrelationrowIndexLoopVar = CorrelationrowIndex;
            try {
                forallChunked(0, PWCUtility.ThreadCount-1, (threadIndex) ->
                {
                    //	Start Code setting partialsum_Correlation
                    double[] TempCorrel = new double[localNcent];
                    int indexlen = PWCUtility.PointsperThread[threadIndex];
                    int beginpoint = PWCUtility.StartPointperThread[threadIndex] - PWCUtility.PointStart_Process;
                    for (int ProcessPointIndex = beginpoint; ProcessPointIndex < indexlen + beginpoint; ProcessPointIndex++)
                    {
                        for (int ClusterIndex = 0; ClusterIndex < localNcent; ClusterIndex++)
                        {
                            TempCorrel[ClusterIndex] = Dist.RunningPWC.Malpha_k_[ProcessPointIndex][ClusterIndex] * Dist.RunningPWC.Balpha_k_[ProcessPointIndex][CorrelationrowIndexLoopVar];
                        }
                            Find_CorrelationRow.addAPoint(threadIndex, TempCorrel);
                    }
                }
               );
            } catch (SuspendableException e) {
                PWCUtility.printAndThrowRuntimeException(e.getMessage());
            }

            //  Form Correlation Row from sum over threads
			Find_CorrelationRow.sumOverThreadsAndMPI();
			for (int ClusterIndex = 0; ClusterIndex < Dist.RunningPWC.Ncent; ClusterIndex++)
			{
				Dist.DistanceCorrelation[CorrelationrowIndex][ClusterIndex] = Find_CorrelationRow.TotalVectorSum[ClusterIndex] / Dist.RunningPWC.C_k_[ClusterIndex];
			}

		} // End Loop over Rows of Correlation Matrix
	} // End Correlation Calculation

	// Return False if current final solution has a problem
	//  Return True if there is no problem or no way of fixing problem
	public static boolean CheckValidSolution() throws MPIException {
		boolean[] RemoveCluster = new boolean[Program.maxNcent];
		int Numbertoosmall = 0;
		double ProbabilitySum = 0.0;
		String documentit = "";
		for (int ClusterIndex = 0; ClusterIndex < Dist.RunningPWC.Ncent; ClusterIndex++)
		{
            RemoveCluster[ClusterIndex] =
                    Dist.RunningPWC.C_k_[ClusterIndex] < (double) Program.minimumclustercount - 0.5;
            RemoveCluster[ClusterIndex] = PWCUtility.synchronizeMPIVariable(RemoveCluster[ClusterIndex]);
			if (RemoveCluster[ClusterIndex])
			{
				++Numbertoosmall;
				documentit += " " + (new Integer(ClusterIndex)).toString() + " # " + String.format("%1$.1f", Dist.RunningPWC.C_k_[ClusterIndex]);
				continue;
			}
			ProbabilitySum += Dist.RunningPWC.P_k_[ClusterIndex];
		}
		if (Numbertoosmall == 0)
		{
			return true;
		}

		//  Need to change solution as one or more clusters too small
		PWCUtility.SALSAPrint(1, (new Integer(Numbertoosmall)).toString() + " Clusters Too Small " + documentit);
		if (!Dist.BestPWC.SolutionSet)
		{
			PWCUtility.SALSAPrint(1, "Solution Not changed as no best solution");
			return true;
		}
		if (Dist.BestPWC.Ncent != Dist.RunningPWC.Ncent)
		{
			PWCUtility.SALSAPrint(1, " Best Solution Not Used to restart as Number of Centers " + Dist.BestPWC.Ncent + " Different from Running Solution with " + Dist.RunningPWC.Ncent);
			return true;
		}

		ClusteringSolution.CopySolution(Dist.BestPWC, Dist.RunningPWC);
		Dist.RunningPWC.OldHammy = 0.0;
		Dist.RunningPWC.PairwiseHammy = 0.0;
		Dist.RunningPWC.Axset = false;
		Dist.RunningPWC.ClustertoSplit = -1;

		for (int clusterindex = 0; clusterindex < Dist.RunningPWC.Ncent; clusterindex++)
		{
			if (!RemoveCluster[clusterindex])
			{
				if (Program.ContinuousClustering)
				{
					Dist.RunningPWC.P_k_[clusterindex] = Dist.RunningPWC.P_k_[clusterindex] / ProbabilitySum;
				}
            }
		}
		int oldNcent = Dist.RunningPWC.Ncent;
		for (int clusterindex = oldNcent - 1; clusterindex >= 0; clusterindex--)
		{
			if (RemoveCluster[clusterindex])
			{
				ClusteringSolution.RemoveCluster(RunningPWC, clusterindex); // This both changes cluster count and shifts up clusters
			}
		}
		Dist.ActualMaxNcent = Dist.RunningPWC.Ncent;
		return false;

	} // End CheckValidSolution

} // End class dist
 // End namespace cluster
