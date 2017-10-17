package edu.indiana.soic.spidal.dapwc;

import com.google.common.base.Joiner;
import com.google.common.base.Strings;
import edu.rice.hj.api.SuspendableException;
import mpi.MPIException;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.charset.Charset;
import java.nio.file.*;
import java.util.Arrays;
import java.util.Date;
import java.util.regex.Pattern;

import static edu.rice.hj.Module1.forallChunked;

public class FindCenters
{
	// Sequence Space Measures
	//  SmallestDistanceMeans Find Sequence that has Minimum Mean Distance to all other sequences in cluster using Smith Waterman Gotoh distance
	//  Bucket Method defaults f= 0.15,0.4,0.75 for buckets 0 1 2
	//  Find value V(f) such that fraction f of all pairwise distances have size < V(f)
	//  Find Center with maximum number of points satisfying distance < V(f) for distances to center

	//  MDS Space
	//  SmallestMDSDistanceMeans Find Sequence that has Minimum Mean Distance to all other sequences in cluster using distances from mapping to 3D
	//  SmallestMDSCoG Take cluster points in 3D. Find geometric center with position gotten by averaging all positions in cluster. Find sequence nearest the "center of gravity"


	//  GroupIndex[GlobalPointIndex] is Cluster (group) number with label starting at StartPosition
	//  NumberofGroups is Number of clusters (families)
	//  PWCUtility.NumberofCenters is Number of Centers of each type to be found
	//  StartPosition is lowest cluster number 0 or 1
	//  CenterOutputFileName is file name of output file containing centers
	//  MDSInputFileName is file containg MDS  data
	//  Ignore group if number of points in group <= PWCUtility.NumberofCenters

	public static double[][] MDSvalues; // Array to store MDS positions if available

    public static void FindGroupCenters(int[] GroupIndex, int NumberofGroups, int StartPosition,
                                        String CenterOutputFileName, String MDSInputFileName,
                                        String LabelInputFileName) throws MPIException {
        if (PWCUtility.NumberofCenters <= 0) {
            return;
        }

        // for Bucket Method
        int NumberofBins = 1000;

        // For total evaluation
        int numberofcategories = 3 + PWCUtility.NumberofBuckets;
        int[] CombinedListIndices = new int[numberofcategories * PWCUtility.NumberofCenters];
        double[] CombinedListRatings = new double[numberofcategories * PWCUtility.NumberofCenters];
        double[] CombinedListRatingsCount = new double[numberofcategories * PWCUtility.NumberofCenters];
        boolean[] CombinedListSourceSeq = new boolean[numberofcategories * PWCUtility.NumberofCenters];
        boolean[] CombinedListSourceMDS = new boolean[numberofcategories * PWCUtility.NumberofCenters];
        double[] TopMeanDistance = new double[PWCUtility.NumberofCenters];
        double[] TopMDSMeanDistance = new double[PWCUtility.NumberofCenters];
        double[] TopMDSCoGDistance = new double[PWCUtility.NumberofCenters];
        double[] TopBucketDistance = new double[PWCUtility.NumberofCenters];

        // For Print Out
        int[] UtilityIndex = new int[numberofcategories * PWCUtility.NumberofCenters];
        String[] PointProperties = new String[numberofcategories * PWCUtility.NumberofCenters];

        double[][] GroupSmallestfromMinMeans =
                new double[NumberofGroups][]; // List means from of centers from min means
        int[][] GroupIndexfromMinMeans =
                new int[NumberofGroups][]; //  List of Global Point indices for centers for min means
        double[] GlobalGroupMean = new double[NumberofGroups]; // Mean distance in groups
        double[] GlobalGroupMax = new double[NumberofGroups]; // Maximum distances in groups
        int[] GroupCount = new int[NumberofGroups]; // Count of points in a group

        // For MDS need to read all points into all processes
        MDSvalues = new double[PWCUtility.PointCount_Global][];
        if (PWCUtility.addMDS > 0) {
            for (int GlobalPointIndex = 0; GlobalPointIndex < PWCUtility.PointCount_Global; GlobalPointIndex++) {
                MDSvalues[GlobalPointIndex] = new double[3];
            }
            ReadMDS(MDSInputFileName, MDSvalues, 0, PWCUtility.PointCount_Global);
        }

        // Read labels if they exist
        final boolean LabelsAvailable = LabelInputFileName.length() > 0;
        String[] SequenceLabels = new String[PWCUtility.PointCount_Global];
        int[] SequenceLengths = new int[PWCUtility.PointCount_Global];
        if (LabelInputFileName.length() > 0) {
            if (PWCUtility.MPI_Rank == 0) {
                ReadLabels(LabelInputFileName, SequenceLabels, SequenceLengths, 0, PWCUtility.PointCount_Global);
                int maxlen = 0;
                for (int looplabel = 0; looplabel < PWCUtility.PointCount_Global; looplabel++) {
                    maxlen = Math.max(maxlen, SequenceLabels[looplabel].length());
                }
                for (int looplabel = 0; looplabel < PWCUtility.PointCount_Global; looplabel++) {
                    for (int addblanks = SequenceLabels[looplabel].length(); addblanks < maxlen; addblanks++) {
                        SequenceLabels[looplabel] += " ";
                    }
                }
            }
            PWCUtility.StartSubTimer(PWCUtility.MPIBROADCASTTiming);
            // Note - MPI Call - Broadcast - int[]
            PWCUtility.mpiOps.broadcast(SequenceLengths, 0);
            PWCUtility.StopSubTimer(PWCUtility.MPIBROADCASTTiming);
        } else {
            PWCUtility.LengthCut1 = -1;
            PWCUtility.LengthCut2 = -1;
        }

        double[][] GroupMDSCoG = new double[NumberofGroups][]; // Mean Position(CoG) of Groups in MDS space
        double[] GlobalGroupMDSMean = new double[NumberofGroups]; // Mean of Group in MDS space

        int[][] GroupIndexfromMDSMinMeans = new int[NumberofGroups][]; // Indices from min mean in MDS space
        double[][] GroupSmallestfromMDSMinMeans = new double[NumberofGroups][]; // min means from min mean in MDS space

        int[][] GroupIndexfromMDSCoG = new int[NumberofGroups][]; // Indices nearest Center of Mass in MDS space
        double[][] GroupSmallestfromMDSCoG = new double[NumberofGroups][]; // Distances from Center of Mass in MDS space

        for (int group = 0; group < NumberofGroups; group++) {
            GroupCount[group] = 0;
            GlobalGroupMean[group] = 0.0;
            GlobalGroupMax[group] = 0.0;
            GroupSmallestfromMinMeans[group] = new double[PWCUtility.NumberofCenters];
            GroupIndexfromMinMeans[group] = new int[PWCUtility.NumberofCenters];

            GlobalGroupMDSMean[group] = 0.0;
            GroupMDSCoG[group] = new double[3];
            GroupIndexfromMDSMinMeans[group] = new int[PWCUtility.NumberofCenters];
            GroupSmallestfromMDSMinMeans[group] = new double[PWCUtility.NumberofCenters];
            GroupIndexfromMDSCoG[group] = new int[PWCUtility.NumberofCenters];
            GroupSmallestfromMDSCoG[group] = new double[PWCUtility.NumberofCenters];
        }

        for (int GlobalPointindex = 0; GlobalPointindex < PWCUtility.PointCount_Global; GlobalPointindex++) {
            int group = GroupIndex[GlobalPointindex] - StartPosition;
            ++GroupCount[group];
        }
        for (int group = 0; group < NumberofGroups; group++) {
            if (GroupCount[group] <= PWCUtility.NumberofCenters) {
                GroupCount[group] = 0;
            }
        }

        // Initialize accumulation data structures for each group
        GlobalReductions.FindDoubleMax[] FindGroupMax = new GlobalReductions.FindDoubleMax[NumberofGroups];
        GlobalReductions.FindMeanSigma[] FindGroupMeansigma = new GlobalReductions.FindMeanSigma[NumberofGroups];
        GlobalReductions.FindMeanSigma[] FindGroupMDSMeansigma = new GlobalReductions.FindMeanSigma[NumberofGroups];
        GlobalReductions.FindMeanSigma[][] FindGroupMDSCoG = new GlobalReductions.FindMeanSigma[NumberofGroups][];

        GlobalReductions.FindManyMinValuewithIndex[] FindGroupOriginalDataCenters_mean =
                new GlobalReductions.FindManyMinValuewithIndex[NumberofGroups];
        GlobalReductions.FindManyMinValuewithIndex[] FindGroupMDSCenters_mean =
                new GlobalReductions.FindManyMinValuewithIndex[NumberofGroups];
        GlobalReductions.FindManyMinValuewithIndex[] FindGroupMDSCenters_CoG =
                new GlobalReductions.FindManyMinValuewithIndex[NumberofGroups];

        for (int group = 0; group < NumberofGroups; group++) {
            if (GroupCount[group] == 0) {
                continue;
            }

            int Countlinks = 0;
            if (LabelsAvailable) {
                for (int GlobalPointIndex = 0; GlobalPointIndex < PWCUtility.PointCount_Global; GlobalPointIndex++) {
                    int group_Point = GroupIndex[GlobalPointIndex] - StartPosition;
                    if (group != group_Point) {
                        continue;
                    }
                    if ((SequenceLengths[GlobalPointIndex] > PWCUtility.LengthCut1) &&
                            (SequenceLengths[GlobalPointIndex] > PWCUtility.LengthCut2)) {
                        ++Countlinks;
                    }

                }
            } else {
                Countlinks = GroupCount[group];
            }
            if (Countlinks < PWCUtility.LinkCountinCenterFinding) {
                GroupCount[group] = 0;
                continue;
            }

            FindGroupMax[group] = new GlobalReductions.FindDoubleMax(PWCUtility.ThreadCount);
            FindGroupMeansigma[group] = new GlobalReductions.FindMeanSigma(PWCUtility.ThreadCount);
            FindGroupMDSMeansigma[group] = new GlobalReductions.FindMeanSigma(PWCUtility.ThreadCount);
            FindGroupMDSCoG[group] = new GlobalReductions.FindMeanSigma[3];
            for (int MDSIndex = 0; MDSIndex < 3; MDSIndex++) {
                FindGroupMDSCoG[group][MDSIndex] = new GlobalReductions.FindMeanSigma(PWCUtility.ThreadCount);
            }

            FindGroupOriginalDataCenters_mean[group] =
                    new GlobalReductions.FindManyMinValuewithIndex(PWCUtility.ThreadCount, PWCUtility.NumberofCenters);
            FindGroupMDSCenters_mean[group] =
                    new GlobalReductions.FindManyMinValuewithIndex(PWCUtility.ThreadCount, PWCUtility.NumberofCenters);
            FindGroupMDSCenters_CoG[group] =
                    new GlobalReductions.FindManyMinValuewithIndex(PWCUtility.ThreadCount, PWCUtility.NumberofCenters);
        }


        //  Loop over points GlobalPointIndex1 to find MDS means which are needed for next step
        if (PWCUtility.addMDS > 0) {
            // Note - parallel for
            try {
                forallChunked(0, PWCUtility.ThreadCount - 1, (threadIndex) -> {
                    int indexlen = PWCUtility.PointsperThread[threadIndex];
                    int beginpoint = PWCUtility.StartPointperThread[threadIndex] - PWCUtility.PointStart_Process;
                    for (int index = beginpoint; index < indexlen + beginpoint; index++) {
                        int GlobalPointIndex1 = index + PWCUtility.PointStart_Process;
                        int group1 = GroupIndex[GlobalPointIndex1] - StartPosition;
                        if ((group1 < 0) || (group1 >= NumberofGroups)) {
                            PWCUtility.printAndThrowRuntimeException(
                                    " Illegal group number " + (new Integer(group1)).toString() + " Point " +
                                            (new Integer(GlobalPointIndex1)).toString());
                        }
                        if (GroupCount[group1] <= 0) {
                            continue;
                        }
                        for (int GlobalPointIndex2 = 0; GlobalPointIndex2 < PWCUtility.PointCount_Global;
                             GlobalPointIndex2++) {
                            if (GlobalPointIndex1 == GlobalPointIndex2) {
                                continue;
                            }
                            int group2 = GroupIndex[GlobalPointIndex2] - StartPosition;
                            if ((group2 < 0) || (group2 >= NumberofGroups)) {
                                PWCUtility.printAndThrowRuntimeException(
                                        " Illegal group number " + (new Integer(group2)).toString() + " Point " +
                                                (new Integer(GlobalPointIndex2)).toString());
                            }
                            if (group1 != group2) {
                                continue;
                            }
                            double tmp = getMDSDistanceValue(GlobalPointIndex1, GlobalPointIndex2);
                            FindGroupMDSMeansigma[group1].addapoint(threadIndex, tmp);
                        }
                        for (int MDSindex = 0; MDSindex < 3; MDSindex++) {
                            FindGroupMDSCoG[group1][MDSindex]
                                    .addapoint(threadIndex, MDSvalues[GlobalPointIndex1][MDSindex]);
                        }
                    }
                });
            } catch (SuspendableException e) {
                PWCUtility.printAndThrowRuntimeException(e.getMessage());
            }

            for (int group = 0; group < NumberofGroups; group++) {
                if (GroupCount[group] == 0) {
                    continue;
                }
                FindGroupMDSMeansigma[group].sumoverthreadsandmpi();
                GlobalGroupMDSMean[group] = FindGroupMDSMeansigma[group].Totalmean;
                for (int MDSindex = 0; MDSindex < 3; MDSindex++) {
                    FindGroupMDSCoG[group][MDSindex].sumoverthreadsandmpi();
                    GroupMDSCoG[group][MDSindex] = FindGroupMDSCoG[group][MDSindex].Totalmean;
                }

            }
        } // End case where MDS values exist

        //  Loop over points GlobalPointIndex1
        // Note - parallel for
        try {
            forallChunked(0, PWCUtility.ThreadCount - 1, (threadIndex) -> {
                int indexlen = PWCUtility.PointsperThread[threadIndex];
                int beginpoint = PWCUtility.StartPointperThread[threadIndex] - PWCUtility.PointStart_Process;
                for (int index = beginpoint; index < indexlen + beginpoint; index++) {
                    int GlobalPointIndex1 = index + PWCUtility.PointStart_Process;
                    int group1 = GroupIndex[GlobalPointIndex1] - StartPosition;
                    if ((group1 < 0) || (group1 >= NumberofGroups)) {
                        PWCUtility.printAndThrowRuntimeException(
                                " Illegal group number " + (new Integer(group1)).toString() + " Point " +
                                        (new Integer(GlobalPointIndex1)).toString());
                    }
                    if (GroupCount[group1] <= 0) {
                        continue;
                    }
                    int LinkCutUsed = PWCUtility.LinkCountinCenterFinding;
                    LinkCutUsed = Math.min(LinkCutUsed, GroupCount[group1] - 20);
                    LinkCutUsed = Math.max(LinkCutUsed, 1);
                    double thispointmean = 0.0;
                    double thispointMDSmean = 0.0;
                    int Countlinks = 0;
                    for (int GlobalPointIndex2 = 0; GlobalPointIndex2 < PWCUtility.PointCount_Global; GlobalPointIndex2++) {
                        if (GlobalPointIndex1 == GlobalPointIndex2) {
                            continue;
                        }
                        int group2 = GroupIndex[GlobalPointIndex2] - StartPosition;
                        if ((group2 < 0) || (group2 >= NumberofGroups)) {
                            PWCUtility.printAndThrowRuntimeException(
                                    " Illegal group number " + (new Integer(group2)).toString() + " Point " +
                                            (new Integer(GlobalPointIndex2)).toString());
                        }
                        if (group1 != group2) {
                            continue;
                        }
                        double tmp = PWCUtility.PointDistances[index*PWCUtility.PointCount_Global+GlobalPointIndex2]* PWCUtility.INV_SHORT_MAX;
                        if (tmp > PWCUtility.MinimumDistanceCut) {
                            if (LabelsAvailable && (SequenceLengths[GlobalPointIndex2] > PWCUtility.LengthCut2)) {
                                ++Countlinks;
                                thispointmean += tmp;
                                FindGroupMax[group1].addapoint(threadIndex, tmp);
                                FindGroupMeansigma[group1].addapoint(threadIndex, tmp);
                            }
                        }
                        if (PWCUtility.addMDS > 0) {
                            tmp = getMDSDistanceValue(GlobalPointIndex1, GlobalPointIndex2);
                            thispointMDSmean += tmp;
                        }
                    }
                    if (Countlinks >= LinkCutUsed) {
                        FindGroupOriginalDataCenters_mean[group1]
                                .addAPoint(threadIndex, GlobalPointIndex1, thispointmean / Countlinks);
                    }
                    if (PWCUtility.addMDS > 0) {
                        double MDSmeandistance = getMDSDistancefromPoint(GlobalPointIndex1, GroupMDSCoG[group1]);
                        System.out.println("Line 315  MDSmeandistance " +  MDSmeandistance);

                        FindGroupMDSCenters_CoG[group1].addAPoint(threadIndex, GlobalPointIndex1, MDSmeandistance);
                        if (Countlinks > 0) {
                            FindGroupMDSCenters_mean[group1]
                                    .addAPoint(threadIndex, GlobalPointIndex1, thispointMDSmean / Countlinks);
                        }
                    }
                }
            });
        } catch (SuspendableException e) {
            PWCUtility.printAndThrowRuntimeException(e.getMessage());
        }

        // Finishup finding centers from sequence distance means and MDS CoG and means
        for (int group = 0; group < NumberofGroups; group++) {
            if (GroupCount[group] <= 0) {
                continue;
            }
            FindGroupMax[group].sumoverthreadsandmpi();
            FindGroupMeansigma[group].sumoverthreadsandmpi();
            GlobalGroupMean[group] = FindGroupMeansigma[group].Totalmean;
            GlobalGroupMax[group] = FindGroupMax[group].TotalMax;

            FindGroupOriginalDataCenters_mean[group].sumOverThreadsAndMPI();
            int TotalFound = (int) (FindGroupOriginalDataCenters_mean[group].TotalNumberofPoints + 0.01);
            PWCUtility.SALSAPrint(0, "Group " + (new Integer(group)).toString() + " Number in Mean " +
                    (new Integer(GroupCount[group])).toString() + " After Cuts " +
                    (new Integer(TotalFound)).toString());

            if (PWCUtility.addMDS > 0) {
                FindGroupMDSCenters_CoG[group].sumOverThreadsAndMPI();
                FindGroupMDSCenters_mean[group].sumOverThreadsAndMPI();
                System.out.println("Line 347  FindGroupMDSCenters_CoG[group] " +  FindGroupMDSCenters_CoG[group]);
            }

            for (int CenterIndex = 0; CenterIndex < PWCUtility.NumberofCenters; CenterIndex++) {
                GroupSmallestfromMinMeans[group][CenterIndex] =
                        FindGroupOriginalDataCenters_mean[group].OrderedMinValue[CenterIndex];
                GroupIndexfromMinMeans[group][CenterIndex] =
                        FindGroupOriginalDataCenters_mean[group].OrderedIndexValue[CenterIndex];

                if (PWCUtility.addMDS > 0) {
                    GroupIndexfromMDSMinMeans[group][CenterIndex] =
                            FindGroupMDSCenters_mean[group].OrderedIndexValue[CenterIndex]; // Indices from min mean
                            // in MDS space
                    GroupSmallestfromMDSMinMeans[group][CenterIndex] =
                            FindGroupMDSCenters_mean[group].OrderedMinValue[CenterIndex]; // min means from min mean
                            // in MDS space

                    GroupIndexfromMDSCoG[group][CenterIndex] =
                            FindGroupMDSCenters_CoG[group].OrderedIndexValue[CenterIndex]; // Indices nearest Center
                            // of Mass in MDS space
                    GroupSmallestfromMDSCoG[group][CenterIndex] =
                            FindGroupMDSCenters_CoG[group].OrderedMinValue[CenterIndex]; // Distances from group CoG
                }
            }
        }

        // Loop over groups -- output Min means MDS and bucket solutions -- also find bucket answers
        for (int group = 0; group < NumberofGroups; group++) {
            if (GroupCount[group] <= 0) {
                continue;
            }

            //  Overall header
            PWCUtility.SALSAPrint(1, "\n++++++++++++++++++\nFind Centers Group=" + (group + StartPosition) + " Count " +
                    (new Integer(GroupCount[group])).toString() + " Mean " +
                    String.format("%1$.4f", GlobalGroupMean[group]) + " MDS Mean " +
                    String.format("%1$.4f", GlobalGroupMDSMean[group]) + " Max " +
                    String.format("%1$.4f", GlobalGroupMax[group]));
            PWCUtility.SALSAPrint(1, "MDS Center of Gravity " + String.format("%1$.4f", GroupMDSCoG[group][0]) + " " +
                    String.format("%1$.4f", GroupMDSCoG[group][1]) + " " +
                    String.format("%1$.4f", GroupMDSCoG[group][2]));

            // Output conventional min Means with weight of 1
            int countindicesfound = 0;
            TopMeanDistance[0] = 0.0;

            //  Output min means from sequence space
            int NumMinMeansFound = 0;
            for (int CenterIndex1 = 0; CenterIndex1 < PWCUtility.NumberofCenters; CenterIndex1++) {
                int GlobalPointIndex1 = GroupIndexfromMinMeans[group][CenterIndex1];
                if (GlobalPointIndex1 < 0) {
                    PWCUtility.SALSAPrint(0, "Group " + (new Integer(group)).toString() +
                            " Center Index for Sequence Min-Mean " + (new Integer(CenterIndex1)).toString() +
                            " Undefined");
                    continue;
                }
                ++NumMinMeansFound;

                // Add into Global Rating List
                int useposition = countindicesfound;
                if (countindicesfound > 0) {
                    for (int loopindexlist = 0; loopindexlist < countindicesfound; loopindexlist++) {
                        if (CombinedListIndices[loopindexlist] != GlobalPointIndex1) {
                            continue;
                        }
                        useposition = loopindexlist;
                        break;
                    }
                }
                if (useposition < countindicesfound) {
                    CombinedListRatingsCount[useposition] += 1.0;
                    CombinedListRatings[useposition] += CenterIndex1;
                    CombinedListSourceSeq[useposition] = true;
                } else {
                    CombinedListIndices[useposition] = GlobalPointIndex1;
                    CombinedListRatingsCount[useposition] = 1.0;
                    CombinedListRatings[useposition] = CenterIndex1;
                    CombinedListSourceSeq[useposition] = true;
                    CombinedListSourceMDS[useposition] = false;
                    ++countindicesfound;
                }

                //  Calculate Averages and list each choice
                double avgmean = 0.0;
                double Avg_CoGcenters = 0.0;
                int owner = PWCParallelism.OwnerforThisPoint(GlobalPointIndex1);
                if (owner == PWCUtility.MPI_Rank) {
                    int procLocalRow = GlobalPointIndex1 - PWCUtility.PointStart_Process;
                    for (int CenterIndex2 = 0; CenterIndex2 < PWCUtility.NumberofCenters; CenterIndex2++) {
                        int GlobalPointIndex2 = GroupIndexfromMinMeans[group][CenterIndex2];
                        if (GlobalPointIndex2 < 0) {
                            PWCUtility.printAndThrowRuntimeException(
                                    "Means Group-2 " + (new Integer(group)).toString() + " Center " +
                                            (new Integer(CenterIndex2)).toString() + " Undefined");
                        }
                        double tmp;
                        if (CenterIndex1 != CenterIndex2) {
                            tmp = PWCUtility.PointDistances[procLocalRow*PWCUtility.PointCount_Global+GlobalPointIndex2]* PWCUtility.INV_SHORT_MAX;
                            if (CenterIndex1 == 0) {
                                TopMeanDistance[CenterIndex2] = tmp;
                            }
                            avgmean = avgmean + tmp;
                        }
                        if (PWCUtility.addMDS > 0) {
                            GlobalPointIndex2 = GroupIndexfromMDSCoG[group][CenterIndex2];
                            tmp = getMDSDistanceValue(GlobalPointIndex1, GlobalPointIndex2);
                            if (CenterIndex2 == 0) {
                                TopMDSCoGDistance[CenterIndex1] = tmp;
                            }
                            Avg_CoGcenters += tmp;
                        }
                    }
                    avgmean = avgmean / (PWCUtility.NumberofCenters - 1);
                    if (PWCUtility.addMDS > 0) {
                        Avg_CoGcenters = Avg_CoGcenters / PWCUtility.NumberofCenters;
                    }
                }
                if (PWCUtility.MPI_Size > 1) {
                    PWCUtility.StartSubTimer(PWCUtility.MPIBROADCASTTiming);
                    // Note - MPI Call - Broadcast - double
                    avgmean = PWCUtility.mpiOps.broadcast(avgmean, owner);
                    if (PWCUtility.addMDS > 0) {
                        // Note - MPI Call - Broadcast - double
                        Avg_CoGcenters = PWCUtility.mpiOps.broadcast(Avg_CoGcenters, owner);
                        // Note - MPI Call - Broadcast - double
                        TopMDSCoGDistance[CenterIndex1] =
                                PWCUtility.mpiOps.broadcast(TopMDSCoGDistance[CenterIndex1], owner);
                    }
                    if (CenterIndex1 == 0) {
                        // Note - MPI Call - Broadcast - double[]
                        PWCUtility.mpiOps.broadcast(TopMeanDistance, owner);
                    }
                    PWCUtility.StopSubTimer(PWCUtility.MPIBROADCASTTiming);
                }

                if (PWCUtility.MPI_Rank == 0) {
                    String Seqlabel = "";
                    if (LabelsAvailable) {
                        Seqlabel = " Sequence=" + SequenceLabels[GlobalPointIndex1] + " Length=" +
                                (new Integer(SequenceLengths[GlobalPointIndex1])).toString();
                    }
                    PointProperties[CenterIndex1] =
                            "Measure=" + String.format("%1$.4f", GroupSmallestfromMinMeans[group][CenterIndex1]) +
                                    " Method=SmallestDistanceMeans Group=" + (group + StartPosition) + Seqlabel +
                                    " Date=\"" + new Date() + "\"";
                    String MDSstring = "";
                    if (PWCUtility.addMDS > 0) {
                        MDSstring = " Avg MDS Distce CoG's " + String.format("%1$.4f", Avg_CoGcenters) + " to top CoG " +
                                String.format("%1$.4f", TopMDSCoGDistance[CenterIndex1]);
                    }
                    PWCUtility.SALSAPrint(1, "Mean " +
                            String.format("%1$.4f", GroupSmallestfromMinMeans[group][CenterIndex1]) + " Index " +
                            PWCUtility.PrintFixedInteger(GlobalPointIndex1, 6) + " Avg Distce-means " +
                            String.format("%1$.3E", avgmean) + " Distce to top mean " +
                            String.format("%1$.3E", TopMeanDistance[CenterIndex1]) + MDSstring + Seqlabel);
                }
            }

            if (PWCUtility.MPI_Rank == 0) {
                if (NumMinMeansFound > 0) {
                    System.arraycopy(GroupIndexfromMinMeans[group], 0, UtilityIndex, 0, NumMinMeansFound);
                    WritePointProperties(CenterOutputFileName, UtilityIndex, PointProperties, NumMinMeansFound, true);
                }
            }

            //  Output MDS results -- min mean method
            if (PWCUtility.addMDS > 0) {
                double scaleMDSratings = 1.0; // Use to scale ratings
                TopMDSMeanDistance[0] = 0.0;
                int NumMDSMeansFound = 0;
                for (int CenterIndex1 = 0; CenterIndex1 < PWCUtility.NumberofCenters; CenterIndex1++) {
                    int GlobalPointIndex1 = GroupIndexfromMDSMinMeans[group][CenterIndex1];
                    if (GlobalPointIndex1 < 0) {
                        PWCUtility.SALSAPrint(0, "MDS Means Group " + (new Integer(group)).toString() + " Center " +
                                (new Integer(CenterIndex1)).toString() + " Undefined");
                        continue;
                    }
                    NumMDSMeansFound++;
                    // Add into Global Rating List
                    int useposition = countindicesfound;
                    if (countindicesfound > 0) {
                        for (int loopindexlist = 0; loopindexlist < countindicesfound; loopindexlist++) {
                            if (CombinedListIndices[loopindexlist] != GlobalPointIndex1) {
                                continue;
                            }
                            useposition = loopindexlist;
                            break;
                        }
                    }
                    if (useposition < countindicesfound) {
                        CombinedListRatingsCount[useposition] += scaleMDSratings;
                        CombinedListRatings[useposition] += scaleMDSratings * CenterIndex1;
                        CombinedListSourceMDS[useposition] = true;
                    } else {
                        CombinedListIndices[useposition] = GlobalPointIndex1;
                        CombinedListRatingsCount[useposition] = scaleMDSratings;
                        CombinedListRatings[useposition] = scaleMDSratings * CenterIndex1;
                        CombinedListSourceSeq[useposition] = false;
                        CombinedListSourceMDS[useposition] = true;
                        ++countindicesfound;
                    }

                    double Avg_MeanMethod = 0.0;
                    double Avg_Internal = 0.0;
                    int owner = PWCParallelism.OwnerforThisPoint(GlobalPointIndex1);
                    if (owner == PWCUtility.MPI_Rank) {
                        int procLocalRow = GlobalPointIndex1 - PWCUtility.PointStart_Process;
                        for (int CenterIndex2 = 0; CenterIndex2 < PWCUtility.NumberofCenters; CenterIndex2++) {
                            int GlobalPointIndex2 = GroupIndexfromMDSMinMeans[group][CenterIndex2];
                            if (GlobalPointIndex2 < 0) {
                                PWCUtility.printAndThrowRuntimeException(
                                        "MDS Means Group-2 " + (new Integer(group)).toString() + " Center " +
                                                (new Integer(CenterIndex2)).toString() + " Undefined");
                            }
                            double tmp;
                            if (CenterIndex1 != CenterIndex2) {
                                tmp = getMDSDistanceValue(GlobalPointIndex1, GlobalPointIndex2);
                                if (CenterIndex1 == 0) {
                                    TopMDSMeanDistance[CenterIndex2] = tmp;
                                }
                                Avg_Internal = Avg_Internal + tmp;
                            }
                            GlobalPointIndex2 = GroupIndexfromMinMeans[group][CenterIndex2];
                            tmp = PWCUtility.PointDistances[procLocalRow*PWCUtility.PointCount_Global+GlobalPointIndex2]* PWCUtility.INV_SHORT_MAX;
                            Avg_MeanMethod = Avg_MeanMethod + tmp;
                            if (CenterIndex2 == 0) {
                                TopMeanDistance[CenterIndex1] = tmp;
                            }
                        }
                    }

                    if (PWCUtility.MPI_Size > 1) {
                        PWCUtility.StartSubTimer(PWCUtility.MPIBROADCASTTiming);
                        // Note - MPI Call - Broadcast - double
                        Avg_Internal = PWCUtility.mpiOps.broadcast(Avg_Internal, owner);
                        // Note - MPI Call - Broadcast - double
                        Avg_MeanMethod = PWCUtility.mpiOps.broadcast(Avg_MeanMethod, owner);
                        // Note - MPI Call - Broadcast - double
                        TopMeanDistance[CenterIndex1] =
                                PWCUtility.mpiOps.broadcast(TopMeanDistance[CenterIndex1], owner);
                        if (CenterIndex1 == 0) {
                            // Note - MPI Call - Broadcast - double[]
                            PWCUtility.mpiOps.broadcast(TopMeanDistance, owner);
                        }
                        PWCUtility.StopSubTimer(PWCUtility.MPIBROADCASTTiming);
                    }
                    Avg_MeanMethod = Avg_MeanMethod / PWCUtility.NumberofCenters;
                    Avg_Internal = Avg_Internal / (PWCUtility.NumberofCenters - 1);

                    if (PWCUtility.MPI_Rank == 0) {
                        String Seqlabel = "";
                        if (LabelsAvailable) {
                            Seqlabel = " Sequence=" + SequenceLabels[GlobalPointIndex1] + " Length=" +
                                    (new Integer(SequenceLengths[GlobalPointIndex1])).toString();
                        }
                        PointProperties[CenterIndex1] =
                                "Measure=" + String.format("%1$.4f", GroupSmallestfromMDSMinMeans[group][CenterIndex1]) +
                                        " Method=SmallestMDSDistanceMeans Group=" + (group + StartPosition) + Seqlabel +
                                        " Date=\"" + new Date() + "\"";
                        PWCUtility.SALSAPrint(1, "MDS Mean " +
                                String.format("%1$.4f", GroupSmallestfromMDSMinMeans[group][CenterIndex1]) + " Index " +
                                PWCUtility.PrintFixedInteger(GlobalPointIndex1, 6) + " Distce-MDSmeans " +
                                String.format("%1$.3E", Avg_Internal) + " Distce to MDS Mean top " +
                                String.format("%1$.3E", TopMDSMeanDistance[CenterIndex1]) + " Distce-OriginalMeans " +
                                String.format("%1$.3E", Avg_MeanMethod) + " Distce Original Mean Top " +
                                String.format("%1$.3E", TopMeanDistance[CenterIndex1]) + Seqlabel);
                    }
                }

                if (PWCUtility.MPI_Rank == 0) {
                    if (NumMDSMeansFound > 0) {
                        System.arraycopy(GroupIndexfromMDSMinMeans[group], 0, UtilityIndex, 0, NumMDSMeansFound);
                        WritePointProperties(CenterOutputFileName, UtilityIndex, PointProperties, NumMDSMeansFound,
                                             true);
                    }
                }

                //  Output MDS results -- CoG method
                TopMDSCoGDistance[0] = 0.0;
                for (int CenterIndex1 = 0; CenterIndex1 < PWCUtility.NumberofCenters; CenterIndex1++) {
                    int GlobalPointIndex1 = GroupIndexfromMDSCoG[group][CenterIndex1];
                    if (GlobalPointIndex1 < 0) {
                        PWCUtility.printAndThrowRuntimeException(
                                "MDS CoG Group-1 " + (new Integer(group)).toString() + " Center " +
                                        (new Integer(CenterIndex1)).toString() + " Undefined");
                    }

                    // Add into Global Rating List
                    int useposition = countindicesfound;
                    if (countindicesfound > 0) {
                        for (int loopindexlist = 0; loopindexlist < countindicesfound; loopindexlist++) {
                            if (CombinedListIndices[loopindexlist] != GlobalPointIndex1) {
                                continue;
                            }
                            useposition = loopindexlist;
                            break;
                        }
                    }
                    if (useposition < countindicesfound) {
                        CombinedListRatingsCount[useposition] += scaleMDSratings;
                        CombinedListRatings[useposition] += scaleMDSratings * CenterIndex1;
                        CombinedListSourceMDS[useposition] = true;
                    } else {
                        CombinedListIndices[useposition] = GlobalPointIndex1;
                        CombinedListRatingsCount[useposition] = scaleMDSratings;
                        CombinedListRatings[useposition] = scaleMDSratings * CenterIndex1;
                        CombinedListSourceSeq[useposition] = false;
                        CombinedListSourceMDS[useposition] = true;
                        ++countindicesfound;
                    }

                    double Avg_MeanMethod = 0.0;
                    int Num_MeanMethod = 0;
                    double Avg_Internal = 0.0;
                    int owner = PWCParallelism.OwnerforThisPoint(GlobalPointIndex1);
                    if (owner == PWCUtility.MPI_Rank) {
                        int procLocalRow = GlobalPointIndex1 - PWCUtility.PointStart_Process;
                        for (int CenterIndex2 = 0; CenterIndex2 < PWCUtility.NumberofCenters; CenterIndex2++) {
                            int GlobalPointIndex2 = GroupIndexfromMDSCoG[group][CenterIndex2];
                            if (GlobalPointIndex2 < 0) {
                                PWCUtility.printAndThrowRuntimeException(
                                        "MDS CoG Group-2 " + (new Integer(group)).toString() + " Center " +
                                                (new Integer(CenterIndex2)).toString() + " Undefined");
                            }
                            double tmp;
                            if (CenterIndex1 != CenterIndex2) {
                                tmp = getMDSDistanceValue(GlobalPointIndex1, GlobalPointIndex2);
                                if (CenterIndex1 == 0) {
                                    TopMDSCoGDistance[CenterIndex2] = tmp;
                                }
                                Avg_Internal = Avg_Internal + tmp;
                            }
                            GlobalPointIndex2 = GroupIndexfromMinMeans[group][CenterIndex2];
                            if (GlobalPointIndex2 < 0) {
                                if (CenterIndex2 == 0) {
                                    TopMeanDistance[CenterIndex1] = 0.0;
                                }
                            } else {
                                tmp = PWCUtility.PointDistances[procLocalRow*PWCUtility.PointCount_Global+GlobalPointIndex2]* PWCUtility.INV_SHORT_MAX;
                                Avg_MeanMethod = Avg_MeanMethod + tmp;
                                ++Num_MeanMethod;
                                if (CenterIndex2 == 0) {
                                    TopMeanDistance[CenterIndex1] = tmp;
                                }
                            }
                        }

                    }

                    if (PWCUtility.MPI_Size > 1) {
                        PWCUtility.StartSubTimer(PWCUtility.MPIBROADCASTTiming);
                        // Note - MPI Call - Broadcast - double
                        Avg_Internal = PWCUtility.mpiOps.broadcast(Avg_Internal, owner);
                        // Note - MPI Call - Broadcast - double
                        Avg_MeanMethod = PWCUtility.mpiOps.broadcast(Avg_MeanMethod, owner);
                        // Note - MPI Call - Broadcast - int
                        Num_MeanMethod = PWCUtility.mpiOps.broadcast(Num_MeanMethod, owner);
                        // Note - MPI Call - Broadcast - double
                        TopMeanDistance[CenterIndex1] =
                                PWCUtility.mpiOps.broadcast(TopMeanDistance[CenterIndex1], owner);
                        if (CenterIndex1 == 0) {
                            // Note - MPI Call - Broadcast - double[]
                            PWCUtility.mpiOps.broadcast(TopMDSCoGDistance, owner);
                        }
                        PWCUtility.StopSubTimer(PWCUtility.MPIBROADCASTTiming);
                    }
                    if (Num_MeanMethod > 0) {
                        Avg_MeanMethod = Avg_MeanMethod / Num_MeanMethod;
                    }
                    Avg_Internal = Avg_Internal / (PWCUtility.NumberofCenters - 1);

                    if (PWCUtility.MPI_Rank == 0) {
                        String Seqlabel = "";
                        if (LabelsAvailable) {
                            Seqlabel = " Sequence=" + SequenceLabels[GlobalPointIndex1] + " Length=" +
                                    (new Integer(SequenceLengths[GlobalPointIndex1])).toString();
                        }
                        PointProperties[CenterIndex1] =
                                "Measure=" + String.format("%1$.4f", GroupSmallestfromMDSCoG[group][CenterIndex1]) +
                                        " Method=SmallestMDSCoG Group=" + (group + StartPosition) + Seqlabel +
                                        " Date=\"" + new Date() + "\"";
                        PWCUtility.SALSAPrint(1, "MDS CoG  " +
                                String.format("%1$.4f", GroupSmallestfromMDSCoG[group][CenterIndex1]) + " Index " +
                                PWCUtility.PrintFixedInteger(GlobalPointIndex1, 6) + " Distce-MDS CoG " +
                                String.format("%1$.3E", Avg_Internal) + " Distce to MDS CoG top " +
                                String.format("%1$.3E", TopMDSCoGDistance[CenterIndex1]) + " Distce-OriginalMeans " +
                                String.format("%1$.3E", Avg_MeanMethod) + " Distce Original Mean Top " +
                                String.format("%1$.3E", TopMeanDistance[CenterIndex1]) + Seqlabel);
                    }
                }

                if (PWCUtility.MPI_Rank == 0) {
                    System.arraycopy(GroupIndexfromMDSCoG[group], 0, UtilityIndex, 0, PWCUtility.NumberofCenters);
                    WritePointProperties(CenterOutputFileName, UtilityIndex, PointProperties,
                                         PWCUtility.NumberofCenters, true);
                }

            } // End Output of MDS

            // Do Buckets *******************************************************************************************
            if (PWCUtility.NumberofBuckets > 0) {
                double[] BucketRadii = new double[PWCUtility.NumberofBuckets];

                //  First Histogram distance values so we can convert bucket fractions into radii
                double fudge = (double) NumberofBins / GlobalGroupMax[group];
                GlobalReductions.FindDoubleArraySum DistanceHistogramBinCounts =
                        new GlobalReductions.FindDoubleArraySum(PWCUtility.ThreadCount, NumberofBins);

                //  Loop over points selecting those in this group
                // Note - parallel for
                final int groupLoopVar = group;
                try {
                    forallChunked(0, PWCUtility.ThreadCount - 1, (threadIndex) -> {
                        int indexlen = PWCUtility.PointsperThread[threadIndex];
                        int beginpoint = PWCUtility.StartPointperThread[threadIndex] - PWCUtility.PointStart_Process;
                        DistanceHistogramBinCounts.startthread(threadIndex);
                        for (int index = beginpoint; index < indexlen + beginpoint; index++) {
                            int GlobalPointIndex1 = index + PWCUtility.PointStart_Process;
                            int group1 = GroupIndex[GlobalPointIndex1] - StartPosition;
                            if (group1 != groupLoopVar) {
                                continue;
                            }
                            for (int GlobalPointIndex2 = 0; GlobalPointIndex2 < PWCUtility.PointCount_Global;
                                 GlobalPointIndex2++) {
                                if (GlobalPointIndex1 == GlobalPointIndex2) {
                                    continue;
                                }
                                int group2 = GroupIndex[GlobalPointIndex2] - StartPosition;
                                if (group2 != groupLoopVar) {
                                    continue;
                                }
                                double tmp = PWCUtility.PointDistances[index*PWCUtility.PointCount_Global+GlobalPointIndex2]* PWCUtility.INV_SHORT_MAX;
                                if (tmp > PWCUtility.MinimumDistanceCut) {
                                    if (LabelsAvailable && (SequenceLengths[GlobalPointIndex2] > PWCUtility.LengthCut2)) {
                                        int itmp = (int) Math.floor(tmp * fudge);
                                        if (itmp >= NumberofBins) {
                                            itmp = NumberofBins - 1;
                                        }
                                        DistanceHistogramBinCounts.addapoint(threadIndex, itmp);
                                    }
                                }
                            }
                        }
                    });
                } catch (SuspendableException e) {
                    PWCUtility.printAndThrowRuntimeException(e.getMessage());
                }

                DistanceHistogramBinCounts.sumoverthreadsandmpi();
                double NumberinHistogram = DistanceHistogramBinCounts.TotalNumberofPoints;

                // Find Bucket Distance Cuts from histogram
                for (int BucketIndex = 0; BucketIndex < PWCUtility.NumberofBuckets; BucketIndex++) {
                    double target = (PWCUtility.BucketFractions[BucketIndex] * NumberinHistogram);
                    double TotalbinCounts = 0;
                    for (int BinIndex = 0; BinIndex < NumberofBins; BinIndex++) {
                        if (TotalbinCounts >= target) {
                            BucketRadii[BucketIndex] = (BinIndex + 0.5) / fudge;
                            break;
                        }
                        TotalbinCounts += DistanceHistogramBinCounts.TotalSum[BinIndex];
                    }
                }

                GlobalReductions.FindManyMinValuewithIndex[] FindCentersbybuckets =
                        new GlobalReductions.FindManyMinValuewithIndex[PWCUtility.NumberofBuckets];
                for (int BucketIndex = 0; BucketIndex < PWCUtility.NumberofBuckets; BucketIndex++) {
                    FindCentersbybuckets[BucketIndex] =
                            new GlobalReductions.FindManyMinValuewithIndex(PWCUtility.ThreadCount,
                                                                           PWCUtility.NumberofCenters);
                }

                //  Loop over points
                // Note - parallel for
                try {
                    forallChunked(0, PWCUtility.ThreadCount - 1, (threadIndex) -> {
                        int indexlen = PWCUtility.PointsperThread[threadIndex];
                        int beginpoint = PWCUtility.StartPointperThread[threadIndex] - PWCUtility.PointStart_Process;
                        for (int index = beginpoint; index < indexlen + beginpoint; index++) {
                            int GlobalPointIndex1 = index + PWCUtility.PointStart_Process;
                            int group1 = GroupIndex[GlobalPointIndex1] - StartPosition;
                            if (groupLoopVar != group1) {
                                continue;
                            }
                            double[] BucketCounts = new double[PWCUtility.NumberofBuckets];
                            for (int BucketIndex = 0; BucketIndex < PWCUtility.NumberofBuckets; BucketIndex++) {
                                BucketCounts[BucketIndex] = 0.0;
                            }
                            for (int GlobalPointIndex2 = 0; GlobalPointIndex2 < PWCUtility.PointCount_Global;
                                 GlobalPointIndex2++) {
                                if (GlobalPointIndex1 == GlobalPointIndex2) {
                                    continue;
                                }
                                int group2 = GroupIndex[GlobalPointIndex2] - StartPosition;
                                if (groupLoopVar != group2) {
                                    continue;
                                }
                                double tmp = PWCUtility.PointDistances[index*PWCUtility.PointCount_Global+GlobalPointIndex2]* PWCUtility.INV_SHORT_MAX;
                                if (tmp > PWCUtility.MinimumDistanceCut) {
                                    if (LabelsAvailable && (SequenceLengths[GlobalPointIndex2] > PWCUtility.LengthCut2)) {
                                        for (int BucketIndex = 0; BucketIndex < PWCUtility.NumberofBuckets; BucketIndex++) {
                                            if (tmp <= BucketRadii[BucketIndex]) {
                                                BucketCounts[BucketIndex] += 1.0;
                                            }
                                        }
                                    }
                                }
                            }
                            for (int BucketIndex = 0; BucketIndex < PWCUtility.NumberofBuckets; BucketIndex++) {
                                if (BucketCounts[BucketIndex] < 0.5) {
                                    continue;
                                }
                                double tmp = (double) GroupCount[groupLoopVar] - BucketCounts[BucketIndex];
                                FindCentersbybuckets[BucketIndex].addAPoint(threadIndex, GlobalPointIndex1, tmp);
                            }
                        }
                    });
                } catch (SuspendableException e) {
                    PWCUtility.printAndThrowRuntimeException(e.getMessage());
                }

                for (int BucketIndex = 0; BucketIndex < PWCUtility.NumberofBuckets; BucketIndex++) {
                    FindCentersbybuckets[BucketIndex].sumOverThreadsAndMPI();
                }

                // Output Bucket results
                for (int BucketIndex = 0; BucketIndex < PWCUtility.NumberofBuckets; BucketIndex++) {
                    PWCUtility.SALSAPrint(1, "\n----------------------\nGroup " + (group + StartPosition) + " Max " +
                            String.format("%1$.4f", GlobalGroupMax[group]) + " Count " +
                            (new Integer(GroupCount[group])).toString() + " Find Buckets CutOff Fraction " +
                            String.format("%1$.3f", PWCUtility.BucketFractions[BucketIndex]) + " Avg Bucket Distce " +
                            String.format("%1$.4f", BucketRadii[BucketIndex]));

                    TopBucketDistance[0] = 0.0;
                    int NumBucketEntriesFound = 0;
                    for (int CenterIndex1 = 0; CenterIndex1 < PWCUtility.NumberofCenters; CenterIndex1++) {
                        int GlobalPointIndex1 = FindCentersbybuckets[BucketIndex].OrderedIndexValue[CenterIndex1];
                        if (GlobalPointIndex1 < 0) {
                            PWCUtility.SALSAPrint(0, "Group " + (new Integer(group)).toString() + " Bucket " +
                                    (new Integer(BucketIndex)).toString() + " Center Index " +
                                    (new Integer(CenterIndex1)).toString() + " Undefined");
                            continue;
                        }
                        ++NumBucketEntriesFound;

                        // Add into Global Rating List
                        int useposition = countindicesfound;
                        if (countindicesfound > 0) {
                            for (int loopindexlist = 0; loopindexlist < countindicesfound; loopindexlist++) {
                                if (CombinedListIndices[loopindexlist] != GlobalPointIndex1) {
                                    continue;
                                }
                                useposition = loopindexlist;
                                break;
                            }
                        }
                        if (useposition < countindicesfound) {
                            CombinedListRatingsCount[useposition] += 1.0 / PWCUtility.NumberofBuckets;
                            CombinedListRatings[useposition] += CenterIndex1 / PWCUtility.NumberofBuckets;
                            CombinedListSourceSeq[useposition] = true;
                        } else {
                            CombinedListIndices[useposition] = GlobalPointIndex1;
                            CombinedListRatingsCount[useposition] = 1.0 / PWCUtility.NumberofBuckets;
                            CombinedListRatings[useposition] = CenterIndex1 / PWCUtility.NumberofBuckets;
                            CombinedListSourceSeq[useposition] = true;
                            CombinedListSourceMDS[useposition] = false;
                            ++countindicesfound;
                        }
                        double Avg_Internal = 0.0;
                        double Avg_MeanMethod = 0.0;
                        int Num_MeanMethod = 0;
                        int owner = PWCParallelism.OwnerforThisPoint(GlobalPointIndex1);
                        if (owner == PWCUtility.MPI_Rank) {
                            int procLocalRow = GlobalPointIndex1 - PWCUtility.PointStart_Process;
                            for (int CenterIndex2 = 0; CenterIndex2 < PWCUtility.NumberofCenters; CenterIndex2++) {
                                int GlobalPointIndex2 =
                                        FindCentersbybuckets[BucketIndex].OrderedIndexValue[CenterIndex2];
                                if (GlobalPointIndex2 < 0) {
                                    PWCUtility.printAndThrowRuntimeException(
                                            "Means Group-2 " + (new Integer(group)).toString() + " Center " +
                                                    (new Integer(CenterIndex2)).toString() + " Undefined");
                                }
                                double tmp;
                                if (CenterIndex1 != CenterIndex2) {
                                    tmp = PWCUtility.PointDistances[procLocalRow*PWCUtility.PointCount_Global+GlobalPointIndex2]* PWCUtility.INV_SHORT_MAX;
                                    if (CenterIndex1 == 0) {
                                        TopBucketDistance[CenterIndex2] = tmp;
                                    }
                                    Avg_Internal = Avg_Internal + tmp;
                                }

                                GlobalPointIndex2 = GroupIndexfromMinMeans[group][CenterIndex2];
                                if (GlobalPointIndex2 < 0) {
                                    if (CenterIndex2 == 0) {
                                        TopMeanDistance[CenterIndex1] = 0.0;
                                    }
                                } else {
                                    tmp = PWCUtility.PointDistances[procLocalRow*PWCUtility.PointCount_Global+GlobalPointIndex2]* PWCUtility.INV_SHORT_MAX;
                                    Avg_MeanMethod = Avg_MeanMethod + tmp;
                                    ++Num_MeanMethod;
                                    if (CenterIndex2 == 0) {
                                        TopMeanDistance[CenterIndex1] = tmp;
                                    }
                                }
                            }
                        }
                        if (PWCUtility.MPI_Size > 1) {
                            PWCUtility.StartSubTimer(PWCUtility.MPIBROADCASTTiming);
                            // Note - MPI Call - Broadcast - double
                            Avg_Internal = PWCUtility.mpiOps.broadcast(Avg_Internal, owner);
                            // Note - MPI Call - Broadcast - double
                            Avg_MeanMethod = PWCUtility.mpiOps.broadcast(Avg_MeanMethod, owner);
                            // Note - MPI Call - Broadcast - int
                            Num_MeanMethod = PWCUtility.mpiOps.broadcast(Num_MeanMethod, owner);
                            if (CenterIndex1 == 0) {
                                // Note - MPI Call - Broadcast - double[]
                                PWCUtility.mpiOps.broadcast(TopBucketDistance, owner);
                            }
                            // Note - MPI Call - Broadcast - double
                            TopMeanDistance[CenterIndex1] =
                                    PWCUtility.mpiOps.broadcast(TopMeanDistance[CenterIndex1], owner);
                            PWCUtility.StopSubTimer(PWCUtility.MPIBROADCASTTiming);
                        }
                        Avg_Internal = Avg_Internal / (PWCUtility.NumberofCenters - 1);
                        if (Num_MeanMethod > 0) {
                            Avg_MeanMethod = Avg_MeanMethod / Num_MeanMethod;
                        }
                        if (PWCUtility.MPI_Rank == 0) {
                            String Seqlabel = "";
                            if (LabelsAvailable) {
                                Seqlabel = " Sequence=" + SequenceLabels[GlobalPointIndex1] + " Length=" +
                                        (new Integer(SequenceLengths[GlobalPointIndex1])).toString();
                            }
                            int NumPointsinBucket = (int) Math.floor((double) GroupCount[group] -
                                                                             FindCentersbybuckets[BucketIndex].OrderedMinValue[CenterIndex1] +
                                                                             0.5);
                            PointProperties[CenterIndex1] = "Measure=" + (new Integer(NumPointsinBucket)).toString() +
                                    " Method=LargestCountsBucket-" + (new Integer(BucketIndex)).toString() + " Group=" +
                                    (group + StartPosition) + Seqlabel + " Date=\"" + new Date() + "\"";
                            PWCUtility.SALSAPrint(1, "Count " + String.format("%1$.2f", NumPointsinBucket*1.0) + " Index " +
                                    PWCUtility.PrintFixedInteger(GlobalPointIndex1, 6) + " Distce-Internal " +
                                    String.format("%1$.3E", Avg_Internal) + " Distce-means " +
                                    String.format("%1$.3E", Avg_MeanMethod) + " Distce to top mean " +
                                    String.format("%1$.3E", TopMeanDistance[CenterIndex1]) + " To Top Bucket " +
                                    String.format("%1$.3E", TopBucketDistance[CenterIndex1]) + Seqlabel);
                        }
                    }

                    // Output
                    if (PWCUtility.MPI_Rank == 0) {
                        if (NumBucketEntriesFound > 0) {
                            System.arraycopy(FindCentersbybuckets[BucketIndex].OrderedIndexValue, 0, UtilityIndex, 0,
                                             NumBucketEntriesFound);
                            WritePointProperties(CenterOutputFileName, UtilityIndex, PointProperties,
                                                 NumBucketEntriesFound, true);
                        }
                    }

                } // End BucketIndex Loop
            }

            //  Global Rating
            if (PWCUtility.MPI_Rank == 0) {
                PWCUtility.SALSAPrint(1, "\n---------------- Best Rated Points");
                double bestrating = -1.0;
                int bestposition = -1;
                for (int loopindexlist = 0; loopindexlist < countindicesfound; loopindexlist++) {
                    CombinedListRatings[loopindexlist] =
                            CombinedListRatings[loopindexlist] / CombinedListRatingsCount[loopindexlist];
                }

                int bestratedpoints = 0;
                for (int loopindexlist1 = 0; loopindexlist1 < countindicesfound; loopindexlist1++) {
                    bestrating = -1.0;
                    bestposition = -1;
                    for (int loopindexlist2 = 0; loopindexlist2 < countindicesfound; loopindexlist2++) {
                        if (CombinedListRatingsCount[loopindexlist2] < -0.5) {
                            continue;
                        }
                        if ((bestposition >= 0) && (CombinedListRatings[loopindexlist2] >= bestrating)) {
                            continue;
                        }
                        bestposition = loopindexlist2;
                        bestrating = CombinedListRatings[loopindexlist2];
                    }
                    if (CombinedListRatingsCount[bestposition] > 0.6) {
                        int GlobalPointIndex = CombinedListIndices[bestposition];
                        if (GlobalPointIndex < 0) {
                            continue;
                        }

                        String Seqlabel = "";
                        if (LabelsAvailable) {
                            Seqlabel = " Sequence=" + SequenceLabels[GlobalPointIndex] + " Length=" +
                                    (new Integer(SequenceLengths[GlobalPointIndex])).toString();
                        }
                        String source = "";
                        if (CombinedListSourceMDS[bestposition]) {
                            if (CombinedListSourceSeq[bestposition]) {
                                source = " Both";
                            } else {
                                source = " MDS ";
                            }
                        } else {
                            source = " Seq ";
                        }

                        UtilityIndex[bestratedpoints] = GlobalPointIndex;
                        PointProperties[bestratedpoints] =
                                "Measure=" + String.format("%1$.2f", CombinedListRatings[bestposition]) + " Count=" +
                                        String.format("%1$.1f", CombinedListRatingsCount[bestposition]) + " Source=" +
                                        source.trim() + " Method=OverallBest" + " Group=" + (group + StartPosition) +
                                        Seqlabel + " Date=\"" + new Date() + "\"";

                        ++bestratedpoints;
                        PWCUtility.SALSAPrint(1, "Index " + PWCUtility.PrintFixedInteger(GlobalPointIndex, 6) +
                                " Rating " + String.format("%1$.2f", CombinedListRatings[bestposition]) + " Count " +
                                String.format("%1$.1f", CombinedListRatingsCount[bestposition]) + source + Seqlabel);
                    }
                    CombinedListRatingsCount[bestposition] = -1.0;
                }
                if (bestratedpoints > 0) {
                    WritePointProperties(CenterOutputFileName, UtilityIndex, PointProperties, bestratedpoints, true);
                }
            } // End Output in rank 0 only of best rated points

        } // End loop over groups
    }

    public static double getMDSDistanceValue(int GlobalPointIndex1, int GlobalPointIndex2) {
        if (GlobalPointIndex1 == GlobalPointIndex2) {
            return 0.0;
        }
        double tmp0 = MDSvalues[GlobalPointIndex1][0] - MDSvalues[GlobalPointIndex2][0];
        double tmp1 = MDSvalues[GlobalPointIndex1][1] - MDSvalues[GlobalPointIndex2][1];
        double tmp2 = MDSvalues[GlobalPointIndex1][2] - MDSvalues[GlobalPointIndex2][2];
        return Math.sqrt(tmp0 * tmp0 + tmp1 * tmp1 + tmp2 * tmp2);
    }

    public static double getMDSDistancefromPoint(int GlobalPointIndex, double[] MDSPosition) {
        double tmp0 = MDSvalues[GlobalPointIndex][0] - MDSPosition[0];
        double tmp1 = MDSvalues[GlobalPointIndex][1] - MDSPosition[1];
        double tmp2 = MDSvalues[GlobalPointIndex][2] - MDSPosition[2];
        return Math.sqrt(tmp0 * tmp0 + tmp1 * tmp1 + tmp2 * tmp2);
    }

    // Read MDS values
    public static void ReadMDS(String fname, double[][] MDSvaluestoberead, int BeginPoint, int PointstoRead) {
        if (Strings.isNullOrEmpty(fname)) {
            PWCUtility.printAndThrowRuntimeException(new IllegalArgumentException(Constants.ERR_EMPTY_FILE_NAME));
        }

        boolean success = false;
        try (BufferedReader br = Files.newBufferedReader(Paths.get(fname), Charset.defaultCharset())) {
            String line;
            int count = 0;
            Pattern pattern = Pattern.compile("[\t ]");
            while ((line = br.readLine()) != null) {
                if (!Strings.isNullOrEmpty(line)) {
                    if (Strings.isNullOrEmpty(line)) {
                        continue; // continue on empty lines - "while" will break on null anyway;
                    }
                    String[] splits = pattern.split(line.trim());
                    int index = Integer.parseInt(splits[0]);
                    if (index < BeginPoint) {
                        continue;
                    }
                    if (index >= BeginPoint + PointstoRead) {
                        break;
                    }
                    MDSvaluestoberead[index - BeginPoint][0] = Double.parseDouble(splits[1]);
                    MDSvaluestoberead[index - BeginPoint][1] = Double.parseDouble(splits[2]);
                    MDSvaluestoberead[index - BeginPoint][2] = Double.parseDouble(splits[3]);
                    count++;
                }
            }
            if (count != PointstoRead) {
                PWCUtility.printAndThrowRuntimeException(
                        "Illegal count on MDS file " + (new Integer(count)).toString() + " " +
                                (new Integer(PointstoRead)).toString());
            }
            success = true;
        } catch (Exception e) {
            PWCUtility.printAndThrowRuntimeException("Failed reading MDS data" + e);
        }
        if (!success) {
            PWCUtility.printAndThrowRuntimeException("MDS File read error " + fname);
        }

    } // End ReadMDS

    // Read Point Labels
    public static void ReadLabels(String fname, String[] SequenceLabels, int[] SequenceLengths, int BeginPoint,
                                  int PointstoRead) {
        if (Strings.isNullOrEmpty(fname)) {
            PWCUtility.printAndThrowRuntimeException(new IllegalArgumentException(Constants.ERR_EMPTY_FILE_NAME));
        }

        boolean success = false;
        String line = "NotSet";
        int count = 0;
        try (BufferedReader br = Files.newBufferedReader(Paths.get(fname), Charset.defaultCharset())) {
            Pattern pattern = Pattern.compile("[\t ]");
            while ((line = br.readLine()) != null) {
                if (Strings.isNullOrEmpty(line)) {
                    continue; // continue on empty lines - "while" will break on null anyway;
                }
                String[] splits = pattern.split(line.trim());
                int index = Integer.parseInt(splits[0]);
                if (index < BeginPoint) {
                    continue;
                }
                if (index >= BeginPoint + PointstoRead) {
                    break;
                }
                if (splits.length == 3) {
                    SequenceLabels[index - BeginPoint] = splits[1];
                    SequenceLengths[index - BeginPoint] = Integer.parseInt(splits[2]);
                } else if (splits.length > 3) {
                    // Assume splits[1] to splits[splits.Length - 2] contribute to the name
                    SequenceLabels[index - BeginPoint] =
                            Joiner.on(" ").join(Arrays.copyOfRange(splits, 1, splits.length - 1));
                    SequenceLengths[index - BeginPoint] = Integer.parseInt(splits[splits.length - 1]);
                } else {
                    // Assume at least two splits where splits[0] is the index and splits[1] is the length
                    SequenceLabels[index - BeginPoint] = "";
                    SequenceLengths[index - BeginPoint] = Integer.parseInt(splits[1]);
                }
                count++;
            }
            if (count != PointstoRead) {
                PWCUtility.printAndThrowRuntimeException(
                        "Illegal count on Label file " + (new Integer(count)).toString() + " " +
                                (new Integer(PointstoRead)).toString());
            }
            success = true;
        } catch (Exception e) {
            PWCUtility.printAndThrowRuntimeException(
                    "Failed reading Label data: count " + count + " Line " + line + "\n" + e);
        }
        if (!success) {
            PWCUtility.printAndThrowRuntimeException(
                    "Label File read error: count " + count + " Line " + line + "\n" + fname);
        }

    } // End ReadSequenceLabels

    // Write General Point results into a file
    public static void WritePointProperties(String file, int[] PointNumbers, String[] labels, int dataPoints,
                                            boolean append) {
        if (PWCUtility.MPI_Rank != 0) {
            return;
        }

        if (Strings.isNullOrEmpty(file)) {
            PWCUtility.printAndThrowRuntimeException(new IllegalArgumentException(Constants.ERR_EMPTY_FILE_NAME));
        }

        Path filePath = Paths.get(file);
        OpenOption mode = append ? StandardOpenOption.APPEND : StandardOpenOption.CREATE;
        try (PrintWriter writer = new PrintWriter(Files.newBufferedWriter(filePath, Charset.defaultCharset(), mode),
                                                  true)) {
            for (int i = 0; i < dataPoints; i++) {
                String outputpoint = PWCUtility.LeftPrintFixedInteger(PointNumbers[i], 6);
                writer.println("PointNumber=" + outputpoint + " " + labels[i]);
            }
            writer.close();
        } catch (IOException e) {
            System.err.format("Failed writing cluster results due to I/O exception: %s%n", e);
        }
    } // End WritePointProperties
} // End class FindCenters