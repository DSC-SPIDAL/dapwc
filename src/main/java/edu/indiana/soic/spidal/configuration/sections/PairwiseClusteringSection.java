package edu.indiana.soic.spidal.configuration.sections;

import java.io.FileInputStream;
import java.io.IOException;
import java.util.Properties;

public class PairwiseClusteringSection {
    public PairwiseClusteringSection(String configurationFilePath) {
        Properties p = new Properties();
        try {
            p.load(new FileInputStream(configurationFilePath));
            ClusterFile = p.getProperty("ClusterFile","cluster.txt");
            DistanceMatrixFile = p.getProperty("DistanceMatrixFile", "distance.bin");
            AddMdsFile = p.getProperty("AddMdsFile","mds.txt");
            ClusterNumberFile = p.getProperty("ClusterNumberFile","cnums.txt");
            CenterPlotFile = p.getProperty("CenterPlotFile","cplot.txt");
            LabelFile = p.getProperty("LabelFile", "labels.txt");
            TimingFile = p.getProperty("TimingFile", "timings.txt");
            SummaryFile = p.getProperty("SummaryFile","summary.txt");

            NumberDataPoints = Integer.parseInt(p.getProperty("NumberDataPoints","-1"));
            ProcessingOption = Integer.parseInt(p.getProperty("ProcessingOption", "0"));
            TransformDimension = Integer.parseInt(p.getProperty("TransformDimension", "4"));

            MaxNcent = Integer.parseInt(p.getProperty("MaxNcent", "8"));
            SplitOrExpandIt = Integer.parseInt(p.getProperty("SplitOrExpandIt", "1"));
            MPIIOStrategy = Integer.parseInt(p.getProperty("MPIIOStrategy", "0"));
            TooSmallToSplit = Integer.parseInt(p.getProperty("TooSmallToSplit", "5"));

            MinEigTest = Double.parseDouble(p.getProperty("MinEigTest", "-0.01"));
            ConvergeIntermediateClusters = Boolean.parseBoolean(p.getProperty("ConvergeIntermediateClusters", "false"));
            WaitIterations = Integer.parseInt(p.getProperty("WaitIterations", "10"));
            EpsiMaxChange = Double.parseDouble(p.getProperty("EpsiMaxChange", "0.001"));

            InitialCoolingFactor = Double.parseDouble(p.getProperty("InitialCoolingFactor", "0.9"));
            FineCoolingFactor = Double.parseDouble(p.getProperty("FineCoolingFactor", "0.99"));
            EigenValueChange = Double.parseDouble(p.getProperty("EigenValueChange", "0.001"));
            EigenVectorChange = Double.parseDouble(p.getProperty("EigenVectorChange", "0.001"));

            IterationAtEnd = Integer.parseInt(p.getProperty("IterationAtEnd", "2000"));
            ConvergenceLoopLimit = Integer.parseInt(p.getProperty("ConvergenceLoopLimit", "2000"));
            FreezingLimit = Double.parseDouble(p.getProperty("FreezingLimit", "0.002"));
            PowerIterationLimit = Integer.parseInt(p.getProperty("PowerIterationLimit", "200"));
            ContinuousClustering = Boolean.parseBoolean(p.getProperty("ContinuousClustering", "false"));

            AddMds = Integer.parseInt(p.getProperty("AddMds", "1"));
            CenterPointsPerCenterTypeInOuput = Integer.parseInt(p.getProperty("CenterPointsPerCenterTypeInOuput", "3"));
            String BucketFractionsString = p.getProperty("BucketFractions", "0.15,0.4,0.75");
            String [] splits = BucketFractionsString.split(",");
            BucketFractions = new double[splits.length];
            for (int i = 0; i < splits.length; ++i){
                BucketFractions[i] = Double.parseDouble(splits[i]);
            }
            NumberOfCenters = Integer.parseInt(p.getProperty("NumberOfCenters", "8"));

            DebugPrintOption = Integer.parseInt(p.getProperty("DebugPrintOption", "1"));
            ConsoleDebugOutput = Boolean.parseBoolean(p.getProperty("ConsoleDebugOutput", "true"));

        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public String ClusterFile;
    public String DistanceMatrixFile;
    public String AddMdsFile;
    public String ClusterNumberFile;
    public String CenterPlotFile;
    public String LabelFile;
    public String TimingFile;
    public String SummaryFile;

    public int NumberDataPoints;
    public int ProcessingOption;
    public int TransformDimension;

    public int MaxNcent;
    public int SplitOrExpandIt;
    public int MPIIOStrategy;
    public int TooSmallToSplit;

    public double MinEigTest;
    public boolean ConvergeIntermediateClusters;
    public int WaitIterations;
    public double EpsiMaxChange;

    public double InitialCoolingFactor;
    public double FineCoolingFactor;
    public double EigenValueChange;
    public double EigenVectorChange;

    public int IterationAtEnd;
    public int ConvergenceLoopLimit;
    public double FreezingLimit;
    public int PowerIterationLimit;
    public boolean ContinuousClustering;

    public int AddMds;
    public double [] BucketFractions;
    public int NumberOfCenters;
    public int CenterPointsPerCenterTypeInOuput;

    public int NodeCount;
    public int ThreadCount;

    public int DebugPrintOption;
    public boolean ConsoleDebugOutput;

    public int dataTypeSize = 2; // 2 for short
    public boolean isBigEndian = false; // true for Java style binary data and false for C# style binary data

    public String getClusterFile() {
        return ClusterFile;
    }

    public void setClusterFile(String clusterFile) {
        ClusterFile = clusterFile;
    }

    public String getDistanceMatrixFile() {
        return DistanceMatrixFile;
    }

    public void setDistanceMatrixFile(String distanceMatrixFile) {
        DistanceMatrixFile = distanceMatrixFile;
    }

    public String getAddMdsFile() {
        return AddMdsFile;
    }

    public void setAddMdsFile(String addMdsFile) {
        AddMdsFile = addMdsFile;
    }

    public String getClusterNumberFile() {
        return ClusterNumberFile;
    }

    public void setClusterNumberFile(String clusterNumberFile) {
        ClusterNumberFile = clusterNumberFile;
    }

    public String getCenterPlotFile() {
        return CenterPlotFile;
    }

    public void setCenterPlotFile(String centerPlotFile) {
        CenterPlotFile = centerPlotFile;
    }

    public String getLabelFile() {
        return LabelFile;
    }

    public void setLabelFile(String labelFile) {
        LabelFile = labelFile;
    }

    public String getTimingFile() {
        return TimingFile;
    }

    public void setTimingFile(String timingFile) {
        TimingFile = timingFile;
    }

    public String getSummaryFile() {
        return SummaryFile;
    }

    public void setSummaryFile(String summaryFile) {
        SummaryFile = summaryFile;
    }

    public int getNumberDataPoints() {
        return NumberDataPoints;
    }

    public void setNumberDataPoints(int numberDataPoints) {
        NumberDataPoints = numberDataPoints;
    }

    public int getProcessingOption() {
        return ProcessingOption;
    }

    public void setProcessingOption(int processingOption) {
        ProcessingOption = processingOption;
    }

    public int getTransformDimension() {
        return TransformDimension;
    }

    public void setTransformDimension(int transformDimension) {
        TransformDimension = transformDimension;
    }

    public int getMaxNcent() {
        return MaxNcent;
    }

    public void setMaxNcent(int maxNcent) {
        MaxNcent = maxNcent;
    }

    public int getSplitOrExpandIt() {
        return SplitOrExpandIt;
    }

    public void setSplitOrExpandIt(int splitOrExpandIt) {
        SplitOrExpandIt = splitOrExpandIt;
    }

    public int getMPIIOStrategy() {
        return MPIIOStrategy;
    }

    public void setMPIIOStrategy(int MPIIOStrategy) {
        this.MPIIOStrategy = MPIIOStrategy;
    }

    public int getTooSmallToSplit() {
        return TooSmallToSplit;
    }

    public void setTooSmallToSplit(int tooSmallToSplit) {
        TooSmallToSplit = tooSmallToSplit;
    }

    public double getMinEigTest() {
        return MinEigTest;
    }

    public void setMinEigTest(double minEigTest) {
        MinEigTest = minEigTest;
    }

    public boolean isConvergeIntermediateClusters() {
        return ConvergeIntermediateClusters;
    }

    public void setConvergeIntermediateClusters(boolean convergeIntermediateClusters) {
        ConvergeIntermediateClusters = convergeIntermediateClusters;
    }

    public int getWaitIterations() {
        return WaitIterations;
    }

    public void setWaitIterations(int waitIterations) {
        WaitIterations = waitIterations;
    }

    public double getEpsiMaxChange() {
        return EpsiMaxChange;
    }

    public void setEpsiMaxChange(double epsiMaxChange) {
        EpsiMaxChange = epsiMaxChange;
    }

    public double getInitialCoolingFactor() {
        return InitialCoolingFactor;
    }

    public void setInitialCoolingFactor(double initialCoolingFactor) {
        InitialCoolingFactor = initialCoolingFactor;
    }

    public double getFineCoolingFactor() {
        return FineCoolingFactor;
    }

    public void setFineCoolingFactor(double fineCoolingFactor) {
        FineCoolingFactor = fineCoolingFactor;
    }

    public double getEigenValueChange() {
        return EigenValueChange;
    }

    public void setEigenValueChange(double eigenValueChange) {
        EigenValueChange = eigenValueChange;
    }

    public double getEigenVectorChange() {
        return EigenVectorChange;
    }

    public void setEigenVectorChange(double eigenVectorChange) {
        EigenVectorChange = eigenVectorChange;
    }

    public int getIterationAtEnd() {
        return IterationAtEnd;
    }

    public void setIterationAtEnd(int iterationAtEnd) {
        IterationAtEnd = iterationAtEnd;
    }

    public int getConvergenceLoopLimit() {
        return ConvergenceLoopLimit;
    }

    public void setConvergenceLoopLimit(int convergenceLoopLimit) {
        ConvergenceLoopLimit = convergenceLoopLimit;
    }

    public double getFreezingLimit() {
        return FreezingLimit;
    }

    public void setFreezingLimit(double freezingLimit) {
        FreezingLimit = freezingLimit;
    }

    public int getPowerIterationLimit() {
        return PowerIterationLimit;
    }

    public void setPowerIterationLimit(int powerIterationLimit) {
        PowerIterationLimit = powerIterationLimit;
    }

    public boolean isContinuousClustering() {
        return ContinuousClustering;
    }

    public void setContinuousClustering(boolean continuousClustering) {
        ContinuousClustering = continuousClustering;
    }

    public int getAddMds() {
        return AddMds;
    }

    public void setAddMds(int addMds) {
        AddMds = addMds;
    }

    public double[] getBucketFractions() {
        return BucketFractions;
    }

    public void setBucketFractions(double[] bucketFractions) {
        BucketFractions = bucketFractions;
    }

    public int getNumberOfCenters() {
        return NumberOfCenters;
    }

    public void setNumberOfCenters(int numberOfCenters) {
        NumberOfCenters = numberOfCenters;
    }

    public int getCenterPointsPerCenterTypeInOuput() {
        return CenterPointsPerCenterTypeInOuput;
    }

    public void setCenterPointsPerCenterTypeInOuput(int centerPointsPerCenterTypeInOuput) {
        CenterPointsPerCenterTypeInOuput = centerPointsPerCenterTypeInOuput;
    }

    public int getNodeCount() {
        return NodeCount;
    }

    public void setNodeCount(int nodeCount) {
        NodeCount = nodeCount;
    }

    public int getThreadCount() {
        return ThreadCount;
    }

    public void setThreadCount(int threadCount) {
        ThreadCount = threadCount;
    }

    public int getDebugPrintOption() {
        return DebugPrintOption;
    }

    public void setDebugPrintOption(int debugPrintOption) {
        DebugPrintOption = debugPrintOption;
    }

    public boolean isConsoleDebugOutput() {
        return ConsoleDebugOutput;
    }

    public void setConsoleDebugOutput(boolean consoleDebugOutput) {
        ConsoleDebugOutput = consoleDebugOutput;
    }

    public int getDataTypeSize() {
        return dataTypeSize;
    }

    public void setDataTypeSize(int dataTypeSize) {
        this.dataTypeSize = dataTypeSize;
    }

    public boolean isBigEndian() {
        return isBigEndian;
    }

    public void setBigEndian(boolean isBigEndian) {
        this.isBigEndian = isBigEndian;
    }
}
