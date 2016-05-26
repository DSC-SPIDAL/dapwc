package edu.indiana.soic.spidal.dapwc;

import edu.indiana.soic.spidal.common.*;
import edu.indiana.soic.spidal.common.Range;
import edu.indiana.soic.spidal.common.RangePartitioner;
import mpi.Intracomm;
import mpi.MPI;
import mpi.MPIException;
import mpi.Op;
import net.openhft.lang.io.ByteBufferBytes;
import net.openhft.lang.io.Bytes;

import java.io.File;
import java.io.IOException;
import java.nio.*;
import java.nio.channels.FileChannel;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.regex.Pattern;
import java.util.stream.IntStream;

public class ParallelOps {
    public static String machineName;
    public static int nodeCount=1;
    public static int threadCount=1;

    public static int nodeId;

    public static Intracomm worldProcsComm;
    public static int worldProcRank;
    public static int worldProcsCount;
    public static int worldProcsPerNode;

    public static Intracomm mmapProcComm;
    // Number of memory mapped groups per process
    public static int mmapsPerNode;
    public static String mmapScratchDir;
    public static int worldProcRankLocalToNode;
    public static int mmapIdLocalToNode;
    public static int mmapProcRank;
    public static int mmapProcsCount;
    public static boolean isMmapLead;
    public static boolean isMmapHead = false;
    public static boolean isMmapTail = false;
    public static int[] mmapProcsWorldRanks;
    public static int mmapLeadWorldRank;
    public static int mmapLeadWorldRankLocalToNode;
    public static int mmapProcsRowCount;

    // mmap leaders form one communicating group and the others (followers)
    // belong to another communicating group.
    public static Intracomm cgProcComm;
    public static int cgProcRank;
    public static int cgProcsCount;
    public static int[] cgProcsMmapRowCounts;
    public static int[] cgProcsMmapXByteExtents;
    public static int[] cgProcsMmapXDisplas;

    public static String parallelPattern;
    public static Range[] procRowRanges;
    public static Range procRowRange;
    public static int procRowStartOffset;
    public static int procRowCount;

    public static Range[] threadRowRanges;
    public static int[] threadRowStartOffsets;
    public static int[] threadRowCounts;

    public static int globalColCount;

    // Buffers for MPI operations
    private static ByteBuffer statBuffer;
    private static DoubleBuffer doubleBuffer;
    private static IntBuffer intBuffer;
    public static LongBuffer threadsAndMPIBuffer;
    public static LongBuffer mpiOnlyBuffer;

    public static String fullXFname;
    public static Bytes fullXBytes;
    public static ByteBuffer fullXByteBuffer;
    public static double[] fullXArray;

    public static Bytes mmapCollPackReadBytes;
    public static ByteBuffer mmapCollPackReadByteBuffer;
    public static Bytes mmapCollPackWriteBytes;
    public static ByteBuffer mmapCollPackWriteByteBuffer;

    public static String mmapCollectiveFileName;
    public static Bytes mmapCollectiveReadBytes;
    public static ByteBuffer mmapCollectiveReadByteBuffer;
    public static Bytes mmapCollectiveWriteBytes;

    public static String ZmmapCollectiveFileName;
    public static Bytes ZmmapCollectiveReadBytes;
    public static ByteBuffer ZmmapCollectiveReadByteBuffer;

    public static MPISecPacket tmpMPISecPacket;

    public static void setupParallelism(String[] args) throws MPIException {
        MPI.Init(args);
        machineName = MPI.getProcessorName();

        /* Allocate basic buffers for communication */
        statBuffer = MPI.newByteBuffer(DoubleStatistics.extent);
        doubleBuffer = MPI.newDoubleBuffer(1);
        intBuffer = MPI.newIntBuffer(1);

        worldProcsComm = MPI.COMM_WORLD; //initializing MPI world communicator
        worldProcRank = worldProcsComm.getRank();
        worldProcsCount = worldProcsComm.getSize();

        /* Create communicating groups */
        worldProcsPerNode = worldProcsCount / nodeCount;
        boolean heterogeneous = (worldProcsPerNode * nodeCount) != worldProcsCount;
        if (heterogeneous) {
            PWCUtility.SALSAPrint(0,"Running in heterogeneous mode");
        }

        /* Logic to identify how many processes are within a node and
        *  the q and r values. These are used to processes to mmap groups
        *  within a node.
        *
        *  Important: the code assumes continues rank distribution
        *  within a node. */
        int[] qr = findQandR();
        int q = qr[0];
        int r = qr[1];

        // Memory mapped groups and communicating groups
        mmapIdLocalToNode =
                worldProcRankLocalToNode < r * (q + 1)
                        ? worldProcRankLocalToNode / (q + 1)
                        : (worldProcRankLocalToNode - r) / q;
        mmapProcsCount = worldProcRankLocalToNode < r*(q+1) ? q+1 : q;


        // Communicator for processes within a  memory map group
        mmapProcComm = worldProcsComm.split((nodeId*mmapsPerNode)+mmapIdLocalToNode, worldProcRank);
        mmapProcRank = mmapProcComm.getRank();

        isMmapLead = mmapProcRank == 0;
        isMmapHead = isMmapLead; // for chain calls
        isMmapTail = (mmapProcRank == mmapProcsCount - 1); // for chain calls
        mmapProcsWorldRanks = new int[mmapProcsCount];
        mmapLeadWorldRankLocalToNode =
                isMmapLead
                        ? worldProcRankLocalToNode
                        : (q * mmapIdLocalToNode + (mmapIdLocalToNode < r
                        ? mmapIdLocalToNode
                        : r));
        mmapLeadWorldRank = worldProcRank - (worldProcRankLocalToNode
                - mmapLeadWorldRankLocalToNode);
        // Assumes continues ranks within a node
        for (int i = 0; i < mmapProcsCount; ++i){
            mmapProcsWorldRanks[i] = mmapLeadWorldRank +i;
        }

        // Leaders talk, their color is 0
        // Followers will get a communicator of color 1,
        // but will make sure they don't talk ha ha :)
        cgProcComm = worldProcsComm.split(isMmapLead ? 0 : 1, worldProcRank);
        cgProcRank = cgProcComm.getRank();
        cgProcsCount = cgProcComm.getSize();

        parallelPattern =
                "---------------------------------------------------------\n"
                        + "Machine:" + machineName + ' ' + threadCount + 'x'
                        + worldProcsPerNode + 'x' + nodeCount;
        if (worldProcRank == 0) {
            PWCUtility.SALSAPrint(0, parallelPattern);
        }
    }

    private static int[] findQandR() throws MPIException {
        int q,r;
        String str = worldProcRank+ "@" +machineName +'#';
        intBuffer.put(0, str.length());
        worldProcsComm.allReduce(intBuffer, 1, MPI.INT, MPI.MAX);
        int maxLength = intBuffer.get(0);
        CharBuffer buffer = MPI.newCharBuffer(maxLength*worldProcsCount);
        buffer.position(maxLength*worldProcRank);
        buffer.put(str);
        for (int i = str.length(); i < maxLength; ++i){
            buffer.put(i, '~');
        }

        worldProcsComm.allGather(buffer, maxLength, MPI.CHAR);
        buffer.position(0);
        Pattern nodeSep = Pattern.compile("#~*");
        Pattern nameSep = Pattern.compile("@");
        String[] nodeSplits = nodeSep.split(buffer.toString());
        HashMap<String, Integer> nodeToProcCount = new HashMap<>();
        HashMap<Integer, String> rankToNode = new HashMap<>();
        String node;
        int rank;
        String[] splits;
        for(String s: nodeSplits){
            splits = nameSep.split(s);
            rank = Integer.parseInt(splits[0].trim());
            node = splits[1].trim();
            if (nodeToProcCount.containsKey(node)){
                nodeToProcCount.put(node, nodeToProcCount.get(node)+1);
            } else {
                nodeToProcCount.put(node, 1);
            }
            rankToNode.put(rank, node);
        }

        // The following logic assumes MPI ranks are continuous within a node
        String myNode = rankToNode.get(worldProcRank);
        HashSet<String> visited = new HashSet<>();
        int rankOffset=0;
        nodeId = 0;
        for (int i = 0; i < worldProcRank; ++i){
            node = rankToNode.get(i);
            if (visited.contains(node)) continue;
            visited.add(node);
            if (node.equals(myNode)) break;
            ++nodeId;
            rankOffset += nodeToProcCount.get(node);
        }
        worldProcRankLocalToNode = worldProcRank - rankOffset;
        final int procCountOnMyNode = nodeToProcCount.get(myNode);
        q = procCountOnMyNode / mmapsPerNode;
        r = procCountOnMyNode % mmapsPerNode;

        return new int[]{q,r};
    }

    public static void tempBreak() throws MPIException {
        if (worldProcRank ==0){
            MPI.Finalize();
            System.exit(0);
        } else {
            MPI.Finalize();
            System.exit(0);
        }
    }

    public static void tearDownParallelism() throws MPIException {
        // End MPI
        MPI.Finalize();

    }

    public static void setParallelDecomposition(int globalRowCount, int targetDimension)
            throws IOException, MPIException {
        //	First divide points among processes
        procRowRanges = edu.indiana.soic.spidal.common.RangePartitioner.partition(globalRowCount,
                worldProcsCount);
        Range rowRange = procRowRanges[worldProcRank]; // The range of points for this process

        procRowRange = rowRange;
        procRowStartOffset = rowRange.getStartIndex();
        procRowCount = rowRange.getLength();
        globalColCount = globalRowCount;

        // Next partition points per process among threads
        threadRowRanges = RangePartitioner.partition(procRowCount, threadCount);
        threadRowCounts = new int[threadCount];
        threadRowStartOffsets = new int[threadCount];
        IntStream.range(0, threadCount)
                .parallel()
                .forEach(threadIdx -> {
                    Range threadRowRange = threadRowRanges[threadIdx];
                    threadRowCounts[threadIdx] =
                            threadRowRange.getLength();
                    threadRowStartOffsets[threadIdx] =
                            threadRowRange.getStartIndex();
                });

        // Allocate timing buffers
        mpiOnlyBuffer = MPI.newLongBuffer(worldProcsCount);
        threadsAndMPIBuffer = MPI.newLongBuffer(worldProcsCount * threadCount);

        cgProcsMmapRowCounts = new int[cgProcsCount];
        cgProcsMmapXByteExtents = new int[cgProcsCount];
        cgProcsMmapXDisplas = new int[cgProcsCount];

        mmapProcsRowCount = IntStream.range(mmapLeadWorldRank,
                mmapLeadWorldRank + mmapProcsCount)
                .map(i -> procRowRanges[i].getLength())
                .sum();
        if (isMmapLead){
            cgProcsMmapRowCounts[cgProcRank] = mmapProcsRowCount;
            cgProcComm.allGather(cgProcsMmapRowCounts, 1, MPI.INT);
            for (int i = 0; i < cgProcsCount; ++i){
                cgProcsMmapXByteExtents[i] = cgProcsMmapRowCounts[i] * targetDimension * Double.BYTES;
            }

            cgProcsMmapXDisplas[0] = 0;
            System.arraycopy(cgProcsMmapXByteExtents, 0, cgProcsMmapXDisplas, 1, cgProcsCount - 1);
            Arrays.parallelPrefix(cgProcsMmapXDisplas, (m, n) -> m + n);
        }

        boolean status = new File(mmapScratchDir).mkdirs();

        String mmapCollPackReadFileName = machineName + ".mmapId." + mmapIdLocalToNode + ".mmapCollPackRead.bin";
        String mmapCollPackWriteFileName = machineName + ".mmapId." + mmapIdLocalToNode + ".mmapCollPackWrite.bin";
        try (FileChannel mmapCollPackReadFc = FileChannel
                .open(Paths.get(mmapScratchDir, mmapCollPackReadFileName),
                        StandardOpenOption.CREATE, StandardOpenOption.READ,
                        StandardOpenOption.WRITE);
             FileChannel mmapCollPackWriteFc = FileChannel
                     .open(Paths.get(mmapScratchDir, mmapCollPackWriteFileName),
                             StandardOpenOption.CREATE, StandardOpenOption.READ,
                             StandardOpenOption.WRITE)) {

            int packetSize = 2 * Integer.BYTES +
                    2 * Program.maxNcent * PWCUtility.PointCount_Largest *
                            Double.BYTES;
            int mmapCollPackReadByteExtent = packetSize * (worldProcsCount);
            int mmapCollPackWriteByteExtent = packetSize * (mmapProcsCount);

            mmapCollPackReadBytes = ByteBufferBytes.wrap(mmapCollPackReadFc.map(
                    FileChannel.MapMode.READ_WRITE, 0L,
                    mmapCollPackReadByteExtent));
            mmapCollPackReadByteBuffer = mmapCollPackReadBytes.sliceAsByteBuffer(

                    mmapCollPackReadByteBuffer);

            mmapCollPackWriteBytes = ByteBufferBytes.wrap(mmapCollPackWriteFc.map(
                    FileChannel.MapMode.READ_WRITE, 0L,
                    mmapCollPackWriteByteExtent));
            mmapCollPackWriteByteBuffer = mmapCollPackWriteBytes.sliceAsByteBuffer(

                    mmapCollPackWriteByteBuffer);

            // TODO - this might not be necessary
            if (isMmapLead){
                for (int i = 0; i < mmapCollPackReadByteExtent; ++i) {
                    mmapCollPackReadBytes.writeByte(i, 0);
                }
                for (int i = 0; i < mmapCollPackWriteByteExtent; ++i) {
                    mmapCollPackWriteBytes.writeByte(i, 0);
                }
            }
        }

        fullXFname = machineName + ".mmapId." + mmapIdLocalToNode +".fullX.bin";
        try (FileChannel fullXFc = FileChannel.open(Paths.get(mmapScratchDir,fullXFname),
                     StandardOpenOption.CREATE,StandardOpenOption.WRITE,StandardOpenOption.READ)) {

            int fullXByteExtent = Program.maxNcent * PWCUtility.PointCount_Largest * Double.BYTES * worldProcsCount;
            fullXArray = new double[fullXByteExtent / Double.BYTES];
            fullXBytes = ByteBufferBytes.wrap(fullXFc.map(FileChannel.MapMode.READ_WRITE, 0L, fullXByteExtent));
            fullXByteBuffer = fullXBytes.sliceAsByteBuffer(fullXByteBuffer);

            if (isMmapLead){
                for (int i = 0; i < fullXByteExtent; ++i)
                    fullXBytes.writeByte(i,0);
            }
        }

        /* Allocate memory maps for collective communications like AllReduce and Broadcast */
        mmapCollectiveFileName = machineName + ".mmapId." + mmapIdLocalToNode + ".mmapCollective.bin";
        try (FileChannel mmapCollectiveFc = FileChannel
                .open(Paths.get(mmapScratchDir, mmapCollectiveFileName),
                        StandardOpenOption.CREATE, StandardOpenOption.READ,
                        StandardOpenOption.WRITE)) {

            // See SharedMemoryCommunicatioNotes for more info on the number 1000
            int mmapAllReduceChunkSizeInBytes = Math.max(Program.maxNcent, 1000)*Double.BYTES;
            int mmapCollectiveReadByteExtent = Math.max(
                    mmapProcsCount * mmapAllReduceChunkSizeInBytes,
                    (Math.max(
                            Program.maxNcent * Double.BYTES,
                            globalColCount*Integer.BYTES)));

            long mmapCollectiveReadByteOffset = 0L;

            mmapCollectiveReadBytes = ByteBufferBytes.wrap(mmapCollectiveFc.map(
                    FileChannel.MapMode.READ_WRITE, mmapCollectiveReadByteOffset,
                    mmapCollectiveReadByteExtent));
            mmapCollectiveReadByteBuffer = mmapCollectiveReadBytes.sliceAsByteBuffer(
                    mmapCollectiveReadByteBuffer);

            mmapCollectiveReadBytes.position(0);
            mmapCollectiveWriteBytes = mmapCollectiveReadBytes.slice(0L,
                    mmapCollectiveReadByteExtent);

            if (isMmapLead){
                for (int i = 0; i < mmapCollectiveReadByteExtent; ++i)
                    mmapCollectiveReadBytes.writeByte(i,0);
            }
        }


        ZmmapCollectiveFileName = machineName + ".mmapId." + mmapIdLocalToNode + ".ZmmapCollective.bin";
        try (FileChannel ZmmapCollectiveFc = FileChannel
                .open(Paths.get(mmapScratchDir, ZmmapCollectiveFileName),
                        StandardOpenOption.CREATE, StandardOpenOption.READ,
                        StandardOpenOption.WRITE)) {

            int chunkSize = 2 * Integer.BYTES +
                    Program.maxNcent * PWCUtility.PointCount_Largest *
                            Double.BYTES;
            int ZmmapCollectiveReadByteExtent = chunkSize * (worldProcsCount);

            long mmapCollectiveReadByteOffset = 0L;

            ZmmapCollectiveReadBytes = ByteBufferBytes.wrap(ZmmapCollectiveFc.map(
                    FileChannel.MapMode.READ_WRITE, mmapCollectiveReadByteOffset,
                    ZmmapCollectiveReadByteExtent));
            ZmmapCollectiveReadByteBuffer = ZmmapCollectiveReadBytes.sliceAsByteBuffer(
                    ZmmapCollectiveReadByteBuffer);

            ZmmapCollectiveReadBytes.position(0);

            /*if (isMmapLead){
                for (int i = 0; i < ZmmapCollectiveReadByteExtent; ++i)
                    ZmmapCollectiveReadBytes.writeByte(i,0);
            }*/
        }

        // functional barrier
        worldProcsComm.barrier();
    }

    public static double[] allGather(double[] array) throws MPIException{
        int offset = array.length*Double.BYTES*mmapProcRank;
        copyToBuffer(array, fullXBytes, array.length, offset);
        worldProcsComm.barrier();
        if (isMmapLead){
            cgProcComm.allGather(fullXByteBuffer, array.length*mmapProcsCount, MPI.DOUBLE);
        }
        worldProcsComm.barrier();
        copyFromBuffer(fullXArray, fullXBytes, array.length*worldProcsCount, 0);
        return fullXArray;
    }

    /*public static void allGather(MPISecPacket packet, MPISecPacket[] packets) throws MPIException {
        int offset = packet.getExtent() * mmapProcRank;
        packet.copyTo(offset, mmapCollPackReadBytes);
        worldProcsComm.barrier();

        if(isMmapLead){
            cgProcComm.allGather(mmapCollPackReadByteBuffer, packet.getExtent()*mmapProcsCount, MPI.BYTE);
        }
        worldProcsComm.barrier();

        for (int i = 0; i < worldProcsCount; ++i){
            packets[i].copyFrom(i*packet.getExtent(), packet.getArrayLength(), mmapCollPackReadBytes);
        }
    }*/

    // TODO - debugs
    public static void allGather(MPISecPacket packet, MPISecPacket[] packets) throws MPIException {
        int offset = packet.getExtent() * mmapProcRank;
        packet.copyTo(offset, mmapCollPackWriteBytes);
        worldProcsComm.barrier();

        if(isMmapLead){
            cgProcComm.allGather(mmapCollPackWriteByteBuffer, packet.getExtent()*mmapProcsCount, MPI.BYTE, mmapCollPackReadByteBuffer, packet.getExtent()*mmapProcsCount, MPI.BYTE);
        }
        worldProcsComm.barrier();

        for (int i = 0; i < worldProcsCount; ++i){
            packets[i].copyFrom(i*packet.getExtent(), packet.getArrayLength(), mmapCollPackReadBytes);
        }

        /*cgProcComm.allGather(packet.getBuffer(), packet.getExtent(), MPI
                .BYTE, mmapCollPackReadByteBuffer, packet.getExtent(), MPI.BYTE);
        for (int i = 0; i < worldProcsCount; ++i){
            packets[i].copyFrom(i*packet.getExtent(), packet.getArrayLength()
                    , mmapCollPackReadByteBuffer);
        }*/

    }

    public static void allGather(MPIPacket packet, MPIPacket[] packets) throws MPIException {
        int offset = packet.getExtent() * mmapProcRank;
        /*packet.copyTo(offset, ZmmapCollectiveReadBytes);*/
        ZmmapCollectiveReadBytes.writeInt(2*Integer.BYTES*mmapProcRank, worldProcRank);
        ZmmapCollectiveReadBytes.writeInt(2*Integer.BYTES*mmapProcRank+Integer.BYTES, 53);

        worldProcsComm.barrier();

        // TODO - debugs
        /*if (worldProcRank == 176) {
            for (int i = 0; i < mmapProcsCount; ++i) {

*//*                packets[i].copyFrom(i *
                        packet.getExtent(), packet.getArrayLength(), ZmmapCollectiveReadByteBuffer);*//*

                System.out.println("++++ mmapProcsCount " + mmapProcsCount + " r " + ZmmapCollectiveReadBytes.readInt(i*packet.getExtent())
                        + " v " + ZmmapCollectiveReadBytes.readInt(i*packet.getExtent()+Integer.BYTES));
            }
        }*/
        worldProcsComm.barrier();

        if(isMmapLead){
            cgProcComm.allGather(ZmmapCollectiveReadByteBuffer, 2*Integer.BYTES*mmapProcsCount, MPI.BYTE);
        }
        worldProcsComm.barrier();

        // TODO - debugs
        if (worldProcRank == 176) {
            for (int i = 0; i < worldProcsCount; ++i) {

/*                packets[i].copyFrom(i *
                        packet.getExtent(), packet.getArrayLength(), ZmmapCollectiveReadByteBuffer);*/

               /* System.out.println("++++  r " + ZmmapCollectiveReadBytes.readInt(2*Integer.BYTES*i)
                        + " v " + ZmmapCollectiveReadBytes.readInt(2*Integer.BYTES*i+Integer.BYTES));*/
            }
//            System.out.println("DONE");
        }

        worldProcsComm.barrier();


        /*for (int i = 0; i < worldProcsCount; ++i){
            packets[i].copyFrom(i*packet.getExtent(), packet.getArrayLength(), ZmmapCollectiveReadBytes);
            // TODO - debugs
            if (worldProcRank == 176){
                System.out.println("**** number of points for " + i + " " + packets[i].getNumberOfPoints() + " frombuff " + ZmmapCollectiveReadBytes.readInt(i*packet.getExtent()+Integer.BYTES) + " i was sending " + packet.getNumberOfPoints());
            }
        }*/
        worldProcsComm.barrier();
    }

    private static void copyToBuffer(double[] array, Bytes to, long size, long offset) {
        for (int i = 0; i < size; ++i){
            to.writeDouble(offset+(i*Double.BYTES), array[i]);
        }
    }

    private static void copyFromBuffer(double[] array, Bytes from, long length, long offset) {
        for (int i = 0; i < length; ++i){
            array[i] = from.readDouble(offset+(i*Double.BYTES));
        }
    }

    private static void copyFromBuffer(double[] array, ByteBuffer from, int size, int offset) {
        for (int i = 0; i < size; ++i){
            array[i] = from.getDouble(offset+(i*Double.BYTES));
        }
    }

    private static void printInOrder(String str) throws MPIException {
        for (int i = 0; i < worldProcsCount; ++i){
            intBuffer.put(0,i);
            worldProcsComm.bcast(intBuffer, 1, MPI.INT, 0);
            int r = intBuffer.get(0);
            if (r == worldProcRank){
                System.out.println(str);
            }
        }
    }


    public static void broadcast(int[] values, int root) throws MPIException {
        int mmapLeaderCgProcCommRankOfRoot = 0;
        if (isMmapLead){
            // Let's find the cgProcComm rank of root's mmap leader
            mmapLeaderCgProcCommRankOfRoot = isRankWithinMmap(root) ? cgProcRank : 0;
            intBuffer.put(0, mmapLeaderCgProcCommRankOfRoot);
            cgProcComm.allReduce(intBuffer, 1, MPI.INT, MPI.SUM);
            mmapLeaderCgProcCommRankOfRoot = intBuffer.get(0);
        }

        if (root == worldProcRank){
            mmapCollectiveWriteBytes.position(0);
            for (int i = 0; i < values.length; ++i){
                mmapCollectiveWriteBytes.writeInt(i*Integer.BYTES, values[i]);
            }
        }
        worldProcsComm.barrier();

        if (ParallelOps.isMmapLead){
            cgProcComm.bcast(mmapCollectiveReadByteBuffer, values.length, MPI.INT, mmapLeaderCgProcCommRankOfRoot);
        }

        worldProcsComm.barrier();

        if (root != worldProcRank){
            for (int i = 0; i < values.length; ++i){
                values[i] = mmapCollectiveReadBytes.readInt(i*Integer.BYTES);
            }
        }
        /*worldProcsComm.barrier();*/
    }

    public static void broadcast(double[] values, int root) throws MPIException {
        int mmapLeaderCgProcCommRankOfRoot = 0;
        if (isMmapLead){
            // Let's find the cgProcComm rank of root's mmap leader
            mmapLeaderCgProcCommRankOfRoot = isRankWithinMmap(root) ? cgProcRank : 0;
            intBuffer.put(0, mmapLeaderCgProcCommRankOfRoot);
            cgProcComm.allReduce(intBuffer, 1, MPI.INT, MPI.SUM);
            mmapLeaderCgProcCommRankOfRoot = intBuffer.get(0);
        }

        if (root == worldProcRank){
            mmapCollectiveWriteBytes.position(0);
            for (int i = 0; i < values.length; ++i){
                mmapCollectiveWriteBytes.writeDouble(i*Double.BYTES, values[i]);
            }
        }
        worldProcsComm.barrier();

        if (ParallelOps.isMmapLead){
            cgProcComm.bcast(mmapCollectiveReadByteBuffer, values.length, MPI.DOUBLE, mmapLeaderCgProcCommRankOfRoot);
        }

        worldProcsComm.barrier();

        if (root != worldProcRank){
            for (int i = 0; i < values.length; ++i){
                values[i] = mmapCollectiveReadBytes.readDouble(i*Double.BYTES);
            }
        }
        /*worldProcsComm.barrier();*/
    }

    private static boolean isRankWithinMmap(int rank){
        return (mmapLeadWorldRank <= rank && rank <= (mmapLeadWorldRank+mmapProcsCount));
    }

    public static void allReduceSum(int[] values) throws MPIException {
        int idx;
        mmapCollectiveWriteBytes.position(0);
        for (int i = 0; i < values.length; ++i){
            idx = (i*mmapProcsCount)+mmapProcRank;
            mmapCollectiveWriteBytes.writeInt(idx*Integer.BYTES, values[i]);
        }
        // Important barrier here - as we need to make sure writes are done
        // to the mmap file.
        // It's sufficient to wait on ParallelOps.mmapProcComm,
        // but it's cleaner for timings if we wait on the whole world
        worldProcsComm.barrier();
        if (ParallelOps.isMmapLead) {
            // Node local reduction using shared memory maps
            int sum;
            int pos;
            for (int i = 0; i < values.length; ++i){
                sum = 0;
                pos = i*mmapProcsCount*Integer.BYTES;
                for (int j = 0; j < mmapProcsCount; ++j){
                    ParallelOps.mmapCollectiveReadBytes.position(pos);
                    sum += mmapCollectiveReadBytes.readInt();
                    pos += Integer.BYTES;
                }
                mmapCollectiveWriteBytes.writeInt(i*Integer.BYTES, sum);
            }

            // Leaders participate in MPI AllReduce
            cgProcComm.allReduce(mmapCollectiveReadByteBuffer, values.length, MPI.INT,MPI.SUM);
        }

        ParallelOps.worldProcsComm.barrier();
        ParallelOps.mmapCollectiveReadBytes.position(0);
        for (int i = 0; i < values.length; ++i){
            values[i] = ParallelOps.mmapCollectiveReadBytes.readInt();
        }
    }

    public static void allReduceSum(double[] values) throws MPIException {
        int idx;
        mmapCollectiveWriteBytes.position(0);
        for (int i = 0; i < values.length; ++i){
            idx = (i*mmapProcsCount)+mmapProcRank;
            mmapCollectiveWriteBytes.writeDouble(idx*Double.BYTES, values[i]);
        }
        // Important barrier here - as we need to make sure writes are done
        // to the mmap file.
        // It's sufficient to wait on ParallelOps.mmapProcComm,
        // but it's cleaner for timings if we wait on the whole world
        worldProcsComm.barrier();
        if (ParallelOps.isMmapLead) {
            // Node local reduction using shared memory maps
            double sum;
            int pos;
            for (int i = 0; i < values.length; ++i){
                sum = 0.0;
                pos = i*mmapProcsCount*Double.BYTES;
                for (int j = 0; j < mmapProcsCount; ++j){
                    ParallelOps.mmapCollectiveReadBytes.position(pos);
                    sum += mmapCollectiveReadBytes.readDouble();
                    pos += Double.BYTES;
                }
                mmapCollectiveWriteBytes.writeDouble(i*Double.BYTES, sum);
            }

            // Leaders participate in MPI AllReduce
            cgProcComm.allReduce(mmapCollectiveReadByteBuffer, values.length, MPI.DOUBLE,MPI.SUM);
        }

        ParallelOps.worldProcsComm.barrier();
        ParallelOps.mmapCollectiveReadBytes.position(0);
        for (int i = 0; i < values.length; ++i){
            values[i] = ParallelOps.mmapCollectiveReadBytes.readDouble();
        }
    }

    public static DoubleStatistics allReduce(DoubleStatistics stat)
            throws MPIException {
        stat.addToBuffer(statBuffer, 0);
        worldProcsComm.allReduce(statBuffer, DoubleStatistics.extent, MPI.BYTE,
                DoubleStatistics.reduceSummaries());
        return DoubleStatistics.getFromBuffer(statBuffer, 0);
    }

    public static double allReduce(double value) throws MPIException{
        doubleBuffer.put(0, value);
        worldProcsComm.allReduce(doubleBuffer, 1, MPI.DOUBLE, MPI.SUM);
        return doubleBuffer.get(0);
    }

    public static int allReduce(int value) throws MPIException{
        intBuffer.put(0, value);
        worldProcsComm.allReduce(intBuffer, 1, MPI.INT, MPI.SUM);
        return intBuffer.get(0);
    }

    public static void partialSAllReduce(Op op) throws MPIException{
        cgProcComm.allReduce(mmapCollectiveReadByteBuffer, 1, MPI.DOUBLE,op);
    }

    public static void broadcast(ByteBuffer buffer, int extent, int root)
            throws MPIException {
        worldProcsComm.bcast(buffer, extent, MPI.BYTE, root);
    }

    public static void gather(LongBuffer buffer, int count, int root)
            throws MPIException {
        worldProcsComm.gather(buffer, count, MPI.LONG, root);
    }
}

