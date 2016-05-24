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
import java.io.RandomAccessFile;
import java.nio.*;
import java.nio.channels.FileChannel;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
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

    public static Bytes mmapXReadBytes;
    public static ByteBuffer mmapXReadByteBuffer;
    public static Bytes mmapXInterSendBytes;
    public static ByteBuffer mmapXInterSendByteBuffer;
    public static Bytes mmapXInterRecvBytes;
    public static ByteBuffer mmapXInterRecvByteBuffer;
    public static Bytes mmapXWriteBytes;
    public static Bytes sendLock;
    public static Bytes recvLock;
    public static Bytes fullXBytes;
    public static ByteBuffer fullXByteBuffer;

    private static int LOCK = 0;
    private static int FLAG = Long.BYTES;

    public static Bytes mmapCollectiveReadBytes;
    public static ByteBuffer mmapCollectiveReadByteBuffer;
    public static Bytes mmapCollectiveWriteBytes;
    public static int mmapAllReduceChunkSizeInBytes;
    public static long mmapCollectiveWriteByteOffset;

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

        final String mmapXFname = machineName + ".mmapId." + mmapIdLocalToNode + ".mmapX.bin";
        /* Note, my send lock file name should match my successors recv lock file name */
        final String sendLockFname = machineName + ".mmapId." + mmapIdLocalToNode + ".lock." + worldProcRank + ".bin";
        final String recvLockFname = machineName + ".mmapId." + mmapIdLocalToNode + ".lock." + (worldProcRank != 0 ? worldProcRank-1 : worldProcsCount-1) +".bin";
        final String fullXFname = machineName + ".mmapId." + mmapIdLocalToNode +".fullX.bin";
        try (FileChannel mmapXFc = FileChannel.open(Paths.get(mmapScratchDir,
                mmapXFname),
                StandardOpenOption
                        .CREATE,
                StandardOpenOption.READ,
                StandardOpenOption
                        .WRITE);
             FileChannel fullXFc = FileChannel.open(Paths.get(mmapScratchDir,
                     fullXFname),
                     StandardOpenOption.CREATE,StandardOpenOption.WRITE,StandardOpenOption.READ)) {

            /* we need buffer space only for (mmapProcsCount - 1) to do internal
             * communication. Then we need 2 extra to do MPI send/receive at the
             * boundaries. These extra 2 are needed only when internal data is
             * represented as arrays and not direct buffers like in sendrecv<double[]>
             *
             * We'll use last chunk to do inter-process receive (of head mmap)
             * and the one before the last to do inter-process send (of tail mmap)*/

            int chunkSize = 2 * Integer.BYTES +
                    2 * Program.maxNcent * PWCUtility.PointCount_Largest *
                            Double.BYTES;
            int mmapXReadByteExtent = chunkSize * (mmapProcsCount + 1);

            long mmapXReadByteOffset = 0L;
            long
                    mmapXWriteByteOffset = 0L;
            int fullXByteExtent = globalRowCount * targetDimension * Double.BYTES;
            long fullXByteOffset = 0L;

            mmapXReadBytes = ByteBufferBytes.wrap(mmapXFc.map(
                    FileChannel.MapMode.READ_WRITE, mmapXReadByteOffset,
                    mmapXReadByteExtent));
            mmapXReadByteBuffer = mmapXReadBytes.sliceAsByteBuffer(
                    mmapXReadByteBuffer);

            mmapXReadBytes.position(0);
            mmapXWriteBytes = mmapXReadBytes.slice(mmapXWriteByteOffset,
                    mmapXReadByteExtent);

            /* Send recv buffers for mmap tail and mmap head */
            if (isMmapTail) {
                mmapXInterSendBytes = mmapXReadBytes.slice(mmapProcRank * chunkSize, chunkSize);
                mmapXInterSendByteBuffer = mmapXInterSendBytes.sliceAsByteBuffer(mmapXInterSendByteBuffer);
            }
            if (isMmapHead) {
                mmapXInterRecvBytes = mmapXReadBytes.slice(mmapProcsCount * chunkSize, chunkSize);
                mmapXInterRecvByteBuffer = mmapXInterRecvBytes.sliceAsByteBuffer(mmapXInterRecvByteBuffer);
            }


            /* Send receive locks */
            if (!isMmapTail){
                MappedByteBuffer mbb = createMMapLockBuffer(sendLockFname);
                sendLock = ByteBufferBytes.wrap(mbb);
            }

            if (!isMmapHead){
                MappedByteBuffer mbb = createMMapLockBuffer(recvLockFname);
                recvLock = ByteBufferBytes.wrap(mbb);
            }

            fullXBytes = ByteBufferBytes.wrap(fullXFc.map(FileChannel.MapMode
                            .READ_WRITE,
                    fullXByteOffset,
                    fullXByteExtent));
            fullXByteBuffer = fullXBytes.sliceAsByteBuffer(fullXByteBuffer);
        }

        /* Allocate memory maps for collective communications like AllReduce and Broadcast */
        final String mmapCollectiveFileName = machineName + ".mmapId." + mmapIdLocalToNode + ".mmapCollective.bin";
        try (FileChannel mmapCollectiveFc = FileChannel
                .open(Paths.get(mmapScratchDir, mmapCollectiveFileName),
                        StandardOpenOption.CREATE, StandardOpenOption.READ,
                        StandardOpenOption.WRITE)) {

            // See SharedMemoryCommunicatioNotes for more info on the number 1000
            mmapAllReduceChunkSizeInBytes = Math.max(Program.maxNcent, 1000)*Double.BYTES;
            int mmapCollectiveReadByteExtent = Math.max(
                    mmapProcsCount * mmapAllReduceChunkSizeInBytes,
                    (Math.max(
                            Program.maxNcent * Double.BYTES,
                            globalColCount*Integer.BYTES)));
            long mmapCollectiveReadByteOffset = 0L;

            mmapCollectiveWriteByteOffset = 0L;

            mmapCollectiveReadBytes = ByteBufferBytes.wrap(mmapCollectiveFc.map(
                    FileChannel.MapMode.READ_WRITE, mmapCollectiveReadByteOffset,
                    mmapCollectiveReadByteExtent));
            mmapCollectiveReadByteBuffer = mmapCollectiveReadBytes.sliceAsByteBuffer(
                    mmapCollectiveReadByteBuffer);

            mmapCollectiveReadBytes.position(0);
            mmapCollectiveWriteBytes = mmapCollectiveReadBytes.slice(mmapCollectiveWriteByteOffset,
                    mmapCollectiveReadByteExtent);
        }
    }

    private static MappedByteBuffer createMMapLockBuffer(String lockFname) throws IOException {
        File lockFile = new File(mmapScratchDir, lockFname);
        FileChannel fc = new RandomAccessFile(lockFile, "rw").getChannel();
        return fc.map(FileChannel.MapMode.READ_WRITE, 0, 64);
    }

    public static Iterator<MPISecPacket> allGather(MPISecPacket packet) throws MPIException {
        int offset = packet.getExtent() * mmapProcRank;
        packet.copyTo(offset, mmapXWriteBytes);

        worldProcsComm.barrier();
        if(isMmapLead){
            cgProcComm.allGather(mmapXReadByteBuffer, packet.getExtent()*mmapProcsCount, MPI.BYTE);
        }
        worldProcsComm.barrier();

        return new Iterator<MPISecPacket>() {
            int idx = 0;
            int extent = packet.getExtent();
            @Override
            public boolean hasNext() {
                mmapXReadBytes.position(idx*extent);
                return mmapXReadBytes.remaining() > extent;
            }

            @Override
            public MPISecPacket next() {
                /*packet.mapAt(idx*extent, extent, mmapXReadByteBuffer);*/
                MPISecPacket p = new MPISecPacket(extent);
                p.copyFrom(idx*extent, mmapXReadBytes);
                ++idx;
                return p;
            }
        };
    }

    /* to rank is my successor, from rank is my predecessor */
    public static void sendRecvPipeLine(MPISecPacket send, int to, int sendTag, MPISecPacket recv, int from, int recvTag) throws MPIException, InterruptedException {
        int extent = send.getExtent();
        if (extent != recv.getExtent()){
            PWCUtility.printAndThrowRuntimeException("Send and recv extents should match");
        }

        if (!isMmapTail) {
            sendLock.busyLockLong(LOCK);
            int offset = extent * mmapProcRank;
            send.copyTo(offset, mmapXWriteBytes);
            sendLock.writeBoolean(FLAG, true);
            sendLock.unlockLong(LOCK);
        }

        if (isMmapHead){
            /* mmap heads receive from tails of the previous mamps (or last mmap)*/
            worldProcsComm.recv(recv.getBuffer(), extent, MPI.BYTE, from, recvTag);
        } else if (isMmapTail) {
            worldProcsComm.send(send.getBuffer(), extent, MPI.BYTE, to, sendTag);
        }

        if (!isMmapHead){
            boolean dataReady = false;
            while (!dataReady) {
                recvLock.busyLockLong(LOCK);
                dataReady = recvLock.readBoolean(FLAG);
                if (dataReady){
                    /* Assumes receives are from previous ranks, i.e. the nature of pipeline */
                    int offset = extent*(mmapProcRank - 1);
                    recv.copyFrom(offset, mmapXWriteBytes);
                    recvLock.writeBoolean(FLAG, false);
                }
                recvLock.unlockLong(LOCK);
            }
        }
        /* Important functional barrier for correctness */
        worldProcsComm.barrier();
    }

    /* to rank is my successor, from rank is my predecessor */
    public static void sendRecvPipeLine(MPIPacket send, int to, int sendTag, MPIPacket recv, int from, int recvTag) throws MPIException, InterruptedException {
        int extent = send.getExtent();
        if (extent != recv.getExtent()){
            PWCUtility.printAndThrowRuntimeException("Send and recv extents should match");
        }

        if (!isMmapTail) {
            sendLock.busyLockLong(LOCK);
            int offset = extent * mmapProcRank;
            send.copyTo(offset, mmapXWriteBytes);
            sendLock.writeBoolean(FLAG, true);
            sendLock.unlockLong(LOCK);
        }

        if (isMmapHead){
            /* mmap heads receive from tails of the previous mamps (or last mmap)*/
            worldProcsComm.recv(recv.getBuffer(), extent, MPI.BYTE, from, recvTag);
        } else if (isMmapTail) {
            worldProcsComm.send(send.getBuffer(), extent, MPI.BYTE, to, sendTag);
        }

        if (!isMmapHead){
            boolean dataReady = false;
            while (!dataReady) {
                recvLock.busyLockLong(LOCK);
                dataReady = recvLock.readBoolean(FLAG);
                if (dataReady){
                    /* Assumes receives are from previous ranks, i.e. the nature of pipeline */
                    int offset = extent*(mmapProcRank - 1);
                    recv.copyFrom(offset, mmapXWriteBytes);
                    recvLock.writeBoolean(FLAG, false);
                }
                recvLock.unlockLong(LOCK);
            }
        }
        /* Important functional barrier for correctness */
        worldProcsComm.barrier();
    }

    /* to rank is my successor, from rank is my predecessor */
    public static void sendRecvPipeLine(double[] send, int to, int sendTag, double[] recv, int from, int recvTag) throws MPIException, InterruptedException {
        int size = send.length;
        int extent = size*Double.BYTES;

        if (extent != recv.length*Double.BYTES){
            PWCUtility.printAndThrowRuntimeException("Send and recv extents should match");
        }

        long offset = extent * mmapProcRank;
        if (!isMmapTail) {
            sendLock.busyLockLong(LOCK);
            copyToBuffer(send, mmapXWriteBytes, size, offset);
            sendLock.writeBoolean(FLAG, true);
            sendLock.unlockLong(LOCK);
        } else {
            // If tail, no need of locks
            copyToBuffer(send, mmapXInterSendBytes, size, 0);
        }

        if (isMmapHead){
            /* mmap heads receive from tails of the previous mamps (or last mmap)*/
            worldProcsComm.recv(mmapXInterRecvByteBuffer, extent, MPI.BYTE, from, recvTag);
        } else if (isMmapTail) {
            worldProcsComm.send(mmapXInterSendByteBuffer, extent, MPI.BYTE, to, sendTag);
        }

        if (!isMmapHead){
            boolean dataReady = false;
            while (!dataReady) {
                recvLock.busyLockLong(LOCK);
                dataReady = recvLock.readBoolean(FLAG);
                if (dataReady){
                    /* Assumes receives are from previous ranks, i.e. the nature of pipeline */
                    offset = extent*(mmapProcRank - 1);
                    copyFromBuffer(recv, mmapXWriteBytes, size, offset);
                    recvLock.writeBoolean(FLAG, false);
                }
                recvLock.unlockLong(LOCK);
            }
        } else {
            // If head, no need of locks
            copyFromBuffer(recv, mmapXInterRecvBytes, size, 0);
        }
        /* Important functional barrier for correctness */
        worldProcsComm.barrier();
    }

    private static void copyToBuffer(double[] array, Bytes to, long size, long offset) {
        for (int i = 0; i < size; ++i){
            to.writeDouble(offset+(i*Double.BYTES), array[i]);
        }
    }

    private static void copyFromBuffer(double[] array, Bytes from, long size, long offset) {
        for (int i = 0; i < size; ++i){
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

    public static void partialXAllGather() throws MPIException {
        cgProcComm.allGatherv(mmapXReadByteBuffer,
                cgProcsMmapXByteExtents[cgProcRank], MPI.BYTE,
                fullXByteBuffer, cgProcsMmapXByteExtents,
                cgProcsMmapXDisplas, MPI.BYTE);
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

