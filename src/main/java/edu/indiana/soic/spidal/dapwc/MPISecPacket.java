package edu.indiana.soic.spidal.dapwc;

import mpi.MPI;
import net.openhft.lang.io.Bytes;

import java.io.Serializable;
import java.nio.ByteBuffer;
import java.nio.IntBuffer;

public class MPISecPacket implements Serializable
{
    private int firstPointOffset = 0;
    private int numberOfPointsOffset = Integer.BYTES;
    private int mArrayOffset = 2*Integer.BYTES;
    private int extent;
    private int arrayLength;
    private int bArrayOffset;
    private ByteBuffer buffer;

    public void mapAt(int offset, int length, ByteBuffer buffer){
        firstPointOffset = offset;
        numberOfPointsOffset = firstPointOffset+Integer.BYTES;
        mArrayOffset = firstPointOffset+2* Integer.BYTES;
        bArrayOffset = mArrayOffset+length*Double.BYTES;
        extent = 2*length*Double.BYTES + 2*Integer.BYTES;
        arrayLength = length;
        this.buffer = buffer;
    }

    public void copyFrom(int offset, ByteBuffer buffer){
        copyFrom(offset, arrayLength, buffer);
    }

    public void copyFrom(int offset, int length, ByteBuffer buffer){
        if (length != this.arrayLength){
            throw new RuntimeException("Array lengths should be equal!");
        }

        buffer.position(offset);
        this.buffer.putInt(firstPointOffset, buffer.getInt(offset+firstPointOffset));
        this.buffer.putInt(numberOfPointsOffset, buffer.getInt(offset+numberOfPointsOffset));
        for (int i = 0; i < 2*length; ++i){
            this.buffer.putDouble(mArrayOffset+i*Double.BYTES, buffer.getDouble(offset+mArrayOffset+i*Double.BYTES));
        }
    }

    public void copyTo(int offset, Bytes buffer){
        buffer.writeInt(offset+firstPointOffset, this.buffer.getInt(firstPointOffset));
        buffer.writeInt(offset+numberOfPointsOffset, this.buffer.getInt(numberOfPointsOffset));
        for (int i = 0; i < 2*arrayLength; ++i){
            try {
                buffer.writeDouble(offset+mArrayOffset+i*Double.BYTES, this.buffer.getDouble(mArrayOffset+i*Double.BYTES));
            } catch (Exception e) {
                if (ParallelOps.worldProcRank == 0) {
                    System.out.println("Rank: " + ParallelOps.worldProcRank +
                            " current write offset: " +
                            (offset + mArrayOffset + i * Double.BYTES) +
                            " buffersize: " + buffer.position(0).remaining());
                    e.printStackTrace();
                }
                throw e;
            }
        }
    }

    public MPISecPacket(int length){
        this(length, null);
    }

	private MPISecPacket(int length, ByteBuffer buffer) {
        extent = 2*length*Double.BYTES + 2*Integer.BYTES; // mArray + bArray + firstPoint + numberOfPoints
        arrayLength = length;
        bArrayOffset = mArrayOffset + length*Double.BYTES;
        this.buffer = (buffer == null) ? MPI.newByteBuffer(extent) : buffer;
	}

	public final void Clear()
	{
        setFirstPoint(0);
        setNumberOfPoints(0);
        for (int i = 0; i < arrayLength; ++i){
            setMArrayDoubleAt(i, 0);
        }
	}

	public static void memberCopy(MPISecPacket from, MPISecPacket to)
	{
        to.buffer = from.buffer;
        to.arrayLength = from.arrayLength;
        to.bArrayOffset = from.bArrayOffset;
        to.extent = from.extent;
	}

    public static MPISecPacket loadMPISecPacket(ByteBuffer buffer, int extent){
        return new MPISecPacket((extent - 2*Integer.BYTES)/(2*Double.BYTES), buffer);
    }

    public void setMArrayDoubleAt(int idx, double value){
        this.buffer.putDouble(mArrayOffset + idx * Double.BYTES, value);
    }

    public double getMArrayDoubleAt(int idx){
        return this.buffer.getDouble(mArrayOffset + idx * Double.BYTES);
    }

    public void setBArrayDoubleAt(int idx, double value){
        this.buffer.putDouble(bArrayOffset + idx * Double.BYTES, value);
    }

    public double getBArrayDoubleAt(int idx){
        return this.buffer.getDouble(bArrayOffset + idx * Double.BYTES);
    }

    public int getFirstPoint(){
        return buffer.getInt(firstPointOffset);
    }

    public void setFirstPoint(int firstPoint){
        buffer.putInt(firstPointOffset, firstPoint);
    }

    public int getNumberOfPoints(){
        return buffer.getInt(numberOfPointsOffset);
    }

    public void setNumberOfPoints(int numberOfPoints){
        buffer.putInt(numberOfPointsOffset, numberOfPoints);
    }

    public int getMArrayLength() {
        return arrayLength;
    }
    public int getBArrayLength() {
        return arrayLength;
    }

    public int getExtent() {
        return extent;
    }

    public ByteBuffer getBuffer() {
        return buffer;
    }

    @Override
    public String toString() {
        return "MPISecPacket{" +
                "extent=" + extent +
                ", arrayLength=" + arrayLength +
                ", bArrayOffset=" + bArrayOffset +
                ", buffer.capacity=" + buffer.capacity() +
                '}';
    }
}