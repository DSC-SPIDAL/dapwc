package edu.indiana.soic.spidal.dapwc;

import mpi.MPI;

import java.io.Serializable;
import java.nio.ByteBuffer;

public class MPISecPacket implements Serializable
{
    private int firstPointOffset = 0;
    private int numberOfPointsOffset = Integer.BYTES;
    private int mArrayOffset = 2*Integer.BYTES;
    private int extent;
    private int arrayLength;
    private int bArrayOffset;
    private ByteBuffer buffer;

    public void mapAt(int offset, int length){
        firstPointOffset = offset;
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