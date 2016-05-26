package edu.indiana.soic.spidal.dapwc;

import mpi.MPI;
import net.openhft.lang.io.Bytes;

import java.nio.ByteBuffer;

/**
 Simple MPI Serializable packet for EM Loop
*/
public abstract class MPIPacket
{
    final int shift;
    private final int extent;
    private final int mArrayLength;
    ByteBuffer buffer;
    private static int firstPointOffset = 0;
    private static int numberOfPointsOffset = Integer.BYTES;
    private static int mArrayOffset = 2*Integer.BYTES;


    private MPIPacket(int mArrayLength, int extent, int shift, ByteBuffer buffer)
    {
        this.extent = extent;
        this.shift = shift;
        this.mArrayLength = mArrayLength;
        this.buffer =  (buffer == null) ? MPI.newByteBuffer(extent) : buffer; // new ByteBuffers are automatically initialized to zero
    }

    public static MPIPacket newIntegerPacket(int mArrayLength){
        return newIntegerPacket(mArrayLength, null);
    }

    private static MPIPacket newIntegerPacket(int mArrayLength, ByteBuffer buffer){
        return new MPIPacket(mArrayLength,(mArrayLength+2) * Integer.BYTES, Integer.BYTES, buffer){
            public int getMArrayIntAt(int idx){
                return this.buffer.getInt(mArrayOffset+idx*this.shift);
            }

            public void setMArrayIntAt(int idx, int value){
                this.buffer.putInt(mArrayOffset+idx*this.shift, value);
            }

            public void copyTo(int offset, Bytes buffer){
                buffer.writeInt(offset+firstPointOffset, this.buffer.getInt(firstPointOffset));
                buffer.writeInt(offset+numberOfPointsOffset, this.buffer.getInt(numberOfPointsOffset));
                for (int i = 0; i < mArrayLength; ++i){
                    buffer.writeInt(offset+mArrayOffset+i*Integer.BYTES, this.buffer.getInt(mArrayOffset+i*Integer.BYTES));
                }
            }

            public void copyFrom(int offset, int length, Bytes buffer) {
                buffer.position(offset);
                this.buffer.putInt(firstPointOffset, buffer.readInt(
                        offset + firstPointOffset));
                this.buffer.putInt(numberOfPointsOffset, buffer.readInt(
                        offset + numberOfPointsOffset));
                for (int i = 0; i < length; ++i) {
                    this.buffer.putInt(
                            mArrayOffset + i * Integer.BYTES, buffer.readInt(
                                    offset + mArrayOffset + i * Integer.BYTES));
                }
            }

            @Override
            public void Clear() {
                setFirstPoint(0);
                setNumberOfPoints(0);
                for (int i = 0; i < mArrayLength; ++i){
                    setMArrayIntAt(i, 0);
                }
            }
        };
    }

    public static MPIPacket loadIntegerPacket(ByteBuffer buffer){
        return newIntegerPacket(((buffer.limit()/Integer.BYTES) - 2), buffer); // buffer.limit() will be a multiplier of Interger.BYTES for this case
    }


    public static MPIPacket newDoublePacket(int mArrayLength){
        return newDoublePacket(mArrayLength, null);
    }

    private static MPIPacket newDoublePacket(int mArrayLength, ByteBuffer buffer){
        return new MPIPacket(mArrayLength, 2 * Integer.BYTES + mArrayLength * Double.BYTES, Double.BYTES, buffer){
            public double getMArrayDoubleAt(int idx){
                return this.buffer.getDouble(mArrayOffset + idx * this.shift);
            }

            public void setMArrayDoubleAt(int idx, double value){
                this.buffer.putDouble(mArrayOffset + idx * this.shift, value);
            }

            public void copyTo(int offset, Bytes buffer){
                buffer.writeInt(offset+firstPointOffset, this.buffer.getInt(firstPointOffset));
                buffer.writeInt(offset+numberOfPointsOffset, this.buffer.getInt(numberOfPointsOffset));
                for (int i = 0; i < mArrayLength; ++i){
                    buffer.writeDouble(offset+mArrayOffset+i*Double.BYTES, this.buffer.getDouble(mArrayOffset+i*Double.BYTES));
                }
            }

            public void copyFrom(int offset, int length, Bytes buffer){
                buffer.position(offset);
                this.buffer.putInt(firstPointOffset, buffer.readInt(
                        offset + firstPointOffset));
                this.buffer.putInt(numberOfPointsOffset, buffer.readInt(
                        offset + numberOfPointsOffset));
                for (int i = 0; i < length; ++i) {
                    this.buffer.putDouble(
                            mArrayOffset + i * Double.BYTES, buffer.readDouble(
                                    offset + mArrayOffset + i * Double.BYTES));
                }
            }

            @Override
            public void Clear() {
                setFirstPoint(0);
                setNumberOfPoints(0);
                for (int i = 0; i < mArrayLength; ++i){
                    setMArrayDoubleAt(i, 0.0);
                }
            }
        };
    }

    public static MPIPacket loadDoublePacket(ByteBuffer buffer){
        return newDoublePacket(((buffer.limit() - 2*Integer.BYTES)/Double.BYTES), buffer);
    }

    public int getMArrayIntAt(int idx){
        throw new UnsupportedOperationException();
    }

    public void setMArrayIntAt(int idx, int value){
        throw new UnsupportedOperationException();
    }

    public double getMArrayDoubleAt(int idx){
        throw new UnsupportedOperationException();
    }

    public void setMArrayDoubleAt(int idx, double value){
        throw new UnsupportedOperationException();
    }

    public void copyTo(int offset, Bytes to){
        throw new UnsupportedOperationException();
    }

    public void copyFrom(int offset, int length,  Bytes from){
        throw new UnsupportedOperationException();
    }

    public ByteBuffer getBuffer() {
        return buffer;
    }

    public int getExtent() {
        return extent;
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

    public int getArrayLength() {
        return mArrayLength;
    }

    public abstract void Clear();


    public enum Type {Integer, Double}
}