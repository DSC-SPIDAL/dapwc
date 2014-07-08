package edu.indiana.soic.spidal.dapwc;

import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.MappedByteBuffer;
import java.nio.channels.FileChannel;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;

public class Matrix {
    private int globalColCount;
    private Range rowRange;
    private ByteBuffer data;

    public Matrix(int globalColCount, Range rowRange, ByteBuffer data) {
        this.globalColCount = globalColCount;
        this.rowRange = rowRange;
        this.data = data;
    }

    public double getDistance(int globalRow, int globalCol){
        int pos = (globalRow - rowRange.getStartIndex()) * globalColCount + globalCol; // element position - not the byte position
        // TODO - fix data type to things other than short
        return data.getShort(pos*2) / (Short.MAX_VALUE*1.0); // pos*2 is the byte position
    }

    public static Matrix readRowRange(String fname, Range rows, int globalColCount, int dataTypeSize, ByteOrder
            endianness){
        try (FileChannel fc = (FileChannel) Files.newByteChannel(Paths.get(fname), StandardOpenOption.READ)) {
            long pos = ((long)rows.getStartIndex())*globalColCount* dataTypeSize;
            MappedByteBuffer mappedBytes = fc.map(FileChannel.MapMode.READ_ONLY, pos,
                                                  rows.getLength() * globalColCount * dataTypeSize);
            mappedBytes.order(endianness);
            return  new Matrix(globalColCount,rows, mappedBytes);
        } catch (IOException e) {
            e.printStackTrace();
        }
        return null;
    }
}
