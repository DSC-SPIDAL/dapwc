package edu.indiana.soic.spidal.dapwc;

import com.google.common.io.LittleEndianDataInputStream;

import java.io.DataInput;
import java.io.DataInputStream;
import java.io.FileInputStream;
import java.io.IOException;
import java.nio.ByteOrder;
import java.nio.MappedByteBuffer;
import java.nio.ShortBuffer;
import java.nio.channels.FileChannel;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;

public abstract class Matrix {
    public static short[][] readRowRange(String fname, int rowStartForProc, int rowCountForProc, int globalColCount, int dataTypeSize, ByteOrder endianness, boolean mmap){
        if (mmap) {
            try (FileChannel fc = (FileChannel) Files.newByteChannel(Paths.get(fname), StandardOpenOption.READ)) {
                long pos = ((long) rowStartForProc) * globalColCount * dataTypeSize; // byte position
                long length = ((long) rowCountForProc) * globalColCount * dataTypeSize; // number of bytes to map
                if (length > Integer.MAX_VALUE){
                    throw new RuntimeException("Error - mapping a region larger than Integer.MAX_VALUE is not supported yet. Try increasing the number of processes.");
                }

                // Once mapped as a ShortBuffer we only need to give position and lengths in terms of shorts not bytes
                ShortBuffer mappedShorts = fc.map(FileChannel.MapMode.READ_ONLY, pos, length).order(endianness).asShortBuffer();
                short [][] distances = new short[rowCountForProc][];
                int shortPos = 0; // short position
                for (int rowCount = 0; rowCount < rowCountForProc; ++rowCount) {
                    distances[rowCount] = new short[globalColCount];
                    mappedShorts.position(shortPos);
                    mappedShorts.get(distances[rowCount]);
                    shortPos += globalColCount;
                }
                return distances;
            } catch (IOException e) {
                e.printStackTrace();
            }
        } else {
            throw new UnsupportedOperationException("Error - other mapping types are not supported. Use memory mapped I/O.");
        }
        return null;
    }
}
