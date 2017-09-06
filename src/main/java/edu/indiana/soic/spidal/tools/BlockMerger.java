package edu.indiana.soic.spidal.tools;

import edu.indiana.soic.spidal.configuration.sections.BlockMergerSection;

import java.io.FileOutputStream;
import java.io.IOException;
import java.nio.Buffer;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.ShortBuffer;
import java.nio.channels.FileChannel;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.ArrayList;

/**
 * Created by pulasthi on 3/18/17.
 * Used to merge blocks from the SmithWaterman algorithm
 */
public class BlockMerger {
    public static void main(String[] args) {
        BlockMergerSection section = new BlockMergerSection(args[0]);
        try {
            FileChannel fc;
            Buffer buffer;
            FileChannel out = new FileOutputStream(section.outFile).getChannel();
            for (int i = 0; i < section.blockCount; i++) {
                ArrayList<short[]> columns =  new ArrayList<short[]>();
                if(section.skipRanks.contains(i)){
                    System.out.printf("Skipping row block %d", i);
                    continue;
                }

                for (int j = 0; j < section.blockCount; j++) {
                    if(section.skipRanks.contains(i)){
                        System.out.printf("Skipping column block %d", j);
                        columns.add(null);
                        continue;
                    }
                    String fileName = getBlockName(section.nameFormat,i,j);
                    Path path = Paths.get(section.blockDir,fileName);
                    if(Files.exists(path)){
                        fc = (FileChannel)Files.newByteChannel(path,StandardOpenOption.READ);
                    }else{
                        // If direct block (i,j) does not exist then transpose (j,i) must exist.
                        fileName = getBlockName(section.nameFormat,i,j);
                        path = Paths.get(section.blockDir,fileName);
                        if(!Files.exists(path)){
                            System.out.println("Unable to find file" + fileName);
                            throw new IOException("Unable to find file" + i + "," + j);
                        }
                        fc = (FileChannel)Files.newByteChannel(path,StandardOpenOption.READ);

                    }

                    ByteBuffer byteBuffer = ByteBuffer.allocate((int)fc.size());
                    if(section.isBigEndian){
                        byteBuffer.order(ByteOrder.BIG_ENDIAN);
                    }else{
                        byteBuffer.order(ByteOrder.LITTLE_ENDIAN);
                    }

                    fc.read(byteBuffer);
                    byteBuffer.flip();
                    out.write(byteBuffer);
                }

                System.out.printf("Done row %d",i);
                System.out.println("Writing column blocks for row: " + i + " ... ");
            }
            out.close();

        }catch (IOException e) {
            e.printStackTrace();
        }

    }

    private static String getBlockName(String nameFormat, int i, int j) {
        String blockName = nameFormat.replace("{1}",""+i);
        blockName = blockName.replace("{2}",""+j);
        return blockName;
    }

}
