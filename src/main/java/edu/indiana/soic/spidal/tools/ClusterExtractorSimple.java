package edu.indiana.soic.spidal.tools;

import java.io.*;
import java.nio.Buffer;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.ShortBuffer;
import java.nio.channels.FileChannel;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.HashMap;
import java.util.List;
import java.util.Properties;

/**
 * Created by pulasthi on 4/1/17.
 */
public class ClusterExtractorSimple {
    public static void main(String[] args) {
        String clusterFile = args[0];
        String distFile = args[1];
        String outDir = args[2];
        String[] clusters = args[3].split(",");
        String[] newclusters = args[4].split(",");
        boolean isBigendian = Boolean.valueOf(args[5]);
        int numPoints = 0;
        HashMap<Integer, List<Integer>> clusterPoints = new HashMap<Integer, List<Integer>>();
        HashMap<Integer, Integer> clusterNumMap = new HashMap<Integer, Integer>();
        for (int i = 0; i < clusters.length; i++) {
            clusterNumMap.put(Integer.valueOf(clusters[i]), Integer.valueOf(newclusters[i]));

        }

        try(BufferedReader br = Files.newBufferedReader(Paths.get(clusterFile))){
            String line;
            int clusterNum;
            int dataPoint;
            while ((line = br.readLine()) != null){
                if(line == "") break;

                String[] splits = line.split("\\s+");
                if(splits.length != 2) break;
                numPoints++;
                clusterNum = Integer.valueOf(splits[1]);

                dataPoint = Integer.valueOf(splits[0]);
                clusterPoints.get(clusterNum).add(dataPoint);
            }
        }catch (IOException e){
            e.printStackTrace();
        }

        //reading from dist file and writing to the new cluster files
        int newClusterNumbers = 0;
        for (String cluster : clusters) {

            int clusterNum = Integer.valueOf(cluster);

            if(!clusterNumMap.containsKey(clusterNum)) continue;

            Path filePath = Paths.get(outDir,""+clusterNum + ".bin");
            try {
                FileChannel fc = (FileChannel) Files
                        .newByteChannel(Paths.get(distFile), StandardOpenOption.READ);
                File file = new File(filePath.toString());
                if(!file.exists()) {
                    file.createNewFile();
                }else{
                    file.delete();
                    file.createNewFile();
                }
                FileChannel fcout = (FileChannel) Files
                        .newByteChannel(filePath, StandardOpenOption.APPEND);
                for (Integer row : clusterPoints.get(clusterNum)) {
                    long offset = row*numPoints*2;
                    ByteBuffer byteBufferRow = ByteBuffer.allocate(numPoints*2);
                    if(isBigendian){
                        byteBufferRow.order(ByteOrder.BIG_ENDIAN);
                    }else{
                        byteBufferRow.order(ByteOrder.LITTLE_ENDIAN);
                    }
                    fc.read(byteBufferRow,offset);
                    byteBufferRow.flip();
                    Buffer buffer = byteBufferRow.asShortBuffer();
                    short[] shortArray = new short[numPoints];
                    ((ShortBuffer)buffer).get(shortArray);

                    //Ouput section of row
                    ByteBuffer outBuffer = ByteBuffer.allocate(clusterPoints.get(clusterNum).size()*2);
                    if(isBigendian){
                        outBuffer.order(ByteOrder.BIG_ENDIAN);
                    }else{
                        outBuffer.order(ByteOrder.LITTLE_ENDIAN);
                    }

                    for (Integer integer : clusterPoints.get(clusterNum)) {
                        outBuffer.putShort(shortArray[integer]);
                    }
                    outBuffer.flip();
                    fcout.write(outBuffer);
                }
                fcout.close();

                //Generate config files for each cluster
                Path outClusterPath = Paths.get(outDir,"cluster_" + clusterNum + ".txt");
                Properties template = new Properties();
                InputStream in = ClusterExtractor.class.getResourceAsStream("/dapwc_config_template.properties");
                template.load(in);
                template.setProperty("DistanceMatrixFile",filePath.toString());
                template.setProperty("ClusterFile",outClusterPath.toString());
                template.setProperty("NumberDataPoints",""+clusterPoints.get(clusterNum).size());
                template.setProperty("MaxNcent",clusterNumMap.get(clusterNum).toString());

                Path confFilePath = Paths.get(outDir,"config_" + clusterNum + ".properties");
                File fileconf = new File(confFilePath.toString());
                if(!fileconf.exists()){
                    fileconf.createNewFile();
                }

                template.store(new FileOutputStream(confFilePath.toString(),false),null);


            }catch (IOException e){
                e.printStackTrace();
            }



        }


    }
}
