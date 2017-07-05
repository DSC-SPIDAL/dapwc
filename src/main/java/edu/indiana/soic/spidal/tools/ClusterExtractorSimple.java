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
import java.util.ArrayList;
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
        // MDS file is optional and will be used only if specified
        String mdsFile = "";
        double[][] mdsPoints = new double[0][];
        boolean writeMDS = false;

        int numPoints = 0;
        String subclusters = "";
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
                if(clusterPoints.containsKey(clusterNum)){
                    clusterPoints.get(clusterNum).add(dataPoint);
                }else{
                    clusterPoints.put(clusterNum, new ArrayList<Integer>());
                    clusterPoints.get(clusterNum).add(dataPoint);
                }
            }
        }catch (IOException e){
            e.printStackTrace();
        }

        if(args.length == 7){
            mdsFile = args[6];
            mdsPoints = new double[numPoints][3];
            writeMDS = true;

            try {
                BufferedReader points = Files.newBufferedReader(Paths.get(mdsFile));

                String line;
                int count = 0;
                while ((line = points.readLine()) != null){
                    String pointssplit[] = line.split("\\s+");
                    mdsPoints[count][0] = Double.parseDouble(pointssplit[1]);
                    mdsPoints[count][1] = Double.parseDouble(pointssplit[2]);
                    mdsPoints[count][2] = Double.parseDouble(pointssplit[3]);
                    count++;
                }
            }catch (IOException e){
                e.printStackTrace();
            }
        }

        //reading from dist file and writing to the new cluster files
        int newClusterNumbers = 0;
        for (String cluster : clusters) {

            int clusterNum = Integer.valueOf(cluster);

            if(!clusterNumMap.containsKey(clusterNum)) continue;

            subclusters += clusterNum + "," + clusterNumMap.get(clusterNum) + "," + clusterNumMap.get(clusterNum) + "|";
            Path filePath = Paths.get(outDir,clusterNum + ".bin");
            Path filePathMDS = Paths.get(outDir,clusterNum + "_mds.txt");
            try {
                FileChannel fc = (FileChannel) Files
                        .newByteChannel(Paths.get(distFile), StandardOpenOption.READ);
                File file = new File(filePath.toString());
                File fileMDS;
                PrintWriter outMDS = null;

                if(writeMDS){
                    fileMDS = new File(filePathMDS.toString());
                    if(!file.exists()) {
                        fileMDS.createNewFile();
                    }else{
                        fileMDS.delete();
                        fileMDS.createNewFile();
                    }
                    outMDS = new PrintWriter(fileMDS);

                }

                if(!file.exists()) {
                    file.createNewFile();
                }else{
                    file.delete();
                    file.createNewFile();
                }

                FileChannel fcout = (FileChannel) Files
                        .newByteChannel(filePath, StandardOpenOption.APPEND);
                int tempcount = 0;
                for (Integer row : clusterPoints.get(clusterNum)) {
                    if(writeMDS) outMDS.println(tempcount + "\t" + mdsPoints[row][0] + "\t" + mdsPoints[row][1] + "\t" + mdsPoints[row][2] + "\t" + "1");
                    System.out.println("Row :" + row );
                    System.out.println("numPoints :" + numPoints );
                    long offset = row.longValue()*(long)numPoints*2;
                    System.out.println("Offset :" + offset );
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
                    tempcount++;
                }
                outMDS.flush();
                outMDS.close();
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

        //Genereate data for next round collate
        Properties template = new Properties();
        template.setProperty("numClusters",""+clusterPoints.size());
        template.setProperty("subClusters","\"" + subclusters + "\"");

        Path dataFilePath = Paths.get(outDir,"conf_shared.prop");
        try{
            File fileconf = new File(dataFilePath.toString());
            if(!fileconf.exists()){
                fileconf.createNewFile();
            }
            template.store(new FileOutputStream(dataFilePath.toString(),false),null);

        }catch (IOException e){
            e.printStackTrace();
        }

        //PrintWriter dataout = new PrintWriter(Paths.get(outDir,))


    }
}
