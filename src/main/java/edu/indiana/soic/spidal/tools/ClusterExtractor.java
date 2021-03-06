package edu.indiana.soic.spidal.tools;

import edu.indiana.soic.spidal.configuration.sections.ClusterExtractorSection;
import javafx.collections.transformation.SortedList;

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
 * Created by pulasthi on 3/19/17.
 */
public class ClusterExtractor {
    public static void main(String[] args) {
        ClusterExtractorSection section = new ClusterExtractorSection(args[0]);
        HashMap<Integer, List<Integer>> clusterPoints = new HashMap<Integer, List<Integer>>();
        HashMap<Integer, Integer> joinClusters = new HashMap<Integer, Integer>();
        HashMap<Integer, Integer> newClusterNumMap = new HashMap<Integer, Integer>();

        String[] clusters = section.clusters.split(",");
        String[] newclusters = section.newclusters_percluster.split(",");

        for (String cluster : clusters) {
            int clust = Integer.valueOf(cluster);
            if(!clusterPoints.containsKey(clust)){
                clusterPoints.put(clust, new ArrayList<Integer>());
            }
        }

        String[] joinpatterns = section.joins.split("\\|");
        for (String joinpattern : joinpatterns) {
            String[] clus = joinpattern.split("\\+");
            if(clus.length <= 1) continue;
            int main = Integer.valueOf(clus[0]);
            for (int i = 1; i < clus.length; i++) {
                int other = Integer.valueOf(clus[i]);
                joinClusters.put(other, main);
            }
        }

        try(BufferedReader br = Files.newBufferedReader(Paths.get(section.clusterFile))){
            String line;
            int clusterNum;
            int dataPoint;
            while ((line = br.readLine()) != null){
                if(line == "") break;

                String[] splits = line.split("\\s+");
                if(splits.length != 2) break;

                clusterNum = Integer.valueOf(splits[1]);
                if(joinClusters.containsKey(clusterNum)){
                    clusterNum = joinClusters.get(clusterNum);
                }
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
            if(joinClusters.containsKey(clusterNum)){
                clusterNum = joinClusters.get(clusterNum);
            }
            if(!newClusterNumMap.containsKey(clusterNum)){
                newClusterNumMap.put(clusterNum,newClusterNumbers++);
            }
            int newClusterNum = newClusterNumMap.get(clusterNum);
            Path filePath = Paths.get(section.outDir,""+newClusterNum + ".bin");
            try {
                FileChannel fc = (FileChannel) Files
                        .newByteChannel(Paths.get(section.distFile), StandardOpenOption.READ);
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
                    long offset = row*section.numPoints*2;
                    ByteBuffer byteBufferRow = ByteBuffer.allocate(section.numPoints*2);
                    if(section.isBigEndian){
                        byteBufferRow.order(ByteOrder.BIG_ENDIAN);
                    }else{
                        byteBufferRow.order(ByteOrder.LITTLE_ENDIAN);
                    }
                    fc.read(byteBufferRow,offset);
                    byteBufferRow.flip();
                    Buffer buffer = byteBufferRow.asShortBuffer();
                    short[] shortArray = new short[section.numPoints];
                    ((ShortBuffer)buffer).get(shortArray);

                    //Ouput section of row
                    ByteBuffer outBuffer = ByteBuffer.allocate(clusterPoints.get(clusterNum).size()*2);
                    if(section.isBigEndian){
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
                Path outClusterPath = Paths.get(section.outDir,"cluster_" + newClusterNum + ".txt");
                Properties template = new Properties();
                InputStream in = ClusterExtractor.class.getResourceAsStream("/dapwc_config_template.properties");
                template.load(in);
                template.setProperty("DistanceMatrixFile",filePath.toString());
                template.setProperty("ClusterFile",outClusterPath.toString());
                template.setProperty("NumberDataPoints",""+clusterPoints.get(clusterNum).size());
                template.setProperty("MaxNcent",newclusters[newClusterNum]);

                Path confFilePath = Paths.get(section.outDir,"config_" + newClusterNum + ".properties");
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
