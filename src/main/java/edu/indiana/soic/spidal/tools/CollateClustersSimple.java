package edu.indiana.soic.spidal.tools;

import edu.indiana.soic.spidal.configuration.sections.CollateClustersSection;

import java.io.*;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.HashMap;

/**
 * Created by pulasthi on 3/19/17.
 */
public class CollateClustersSimple {
    static final int dustClusterId = 100000; // this needs to be same with the value in ClusterOutlierExtractor
    public static void main(String[] args) {
        String pointsFile = args[0];
        String clusterDirPattern = args[1];
        String clusterFilePattern = args[2];
        String outDir = args[3];
        int numClusters = Integer.valueOf(args[4]);
        String clusterPattern = args[5];

        HashMap<Integer, String> filenameMappings = new HashMap<Integer, String>();
        HashMap<Integer, Integer> newclusterstartpoints = new HashMap<Integer, Integer>();
        HashMap<Integer, Integer> newclusterlengths = new HashMap<Integer, Integer>();
        HashMap<Integer, Integer> origclusterlengths = new HashMap<Integer, Integer>();
        HashMap<Integer, HashMap<Integer,Integer>> finalMappings = new HashMap<Integer, HashMap<Integer,Integer>>();


        String[] subclustersplits = clusterPattern.split("\\|");
        for (String subclustersplit : subclustersplits) {
            String[] vals = subclustersplit.split(",");
            if(vals.length < 1) continue;

            int clusterNum = Integer.valueOf(vals[0]);
            int sugestedClusters = Integer.valueOf(vals[1]);
            int actualClusters = Integer.valueOf(vals[2]);
            String filename = clusterFilePattern.replace("{0}",""+sugestedClusters);
            filename = filename.replace("{1}",""+actualClusters);
            filenameMappings.put(clusterNum,filename);
            origclusterlengths.put(clusterNum,actualClusters);

            newclusterlengths.put(clusterNum,actualClusters);

        }
        for (int i = 0; i < numClusters; i++) {
            if(!newclusterlengths.containsKey(i)){
                newclusterlengths.put(i,1);
                origclusterlengths.put(i,1);
            }
        }

        int currentCount = 0;
        int totalclusterCount = 0;
        for (int i = 0; i < numClusters; i++) {
            newclusterstartpoints.put(i,currentCount);

            HashMap<Integer,Integer> curr = new HashMap<Integer,Integer>();
            for (int j = 0; j < origclusterlengths.get(i); j++) {
                curr.put(j,totalclusterCount++);
            }
            finalMappings.put(i,curr);
            currentCount += newclusterlengths.get(i);

        }
        try{
            BufferedReader original = Files.newBufferedReader(Paths.get(pointsFile));
            HashMap<Integer,BufferedReader> clusterFiles = new HashMap<Integer,BufferedReader>();
            HashMap<Integer,PrintWriter> clusterFilesplots = new HashMap<Integer,PrintWriter>();
            HashMap<Integer,Integer> clusterFilesplotsCount = new HashMap<Integer,Integer>();
            String outFile = Paths.get(outDir,"plot.webplotviz").toString();
            String outFileClust = Paths.get(outDir,"cluster.txt").toString();
            PrintWriter printWriter = new PrintWriter(new FileWriter(outFile));
            PrintWriter printWriterout = new PrintWriter(new FileWriter(outFileClust));
            for (int i = 0; i < numClusters; i++) {
                if(filenameMappings.containsKey(i)){
                    BufferedReader temp = Files.newBufferedReader(Paths.get(clusterDirPattern.replace("{0}",""+i),filenameMappings.get(i)));
                    PrintWriter tempout = new PrintWriter(Paths.get(outDir,"plot_" + i + ".webplotviz").toString());
                    clusterFiles.put(i,temp);
                    clusterFilesplots.put(i,tempout);
                    clusterFilesplotsCount.put(i,0);
                }
            }

            String origLine;
            while ((origLine = original.readLine()) != null){
                if(origLine == ""){
                    System.out.println("Hit error while reading file");
                    break;
                }

                String[] splits = origLine.split("\\s+");
                int index = Integer.valueOf(splits[0]);
                int oriCluster = Integer.valueOf(splits[4]);
                int newClusterId;
                String line;
                if(oriCluster != dustClusterId ) {
                    if (filenameMappings.containsKey(oriCluster)) {
                        line = clusterFiles.get(oriCluster).readLine();
                        String[] data = line.split("\\s+");
                        int subClusterNum = Integer.valueOf(data[1]);
                        newClusterId = finalMappings.get(oriCluster).get(subClusterNum);

                        clusterFilesplots.get(oriCluster).println(clusterFilesplotsCount.get(oriCluster) + " " + splits[1] + " " + splits[2] + " " + splits[3] + " " + newClusterId + " " + newClusterId);
                        clusterFilesplotsCount.put(oriCluster, clusterFilesplotsCount.get(oriCluster) + 1);
                    } else {
                        newClusterId = finalMappings.get(oriCluster).get(0);
                    }

                    printWriter.println(index + " " + splits[1] + " " + splits[2] + " " + splits[3] + " " + newClusterId + " " + newClusterId);
                    printWriterout.println(index + "\t" + newClusterId);
                }else{
                    printWriter.println(index + " " + splits[1] + " " + splits[2] + " " + splits[3] + " " + dustClusterId + " " + "Dust");
                    printWriterout.println(index + "\t" + dustClusterId);
                }

            }
            printWriter.flush();
            printWriter.close();
            printWriterout.flush();
            printWriterout.close();
            original.close();

            for (PrintWriter writer : clusterFilesplots.values()) {
                writer.flush();
                writer.close();
            }
        }catch (IOException e){
            e.printStackTrace();
        }


    }
}
