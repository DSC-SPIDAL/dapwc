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
public class CollateClusters {
    public static void main(String[] args) {
        CollateClustersSection section = new CollateClustersSection(args[0]);

        HashMap<Integer, Integer> newclustermappings = new HashMap<Integer, Integer>();
        HashMap<Integer, String> filenameMappings = new HashMap<Integer, String>();
        HashMap<Integer, Integer> newclusterstartpoints = new HashMap<Integer, Integer>();
        HashMap<Integer, Integer> newclusterlengths = new HashMap<Integer, Integer>();
        HashMap<Integer, Integer> origclusterlengths = new HashMap<Integer, Integer>();
        HashMap<Integer, HashMap<Integer,Integer>> joinMappings = new HashMap<Integer, HashMap<Integer,Integer>>();
        HashMap<Integer, HashMap<Integer,Integer>> finalMappings = new HashMap<Integer, HashMap<Integer,Integer>>();

        //Read join info and populate map
        String[] joins = section.joinPattern.split("\\|");
        for (String join : joins) {
            String[] splits = join.split(",");
            if(splits.length <= 1) continue;
            int clusterNum = Integer.valueOf(splits[0]);
            joinMappings.put(clusterNum,new HashMap<Integer,Integer>());
            for (int i = 1; i < splits.length; i++) {
                String[] splitjoins = splits[i].split("\\+");
                if(splitjoins.length <= 1) continue;
                int main = Integer.valueOf(splitjoins[0]);
                for (int j = 1; j < splitjoins.length; j++) {
                    int splitjoin = Integer.valueOf(splitjoins[j]);
                    joinMappings.get(clusterNum).put(splitjoin,main);
                }
            }

        }

        String[] subclustersplits = section.clusters.split("\\|");
        for (String subclustersplit : subclustersplits) {
            String[] vals = subclustersplit.split(",");
            if(vals.length < 1) continue;

            int clusterNum = Integer.valueOf(vals[0]);
            int sugestedClusters = Integer.valueOf(vals[1]);
            int actualClusters = Integer.valueOf(vals[2]);
            String filename = section.clusterFilePattern.replace("{0}",""+sugestedClusters);
            filename = filename.replace("{1}",""+actualClusters);
            filenameMappings.put(clusterNum,filename);
            origclusterlengths.put(clusterNum,actualClusters);
            if(joinMappings.containsKey(clusterNum)){
                actualClusters = actualClusters - joinMappings.get(clusterNum).size();
            }

            newclusterlengths.put(clusterNum,actualClusters);

        }
        for (int i = 0; i < section.numClusters; i++) {
            if(!newclusterlengths.containsKey(i)){
                newclusterlengths.put(i,1);
                origclusterlengths.put(i,1);
            }
        }

        int currentCount = 0;
        int totalclusterCount = 0;
        for (int i = 0; i < section.numClusters; i++) {
            newclusterstartpoints.put(i,currentCount);

            HashMap<Integer,Integer> curr = new HashMap<Integer,Integer>();
            for (int j = 0; j < origclusterlengths.get(i); j++) {
                if(joinMappings.containsKey(i)){
                    if(joinMappings.get(i).containsKey(j)){
                        int mappedint = joinMappings.get(i).get(j);
                        curr.put(j,curr.get(mappedint));
                    }else{
                        curr.put(j,totalclusterCount++);
                    }
                }else{
                    curr.put(j,totalclusterCount++);
                }
            }
            finalMappings.put(i,curr);
            currentCount += newclusterlengths.get(i);

        }
        try{
            BufferedReader original = Files.newBufferedReader(Paths.get(section.pointsDir,section.pointFilePattern));
            HashMap<Integer,BufferedReader> clusterFiles = new HashMap<Integer,BufferedReader>();
            String outFile = Paths.get(section.outDir,section.outFilePattern).toString();
            PrintWriter printWriter = new PrintWriter(new FileWriter(outFile));
            for (int i = 0; i < section.numClusters; i++) {
                if(filenameMappings.containsKey(i)){
                    BufferedReader temp = Files.newBufferedReader(Paths.get(section.clustersDir,section.clusterDirpattern.replace("{0}",""+i),filenameMappings.get(i)));
                    clusterFiles.put(i,temp);
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
                if(filenameMappings.containsKey(oriCluster)){
                    line = clusterFiles.get(oriCluster).readLine();
                    String[] data = line.split("\\s+");
                    int subClusterNum = Integer.valueOf(data[1]);
                    newClusterId = finalMappings.get(oriCluster).get(subClusterNum);
                }else{
                    newClusterId = finalMappings.get(oriCluster).get(0);
                }
                printWriter.println(index + " " + splits[1] + " " + splits[2] + " " + splits[3] + " " + newClusterId + " " + newClusterId);
            }
            printWriter.flush();
            printWriter.close();
            original.close();
        }catch (IOException e){
            e.printStackTrace();
        }


    }
}
