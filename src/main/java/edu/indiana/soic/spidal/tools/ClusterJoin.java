package edu.indiana.soic.spidal.tools;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.HashMap;
import java.util.HashSet;

/**
 * Created by pulasthi on 3/31/17.
 */
public class ClusterJoin {
    static final int dustClusterId = 100000; // this needs to be same with the value in ClusterOutlierExtractor

    public static void main(String[] args) {
        String inputFileName = args[0];
        String outFileName = args[1];
        String[] joins = args[2].split("\\|");
        HashMap<Integer, Integer> newClusterNumMap = new HashMap<Integer, Integer>();
        HashMap<Integer, Integer> joinClusterNumMap = new HashMap<Integer, Integer>();
        HashSet<Integer> originalCluster = new HashSet<Integer>();

        try(BufferedReader br = Files.newBufferedReader(Paths.get(inputFileName))){

            String line;
            while ((line = br.readLine()) != null){
                String splits[] = line.split("\\s+");
                originalCluster.add(Integer.valueOf(splits[1]));
            }
        }catch (IOException e) {
            e.printStackTrace();
        }

        for (String joinpattern : joins) {
            String[] clus = joinpattern.split("\\+");
            if(clus.length <= 1) continue;
            int main = Integer.valueOf(clus[0]);
            for (int i = 1; i < clus.length; i++) {
                int other = Integer.valueOf(clus[i]);
                joinClusterNumMap.put(other, main);
            }
        }

        int newClusterNum = 0;
        for (int i = 0; i < originalCluster.size(); i++) {
            int tempcluster = newClusterNum;
            if(joinClusterNumMap.containsKey(i)){
                tempcluster = newClusterNumMap.get(joinClusterNumMap.get(i));
            }else{
                newClusterNum++;
            }
            newClusterNumMap.put(i,tempcluster);

        }

        try(BufferedReader br = Files.newBufferedReader(Paths.get(inputFileName))){
            PrintWriter pr = new PrintWriter(Paths.get(outFileName).toString());
            String line;
            while ((line = br.readLine()) != null){
                String splits[] = line.split("\\s+");
                int cluster = Integer.valueOf(splits[1]);
                if(cluster == dustClusterId){
                    pr.println(splits[0] + "\t" + dustClusterId);
                }else{
                    int newCluster = newClusterNumMap.get(Integer.valueOf(splits[1]));
                    pr.println(splits[0] + "\t" + newCluster);
                }

            }
            pr.flush();
            pr.close();
        }catch (IOException e) {
            e.printStackTrace();
        }
    }
}
