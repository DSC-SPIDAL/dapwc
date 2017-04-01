package edu.indiana.soic.spidal.tools;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.Paths;

/**
 * Created by pulasthi on 3/31/17.
 */
public class CreatePlotSimple {
    public static void main(String[] args) {
        String pointsFile = args[0];
        String clusterFile = args[1];
        String outFile = args[2];


        try {
            BufferedReader points = Files.newBufferedReader(Paths.get(pointsFile));
            BufferedReader cluster = Files.newBufferedReader(Paths.get(clusterFile));
            PrintWriter out = new PrintWriter(Paths.get(outFile).toString());
            String line;
            String lineclus;
            while ((line = points.readLine()) != null){
                lineclus = cluster.readLine();
                String pointssplit[] = line.split("\\s+");
                String clustersplit[] = lineclus.split("\\s+");
                out.println(pointssplit[0] + " " + pointssplit[1] + " " + pointssplit[2] + " " + pointssplit[3] + " " + clustersplit[1] + " " + clustersplit[1]);
            }

            out.flush();
            out.close();
        }catch (IOException e){
            e.printStackTrace();
        }
    }
}
