package edu.indiana.soic.spidal.configuration.sections;

import java.io.FileInputStream;
import java.io.IOException;
import java.util.HashSet;
import java.util.Properties;

/**
 * Created by pulasthi on 3/19/17.
 */
public class CollateClustersSection {
    public String pointsDir;
    public String pointFilePattern;
    public String clustersDir;
    public String clusterDirpattern;
    public String clusterFilePattern;
    public String outDir;
    public String outFilePattern;
    public String outPlotFilePattern;
    public String outLablesFilePattern;
    public String outsepereateFilesPattern;
    public int numClusters;
    public String joinPattern;
    public String clusters;

    public CollateClustersSection(String configurationFilePath){
        Properties p = new Properties();
        try {
            p.load(new FileInputStream(configurationFilePath));
            pointsDir = getProperty(p,"pointsDir","/home/");
            pointFilePattern =getProperty(p,"pointFilePattern","points.txt");
            clustersDir = getProperty(p,"clustersDir","/home/");
            clusterDirpattern =getProperty(p,"clusterDirpattern","{0}");
            clusterFilePattern =getProperty(p,"clusterFilePattern","clusters.txt");
            outDir = getProperty(p,"outDir","/home/");
            outFilePattern = getProperty(p,"outFilePattern","out.txt");
            outPlotFilePattern = getProperty(p, "outPlotFilePattern", "outplot.txt");
            outLablesFilePattern = getProperty(p, "outLablesFilePattern", "outlabel.txt");
            outsepereateFilesPattern = getProperty(p,"outsepereateFilesPattern","sepereated.txt");
            numClusters = Integer.valueOf(getProperty(p,"numClusters","0"));
            joinPattern = getProperty(p,"join","");
            clusters = getProperty(p, "clusters", "");
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    private static String getProperty(Properties p, String name, String def) {
        String val = System.getProperty(name);
        if (val == null) {
            if (def != null) {
                val = p.getProperty(name, def);
            } else {
                val = p.getProperty(name);
            }
        }
        return val;
    }

    private static HashSet<Integer> getHashSet(String ranks){
        HashSet<Integer> set = new HashSet<Integer>();
        String[] skips = ranks.split(",");
        for (String skip : skips) {
            set.add(Integer.valueOf(skip));
        }
        return set;
    }
}


