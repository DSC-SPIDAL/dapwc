package edu.indiana.soic.spidal.configuration.sections;

import edu.indiana.soic.spidal.tools.ClusterExtractor;

import java.io.FileInputStream;
import java.io.IOException;
import java.util.HashSet;
import java.util.Properties;
import java.util.StringJoiner;

/**
 * Created by pulasthi on 3/19/17.
 */
public class ClusterExtractorSection {
    public String clusterFile;
    public String distFile;
    public String outDir;
    public int numPoints;
    public String clusters;
    public String newclusters_percluster;
    public boolean isBigEndian = false;
    public String runLine;
    public String joins;
    public ClusterExtractorSection(String configurationFilePath){
        Properties p = new Properties();
        try {
            p.load(new FileInputStream(configurationFilePath));
            clusterFile = getProperty(p,"clusterFile","cluster.txt");
            distFile = getProperty(p,"distFile","in.bin");
            numPoints = Integer.valueOf(getProperty(p,"numPoints","-1"));
            outDir = getProperty(p,"outDir","outfile.bin");
            isBigEndian = Boolean.valueOf(getProperty(p,"isBigEndian","false"));
            clusters = getProperty(p,"clusters","None");
            newclusters_percluster = getProperty(p,"newclusters","1");
            runLine = getProperty(p,"runLine","");
            joins = getProperty(p,"joins","");
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
