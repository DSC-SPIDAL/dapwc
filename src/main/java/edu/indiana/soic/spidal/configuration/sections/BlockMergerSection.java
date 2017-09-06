package edu.indiana.soic.spidal.configuration.sections;

import java.io.FileInputStream;
import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Properties;
import java.util.function.BooleanSupplier;

/**
 * Created by pulasthi on 3/18/17.
 */
public class BlockMergerSection {
    public String blockDir;
    public String nameFormat;
    public int blockCount;
    public String outFile;
    public HashSet<Integer> skipRanks;
    public boolean isBigEndian;

    public BlockMergerSection(String configurationFilePath) {
        Properties p = new Properties();
        try {
            p.load(new FileInputStream(configurationFilePath));
            blockDir = getProperty(p,"blockDir","/home/blockDir");
            nameFormat = getProperty(p,"nameFormat","{1}-{2}.bin");
            blockCount = Integer.valueOf(getProperty(p,"blockCount","1"));
            outFile = getProperty(p,"outDir","/home/outfile.bin");
            String skipRanksStr = getProperty(p,"skipRanks","");
            skipRanks = getHashSet(skipRanksStr);
            isBigEndian = Boolean.valueOf(getProperty(p,"isBigEndian","false"));
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
        if(ranks.equals("")) return set;
        String[] skips = ranks.split(",");
        for (String skip : skips) {
            set.add(Integer.valueOf(skip));
        }
        return set;
    }

}
