package edu.indiana.soic.spidal.configuration;

import edu.indiana.soic.spidal.configuration.sections.PairwiseClusteringSection;

public class ConfigurationMgr {
    private String configurationFilePath;
    public PairwiseClusteringSection pairwiseClusteringSection;

    public ConfigurationMgr(String configurationFilePath) {
        this.configurationFilePath = configurationFilePath;
        pairwiseClusteringSection = new PairwiseClusteringSection(configurationFilePath);


    }

    public static ConfigurationMgr LoadConfiguration(String configurationFilePath){
        // TODO - Fix configuration management
        return new ConfigurationMgr(configurationFilePath);
    }
}
