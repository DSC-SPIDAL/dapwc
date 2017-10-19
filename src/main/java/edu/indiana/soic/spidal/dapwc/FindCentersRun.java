package edu.indiana.soic.spidal.dapwc;

import com.google.common.base.Optional;
import edu.indiana.soic.spidal.configuration.sections.PairwiseClusteringSection;
import mpi.MPIException;
import org.apache.commons.cli.*;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.regex.Pattern;

/**
 * Created by pulasthi on 10/17/17.
 */
public class FindCentersRun {
    private static Options programOptions = new Options();
    static {
        programOptions.addOption(String.valueOf(Constants.CMD_OPTION_SHORT_C),Constants.CMD_OPTION_LONG_C, true,
                Constants.CMD_OPTION_DESCRIPTION_C);
        programOptions.addOption(String.valueOf(Constants.CMD_OPTION_SHORT_N),Constants.CMD_OPTION_LONG_N,true,
                Constants.CMD_OPTION_DESCRIPTION_N);
        programOptions.addOption(String.valueOf(Constants.CMD_OPTION_SHORT_T),Constants.CMD_OPTION_LONG_T,true,
                Constants.CMD_OPTION_DESCRIPTION_T);
    }

    public static void main(String[] args) throws MPIException {
        Optional<CommandLine> parserResult = parseCommandLineArguments(args, programOptions);
        if (!parserResult.isPresent()){
            System.out.println(Constants.ERR_PROGRAM_ARGUMENTS_PARSING_FAILED);
            new HelpFormatter().printHelp(Constants.PROGRAM_NAME, programOptions);
            return;
        }

        CommandLine cmd = parserResult.get();
        if (!(cmd.hasOption(Constants.CMD_OPTION_LONG_C) &&
                cmd.hasOption(Constants.CMD_OPTION_LONG_N) &&
                cmd.hasOption(Constants.CMD_OPTION_LONG_T))){
            System.out.println(Constants.ERR_INVALID_PROGRAM_ARGUMENTS);
            new HelpFormatter().printHelp(Constants.PROGRAM_NAME, programOptions);
            return;
        }

        //  Read Metadata using this as source of other metadata
        Program.ReadControlFile(cmd);

        //  Set up Parallelism
        //  Set up MPI and threads parallelism
        try {
            PWCParallelism.SetupParallelism(args);
        } catch (MPIException e) {
            PWCUtility.printException(e);
            return; // End program on error
        }

        // Set up Decomposition of USED points
        PWCParallelism.SetParallelDecomposition();

        // Initial Processing Complete
        PWCUtility.mpiOps.barrier(); // Make certain all processes have processed original data before writing updated
        //  read data into memory
        PWCParallelism.ReadDataFromFile(Program.config.DistanceMatrixFile);
        Program.ClusterAssignments = new int[PWCUtility.PointCount_Global]; // Set whenever clusters output
        readClusterAssignments(Program.config.AddDapwcFile);
        String file = "CenterFile-M" + (new Integer(Program.maxNcent)).toString() + "-C" + (new Integer(Program.maxNcent)).toString() + ".txt";
        String directory1 = (new File(Program.config.AddDapwcFile)).getParent();
        String CenterFileName = Paths.get(directory1, file).toString();
        String place = (new File(Program.config.DistanceMatrixFile)).getParent();
        String MDSFileName =PWCUtility.addMDSfile;
        String FullLabelFileName = PWCUtility.Labelfile;
        FindCenters.FindGroupCenters(Program.ClusterAssignments, Program.maxNcent, 0, CenterFileName, MDSFileName, FullLabelFileName);

        PWCParallelism.TearDownParallelism(); //  Finalize MPI

    }

    private static void readClusterAssignments(String addDapwcFile) {
        int numPointsRead = 0;
        Pattern pattern = Pattern.compile("\\s+");
        try(BufferedReader reader = Files.newBufferedReader(Paths.get(addDapwcFile))){
            String line;
            while ((line = reader.readLine()) != null){
                String[] splits = pattern.split(line.trim());
                int index = Integer.parseInt(splits[0]);
                Program.ClusterAssignments[index] = Integer.parseInt(splits[1]);
                numPointsRead++;
            }
            if(numPointsRead != Program.config.NumberDataPoints){
                PWCUtility.printAndThrowRuntimeException(
                        "Illegal count on MDS file " + numPointsRead + " " +
                                Program.config.NumberDataPoints);
            }
        }catch (IOException e){
            PWCUtility.printAndThrowRuntimeException("Error While reading PWC results file " + e.getMessage());
        }
    }


    private static Optional<CommandLine> parseCommandLineArguments(String [] args, Options opts){

        CommandLineParser optParser = new GnuParser();

        try {
            return Optional.fromNullable(optParser.parse(opts, args));
        } catch (ParseException e) {
            System.out.println(e);
        }
        return Optional.fromNullable(null);
    }
}
