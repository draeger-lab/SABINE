package main;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.StringTokenizer;

public class SABINE_Validator {

	public static final String validatorInputFile = "input/test.tf";
	public static final String validatorOutputFile = "input/test.tf.out";
	public static final String validatorOutputFile64 = "input/test.tf.out64";
	public static final String validatorBaseDir = "tmp/validator_tmp/";
	
	public void verifyInstallation() {
		
		// set parameters
		FBPPredictor.maxPFMsimilarity = Double.MAX_VALUE;
		FBPPredictor predictor = new FBPPredictor();
		predictor.dynamic_threshold = true;
		predictor.best_match_threshold = 0.95;
		predictor.max_number_of_best_matches = 5;
		predictor.outlier_filter_threshold = 0.5;
		
		
		// generate base directory and subdirectories
		File base_dir_path = new File(validatorBaseDir); 
		if (!base_dir_path.exists() && !base_dir_path.mkdir()) {
			System.out.println("\nInvalid base directory. Aborting.");
			System.out.println("Base directory: " + validatorBaseDir + "\n");
			System.exit(0);
		}
		SABINE_Main dir_creator = new SABINE_Main();
		dir_creator.createTempDirectories(validatorBaseDir);
		
		// run SABINE on generated input file
		System.out.print("Checking installation of SABINE...");
		predictor.silent = true;
		predictor.predictFBP(validatorInputFile, validatorBaseDir, FBPPredictor.public_trainingset, FBPPredictor.defaultModelDir);
		System.out.println("done.\n");
		// read output
		boolean valid = compareOutfiles(validatorOutputFile, validatorBaseDir + "prediction.out");
		boolean valid64 = compareOutfiles(validatorOutputFile64, validatorBaseDir + "prediction.out");
		
		
		if (valid || valid64) {
			System.out.println("The installation was successfully validated.");
			System.out.println("Type \"sabine.sh --gui\" to start the graphical interface or use the " +
							   "command \"sabine.sh --help\" to view a list of the provided command line options.");
		}
		else {
			System.out.println("The installation of the tool is invalid.");
			System.out.println("Please use the webservice version of SABINE (http://webservices.cs.uni-tuebingen.de/?tool=sabine).");
		}
		System.out.println("See the \"readme.txt\" or visit the SABINE website (http://www.cogsys.cs.uni-tuebingen.de/software/SABINE/) for a " +
		                   "comprehensive documentation of the tool.");
	}
	
	public boolean compareOutfiles(String reffile, String outfile) {
		
		try {
			BufferedReader br1 = new BufferedReader(new FileReader(new File(reffile)));
			BufferedReader br2 = new BufferedReader(new FileReader(new File(outfile)));
			
			br1.readLine(); 	// NA
			br1.readLine();		// XX
			
			StringTokenizer strtok1 = null;
			StringTokenizer strtok2 = null;
			
			String line1 = br1.readLine();
			String line2 = br2.readLine();
			
			// compare computed matrix to reference matrix
			while (line1 != null && !line1.startsWith("//")) {
				
				if (line2 == null) {
					return false;
				}
				
				// line in computed output file differs from reference output file
				if (!line1.equals(line2)) {
					
					// test if tokens are equal
					boolean equalLine = true;
					strtok1 = new StringTokenizer(line1);
					strtok2 = new StringTokenizer(line2);
					while(strtok1.hasMoreTokens() || strtok2.hasMoreTokens()) {
						if (!strtok1.nextToken().equals(strtok2.nextToken())) {
							equalLine = false;
							break;
						}
					}
					return equalLine;
				}
				line1 = br1.readLine();
				line2 = br2.readLine();
			}
		} catch (IOException e) {
			e.printStackTrace();
		}
		return true;
	}
	
	public static void main(String[] args) {
		
		SABINE_Validator validator = new SABINE_Validator();
		validator.verifyInstallation();
	}
	
}
