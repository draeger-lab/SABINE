import gui.SequenceLogo;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.StringTokenizer;

import org.biojava.bio.gui.DistributionLogo;


public class SABINE_Validator {

	
	public void verifyInstallation() {
		
		// set parameters
		FBPPredictor predictor = new FBPPredictor();
		predictor.dynamic_threshold = true;
		predictor.best_match_threshold = 0.95;
		predictor.max_number_of_best_matches = 5;
		predictor.outlier_filter_threshold = 0.5;
		
		// generate base directory and subdirectories
		String base_dir = "tmp/validator_tmp/";
		File base_dir_path = new File(base_dir);
		if (!base_dir_path.exists() && !base_dir_path.mkdir()) {
			System.out.println("\nInvalid base directory. Aborting.");
			System.out.println("Base directory: " + base_dir + "\n");
			System.exit(0);
		}
		SABINE_Caller dir_creator = new SABINE_Caller();
		dir_creator.createTempDirectories(base_dir);
		
		// run SABINE on generated input file
		System.out.print("Checking installation of SABINE...");
		predictor.silent = true;
		predictor.predictFBP("input/test.tf", base_dir, "trainingsets_public/", null);
		System.out.println("done.\n");
		// read output
		boolean valid = compareOutfiles("input/test.tf.out", base_dir + "prediction.out");
		
		if (valid) {
			System.out.println("The installation was successfully validated.");
			System.out.println("Type \"sabine.sh --gui\" to start the graphical interface or use the " +
							   "command \"sabine.sh --help\" to view a list of the provided command line options.");
		}
		else {
			System.out.println("The installation of the tool is invalid.");
			System.out.println("Please use the webservice version of SABINE (http://webservices.cs.uni-tuebingen.de/?tool=sabine).");
		}
		System.out.println("See the \"readme.txt\" or visit the SABINE website (http://www.ra.cs.uni-tuebingen.de/software/SABINE/) for a " +
		                   "comprehensive documentation of the tool.");
	}
	
	public boolean compareOutfiles(String reffile, String outfile) {
		
		try {
			BufferedReader br1 = new BufferedReader(new FileReader(new File(reffile)));
			BufferedReader br2 = new BufferedReader(new FileReader(new File(outfile)));
			
			br1.readLine(); 	// NA
			br1.readLine();		// XX
			
			String line1 = br1.readLine();
			String line2 = br2.readLine();
			
			// compare computed matrix to reference matrix
			while (line1 != null && !line1.startsWith("//")) {
				
				// computed output file differs from reference output file
				if (line2 == null || !line1.equals(line2)) {
					System.out.println("Reference: " + line1);
					System.out.println("Result:    " + line2);
					return false;
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
