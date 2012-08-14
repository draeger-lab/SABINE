/*
    SABINE predicts binding specificities of transcription factors.
    Copyright (C) 2009 ZBIT, University of Tübingen, Johannes Eichner

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

package model;


import help.FeatureFileGenerator;
import help.FileCopier;
import help.LibSVMFeatureScaler;
import help.RawDataPreprocessor;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.ArrayList;
import java.util.Calendar;
import java.util.StringTokenizer;

import main.FBPPredictor;

import core.DomainFeatureCalculator;
import core.FeatureRepairer;
import core.SVMPairwiseFeatureCalculator;
import core.SequenceFeatureCalculator;
import core.SpeciesFeatureCalculator;

public class ModelGenerator {
	
	String class_id = "";
	String train_dir = "/trainingsets";
	String base_dir = null;
	boolean silent = false;
	
	public ModelGenerator(String id) {
	  	class_id = id;
	}

	public ArrayList<String> computeFeatures(String feature_dir, String temp_dir) {
		
		/*
		 *  read TFs in training set
		 */
		
		DomainFeatureCalculator domaincalculator = new DomainFeatureCalculator();
		SequenceFeatureCalculator sequencecalculator = new SequenceFeatureCalculator();
		SpeciesFeatureCalculator speciescalculator = new SpeciesFeatureCalculator();
		SVMPairwiseFeatureCalculator svmpairwisecalculator = new SVMPairwiseFeatureCalculator();
		
		domaincalculator.silent = true;
		domaincalculator.basedir = temp_dir;
		sequencecalculator.silent = true;
		sequencecalculator.basedir = temp_dir;

		int num_entries = domaincalculator.getTrainingSetSize(class_id, train_dir);
	    boolean[] irrelevantPairs = new boolean[num_entries]; 
		  
	    domaincalculator.parseRelevantDomains(irrelevantPairs, class_id, train_dir);
	    sequencecalculator.parseRelevantDomainsAndSequences(irrelevantPairs, class_id, train_dir);
	    sequencecalculator.parseRelevantSecondaryStructures(irrelevantPairs, class_id, train_dir);
		
	    /*
	     *  find relevant pairs of TFs 
	     */
	    
	    
	    // initialize relevant pairs 
	    ArrayList<String> relevant_pairs = new ArrayList<String>();
	    for (int i=0; i<num_entries-1; i++) {
	    	for (int j=i+1; j<num_entries; j++) { 
	    		
	    		relevant_pairs.add(i + "_" + j);
	    	}
	    }
	    System.out.println("  Filtering transcription factor pairs.");
	    int num_TFpairs = relevant_pairs.size();
	    
	    /*
	     *  find TF pairs with BLOSUM score > 0.3
	     */
	    
	    domaincalculator.calculateAllDomainFeatures("SMBasedSimilarity", null, FBPPredictor.matrix_dir + "BLOSUM_62.dat"   , feature_dir + "all_pairs_BLOSUM62_score.out", relevant_pairs);
	    relevant_pairs = findRelevantPairs(relevant_pairs, feature_dir + "all_pairs_BLOSUM62_score.out", false);
	    System.out.println("    " + relevant_pairs.size() + " / " + num_TFpairs + " factor pairs have sufficient domain sequence similarity.");
	    num_TFpairs = relevant_pairs.size();
	    
	    /*
	     *  find TF pairs with MoSta score < 1
	     */
	    
	    /*
	    LabelFileGenerator labelgenerator = new LabelFileGenerator();
		labelgenerator.calculateLabelFile("trainingsets/trainingset_" + class_id + ".fbps", feature_dir + "similar_pairs_MoSta_score.out", relevant_pairs);
	    relevant_pairs = findRelevantPairs(relevant_pairs, feature_dir + "similar_pairs_MoSta_score.out", true);
	    System.out.println("    " + relevant_pairs.size() + " / " + num_TFpairs + " factor pairs with different matrices found (MoSta score < 1).\n");
	    */
	    
	    relevant_pairs = removeRedundantFactors(relevant_pairs, train_dir + "trainingset_" + class_id + ".fbps");
	    System.out.println("    " + relevant_pairs.size() + " / " + num_TFpairs + " factor pairs have different PFMs.\n");
	    
	    System.out.println("  Computing features.\n" +
	    				   "     (directory: " + feature_dir + ")\n");
	    /*
	    domaincalculator.calculateAllDomainFeatures("SMBasedSimilarity", null, FBPPredictor.matrix_dir + "BLOSUM_62.dat"   , feature_dir + "domain_scores_BLOSUM_62.out", relevant_pairs);
	    domaincalculator.calculateAllDomainFeatures("SMBasedSimilarity", null, FBPPredictor.matrix_dir + "PAM_080.dat"     , feature_dir + "domain_scores_PAM_080.out", relevant_pairs);
		
	    domaincalculator.calculateAllDomainFeatures("SMBasedSimilarity", null, FBPPredictor.matrix_dir + "PAM_010.dat"     , feature_dir + "domain_scores_PAM_010.out", relevant_pairs);
		domaincalculator.calculateAllDomainFeatures("SMBasedSimilarity", null, FBPPredictor.matrix_dir + "WEIL970101.dat"  , feature_dir + "domain_scores_WEIL970101.out", relevant_pairs);
		domaincalculator.calculateAllDomainFeatures("SMBasedSimilarity", null, FBPPredictor.matrix_dir + "MEHP950101.dat"  , feature_dir + "domain_scores_MEHP950101.out", relevant_pairs);
		domaincalculator.calculateAllDomainFeatures("SMBasedSimilarity", null, FBPPredictor.matrix_dir + "MEHP950102.dat"  , feature_dir + "domain_scores_MEHP950102.out", relevant_pairs);
		domaincalculator.calculateAllDomainFeatures("SMBasedSimilarity", null, FBPPredictor.matrix_dir + "MEHP950103.dat"  , feature_dir + "domain_scores_MEHP950103.out", relevant_pairs);
		domaincalculator.calculateAllDomainFeatures("SMBasedSimilarity", null, FBPPredictor.matrix_dir + "LUTR910102.dat"  , feature_dir + "domain_scores_LUTR910102.out", relevant_pairs);
		domaincalculator.calculateAllDomainFeatures("SMBasedSimilarity", null, FBPPredictor.matrix_dir + "NIEK910102.dat"  , feature_dir + "domain_scores_NIEK910102.out", relevant_pairs);
		domaincalculator.calculateAllDomainFeatures("SMBasedSimilarity", null, FBPPredictor.matrix_dir + "RISJ880101.dat"  , feature_dir + "domain_scores_RISJ880101.out", relevant_pairs);
		domaincalculator.calculateAllDomainFeatures("SMBasedSimilarity", null, FBPPredictor.matrix_dir + "MIYS930101.dat"  , feature_dir + "domain_scores_MIYS930101.out", relevant_pairs);
		domaincalculator.calculateAllDomainFeatures("SMBasedSimilarity", null, FBPPredictor.matrix_dir + "MIYT790101.dat"  , feature_dir + "domain_scores_MIYT790101.out", relevant_pairs);
		
		domaincalculator.calculateAllDomainFeatures("SequenceIdentity" , null, FBPPredictor.matrix_dir + "BLOSUM_62.dat"   , feature_dir + "domain_scores_BLOSUM_62_si.out", relevant_pairs);
	    
		domaincalculator.calculateAllDomainFeatures("SMBasedIdentity" , new String[] {"1.0"} , FBPPredictor.matrix_dir + "BLOSUM_62.dat"   , feature_dir + "domain_scores_BLOSUM_62_t=1.0.out", relevant_pairs);
		domaincalculator.calculateAllDomainFeatures("SMBasedIdentity" , new String[] {"3.0"} , FBPPredictor.matrix_dir + "BLOSUM_62.dat"   , feature_dir + "domain_scores_BLOSUM_62_t=3.0.out", relevant_pairs);
		domaincalculator.calculateAllDomainFeatures("SMBasedIdentity" , new String[] {"5.0"} , FBPPredictor.matrix_dir + "BLOSUM_62.dat"   , feature_dir + "domain_scores_BLOSUM_62_t=5.0.out", relevant_pairs);
		
		domaincalculator.calculateAllDomainFeatures("LocalAlignmentKernel", new String[] {"GCBopt.dat"} 		, null, feature_dir + "domain_scores_lak_GCBopt.out", relevant_pairs);
		domaincalculator.calculateAllDomainFeatures("LocalAlignmentKernel", new String[] {"JTTopt.dat"} 		, null, feature_dir + "domain_scores_lak_JTTopt.out", relevant_pairs);
		domaincalculator.calculateAllDomainFeatures("LocalAlignmentKernel", new String[] {"BLOSUM_62opt.dat"} 	, null, feature_dir + "domain_scores_lak_BLOSUM_62opt.out", relevant_pairs);
		domaincalculator.calculateAllDomainFeatures("LocalAlignmentKernel", new String[] {"PAM_250opt.dat"} 	, null, feature_dir + "domain_scores_lak_PAM_250opt.out", relevant_pairs);
		domaincalculator.calculateAllDomainFeatures("LocalAlignmentKernel", new String[] {"blosum62.dat"} 		, null, feature_dir + "domain_scores_lak_blosum62.out", relevant_pairs);
		
		domaincalculator.calculateAllDomainFeatures("MismatchKernel", new String[] {"4", "1"}, null, feature_dir + "domain_scores_mmk_4_1.out", relevant_pairs);
		domaincalculator.calculateAllDomainFeatures("MismatchKernel", new String[] {"5", "1"}, null, feature_dir + "domain_scores_mmk_5_1.out", relevant_pairs);
		domaincalculator.calculateAllDomainFeatures("MismatchKernel", new String[] {"6", "1"}, null, feature_dir + "domain_scores_mmk_6_1.out", relevant_pairs);
		
		sequencecalculator.calculateAllSequenceFeatures("SecondaryStructure", null, FBPPredictor.matrix_dir + "BLOSUM_62.dat", feature_dir + "domain_scores_secstr_blo62.out", relevant_pairs);

		sequencecalculator.calculateAllSequenceFeatures("Environments", new String[] {"25"}, FBPPredictor.matrix_dir + "BLOSUM_62.dat", feature_dir + "domain_scores_env_25_BLOSUM_62.out", relevant_pairs);
		sequencecalculator.calculateAllSequenceFeatures("Environments", new String[] {"50"}, FBPPredictor.matrix_dir + "BLOSUM_62.dat", feature_dir + "domain_scores_env_50_BLOSUM_62.out", relevant_pairs);
		*/
		speciescalculator.calculateAllPhylogeneticDistances(class_id, irrelevantPairs, train_dir + "new_phylogenetic_distances.out", feature_dir + "domain_scores_phyl_dist.out", relevant_pairs, train_dir);
		
		svmpairwisecalculator.calculateAllSVMPairwiseScores(train_dir + "trainingset_" + class_id + ".blo62", feature_dir + "domain_scores_svm_pairwise_BLOSUM_62.out", relevant_pairs);
		svmpairwisecalculator.calculateAllSVMPairwiseScores(train_dir + "trainingset_" + class_id + ".pam80", feature_dir + "domain_scores_svm_pairwise_PAM_080.out", relevant_pairs);
		
		return relevant_pairs;
	}
	
	public void repairFeatures(String feature_dir) {
		FeatureRepairer repairer = new FeatureRepairer();
		repairer.silent = true;
		repairer.repairAllFeatureFiles(feature_dir.substring(0,feature_dir.length()-1));
	}
	
	public void computeLabels(String label_dir, String temp_dir, ArrayList<String> relevant_pairs) {
		LabelFileGenerator labelfilegenerator = new LabelFileGenerator();
		labelfilegenerator.basedir = temp_dir;
		labelfilegenerator.calculateLabelFile(train_dir + "trainingset_" + class_id + ".fbps", label_dir + "mosta_labels_" + class_id + ".out", relevant_pairs);
	}
	
	public void generateSVMinput(String feature_dir, String label_dir, String libsvm_dir) {
		FeatureFileGenerator libsvmfilegenerator = new FeatureFileGenerator();
		libsvmfilegenerator.silent = true;
		libsvmfilegenerator.generateLibSVMTrainingSet(feature_dir.substring(0, feature_dir.length()-1), 
													  train_dir, 
													  class_id,
													  label_dir + "mosta_labels_" + class_id + ".out", 
													  libsvm_dir + "trainingset_unscaled.lp." + class_id + ".att");
		
		// copy unscaled libsvm file to training directory
		FileCopier.copy(libsvm_dir + "trainingset_unscaled.lp." + class_id + ".att", train_dir + "feature_scaling_orientation_" + class_id + ".out");
	
		// scale libsvm file
		LibSVMFeatureScaler scaler = new LibSVMFeatureScaler();
		scaler.scaleFeatureFile(libsvm_dir + "trainingset_unscaled.lp." + class_id + ".att", 
								libsvm_dir + "trainingset_unscaled.lp." + class_id + ".att",
								libsvm_dir + "trainingset.lp." + class_id + ".att");
	}
	
	public void runCrossvalidation(String infile) {
						
					
			/*
			 *  write shell script
			 */
		
			// grid range and stepping for svm-param c
				  
				int c_start = -2;
				int c_stop  =  6;
				int c_incr  =  2;
					
					
			// grid range and stepping for kernel-param gamma
					  
				int g_start = -10;
				int g_stop  =  2;
				int g_incr  =  2;	
					
					
			// grid range and stepping for svm-param epsilon
					
				int e_start = -8;
				int e_stop  = -6;
				int e_incr  =  2;
					
					
			// number of folds for cross-validation
					
				int folds   =  5;
					
	
			// number of multiruns for smoothing of cv-results
					
				int runs 	=  1;
				
				String sabdir = System.getProperty("user.dir");
				sabdir = sabdir.substring(sabdir.lastIndexOf("/") + 1) + "/";
				
				String cmdString =	  //"/usr/lib/jvm/java-1.5.0-sun-1.5.0.14/bin/"
									  "java -cp libsvm.jar:"
					          	 	+ "lib/cli/commons-cli-1.1.jar:"
					          	 	+ "lib/math/commons-math-1.2.jar "
					          	 	+ "main/CVToolRunner -g "
									+ "-f ../" + sabdir + infile + " "
									+ "-cs " + c_start + " "
									+ "-ce " + c_stop  + " "
									+ "-ci " + c_incr  + " "
									+ "-gs " + g_start + " "
									+ "-ge " + g_stop  + " "
									+ "-gi " + g_incr  + " "
									+ "-es " + e_start + " "
									+ "-ee " + e_stop  + " "
									+ "-ei " + e_incr  + " "
									+ "-k "  + folds   + " "
									+ "-r "  + runs;
				
				
				try {
				
					String line = null;
				
					BufferedWriter bw = new BufferedWriter(new FileWriter(new File("runCV.sh")));
					
					// write script
					bw.write("#!/bin/bash\n\n" +
							"cd ../CVTool/\n" + 
							cmdString + "\n" +
							"mv " + infile.substring(infile.lastIndexOf("/")+1) + ".dat ../" + sabdir + infile + ".dat\n" +
							"cd ../"+ sabdir + "\n");
					
					
					bw.flush();
					bw.close();
				
				/*
				 *  run shell script
				 */	
				
					Process proc = Runtime.getRuntime().exec("sh runCV.sh");
					
					BufferedReader input = new BufferedReader(new InputStreamReader(proc.getErrorStream()));
				     
					while ((line = input.readLine()) != null) {
						
				        System.out.println(line);
				    
					}
					
					input.close();
					proc.waitFor();
					
				}
				catch(IOException e) {
					System.out.println("IOException while executing Script!");
					e.printStackTrace();
				} 
				catch (InterruptedException e) {
					System.out.println("InterruptedException while executing Shell Script!");
				}
		}
	
	
	/*
	 *  remove TF pairs with BLOSUM score < 0.3 or MoSta score = 1 
	 */
	
	public void filterSVMinput(String infile, String outfile) {
		
		try {
			
			/*
			 *  identify relevant TF pairs
			 */
			
			BufferedReader br = new BufferedReader(new FileReader(new File(infile)));
			StringTokenizer strtok;
			String line, label, blosum_score;
			ArrayList<String> relevant_lines = new ArrayList<String>();
			
			while((line = br.readLine()) != null) {
				
				strtok = new StringTokenizer(line);
				label = strtok.nextToken();
				blosum_score = strtok.nextToken();
				
				if (Double.parseDouble(label) < 1 && Double.parseDouble(blosum_score.substring(2)) > 0.3) {
					relevant_lines.add(line);
				} 
			}
			br.close();
						
			/*
			 *  write filtered libsvm file
			 */
			
			BufferedWriter bw = new BufferedWriter(new FileWriter(new File(outfile)));
			
			for (int i=0; i<relevant_lines.size(); i++) {
				bw.write(relevant_lines.get(i) + "\n");
			}
			bw.flush();
			bw.close();
		}
		catch (IOException e) {
			System.out.println(e.getMessage());
			System.out.println("Exception occured while filtering LibSVM file.");
		}
	}
	
	
	public double[] parseCVResults(String infile, String outfile) {
		
		String line = null;
		
		StringTokenizer strtok = null;
		
		double c = 0.0;
		double g = 0.0;
		double e = 0.0;
		
		double mse 			= 0.0;
		double mse_stdev 	= 0.0;
		double aae 			= 0.0;
		double aae_stdev 	= 0.0;
		double pcc 			= 0.0;
		double pcc_stdev 	= 0.0;
		
		
		double best_c = 0.0;
		double best_g = 0.0;
		double best_e = 0.0;
		
		double best_mse 		= Double.POSITIVE_INFINITY;
		double best_mse_stdev 	= 0.0;
		double best_aae 		= Double.POSITIVE_INFINITY;
		double best_aae_stdev 	= 0.0;
		double best_pcc 		= Double.NEGATIVE_INFINITY;
		double best_pcc_stdev 	= 0.0;
		
		
		try {
			
			
			/*
			 * go over all param combinations
			 */
			
			BufferedReader br = new BufferedReader(new FileReader(new File(infile)));
			
			while((line = br.readLine()) != null && line.startsWith("/")) {
				
				strtok = new StringTokenizer(line);
				
				strtok.nextToken();										// input filename
				
				c = Double.parseDouble(strtok.nextToken());
				g = Double.parseDouble(strtok.nextToken());
				e = Double.parseDouble(strtok.nextToken());
				
				mse 		= Double.parseDouble(strtok.nextToken());
				mse_stdev 	= Double.parseDouble(strtok.nextToken());
				aae 		= Double.parseDouble(strtok.nextToken());
				aae_stdev 	= Double.parseDouble(strtok.nextToken());
				pcc 		= Double.parseDouble(strtok.nextToken());
				pcc_stdev 	= Double.parseDouble(strtok.nextToken());
				
				
				if(aae < best_aae) {
					
					
					best_c = c;
					best_g = g;
					best_e = e;
					
					best_mse 		= mse;
					best_mse_stdev 	= mse_stdev;
					best_aae 		= aae;
					best_aae_stdev 	= aae_stdev;
					best_pcc 		= pcc;
					best_pcc_stdev 	= pcc_stdev;
					
					
				}
				
			}
			
			BufferedWriter bw = new BufferedWriter(new FileWriter(new File(outfile)));
					
			bw.write("AAE = " + best_aae + "\tQ^2 = " + best_pcc + "\tMSE = " + best_mse + "\tfor (c,g,e) = (" + best_c + "," + best_g + "," + best_e + ")\n");
			bw.flush();
			bw.close();
			
			if(!silent) {
				
				DecimalFormat fmt = new DecimalFormat();
				DecimalFormatSymbols symbs = fmt.getDecimalFormatSymbols();
		        symbs.setDecimalSeparator('.');
		        fmt.setDecimalFormatSymbols(symbs);
		        fmt.setMaximumFractionDigits(2);
		        fmt.setMinimumFractionDigits(2);
				
				System.out.println("    MSE       = " + fmt.format(best_mse));
				//System.out.println("    std(MSE)  = " + fmt.format(best_mse_stdev));
				System.out.println("    AAE       = " + fmt.format(best_aae));
				//System.out.println("    std(AAE)  = " + fmt.format(best_aae_stdev));
				System.out.println("    PCC       = " + fmt.format(best_pcc));
				//System.out.println("    std(PCC)  = " + fmt.format(best_pcc_stdev));
				System.out.println();
			
			}
		}
		catch(IOException ioe) {
			System.out.println("IOException occurred while parsing results of this festure.");
			System.out.println(ioe.getMessage());
		}
		
		return new double[] {best_c , best_g , best_e, best_aae, best_pcc , best_mse};
	}
	
	
	
	public void trainSVM(String infile, String outfile, double[] params) {
		
		
		
        if (! silent) {
        	
        	DecimalFormat fmt = new DecimalFormat();
    		DecimalFormatSymbols symbs = fmt.getDecimalFormatSymbols();
            symbs.setDecimalSeparator('.');
            fmt.setDecimalFormatSymbols(symbs);
            fmt.setMaximumFractionDigits(2);
            fmt.setMinimumFractionDigits(2);
        
            System.out.println("\nTraining SVR with params ( c , g , e ) = ( " + 
				 fmt.format(params[0]) + " , " + 
				 fmt.format(params[1]) + " , " +  
				 fmt.format(params[2]) + " ) .\n");
		
		
            System.out.println("  Training Set : " + infile.substring(infile.lastIndexOf("/") + 1));
            System.out.println("  Model File   : " + outfile.substring(outfile.lastIndexOf("/") + 1) + "\n");
		
            System.out.println("Expected Performance : \n");
            System.out.println("  AAE = " +  fmt.format(params[3]));
            System.out.println("  MSE = " +  fmt.format(params[5]));
            System.out.println("  Q^2 = " +  fmt.format(params[4]) + "\n");
        }
            
		String[] args = new String[12];
		
		args[0] = "-s";		
		args[1] = "3";
		args[2] = "-t";
		args[3] = "2";
		args[4] = "-g";
		args[5] = "" + params[1];
		args[6] = "-c";
		args[7] = "" + params[0];
		args[8] = "-p";
		args[9] = "" + params[2];
		args[10] = infile;
		args[11] = outfile;
		
		try {
			libsvm.svm_train.main(args);
		}
		catch(IOException ioe) {
			System.out.println(ioe.getMessage());
		}
	}
	
	
	public ArrayList<String> findRelevantPairs(ArrayList<String> relevant_pairs, String infile, boolean filterLabels) {
		
		String line;
		StringTokenizer strtok;
		double curr_score;
		boolean[] rel_positions = new boolean[relevant_pairs.size()];
		int entry_counter = 0;
		
		try {
			
			BufferedReader br = new BufferedReader(new FileReader(new File(infile)));
			
			while ((line = br.readLine()) != null) {
				
				strtok = new StringTokenizer(line, ":");
				strtok.nextToken();
				
				curr_score = Double.parseDouble(strtok.nextToken().trim());
				
				if      (!filterLabels && curr_score > 0.3) rel_positions[entry_counter] = true;
				else if ( filterLabels && curr_score < 1.0) rel_positions[entry_counter] = true;
				entry_counter++;
			}
			br.close();
			
			for (int i=rel_positions.length-1; i>=0; i--) {
				if (! rel_positions[i]) relevant_pairs.remove(i);
			}
		}
		
		catch(IOException ioe) {
			System.out.println(ioe.getMessage());
		}
		
		return relevant_pairs;
	}
	
	/*
	 * removes factor pairs with identical matrices
	 */
	
	public ArrayList<String> removeRedundantFactors(ArrayList<String> relevant_pairs, String infile) {
		
		LabelFileGenerator labelgenerator = new LabelFileGenerator();
		ArrayList<ArrayList<String>> all_matrices = labelgenerator.getAllFBPs(infile);
		boolean[] rel_positions = new boolean[relevant_pairs.size()];
		int entry_counter = 0;
		
		for (int i=0; i<all_matrices.size()-1; i++) {
			for (int j=i+1; j<all_matrices.size(); j++) {
				
				if (! relevant_pairs.contains(i + "_" + j)) continue;
				if (! all_matrices.get(i).equals(all_matrices.get(j))) rel_positions[entry_counter] = true;
				entry_counter++;
			}
		}
		
		for (int i=rel_positions.length-1; i>=0; i--) {
			if (! rel_positions[i]) relevant_pairs.remove(i);
		}
		
		return relevant_pairs;
	}
	
	public void generateTrainingSet(String input_file) {
		RawDataPreprocessor trainset_generator = new RawDataPreprocessor();
		trainset_generator.runPreprocessor(input_file, class_id, train_dir);
	}
	
	public void generateModelFile() {
		
		/*
		 *  set base directory
		 */
		
		String curr_dir = System.getProperty("user.dir") + "/";
		
		if (base_dir == null) {
		
			DecimalFormat fmt = new DecimalFormat();
			fmt.setMaximumIntegerDigits(2);
			fmt.setMinimumIntegerDigits(2);
		
			// get current time and date
			Calendar cal = Calendar.getInstance ();
			String curr_date = (fmt.format(cal.get(Calendar.DAY_OF_MONTH)) + "." + 
								fmt.format((cal.get(Calendar.MONTH) + 1)) + "." + 
								cal.get(Calendar.YEAR));
			String curr_time = (fmt.format(cal.get(Calendar.HOUR_OF_DAY)) + ":" +
								fmt.format(cal.get(Calendar.MINUTE)));
		
		
			base_dir = curr_dir + "svr_models/" + curr_date + "_" + curr_time + "/";
			
			if (! new File(base_dir).mkdir()) {
				System.out.println("\nInvalid base directory. Aborting.");
				System.out.println("Base directory: " + base_dir + "\n");
				System.exit(0);
			}
		}
		
		if (! base_dir.endsWith("/")) base_dir += "/";
		
		String feature_dir = base_dir + "feature_files/";
		String label_dir = base_dir + "label_files/";
		String libsvm_dir = base_dir + "libsvm_files/";
		String temp_dir = base_dir + "temp/";
		
		String mmkernel_dir = temp_dir + "mismatchkernel/";
		String matlign_dir = temp_dir + "matlign/";
		String mosta_dir = temp_dir + "mosta/";
		String stamp_dir = temp_dir + "stamp/";
		String psipred_dir = temp_dir + "psipred/";
		
		if ((! new File(feature_dir).exists() && ! new File(feature_dir).mkdir()) | 
			(! new File(label_dir).exists() && ! new File(label_dir).mkdir()) |
			(! new File(libsvm_dir).exists() && ! new File(libsvm_dir).mkdir()) |
			(! new File(temp_dir).exists() && ! new File(temp_dir).mkdir()) |
			
			(! new File(mmkernel_dir).exists() && ! new File(mmkernel_dir).mkdir()) |
			(! new File(matlign_dir).exists() && ! new File(matlign_dir).mkdir()) | 
			(! new File(mosta_dir).exists() && ! new File(mosta_dir).mkdir()) | 
			(! new File(stamp_dir).exists() && ! new File(stamp_dir).mkdir()) | 
			(! new File(psipred_dir).exists() && ! new File(psipred_dir).mkdir())) {
			
			System.out.println("\nInvalid base directory. Aborting.");
			System.out.println("Base directory: " + base_dir + "\n");
			System.exit(0);
		}

		
		/*
		 *  compute feature files
		 */
		ArrayList<String> relevant_pairs = computeFeatures(feature_dir, temp_dir);
		repairFeatures(feature_dir);
		
		/*
		 *  compute label file containing normalized MoSta scores of all TF pairs
		 */
		System.out.	println("  Computing labels.\n" + 
				   		   "     (directory: " + label_dir + ")\n");
		computeLabels(label_dir, temp_dir, relevant_pairs);
		
		System.out.println("  Writing LibSVM input file.\n" + 
		   		   		   "     (file: " + libsvm_dir + "trainingset.lp." + class_id + ".att" + ")\n");
		generateSVMinput(feature_dir, label_dir, libsvm_dir);
		
		System.out.println("  Performing cross-validation.\n");
		
		String pwd = System.getProperty("user.dir");
		String rel_libsvm_dir = libsvm_dir;
		if (libsvm_dir.startsWith("/") && libsvm_dir.startsWith(pwd)) {
			rel_libsvm_dir = libsvm_dir.substring(pwd.length());
		}
		
		runCrossvalidation(rel_libsvm_dir + "trainingset.lp." + class_id + ".att");
		
		
		double[] params = parseCVResults(libsvm_dir + "trainingset.lp." + class_id + ".att.dat",
										 libsvm_dir + "trainingset.lp." + class_id + ".par");
		
		System.out.println("  Training SVR model.\n" +
						   "     (file: " + base_dir + "mosta.lp." + class_id + ".model" + ")\n");
		
		this.silent = true;
		trainSVM(libsvm_dir + "trainingset.lp." + class_id + ".att", 
				 base_dir + "mosta.lp." + class_id + ".model", params);
		
		System.out.print("\r");
		
	}
	
	public static void main(String[] args) {
		
		ModelGenerator svr_model = new ModelGenerator("class0");
		int superclass = 0;
		
		/*
		 * parse options
		 */
		
		if (args.length == 0 || (args.length == 1 && (args[0].equals("-help") || args[0].equals("--help")))) {
			svr_model.usage();
		}
		
		String input_file = args[0];
		
		for(int i=1; i<args.length-1; i+=2) {
			
			if(args[i].equals("-s")) { superclass		   	= Integer.parseInt(args[i+1]); 	continue; }
			if(args[i].equals("-m")) { svr_model.base_dir	= args[i+1]; 					continue; }
			if(args[i].equals("-t")) { svr_model.train_dir 	= args[i+1]; 					continue; }
			
			if( !args[i].equals("-s") && !args[i].equals("-t") && !args[i].equals("-m") ) {	
				System.out.println("\n  Invalid argument: " + args[i]);
				svr_model.usage();
			}
		}
		
		
		System.out.println();
		System.out.println("  ---------------------------------------------------------------------");
		System.out.println("  Model Generator for SABINE - StandAlone BINding specificity Estimator");
		System.out.println("  ---------------------------------------------------------------------");
		System.out.println("\n");
		System.out.println("  Superclass: " + superclass);
		System.out.println("\n");
		
		
		if (! svr_model.train_dir.endsWith("/")) {
			svr_model.train_dir += "/";
		}
		
		svr_model.class_id = "class" + superclass;
		//svr_model.generateTrainingSet(input_file);
		svr_model.generateModelFile();
	}
	
	private void usage() {
		
		System.out.println();
		System.out.println("  ---------------------------------------------------------------------");
		System.out.println("  Model Generator for SABINE - StandAlone BINding specificity Estimator");
		System.out.println("  ---------------------------------------------------------------------");
		System.out.println("\n");
		System.out.println("  Usage   : modelgen <input_file> [OPTIONS]\n");
		System.out.println("  OPTIONS : -s <superclass>       (TRANSFAC superclass)                    default = 0");
		System.out.println("            -t <train_directory>  (directory to save the training set)");
		System.out.println("            -m <model_directory>  (directory to save the model file) \n\n");
		
		System.exit(0);
		
	}

}


