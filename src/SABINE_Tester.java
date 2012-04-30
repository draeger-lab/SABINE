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


import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.ArrayList;
import java.util.Random;

import libsvm.LabeledTF;

import extension.PredictionEvaluator;
import extension.SABINEInputOutputParser;
import extension.SABINEParser;
import extension.SABINETrainingsetParser;

public class SABINE_Tester {
	
	SABINEParser[] factors;
	ArrayList<ArrayList<Integer>> samples = new ArrayList<ArrayList<Integer>>();
	
	int superclass = 0;
	int[] other_classes = new int[4];
	int samplesize = 20;
	boolean matrices_available = true;
	
	
	public SABINE_Tester() {
		
		SABINEParser class0_factors;
		SABINEParser class1_factors;
		SABINEParser class2_factors;
		SABINEParser class3_factors;
		SABINEParser class4_factors ;
		
		class0_factors = new SABINETrainingsetParser();
		class1_factors = new SABINETrainingsetParser();
		class2_factors = new SABINETrainingsetParser();
		class3_factors = new SABINETrainingsetParser();
		class4_factors = new SABINETrainingsetParser();
		
		class0_factors.parseAll(0, null, null);
		class1_factors.parseAll(1, null, null);
		class2_factors.parseAll(2, null, null);
		class3_factors.parseAll(3, null, null);
		class4_factors.parseAll(4, null, null);
		
		factors = new SABINEParser[] { class0_factors , class1_factors , class2_factors , 
									   class3_factors , class4_factors };
	}
	
	public SABINE_Tester(String inputfile) {
		
		SABINEParser class0_factors;
		SABINEParser class1_factors;
		SABINEParser class2_factors;
		SABINEParser class3_factors;
		SABINEParser class4_factors ;
		
		class0_factors = new SABINEInputOutputParser();
		class1_factors = new SABINEInputOutputParser();
		class2_factors = new SABINEInputOutputParser();
		class3_factors = new SABINEInputOutputParser();
		class4_factors = new SABINEInputOutputParser();
		
		class0_factors.parseAll(0, inputfile, null);
		class1_factors.parseAll(1, inputfile, null);
		class2_factors.parseAll(2, inputfile, null);
		class3_factors.parseAll(3, inputfile, null);
		class4_factors.parseAll(4, inputfile, null);
		
		factors = new SABINEParser[] { class0_factors , class1_factors , class2_factors , 
									   class3_factors , class4_factors };
		
		matrices_available = false;
	}
	
	
	public SABINE_Tester(String inputfile, String labelfile) {
		
		SABINEParser class0_factors;
		SABINEParser class1_factors;
		SABINEParser class2_factors;
		SABINEParser class3_factors;
		SABINEParser class4_factors ;
		
		class0_factors = new SABINEInputOutputParser();
		class1_factors = new SABINEInputOutputParser();
		class2_factors = new SABINEInputOutputParser();
		class3_factors = new SABINEInputOutputParser();
		class4_factors = new SABINEInputOutputParser();
		
		class0_factors.parseAll(0, inputfile, labelfile);
		class1_factors.parseAll(1, inputfile, labelfile);
		class2_factors.parseAll(2, inputfile, labelfile);
		class3_factors.parseAll(3, inputfile, labelfile);
		class4_factors.parseAll(4, inputfile, labelfile);
		
		factors = new SABINEParser[] { class0_factors , class1_factors , class2_factors , 
									   class3_factors , class4_factors };
	}
	
	public void drawSamples() {
		
		if (factors[superclass].get_tf_names().isEmpty()) {
			System.out.println("Input file contains no factors of class " + superclass + ". Aborting.");
			System.exit(0);
		}
		
		Random rand = new Random();
		rand.setSeed(15);
		
		for (int i=0; i<5; i++) samples.add(new ArrayList<Integer>());
		ArrayList<Integer> curr_indices;
		
		ArrayList<Integer> other_class_list = new ArrayList<Integer>();
		for (int i=0; i<5; i++) { 
			if (i != superclass && !factors[i].get_tf_names().isEmpty()) other_class_list.add(i); 
		}
		other_classes = new int[other_class_list.size()];
		for (int i=0; i<other_class_list.size(); i++) {
			other_classes[i] = other_class_list.get(i);
		}
		
		int other_class;
		
		for (int i=0; i<samplesize; i++) {
			
			// randomly choose factor of given superclass
			curr_indices = samples.get(superclass);
			curr_indices.add(rand.nextInt(factors[superclass].get_tf_names().size()));
			samples.set(superclass, curr_indices);
			
			// randomly choose other superclass
			other_class = other_classes[rand.nextInt(other_classes.length)];
			
			// randomly choose factor of this class
			curr_indices = samples.get(other_class);
			curr_indices.add(rand.nextInt(factors[other_class].get_tf_names().size()));
			samples.set(other_class, curr_indices);
		}
	}
	
	public void writeOutfiles(String outdir) {
	
	  	try {
			BufferedWriter bw_pos = new BufferedWriter(new FileWriter(new File(outdir + "sample_pos.input")));
			
			/*
			 *  write output file for given superclass
			 */
			
			int[] indices = new int[samples.get(superclass).size()];  
			for (int i=0; i<samples.get(superclass).size(); i++) indices[i] = samples.get(superclass).get(i);
			
			factors[superclass].writeInputFile(bw_pos, indices, superclass);
			bw_pos.flush();
			bw_pos.close();
			
			// write labelfile
			ArrayList<String> curr_matrix;
			
			if (matrices_available) {
				
				bw_pos = new BufferedWriter(new FileWriter(new File(outdir + "sample_pos.labels")));
			
				for (int i=0; i<indices.length; i++) {
				
					curr_matrix = factors[superclass].get_matrices().get(indices[i]);
				
					for (int j=0; j<curr_matrix.size(); j++) {
						bw_pos.write("MA  " + curr_matrix.get(j) + "\n");
					}
					bw_pos.write("XX  \n");
				}
				bw_pos.flush();
				bw_pos.close();
			}
			
			/*
			 *  write output file for other classes (with wrong "CL"-entry )
			 */
			
			BufferedWriter bw_neg = new BufferedWriter(new FileWriter(new File(outdir + "sample_neg.input")));
			
			for (int c=0; c<other_classes.length; c++) {
				indices = new int[samples.get(other_classes[c]).size()];  
				
				for (int i=0; i<samples.get(other_classes[c]).size(); i++) {
					indices[i] = samples.get(other_classes[c]).get(i);
				}
				if (indices.length > 0) {
					factors[other_classes[c]].writeInputFile(bw_neg, indices, superclass);
				}
			}
			bw_neg.flush();
			bw_neg.close();
			
			
			// write labelfile
			if (matrices_available) {
				bw_neg = new BufferedWriter(new FileWriter(new File(outdir + "sample_neg.labels")));
			
				for (int c=0; c<other_classes.length; c++) {
					for (int i=0; i<samples.get(other_classes[c]).size(); i++) {
					
						curr_matrix = factors[other_classes[c]].get_matrices().get(samples.get(other_classes[c]).get(i));
					
						for (int j=0; j<curr_matrix.size(); j++) {
							bw_neg.write("MA  " + curr_matrix.get(j) + "\n");
						}	
						bw_neg.write("XX  \n");
					}	
				}
				bw_neg.flush();
				bw_neg.close();
			}
			
		} catch (IOException e) {
			System.out.println(e.getMessage());
			System.out.println("IOException occurred while writing output file for SABINE.");
		}
	}
	
	public void readOutfile(String outfile, String labelfile) {
		
		PredictionEvaluator evaluator = new PredictionEvaluator();
		evaluator.parse_predictions(outfile);
		if (matrices_available) evaluator.parse_labels(labelfile);
		
		DecimalFormat fmt = new DecimalFormat();
		DecimalFormatSymbols symbs = fmt.getDecimalFormatSymbols();
        symbs.setDecimalSeparator('.');
        fmt.setDecimalFormatSymbols(symbs);
        fmt.setMaximumFractionDigits(2);
        fmt.setMinimumFractionDigits(2);
		
		double mean_num_matches = meanNumBestMatches(evaluator.best_matches);
		System.out.println("  Mean number of Matches:    " + fmt.format(mean_num_matches));
		if (matrices_available) {
			evaluator.eval_predictions();
		}
		else {
			int pwm_counter = 0;
			for (int i=0; i<evaluator.predicted_matrices.size(); i++) {
				
				if (! evaluator.predicted_matrices.get(i).isEmpty()) pwm_counter++;
			}
			double pred_rate = ((double) pwm_counter) / evaluator.predicted_matrices.size();
			System.out.println("  Prediction rate:           " + fmt.format(pred_rate) + "\n");
		}
	}
	
	public double meanNumBestMatches(ArrayList<ArrayList<LabeledTF>> best_matches) {
	
		int num_matches = 0;
		for (int i=0; i<best_matches.size(); i++) {
			num_matches += best_matches.get(i).size();
		}
		return (double) num_matches / best_matches.size();
	}	
	
	
	public void testSABINE(String basedir, double threshold, int max_num_matches, double outlier_threshold) {
		
		/*
		 *  draw samples
		 */
		
		drawSamples();
		writeOutfiles(basedir);
		
		/*
		 *  call SABINE
		 */
		System.out.println("  Predicting PWMs for Sample 1 (factors of class " + superclass + ") using model for class " + superclass);
		SABINE_Caller.main(new String[] {basedir + "sample_pos.input", 
	              						"-s", "" + threshold, 
	              						"-m", "" + max_num_matches,
	              						"-o", "" + outlier_threshold,  
	              						"-v", "h" });
		
		System.out.println("\n  Predicting PWMs for Sample 2 (factors of other classes) using model for class " + superclass);
		SABINE_Caller.main(new String[] {basedir + "sample_neg.input", 
										"-s", "" + threshold, 
										"-m", "" + max_num_matches,
										"-o", "" + outlier_threshold,  
										"-v", "h" });	
		System.out.println();
	
		/*
		 *  read and evaluate results 
		 *  (compute prediction rate, mean MoSta score, mean number of best matches)
		 */
		
		System.out.println("  Sample 1:");
		readOutfile(basedir + "sample_pos.output", basedir + "sample_pos.labels");

		System.out.println("  Sample 2:");
		readOutfile(basedir + "sample_neg.output", basedir + "sample_neg.labels");
	}
	
	

	public static void main(String[] args) {
		
		/*
		 *  default options 
		 */
		String inputfile = null;
		String labelfile = null;
		String basedir = "internal/";
		int superclass = 1;
		int samplesize = 5;
		double threshold = 0.95;
		int max_num_matches = 5;
		double outlier_threshold = 0.5;
		
		if (args.length == 1 && (args[0].equals("-help") || args[0].equals("--help"))) {
			usage();
		}
		
		for(int i=0; i<args.length-1; i+=2) {
			
			if(args[i].equals("-i")) { inputfile	 	 = args[i+1]; 						continue; }
			if(args[i].equals("-l")) { labelfile	     = args[i+1]; 						continue; }
			if(args[i].equals("-t")) { superclass 		 = Integer.parseInt(args[i+1]);		continue; }
			if(args[i].equals("-n")) { samplesize 		 = Integer.parseInt(args[i+1]);		continue; }
			if(args[i].equals("-s")) { threshold  	     = Double.parseDouble(args[i+1]); 	continue; }
			if(args[i].equals("-m")) { max_num_matches	 = Integer.parseInt(args[i+1]);		continue; }
			if(args[i].equals("-o")) { outlier_threshold = Double.parseDouble(args[i+1]);	continue; }
			
			if( !args[i].equals("-i") && !args[i].equals("-l") && !args[i].equals("-t")
					&& !args[i].equals("-n") && !args[i].equals("-s") 
					&& !args[i].equals("-m") && !args[i].equals("-o")) {	
				
				System.out.println("\n  Invalid argument: " + args[i]);
				usage();
			}
		}
		
		SABINE_Tester tester;
		if (inputfile != null && labelfile != null) {
			tester = new SABINE_Tester(inputfile, labelfile);
		}
		else if (inputfile != null) {
			tester = new SABINE_Tester(inputfile);
		}
		else {
			tester = new SABINE_Tester();
		}
		tester.superclass = superclass;
		tester.samplesize = samplesize;
		
		System.out.println();
		System.out.println("  ---------------------------------------------------------------------");
		System.out.println("  Statistics Tool for SABINE - StandAlone BINding specificity Estimator");
		System.out.println("  ---------------------------------------------------------------------");
		System.out.println("\n");
		System.out.println("  Superclass:  " + tester.superclass);
		System.out.println("  Sample size: " + tester.samplesize + "\n");
		System.out.println("  SABINE parameters:");
		System.out.println("    Min. best match similarity  :  " + threshold);
		System.out.println("    Max. number best matches    :  " + max_num_matches);
		System.out.println("    Max. outlier-FBP deviation  :  " + outlier_threshold + "\n");
		
		tester.testSABINE(basedir, threshold, max_num_matches, outlier_threshold);
	}
	
	
	
	private static void usage() {
		
		System.out.println();
		System.out.println("  ---------------------------------------------------------------------");
		System.out.println("  Statistics Tool for SABINE - StandAlone BINding specificity Estimator");
		System.out.println("  ---------------------------------------------------------------------");
		System.out.println("\n");
		System.out.println("  Usage   : stattool [OPTIONS]\n");
		System.out.println("  OPTIONS : -i <input_file>           (e.g. testdir/hs_testset.input)");
		System.out.println("            -l <label_file>           (e.g. testdir/hs_testset.labels)");
		System.out.println("            -t <transfac_superclass>  (e.g. 0, 1, 2, 3, 4 )                   default = 1");
		System.out.println("            -n <sample_size>          (number of drawn factors per sample)    default = 5");
		System.out.println("            -s <similarity_threshold> (min. FBP-similarity of best matches)   default = 0.95");
		System.out.println("            -m <max_num_best_matches> (max. number of best matches)           default = 5");		
		System.out.println("            -o <outlier_filter_param> (max. deviation of a single best match) default = 0.5\n\n");
		
		System.exit(0);
		
	}
}

