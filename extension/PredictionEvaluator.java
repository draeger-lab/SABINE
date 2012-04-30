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

package extension;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.ArrayList;
import java.util.StringTokenizer;

import libsvm.LabeledTF;

import optimization.MoStaOptimizer;
import plot.PlotGenerator;

public class PredictionEvaluator {

	public ArrayList<ArrayList<String>> predicted_matrices = new ArrayList<ArrayList<String>>(); 
	public ArrayList<ArrayList<String>> annotated_matrices = new ArrayList<ArrayList<String>>();
	
	public ArrayList<String> tf_names = new ArrayList<String>();
	public ArrayList<ArrayList<LabeledTF>> best_matches = new ArrayList<ArrayList<LabeledTF>>();

	double[] mosta_scores;
	boolean plot_hist = false;
	
	public void parse_predictions(String infile) {
		
		String line;
		ArrayList<LabeledTF> curr_best_matches;
		ArrayList<String> curr_matrix;
		StringTokenizer strtok;
		
		try {
			BufferedReader br = new BufferedReader(new FileReader(new File(infile)));
 			

			while ((line = br.readLine()) != null) {

				// parse factor name
				if (line.startsWith("NA  ")) {
					tf_names.add(line.substring(4).trim());
				}
				
				// parse best matches
				if (line.startsWith("BM  ")) {
					curr_best_matches = new ArrayList<LabeledTF>();
					
					if (! line.startsWith("BM  none")) {
					
						while (line.startsWith("BM  ")) {
								
							strtok = new StringTokenizer(line.substring(4));
							curr_best_matches.add(new LabeledTF(strtok.nextToken().trim(), 
																Double.parseDouble(strtok.nextToken().trim())));
							line = br.readLine();
						}
					}
					best_matches.add(curr_best_matches);
				}
				
				// parse predicted matrix
				if (line.startsWith("MA  ")) {
					
					curr_matrix = new ArrayList<String>();
					
					if (! line.substring(3).trim().equals("none")) {
										
						while (line.startsWith("MA  ")) {
							curr_matrix.add(line.substring(4).trim());
							line = br.readLine();
						}
					}
					predicted_matrices.add(curr_matrix);
				}
			}
		}
		catch(IOException ioe) {
			System.out.println(ioe.getMessage());
			System.out.println("IOException occurred while parsing predictions.");
		}
	}
	
	
	public void parse_labels(String infile) {
		
		String line;
		annotated_matrices = new ArrayList<ArrayList<String>>();
		
		try {
			BufferedReader br = new BufferedReader(new FileReader(new File(infile)));
			
			/*
			 *  read matrices from label file
			 */
			
			while ((line = br.readLine()) != null) {
				
				if (line.startsWith("MA  ")) {
					ArrayList<String> curr_matrix = new ArrayList<String>();
					
					while (line.startsWith("MA  ")) {
						curr_matrix.add(line.substring(4).trim());
						line = br.readLine();
					}
					annotated_matrices.add(curr_matrix);
				}
			}
		}
		catch(IOException ioe) {
			System.out.println(ioe.getMessage());
			System.out.println("IOException occurred while parsing labels.");
		}
	}
	
	public void eval_predictions() {
		
		if (annotated_matrices.isEmpty() || predicted_matrices.isEmpty()) {
			System.out.println("Fatal Error. Unable to evaluate predictions. Global variable \"annotated_matrices\" or \"predicted_matrices\" was not initialized. Aborting.");
			System.exit(0);
		}
		
		MoStaOptimizer optimizer = new MoStaOptimizer();
		double mosta_score_xx, mosta_score_yy, mosta_score_xy;
		mosta_scores = new double[annotated_matrices.size()];
		double sum_score = 0;
		int pwm_counter = 0;
		
		for (int i=0; i<annotated_matrices.size(); i++) {
			
			if (! predicted_matrices.get(i).isEmpty()) {
				mosta_score_xy = optimizer.compareMatrices(predicted_matrices.get(i), annotated_matrices.get(i), null);
				mosta_score_xx = optimizer.compareMatrices(predicted_matrices.get(i), predicted_matrices.get(i), null);
				mosta_score_yy = optimizer.compareMatrices(annotated_matrices.get(i), annotated_matrices.get(i), null);
				
				mosta_scores[i] = mosta_score_xy / Math.sqrt(mosta_score_xx * mosta_score_yy);
				//System.out.println("Current score: " + mosta_scores[i]);
				
				if (mosta_scores[i] == Double.POSITIVE_INFINITY) {
					
					System.out.println(mosta_score_xy + " / (sqrt(" + mosta_score_xx + " * " + mosta_score_yy + ")");
					System.exit(0);
				}
				
				pwm_counter++; 
				sum_score += mosta_scores[i];
			}
		}
		
		double mean_score = sum_score / pwm_counter;
		double pred_rate = ((double) pwm_counter) / annotated_matrices.size();
		
		/*
		 *  print results
		 */
		
		DecimalFormat fmt = new DecimalFormat();
		DecimalFormatSymbols symbs = fmt.getDecimalFormatSymbols();
        symbs.setDecimalSeparator('.');
        fmt.setDecimalFormatSymbols(symbs);
        fmt.setMaximumFractionDigits(2);
        fmt.setMinimumFractionDigits(2);
		
        if (new Double(mean_score).isInfinite() || new Double(mean_score).isNaN())
        	System.out.println("  Mean MoSta score:          -");
        else
        	System.out.println("  Mean MoSta score:          " + fmt.format(mean_score));
	    System.out.println("  Prediction rate:           " + fmt.format(pred_rate) + "\n");
		
		if (plot_hist) {
			PlotGenerator hist_plotter = new PlotGenerator();
			hist_plotter.plotHist(mosta_scores, 30, "Distribution of MoSta scores");
		}
	}
	

	/*
	 *  calculate mean MoSta scores for varying best match thresholds
	 */
	
	public void varyBestMatchThreshold(double bmt, double init, double incr, double term, int max_num_matches, String inputfile) {
		
		ArrayList<ArrayList<LabeledTF>> all_best_matches;
		ArrayList<LabeledTF> curr_best_matches;
		
		ArrayList<ArrayList<String>> best_match_matrices;
		ArrayList<ArrayList<String>> all_merged_matrices = new ArrayList<ArrayList<String>>();
		ArrayList<String> merged_matrix = new ArrayList<String>();
		
		MatrixFileParser matrixfile_parser = new MatrixFileParser();
		//String matrixfile = "trainingsets/trainingset_all_classes.fbps";
		String matrixfile = "trainingsets_procaryotes/trainingset_class0.fbps";
		matrixfile_parser.readMatrices(matrixfile);
		
		PFMFormatConverter pfm_converter = new PFMFormatConverter();
		MoStaOptimizer optimizer = new MoStaOptimizer();
		
		double[] all_mosta_scores, filtered_mosta_scores, mean_mosta_scores, prediction_rates;
		int curr_matrix_counter;
		boolean[] valid_mosta_scores;
		double curr_sum_score, curr_mean_score;
		double mosta_score_xy, mosta_score_xx, mosta_score_yy; 
		
		DecimalFormat fmt = new DecimalFormat();
		DecimalFormatSymbols symbs = fmt.getDecimalFormatSymbols();
        symbs.setDecimalSeparator('.');
        fmt.setDecimalFormatSymbols(symbs);
        fmt.setMaximumFractionDigits(2);
        fmt.setMinimumFractionDigits(2);
		
		// define grid
        ArrayList<Double> grid;
        double[] bmt_grid;
        double curr_val;
        
        if (! new Double(bmt).isNaN()) 
        	bmt_grid = new double[] { bmt };
        else {
        	grid = new ArrayList<Double>();
        	curr_val = init;
        	while (curr_val <= term) {
        		grid.add(curr_val);
        		curr_val += incr;
        	}
        	bmt_grid = new double[grid.size()];
        	for (int i=0; i<grid.size(); i++) bmt_grid[i] = grid.get(i);
        		
        	if (bmt_grid.length < 20) {
        		System.out.print("  Best Match Thresholds:  [ ");
        		for (int i=0; i<bmt_grid.length; i++) System.out.print(fmt.format(bmt_grid[i]) + " ");
        		System.out.println("]\n");
        	}
        }
		
		/*
		 * parse best matches for current threshold
		 */ 
		mean_mosta_scores = new double[bmt_grid.length];
		prediction_rates = new double[bmt_grid.length];
		
		for (int i=0; i<bmt_grid.length; i++) {
			
			System.out.println("  Best Match Threshold: " + fmt.format(bmt_grid[i]) + "\n");
			System.out.println("    Parsing best matches.");
			
			all_best_matches = new ArrayList<ArrayList<LabeledTF>>();
			for (int j=0; j<best_matches.size(); j++) {
				
				curr_best_matches = new ArrayList<LabeledTF>();
				for (int k=0; k<best_matches.get(j).size(); k++) {
					
					if (best_matches.get(j).get(k).getLabel() > bmt_grid[i]) {
						
						curr_best_matches.add(best_matches.get(j).get(k));
					}
				}
				all_best_matches.add(curr_best_matches);
			}
			
			// obtain matrices for current best matches
			System.out.println("    Merging best match matrices.");
			all_merged_matrices = new ArrayList<ArrayList<String>>();
			
			for (int j=0; j<all_best_matches.size(); j++) {
								
				best_match_matrices = new ArrayList<ArrayList<String>>();
				for (int k=0; k<Math.min(all_best_matches.get(j).size(), max_num_matches); k++) {
					best_match_matrices.add(matrixfile_parser.obtainMatrix(all_best_matches.get(j).get(k).getName()));
				}
				
				// filter outliers and merge matrices
				merged_matrix = new ArrayList<String>();
				if (best_match_matrices.size() > 1) {
					merged_matrix = pfm_converter.mergePFMs(best_match_matrices);
				}
				else if (best_match_matrices.size() == 1) {
					merged_matrix = best_match_matrices.get(0);
				}
				all_merged_matrices.add(merged_matrix);
			}
			
			curr_matrix_counter = 0;
			filtered_mosta_scores = new double[curr_matrix_counter];
			
			if (! annotated_matrices.isEmpty()) {
			
				// compare to annotated matrix and compute MoSta score
				System.out.println("    Calculating MoSta scores.\n");
				all_mosta_scores = new double[annotated_matrices.size()];
				valid_mosta_scores = new boolean[annotated_matrices.size()];
				curr_sum_score = 0;

				for (int j=0; j<annotated_matrices.size(); j++) {
					
					//System.out.println("Number of merged matrices: " + all_merged_matrices.get(j).get(0));
					
					if (! all_merged_matrices.get(j).isEmpty()) {
						
						mosta_score_xy = optimizer.compareMatrices(all_merged_matrices.get(j), annotated_matrices.get(j), null);
						mosta_score_xx = optimizer.compareMatrices(all_merged_matrices.get(j), all_merged_matrices.get(j), null);
						mosta_score_yy = optimizer.compareMatrices(annotated_matrices.get(j), annotated_matrices.get(j), null);
						
						all_mosta_scores[j] = mosta_score_xy / Math.sqrt(mosta_score_xx * mosta_score_yy);
						valid_mosta_scores[j] = true;
						
						curr_sum_score += all_mosta_scores[j];
						curr_matrix_counter++;
					}
				}
				int cnt = 0;
				for (int j=0; j<valid_mosta_scores.length; j++) {
					if (valid_mosta_scores[j]) filtered_mosta_scores[cnt++] = all_mosta_scores[j];
				}
				
				curr_mean_score = curr_sum_score / curr_matrix_counter; 
				mean_mosta_scores[i] = curr_mean_score;
				prediction_rates[i] = ((double) curr_matrix_counter) / annotated_matrices.size();
				
				System.out.println("  Mean MoSta score: " + fmt.format(curr_mean_score));
			}
			else {
				for (int j=0; j<all_merged_matrices.size(); j++) {
					if (! all_merged_matrices.get(j).isEmpty()) 
						curr_matrix_counter++;
				}
				prediction_rates[i] = ((double) curr_matrix_counter) / all_merged_matrices.size();
			}
			System.out.println("  Prediction rate:  " + fmt.format(prediction_rates[i]));
			
			
			if (i < bmt_grid.length-1) {
				System.out.println(" ______________________________" + "\n");
			}
			else {
				String outfile = inputfile + "_t=" + fmt.format(bmt_grid[i]) + "_m=" + max_num_matches;
				rewriteOutfile(outfile, all_best_matches, all_merged_matrices);
			}
			
			if (plot_hist && !annotated_matrices.isEmpty()) {
				PlotGenerator hist_plotter = new PlotGenerator();
				hist_plotter.plotHist(filtered_mosta_scores, 30, "Distribution of MoSta scores for Best Match Threshold = " + fmt.format(bmt_grid[i]));
			}
		}

		 
		if (plot_hist && bmt_grid.length >= 3 && !annotated_matrices.isEmpty()) {
			PlotGenerator line_plotter = new PlotGenerator();
			line_plotter.plotLine(bmt_grid, mean_mosta_scores, prediction_rates);
		}
		System.out.println();
	}	
		
	
	
	/*
	 *  write scores to textfile
	 */
	
	public void writeMoStaScores(String outfile) {
		
		try {
			BufferedWriter bw = new BufferedWriter(new FileWriter(new File(outfile)));
			
			for (int i=0; i<annotated_matrices.size(); i++) {	
				if (! predicted_matrices.get(i).isEmpty()) {
					bw.write(tf_names.get(i) + ": " + mosta_scores[i] + "\n");
				}
			}
			bw.flush();
			bw.close();
		}
		catch(IOException ioe) {
			System.out.println(ioe.getMessage());
			System.out.println("IOException occurred writing MoSta scores to text file.");
		}
	}
	
	
	public void writeAnnotatedMatrices(String outfile) {
		
		try {
			BufferedWriter bw = new BufferedWriter(new FileWriter(new File(outfile)));
			
			for (int i=0; i<annotated_matrices.size(); i++) {
				bw.write("NA  " + tf_names.get(i) + "\nXX\n");
				for (int j=0; j<annotated_matrices.get(i).size(); j++) {
					
					bw.write("MA  " + annotated_matrices.get(i).get(j) + "\n");
				}
				bw.write("XX\n");
			}
			bw.flush();
			bw.close();
		
		}
		catch(IOException ioe) {
			System.out.println(ioe.getMessage());
			System.out.println("IOException occurred writing MoSta scores to text file.");
		}
	}
	
	
	public void rewriteOutfile(String outfile, ArrayList<ArrayList<LabeledTF>> all_best_matches, ArrayList<ArrayList<String>> all_merged_matrices) {

		
		try {
			
			/*
			 *  write new oufile
			 */
			
			BufferedWriter bw = new BufferedWriter(new FileWriter(new File(outfile)));
			
			for (int i=0; i<tf_names.size(); i++) {
				
				// write name
				if (! all_best_matches.get(i).isEmpty()) {
					bw.write("NA  " + tf_names.get(i) + "\nXX\n");
					
					// write best matches
					for (int j=0; j<all_best_matches.get(i).size(); j++) {
						bw.write("BM  " + all_best_matches.get(i).get(j).getName() + "  " + 
										+ all_best_matches.get(i).get(j).getLabel() + "\n");
					}
					bw.write("XX\n");
					
					// write matrix
					for (int j=0; j<all_merged_matrices.get(i).size(); j++) {	
						bw.write("MA  " + all_merged_matrices.get(i).get(j) + "\n");
					}
					bw.write("XX\n" +
							 "//\n" +
							 "XX\n");
				}
			}
			bw.flush();
			bw.close();
		}
		catch(IOException ioe) {
			System.out.println(ioe.getMessage());
			System.out.println("IOException occurred writing MoSta scores to text file.");
		}
	}
	
	
	public static void main(String[] args) {
		
		
		PredictionEvaluator evaluator = new PredictionEvaluator();
		
		/*
		 * parse options
		 */
		
		if (args.length < 1) evaluator.usage();
		
		

		String inputfile, threshold, plotoption;
		inputfile = args[0];
		
		String labelfile = null;
		int start_idx = 1;
		if (! (args[1].equals("-t") || args[1].equals("-m") || args[1].equals("-p"))) {
			labelfile = args[1];
			start_idx = 2;
		}
		
		labelfile = args[1];
		threshold = null;
		int max_num_matches = 5; 	// default
		plotoption = "n";  			// default
		
		
		for(int i=start_idx; i<args.length-1; i+=2) {
			
			if(args[i].equals("-t")) { threshold	 	= args[i+1]; 					continue; }
			if(args[i].equals("-m")) { max_num_matches  = Integer.parseInt(args[i+1]); 	continue; }
			if(args[i].equals("-p")) { plotoption  		= args[i+1]; 					continue; }
			
			if( !args[i].equals("-t") && !args[i].equals("-p") && !args[i].equals("-m")) {	
				System.out.println("\n  Invalid argument: " + args[i]);
				evaluator.usage();
			}
		}
		
		if (plotoption.toUpperCase().equals("Y") || plotoption.toUpperCase().equals("YES")) 
			evaluator.plot_hist = true;
		
		System.out.println();
		System.out.println("  ---------------------------------------------------------------------");
		System.out.println("  Evaluation Tool for SABINE - StandAlone BINding specificity Estimator");
		System.out.println("  ---------------------------------------------------------------------");
		System.out.println("\n");
		
		evaluator.parse_predictions(inputfile);
		
		if (start_idx == 2) {
			evaluator.parse_labels(labelfile);
		}
		double bmt, init, incr, term;
		bmt = init = incr = term = Double.NaN;
		
		if (threshold == null && start_idx == 2) {
			evaluator.eval_predictions();
			evaluator.writeMoStaScores(inputfile.replace(".output", ".scores"));
		}
		else {
			if (! threshold.contains(":"))
				bmt = Double.parseDouble(threshold);
			else {
				String[] split = threshold.split(":");
				init = Double.parseDouble(split[0]);
				incr = Double.parseDouble(split[1]);
				term = Double.parseDouble(split[2]);
			}
			evaluator.varyBestMatchThreshold(bmt, init, incr, term, max_num_matches, inputfile);
		}
		
	}
	
	
	private void usage() {
		
		System.out.println();
		System.out.println("  ---------------------------------------------------------------------");
		System.out.println("  Evaluation Tool for SABINE - StandAlone BINding specificity Estimator");
		System.out.println("  ---------------------------------------------------------------------");
		System.out.println("\n");
		System.out.println("  Usage   : evaltool <input_file> <label_file> [OPTIONS]\n");
		System.out.println("  OPTIONS : -t <best_match_threshold>  (e.g. 0.95, 0.5:0.1:1) ");	
		System.out.println("            -m <max_num_best_matches>  (max. number of best matches)      default = 5");		
		System.out.println("            -p <plot_histogram>        (e.g. y (yes), n (no)) \n\n");
		
		System.exit(0);
		
	}

}

