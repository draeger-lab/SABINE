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

import gui.SequenceLogo;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;

import optimization.MoStaOptimizer;

public class PWMSimilarityCalculator {
	
	ArrayList<String[]> pwms = new ArrayList<String[]>();
	ArrayList<ArrayList<String>> stamp_pwms = new ArrayList<ArrayList<String>>();
	String basedir = "tmp/my_basedir/";
	
	public void readPWMs(String infile) {
		
		try {
			 
			 BufferedReader bw = new BufferedReader(new FileReader(new File(infile)));
			 String[] curr_matrix;
			 String line;
			 
			 while ((line = bw.readLine()) != null) {

				 if (line.startsWith("MA")) {
					curr_matrix = new String[4];
					curr_matrix[0] = line.substring(3).trim();
					curr_matrix[1] = bw.readLine().substring(3).trim();
					curr_matrix[2] = bw.readLine().substring(3).trim();
					curr_matrix[3] = bw.readLine().substring(3).trim();
					pwms.add(curr_matrix);
				 }
			 }
			 bw.close();
			 System.out.println(pwms.size() + " matrices parsed.");
		}
		catch(IOException ioe) {
			System.out.println(ioe.getMessage());
			System.out.println("IOException occurred while writing species mapping.");
		}
	}
	
	
	public void computePairwiseMostaScores(String basedir) {
		
		// convert matrices to STAMP format
		PFMFormatConverter pwm_converter = new PFMFormatConverter();
		stamp_pwms = pwm_converter.convertAllTransfacToSTAMP(pwms);
		
		// calculate all pairwise scores
		
		MoStaOptimizer optimizer = new MoStaOptimizer();
		optimizer.basedir = basedir;
		
		double[][] distances = new double[stamp_pwms.size()][stamp_pwms.size()];
		double sum = 0;
		
		for(int i=0; i<stamp_pwms.size(); i++) {
			
			for(int j=i; j<stamp_pwms.size(); j++) {
				
				distances[i][j] = optimizer.compareMatrices(stamp_pwms.get(i), stamp_pwms.get(j), null);
				
			}
			
		}
			
		
		// normalize pairwise scores	
		
		for(int i=0; i<stamp_pwms.size()-1; i++) {
			
			for(int j=i+1; j<stamp_pwms.size(); j++) {
				
				distances[i][j] = distances[i][j] / Math.sqrt(distances[i][i] * distances[j][j]);
				distances[j][i] = distances[i][j];
				
				sum += distances[i][j];
			}
		
		}
		
		for(int i=0; i<stamp_pwms.size(); i++) {
			
			distances[i][i] = 1.0;
			
		}
		
		int num_pairs = (stamp_pwms.size() * stamp_pwms.size()-1)/2;
		
		// compute mean pairwise score
		System.out.println("Mean MoSta score: " + sum/num_pairs);
	}
	
	public static void main(String[] args) {
		
		PWMSimilarityCalculator sim_calc = new PWMSimilarityCalculator();
		
		ArrayList<String> best_matches = new ArrayList<String>();
		ArrayList<String> merged_matrix = new ArrayList<String>();
		ArrayList<String> annotated_matrix = new ArrayList<String>();
    	
		best_matches.add("T01006");
    	best_matches.add("T01005");
    	best_matches.add("T01009");
    	
    	MatrixFileParser matrix_parser = new MatrixFileParser();
    	matrix_parser.readMatrices("/rahome/eichner/projects/sabine/src/SABINE_1.0/trainingsets_public/trainingset_all_classes.fbps");
    	
    	ArrayList<ArrayList<String>> best_match_matrices = new ArrayList<ArrayList<String>>();
    	for (int i=0; i<best_matches.size(); i++) {
    		best_match_matrices.add(matrix_parser.obtainMatrix(best_matches.get(i)));
    	}
    	
    	PFMFormatConverter pfm_converter = new PFMFormatConverter();
    	pfm_converter.basedir = sim_calc.basedir;
    	
		merged_matrix = pfm_converter.mergePFMs(best_match_matrices);
		
    	for (int i=0; i<best_match_matrices.size(); i++) {
    		for (int j=0; j<best_match_matrices.get(i).size(); j++) {
    		
    			System.out.println(best_match_matrices.get(i).get(j));
    		}
    		System.out.println("\n");
    	}
    	System.out.println("_______________\n\n");
		
    	for (int i=0; i<merged_matrix.size(); i++) {
    		System.out.println(merged_matrix.get(i));
    	}
		
		SequenceLogo.printLogo("predicted PFM", merged_matrix);
		
		annotated_matrix = matrix_parser.obtainMatrix("T00505");
		
		MoStaOptimizer optimizer = new MoStaOptimizer();
		optimizer.basedir = sim_calc.basedir;
		double dist_xy = optimizer.compareMatrices(merged_matrix, annotated_matrix, null);
		double dist_xx = optimizer.compareMatrices(merged_matrix, merged_matrix, null);
		double dist_yy = optimizer.compareMatrices(annotated_matrix, annotated_matrix, null);
		
		double dist = dist_xy / Math.sqrt(dist_xx * dist_yy);
		
		System.out.println("Distance: " + dist);
	}

}

