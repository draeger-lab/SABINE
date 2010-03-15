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
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;

import optimization.MoStaOptimizer;

public class PWMSimilarityCalculator {
	
	ArrayList<String[]> pwms = new ArrayList<String[]>();
	ArrayList<ArrayList<String>> stamp_pwms = new ArrayList<ArrayList<String>>();
	
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
		sim_calc.readPWMs("trainingsets_procaryotes/training_factors_procaryotes.txt");
		sim_calc.computePairwiseMostaScores("pwm_similarity/");
	}

}

