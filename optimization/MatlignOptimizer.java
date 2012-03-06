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

package optimization;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.StringTokenizer;



/*
 * 
 * calculates the pairwise MATLIGN-scores of a set of TFs
 * 
 */
public class MatlignOptimizer implements Optimizer {
	
	
	boolean silent = false;
	public String basedir = "";
	
	public void calculateFBPScores(String fbpfile, String[] params, String outfile) {
		
		BufferedReader br = null;
		BufferedWriter bw = null;
		
		try {
			
			br = new BufferedReader(new FileReader(new File(fbpfile)));
			bw = new BufferedWriter(new FileWriter(new File(outfile)));
			
			String line = null;
			
			ArrayList<ArrayList<String>> matrices = new ArrayList<ArrayList<String>>();
			
			ArrayList<String> names = new ArrayList<String>();
			
			
		// parse all FBPs and store them in list "matrices"
			
			while((line = br.readLine()) != null) {
				
				// System.out.println("Line: " + line);
				
				if(!line.startsWith("DE\t")) {
					System.out.println("Parse Error. DE expected. Aborting.");
					System.exit(0);
				}
				
				names.add(line.substring(3));
				
				if(!silent) System.out.println("Parsing FBP of " + line.substring(3) + ".");
					
				ArrayList<String> matrix = new ArrayList<String>();
				
				
			// parse a single FBP	
				
				while(!(line = br.readLine()).startsWith("XX")) {
					
					matrix.add(line);
					
				}
				
				matrices.add(matrix);
					
			}
			
			if(!silent) System.out.println("\ndone.\n");
			
		
			double[] self_scores = new double[matrices.size()];
			
		
		// compute all self scores
		
			for(int i=0; i<matrices.size(); i++) {
				
				if(!silent) System.out.println("Calculating self-score of " + names.get(i) + ".");
				
				self_scores[i] = compareMatrices(matrices.get(i), matrices.get(i), params); 
				
			}
			
			
		// compare all FBPs vs all FBPs
			
			for(int i=0; i<matrices.size()-1; i++) {
				
				if(!silent) System.out.println("Comparing " + names.get(i));
				
				for(int j=i+1; j<matrices.size(); j++) {
					
					if(!silent) System.out.println("  with " + names.get(j) + ".");
					
					bw.write(
							
							 compareMatrices(matrices.get(i), matrices.get(j), params)
							 
							 /
							 
							 Math.sqrt(self_scores[i] * self_scores[j]) 
							 
							 + "\n"
							 
					);
				
					bw.flush();	
				}
				
			}
			
			
			
			br.close();
			bw.flush();
			bw.close();
		}
		catch(IOException ioe) {
			System.out.println(ioe.getMessage());
			System.out.println("IOException occurred while comparing FBPs with MATLIGN.");
		}
		
		
	}
	
	
	/*
	 *  gets two fbp matrices in STAMP format and compares them with Matlign
	 */
	
	public double compareMatrices(ArrayList<String> matrix1, ArrayList<String> matrix2, String[] params) {
		
		MatrixConverter converter = new MatrixConverter();
		
		StringTokenizer strtok = null;
		String a,c,g,t;
		
		ArrayList<String> matrix1forward = new ArrayList<String>();
		ArrayList<String> matrix2forward = new ArrayList<String>();
		
		
	// generate reverse complements of matrices	
		
		ArrayList<String> matrix1reversed = converter.generateReverseComplement(matrix1);	
		ArrayList<String> matrix2reversed = converter.generateReverseComplement(matrix2);
		
		
	// reduce matrices to essential info	
		
		for(int i=0; i<matrix1.size(); i++) {
			
			strtok = new StringTokenizer(matrix1.get(i));
			 
			strtok.nextToken();
			 
			a = strtok.nextToken();
			c = strtok.nextToken();
			g = strtok.nextToken();
			t = strtok.nextToken();
			 
			matrix1forward.add(a + "\t" + c + "\t" + g + "\t" + t);
			
		}
		
		for(int i=0; i<matrix2.size(); i++) {
			
			strtok = new StringTokenizer(matrix2.get(i));
			 
			strtok.nextToken();
			 
			a = strtok.nextToken();
			c = strtok.nextToken();
			g = strtok.nextToken();
			t = strtok.nextToken();
			 
			matrix2forward.add(a + "\t" + c + "\t" + g + "\t" + t);
			
		}
		
		
		// run Matlign 4 times
		
		BufferedWriter bw = null;
		
		double maxScore = Double.NEGATIVE_INFINITY;
		
		double score;
		
		try {
			
			// write matrices to file	
				
			ArrayList<ArrayList<String>> matrices = new ArrayList<ArrayList<String>>();
			
			matrices.add(matrix1forward);
			matrices.add(matrix2forward);
			matrices.add(matrix1reversed);
			matrices.add(matrix2reversed);
			
			int[][] comparisons = new int[][] {{0,1},{0,3},{1,2},{2,3}};
			
			
			for(int k=0; k<comparisons.length; k++) {
				
				bw = new BufferedWriter(new FileWriter(new File(basedir + "matlign/matrices.raw")));
				
				bw.write(">Matrix1\n");
				
				for(int i=0; i<matrices.get(comparisons[k][0]).size(); i++) {
					bw.write(matrices.get(comparisons[k][0]).get(i) + "\n");
				}
				
				bw.write("\n");
				
				bw.write(">Matrix2\n");
				
				for(int i=0; i<matrices.get(comparisons[k][1]).size(); i++) {
					bw.write(matrices.get(comparisons[k][1]).get(i) + "\n");
				}
				
				bw.flush();
				bw.close();
				
				
				if((score = compareMatricesMatlign(params)) > maxScore)
					
					maxScore = score;
			}	
				
		}
		catch(IOException ioe) {
			System.out.println("IOException occured while comparing matrices.");
			System.out.println(ioe.getMessage());
		}
		
		
		return maxScore;
		
	}
	
	
	
	public double compareMatricesMatlign(String[] params) {
		
		double minScore = Double.POSITIVE_INFINITY;
		
		if(params == null) {
			
			params = new String[] {"-4","-4","5","-10","-1","5","235","1","0","1","0","0","0.5","0.5"};
			
			
		}
		
	// execute method	
		
		try {
			
			
			String cmdString = "./Matlign/matlign " +
					
					"-input=" + basedir + "matlign/matrices.raw" +
					
					" -transv=" 	+ params[0]  +
					" -transi="		+ params[1]  +
					" -match=" 		+ params[2]  +
					" -gopen=" 		+ params[3]  +
					" -gext=" 		+ params[4]  +
					" -term=" 		+ params[5]  +
					" -norm=" 		+ params[6]  +
					
					" -out=" + basedir + "matlign/result" 	+
					
					" -spacer=" 	+ params[7]  +
					" -mode=" 		+ params[8]  +
					" -zscore=" 	+ params[9]  +
					" -random=" 	+ params[10] +
					" -pseudo=" 	+ params[11] +
					" -freqat=" 	+ params[12] +
					" -freqcg=" 	+ params[13]
					;

			Process proc = Runtime.getRuntime().exec(cmdString);
			
			proc.getOutputStream().flush();
			proc.getOutputStream().close();
			
			proc.waitFor();
			proc.destroy();
			
		}
		catch(IOException e) {
			System.out.println(e.getMessage());
			System.out.println("IOException while executing Matlign!");
		} 
		catch (InterruptedException e) {
			System.out.println(e.getMessage());
			System.out.println("InterruptedException while executing Matlign!");
		}
		
		
	// extract result
		
		BufferedReader br = null;
		
		try {
				
				br = new BufferedReader(new FileReader(new File(basedir + "matlign/result.mtrx")));
				
				String line = br.readLine();
				
				StringTokenizer strtok = new StringTokenizer(line, "|");
				
				strtok.nextToken();
				minScore = Double.parseDouble(strtok.nextToken());
				
				br.close();

				if (basedir.isEmpty()) {
					
					new File("result.fdrs").delete();
					new File("result.pswm").delete();
					new File("result.tree").delete();
					new File("result.vara").delete();
					new File("result.freq").delete();
					new File("matrices.raw").delete();
					new File("result.mtrx").delete();
				}
				
				br.close();

		}
		
		catch(IOException ioe) {
			System.out.println("IOException occured while parsing result of matrix comparison.");
		}
		
		
		return minScore;
		
	}
	
	public static void main(String[] args) {
		
		MatlignOptimizer optimizer = new MatlignOptimizer();
		
		
		String[] params = new String[] {
				
				"-4",		// mismatch score for transversions				(e.g. A --> C)
				"-4",		// mismatch score for transitions				(e.g. A --> G)
				"5",		// match score
				"-10",		// gap open penalty
				"-1",		// gap extension penalty
				"5",		// maximum lengt of terminal gaps 				(determines minimum PFM-overlap)
				"235",		// combination of distance metrics				(Kendall, Spearman, Pearson, Euclidian, Viterbi)
				"1",		// enable/disable internal spacers
				"0",		// cluster motifs 
				"1",		// enable/disable z-score-transformed scores
				"0",		// perform N permutations and calculate FDR
				"0",		// add pseudocounts to base frequencies			(pseudo in [0,1])
				"0.5",		// A/T content
				"0.5"		// C/G content
				
		};
		
		
		optimizer.calculateFBPScores(
				
				"FBPs_rat.out", 
				
				params, 
				
				"Matlign_scores_rat.out");
		
		
	}
	
	
}

