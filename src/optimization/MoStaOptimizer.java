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
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.StringTokenizer;

import main.FBPPredictor;



/*
 * 
 * calculates the pairwise MOSTA-scores of a set of TFs
 * 
 */
public class MoStaOptimizer implements Optimizer {
	
	
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
			System.out.println("IOException occurred while comparing FBPs with MOSTA.");
		}
		
		
	}
	
	
	/*
	 *  gets two fbp matrices in STAMP format and compares them with MoSta
	 */
	public double compareMatrices(ArrayList<String> matrix1, ArrayList<String> matrix2, String[] params) {
		
		
		StringTokenizer strtok = null;
		String a,c,g,t;
		
		int lineCount = 0;
		String lineID = null;
		
		ArrayList<String> matrix1forward = new ArrayList<String>();
		ArrayList<String> matrix2forward = new ArrayList<String>();
		
	// reduce matrices to essential info	
		
		for(int i=0; i<matrix1.size(); i++) {
			
			strtok = new StringTokenizer(matrix1.get(i));
			 
			strtok.nextToken();
			 
			a = strtok.nextToken();
			c = strtok.nextToken();
			g = strtok.nextToken();
			t = strtok.nextToken();
			
			lineCount = 100 + i + 1;
			lineID = ("" + lineCount).substring(1);
			
			matrix1forward.add(lineID + "\t" + a + "\t" + c + "\t" + g + "\t" + t);
			
		}
		
		for(int i=0; i<matrix2.size(); i++) {
			
			strtok = new StringTokenizer(matrix2.get(i));
			 
			strtok.nextToken();
			 
			a = strtok.nextToken();
			c = strtok.nextToken();
			g = strtok.nextToken();
			t = strtok.nextToken();
			
			lineCount = 100 + i + 1;
			lineID = ("" + lineCount).substring(1);
			
			matrix2forward.add(lineID + "\t" + a + "\t" + c + "\t" + g + "\t" + t);
			
		}
		
	// write matrices to input file for external program
		
		BufferedWriter bw = null;
		
		
		try {
			
			bw = new BufferedWriter(new FileWriter(new File(basedir + "mosta/matrices.raw")));
			
			bw.write("ID Matrix1\n");
			bw.write("P0\tA\tC\tG\tT\n");
			
			for(int i=0; i<matrix1forward.size(); i++) {
				bw.write(matrix1forward.get(i) + "\n");
			}
			
			bw.write("XX\n");
			bw.write("//\n");
			
			bw.write("ID Matrix2\n");
			bw.write("P0\tA\tC\tG\tT\n");
			
			for(int i=0; i<matrix2forward.size(); i++) {
				bw.write(matrix2forward.get(i) + "\n");
			}
			
			bw.write("XX\n");
			bw.write("//\n");
			
			bw.flush();
			bw.close();
		
		}
		catch(IOException ioe) {
			System.out.println("IOException occured while comparing matrices.");
			System.out.println(ioe.getMessage());
		}
		
		return compareMatricesMoSta(params);
		
	}
	
	
	
	
	public double compareMatricesMoSta(String[] params) {
		
		double minScore = Double.POSITIVE_INFINITY;
		
		if(params == null) {
			
			params = new String[] {".4","balanced",".1","0"};
			
		}
		
	// execute method	
		
		try {
			
			
			String cmdString = "./" + FBPPredictor.MoStaDir + "code/sstat " +
					
					params[0] + " "  +
					basedir + "mosta/matrices.raw "  +
					params[1]  + " " +
					params[2];

			Process proc = Runtime.getRuntime().exec(cmdString);
			
			proc.getOutputStream().flush();
			proc.getOutputStream().close();
			
			proc.waitFor();
			
			String aline = null;
			String atoken = null;
			StringTokenizer atok = null;
			
			BufferedReader input = new BufferedReader(new InputStreamReader(proc.getInputStream()));
			
			
			while ((aline = input.readLine()) != null) {
		        
		        
		        if(aline.startsWith("Matrix1")) {
		        	
		        	atok = new StringTokenizer(aline);
		        	
		        	atok.nextToken(); 				// Matrix1
		        	atok.nextToken(); 				// Matrix2
		        	
		        	atoken = atok.nextToken();  	// SMAX
		        	
		        	if(params[3].equals("1")) {
		        		atoken = atok.nextToken();  // SSUM
		        	}
		        	
		        	if (atoken.equals("-inf")) {
		        		
		        		minScore =  -1.0;
		        	
		        	}
		        	
		        	else {
		        		
		        		minScore = Double.parseDouble(atoken);	
		        		
		        	}	
		        	
		        }
				
		    }
			
			input.close();
			
			proc.destroy();
			
		}
		catch(IOException e) {
			System.out.println(e.getMessage());
			System.out.println("IOException while executing Mosta!");
		} 
		catch (InterruptedException e) {
			System.out.println(e.getMessage());
			System.out.println("InterruptedException while executing Mosta!");
		}
		
		return minScore;
		
	}
	
	
	
	public static void main(String[] args) {
		
		MoStaOptimizer optimizer = new MoStaOptimizer();
		
		
		String[] params = new String[] {
				
				".4",			// G/C content of background genome						
				"balanced",		// threshold method						(typeI, typeII, balanced, typeIext)
				".1",			// threshold param						(irrelevant if method == balanced)
				"0"				// output score							(0 == S_max , 1 == S_sum)
				
		};
		
		
		optimizer.calculateFBPScores(
				
				"FBPs_rat.out", 
				
				params, 
				
				"MoSta_scores_rat.out");
		
		
	}
	
	
	
	
}

