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

package core;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.StringTokenizer;

public class SVMPairwiseFeatureCalculator {
	
	
	public void calculateSVMPairwiseScores(String tfname, String orientationfile, String testscoresfile, String trainscoresfile, String outfile) {
		
		BufferedReader br = null;
		BufferedWriter bw = null;
		
		
		ArrayList<String> testscores = new ArrayList<String>();
		ArrayList<String> tfnames = new ArrayList<String>();
		
		try {
			
			
		// read scores of the test tf	
			
			br = new BufferedReader(new FileReader(new File(testscoresfile)));
			
			String line = null;
			
			StringTokenizer strtok = null;
			
			
			
			while((line = br.readLine()) != null) {
				
				strtok = new StringTokenizer(line);
				
				strtok.nextToken();							// test tf
				
				strtok.nextToken(); 						// vs.
				
				tfnames.add(strtok.nextToken()); 			// training tf
				
				strtok.nextToken();							// :
				
				testscores.add(strtok.nextToken());			// score
				
			}
				

			br.close();
			
			
			double[][] trainscores = new double[testscores.size()][testscores.size()];
			
			

		// read scores of the training tfs	
				
			br = new BufferedReader(new FileReader(new File(trainscoresfile)));
				

				
			for(int i=0; i<trainscores.length; i++) {
					
					
				for(int j=0; j<trainscores.length; j++) {
					
					strtok = new StringTokenizer(br.readLine());
						
					strtok.nextToken();							// training tf i
					strtok.nextToken(); 						// vs.
					strtok.nextToken(); 						// training tf j
					strtok.nextToken();							// :
						
					trainscores[i][j] = Double.parseDouble(strtok.nextToken());
					
				}
					
			}
					

			br.close();
				
			
		// calculate dot products only for relevant pairs 
			
			br = new BufferedReader(new FileReader(new File(orientationfile)));
			bw = new BufferedWriter(new FileWriter(new File(outfile)));
			
			int index = 0;
			
			double dot_product = 0.0;
			
			while((line = br.readLine()) != null) {
			
					
				strtok = new StringTokenizer(line);
						
				strtok.nextToken();								// test tf
				strtok.nextToken(); 							// vs.
				
				index = tfnames.indexOf(strtok.nextToken());
				
				strtok.nextToken();								// :
					
				
				dot_product = 0.0;
				
				for(int i=0; i<tfnames.size(); i++) {
					
					dot_product += Double.parseDouble(testscores.get(i)) * trainscores[index][i];
					
				}
				
				bw.write(tfname + " vs. " + tfnames.get(index) + " : " + dot_product + "\n");
					
			}
					

			br.close();
				
			
			
			bw.flush();
			bw.close();
			
		}
		catch(IOException ioe) {
			System.out.println(ioe.getMessage());
			System.out.println("IOException occurred while precalculating SVMPairwise Helpfiles.");
		}
		
	}
	
public void calculateAllSVMPairwiseScores(String trainscoresfile, String outfile, ArrayList<String> relevant_pairs) {
		
		BufferedReader br = null;
		BufferedWriter bw = null;
		
		ArrayList<String> tfnames = new ArrayList<String>();
		
		try {
			
			
		// read names of the training tf	
			
			br = new BufferedReader(new FileReader(new File(trainscoresfile)));
			
			String line = null;
			
			StringTokenizer strtok = null;
			
			boolean first = true;
			String first_name = "";
			
			while((line = br.readLine()) != null) {
				
				strtok = new StringTokenizer(line);
														   
				if (first) {								// first training tf
					first_name = strtok.nextToken();
					first = false;
				}
				else if (! strtok.nextToken().equals(first_name)) break;
				
				strtok.nextToken(); 						// vs.
				
				tfnames.add(strtok.nextToken()); 			// second training tf
				
				strtok.nextToken();							// :
				
				strtok.nextToken();							// score
				
			}
				

			br.close();
			
			
			double[][] trainscores = new double[tfnames.size()][tfnames.size()];
			
			

		// read scores of the training tfs	
				
			br = new BufferedReader(new FileReader(new File(trainscoresfile)));
				

				
			for(int i=0; i<trainscores.length; i++) {
					
					
				for(int j=0; j<trainscores.length; j++) {
					
					strtok = new StringTokenizer(br.readLine());
						
					strtok.nextToken();							// training tf i
					strtok.nextToken(); 						// vs.
					strtok.nextToken(); 						// training tf j
					strtok.nextToken();							// :
						
					trainscores[i][j] = Double.parseDouble(strtok.nextToken());
					
				}
					
			}
					

			br.close();
				
			
		// calculate dot products for all pairs of TFs 
			
			bw = new BufferedWriter(new FileWriter(new File(outfile)));
			
			double dot_product;
				
			for(int i=0; i<tfnames.size()-1; i++) {
				
				for(int j=i+1; j<tfnames.size(); j++) {
					
					if (! relevant_pairs.contains(i + "_" + j)) continue;
					
					dot_product = 0.0;
					
					for(int k=0; k<tfnames.size(); k++) {
					
						dot_product += trainscores[i][k] * trainscores[j][k];
					}
					
					bw.write(tfnames.get(i) + " vs. " + tfnames.get(j) + " : " + dot_product + "\n");
					
				}		
			}
					
			bw.flush();
			bw.close();
			
		}
		catch(IOException ioe) {
			System.out.println(ioe.getMessage());
			System.out.println("IOException occurred while precalculating SVMPairwise Helpfiles.");
		}
		
	}
}

