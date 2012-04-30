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


import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.StringTokenizer;

import optimization.MoStaOptimizer;

public class LabelFileGenerator {

	public String basedir = "";
	
	public ArrayList<ArrayList<String>> getAllFBPs(String infile) { 

		ArrayList<ArrayList<String>> fbps = new ArrayList<ArrayList<String>>();
		ArrayList<String> curr_fbp;
		
		BufferedReader br = null;
		String line = null;
		
		try {
			
			br = new BufferedReader(new FileReader(new File(infile)));
			
			while((line = br.readLine()) != null) {
				
				if (! line.startsWith("DE\t")) {
					System.out.println("Parse Error. Error occured while parsing FBPs. \nLine: " + line 
										+ "\nexpected to start with DE.");
					System.exit(0);
				}
				
				curr_fbp = new ArrayList<String>();
				
				while(! (line = br.readLine()).startsWith("XX")) {
					
					curr_fbp.add(line);
				
				}
				
				fbps.add(curr_fbp);
				
			}
			
			br.close();
				
			
		} catch (IOException e) {
			System.out.println(e.getMessage());
			System.out.println("Exception occured while parsing FBPs.");
		}
	
		return fbps;
	}
	
	
	
	public ArrayList<String> getFactorNames(String infile) {
	
		ArrayList<String> tfnames = new ArrayList<String>();
		
		BufferedReader br = null;
		String line = null;
		StringTokenizer strtok = null;
		
		try {
			
			br = new BufferedReader(new FileReader(new File(infile)));
			
			while((line = br.readLine()) != null) {
				
				if (! line.startsWith("DE\t")) {
					System.out.println("Parse Error. Error occured while parsing FBPs. \nLine: " + line 
										+ "\nexpected to start with DE.");
					System.exit(0);
				}
				
				strtok = new StringTokenizer(line);
				
				strtok.nextToken();                   // DE	
				
				tfnames.add(strtok.nextToken());	  // name
				
				while(! (line = br.readLine()).startsWith("XX"));    // go to next TF  
				
			}
			
			br.close();
				
			
		} catch (IOException e) {
			System.out.println(e.getMessage());
			System.out.println("Exception occured while parsing TF names.");
		}
	
		return tfnames;
		
	}
	
	public void calculateLabelFile(String fbpfile, String outfile, ArrayList<String> relevant_pairs) {
	
		BufferedWriter bw = null;
		MoStaOptimizer optimizer = new MoStaOptimizer();
		optimizer.basedir = basedir;
		double curr_score, curr_score_ii, curr_score_jj, curr_score_ij;
		
		ArrayList<ArrayList<String>> fbps = getAllFBPs(fbpfile);
		ArrayList<String> tfnames = getFactorNames(fbpfile);
		
		try {
			bw = new BufferedWriter(new FileWriter(new File(outfile)));
		
			for(int i=0; i<fbps.size()-1; i++) {
			
				for(int j=i+1; j<fbps.size(); j++) {
				
					if (! relevant_pairs.contains(i + "_" + j)) continue;
					
					curr_score_ij = optimizer.compareMatrices(fbps.get(i), fbps.get(j), null);
					curr_score_ii = optimizer.compareMatrices(fbps.get(i), fbps.get(i), null);
					curr_score_jj = optimizer.compareMatrices(fbps.get(j), fbps.get(j), null);
					
					curr_score = curr_score_ij / Math.sqrt(curr_score_ii * curr_score_jj);
				
					bw.write(tfnames.get(i) + " vs. " + tfnames.get(j) + " : " + curr_score + "\n");
				}
			}
			bw.flush();
			bw.close();
			
		} catch (IOException e) {
			System.out.println(e.getMessage());
			System.out.println("Exception occured while writing label file.");
		}
	}
}

