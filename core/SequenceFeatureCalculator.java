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

import help.RawDataPreprocessor;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.StringTokenizer;

import org.biojava.bio.BioException;


public class SequenceFeatureCalculator {
	
	
	SequenceAligner aligner = null;
	
	SequenceAligner aligner_dna = null;
	
	RawDataPreprocessor helper = new RawDataPreprocessor();
	
	
	public boolean silent = false;
	public boolean predicted_domains = false;
	public String basedir = "internal/";
	
	ArrayList<ArrayList<String>> other_domains = null;
	
	ArrayList<String> other_sequences1 = null;
	ArrayList<String> other_sequences2 = null;
	
	ArrayList<String> other_secstructs1 = null;
	ArrayList<String> other_secstructs2 = null;
	
	ArrayList<String> other_names = null;
	
	
	/*
	 * 
	 * initializes global variables "other_names" and "other_domains"
	 * 
	 */
	
	public void parseRelevantDomainsAndSequences(boolean[] irrelevantPairs, String class_id, String train_dir) {
		
		
		other_domains 	= new ArrayList<ArrayList<String>>();
		
		other_sequences1 = new ArrayList<String>();
		other_sequences2 = new ArrayList<String>();
		
		other_names 	= new ArrayList<String>();
		
		
		String line = null;
			
		int entry_counter = 0;
		 
		BufferedReader br = null;
		
		 
		try {
			 
			 br = new BufferedReader(new FileReader(new File(train_dir + "trainingset_" + class_id + ".rawdata")));
			 
			 
			/*
			 * 
			 * parse all relevant domains and sequences of the training set
			 * 
			 */
				
			while((line = br.readLine()) != null) {
					
				
				if(!line.startsWith("NA  ")) {
					System.out.println("Parse Error. NA expected. Aborting.");
					System.exit(0);
				}
					
					
			// tf pair is irrelevant (i.e. BLOSUM-score < 0.3) ? --> skip this tf entry !	
					
				if(irrelevantPairs[entry_counter++])  {
								
					while(!line.startsWith("//")) { line = br.readLine(); }
								
					line = br.readLine(); 		// XX
								
					continue;
								
				}
							
							
			// tf pair is relevant ! --> parse domain information of this tf entry !
							
				other_names.add(line.substring(4));
					
				if(!silent) System.out.println("Parsing domains and sequences of " + line.substring(3) + ".");
						
				ArrayList<String> domaingroup = new ArrayList<String>();
					
				line = br.readLine(); 		// XX
					
					
			// parse domains
					
				while((line = br.readLine()).startsWith("DO  ")) {
							
					domaingroup.add(line.substring(4));
						
				}
					
				other_domains.add(domaingroup);
					
				line = br.readLine();
					
					
			// parse first sequence	
					
				if(line.startsWith("S1")) {
						
					String seq1 = "";
						
					while(line.startsWith("S1")) {
							
						seq1 += line.substring(4);
							
						line = br.readLine();
							
					}
						
					other_sequences1.add(seq1);
						
					line = br.readLine();
						
				}
					
				else {
						
					other_sequences1.add("NO SEQUENCE.");
						
				}
					
					
			// parse second sequence	
					
				if(line.startsWith("S2")) {
						
					String seq2 = "";
						
					while(line.startsWith("S2")) {
							
						seq2 += line.substring(4);
							
						line = br.readLine();
							
					}
						
					other_sequences2.add(seq2);
						
					line = br.readLine();
						
				}
					
				else {
						
					other_sequences2.add("NO SEQUENCE.");
						
				}
					
					
				while(!line.startsWith("//")) { line = br.readLine(); };
					
				line = br.readLine(); 		// XX
						
			
			}	// proceed to next tf entry
					
					
					
			if(!silent) System.out.println("\ndone.\n");
			
			
			br.close();
			 
		}
		catch(IOException ioe) {
			System.out.println(ioe.getMessage());
			System.out.println("IOException occurred while parsing relevant sequences and domains.");
		}
		 
		 
	}
	
	
	
	/*
	 * 
	 * initializes global variables "other_names" and "other_domains"
	 * 
	 */
	
	public void parseRelevantSecondaryStructures(boolean[] irrelevantPairs, String class_id, String train_dir) {
		
		
		other_secstructs1 = new ArrayList<String>();
		other_secstructs2 = new ArrayList<String>();
		
		
		String line = null;
			
		int entry_counter = 0;
		 
		BufferedReader br = null;
		
		 
		try {
			 
			 br = new BufferedReader(new FileReader(new File(train_dir + "trainingset_" + class_id + ".secstruct")));
			 
			 
			/*
			 * 
			 * parse all relevant domains and sequences of the training set
			 * 
			 */
				
			while((line = br.readLine()) != null) {
					
				
				if(!line.startsWith("NA  ")) {
					System.out.println("Parse Error. DE expected. Aborting.");
					System.exit(0);
				}
					
					
			// tf pair is irrelevant (i.e. BLOSUM-score < 0.3) ? --> skip this tf entry !	
					
				if(irrelevantPairs[entry_counter++])  {
								
					while(!line.startsWith("//")) { line = br.readLine(); }
								
					line = br.readLine(); 		// XX
								
					continue;
								
				}
							
							
			// tf pair is relevant ! --> parse domain information of this tf entry !
							
					
				if(!silent) System.out.println("Parsing domains and sequences of " + line.substring(3) + ".");
						
					
				line = br.readLine(); 		// XX
					
					
			// parse domains
					
				while((line = br.readLine()).startsWith("DO  ")) ;
					
				line = br.readLine();
					
					
			// parse first sequence	
					
				if(line.startsWith("S1")) {
						
					String ss1 = "";
						
					while(line.startsWith("S1")) {
							
						ss1 += line.substring(4);
							
						line = br.readLine();
							
					}
						
					other_secstructs1.add(ss1);
						
					line = br.readLine();
						
				}
					
				else {
						
					other_secstructs1.add("NO SEQUENCE.");
						
				}
					
					
			// parse second sequence	
					
				if(line.startsWith("S2")) {
						
					String ss2 = "";
						
					while(line.startsWith("S2")) {
							
						ss2 += line.substring(4);
							
						line = br.readLine();
							
					}
						
					other_secstructs2.add(ss2);
						
					line = br.readLine();
						
				}
					
				else {
						
					other_secstructs2.add("NO SEQUENCE.");
						
				}
					
					
				while(!line.startsWith("//")) { line = br.readLine(); };
					
				line = br.readLine(); 		// XX
						
			
			}	// proceed to next tf entry
					
					
					
			if(!silent) System.out.println("\ndone.\n");
			
			
			br.close();
			 
		}
		catch(IOException ioe) {
			System.out.println(ioe.getMessage());
			System.out.println("IOException occurred while parsing relevant secondary structures.");
		}
		 
		 
	}
	
	
	public void calculateSequenceFeatureFile(String tfname, ArrayList<String> domains, String sequence1, String sequence2, String featuretype, String[] params, String matrix, String outfile) {
		
		
		if(other_names == null) {	
			System.out.println("Fatal Error. Unable to calculate domain features. Global variable \"other_names\" was not initialized. Aborting.");
			System.exit(0);
		}
		
		if(other_domains == null) {	
			System.out.println("Fatal Error. Unable to calculate domain features. Global variable \"other_domains\" was not initialized. Aborting.");
			System.exit(0);
		}
		
		if(other_sequences1 == null) {	
			System.out.println("Fatal Error. Unable to calculate domain features. Global variable \"other_sequences1\" was not initialized. Aborting.");
			System.exit(0);
		}
		
		if(other_sequences2 == null) {	
			System.out.println("Fatal Error. Unable to calculate domain features. Global variable \"other_sequences2\" was not initialized. Aborting.");
			System.exit(0);
		}
		
		
		BufferedWriter bw = null;		
		
		helper.basedir = basedir;
		
		
		try {
			
			String secstruct1 = null;
			String secstruct2 = null;
			
			
			bw = new BufferedWriter(new FileWriter(new File(outfile)));
			
			if(matrix != null) {
				
				aligner = new SequenceAligner(matrix, "NW");
			
			}
			
			if(featuretype.equals("SecondaryStructure")) {
				
				aligner_dna = new SequenceAligner("DNA", "SS-matrix.dat", "NW");
				
			}
			
			
			double maxScore = 0.0;
			double score = 0.0;	

			
			int[] maxPair = new int[] {0,0} ;
			
			StringTokenizer strtok = null;
			
			
			/*
			 * 
			 * compare domains of input tf vs. all relevant domains of training set
			 * 
			 */ 
			
			ArrayList<String> all_domains = new ArrayList<String>();
			if (predicted_domains) {
				for(int i=0; i<domains.size(); i++) {
					all_domains.add(domains.get(i));
				}
			}
			
	
			if(!silent) System.out.println("Comparing query");
				
			
			for(int j=0; j<other_domains.size(); j++) {
					
					
				if(!silent) System.out.println("  with " + other_names.get(j) + ".");
						
				maxScore = Double.NEGATIVE_INFINITY; 
						
				if (predicted_domains) {
					domains = new ArrayList<String>();
					domains.add(all_domains.get(j));
				}
						
			// determine score-maximizing domain pair	
						
					for(int k=0; k<domains.size(); k++) {
							
						for(int l=0; l<other_domains.get(j).size(); l++) {
								
								
							score = aligner.getSMBasedSimilarity  (new StringTokenizer(domains.get(k)).nextToken(), new StringTokenizer(other_domains.get(j).get(l)).nextToken());
								
								
							if(score > maxScore)  {
									
								maxScore = score;
									
								maxPair[0] = k;
								maxPair[1] = l;
								
							}
								
						}
						
					}
					
					
				// get the SM-Score-maximizing domains	
					
					String dom1 = domains.get(maxPair[0]);
					String dom2 = other_domains.get(j).get(maxPair[1]);
					
					String seq1 = null;
					String seq2 = null;
					
				// store the respective protein sequences in seq1 & seq2
					
					strtok = new StringTokenizer(dom1);
					strtok.nextToken();
					
					if(strtok.nextToken().equals("1")) {
						seq1 = sequence1;
					}
					else {
						seq1 = sequence2;
					}
					
					strtok = new StringTokenizer(dom2);
					strtok.nextToken();
					
					if(strtok.nextToken().equals("1")) {
						seq2 = other_sequences1.get(j);
					}
					else {
						seq2 = other_sequences2.get(j);
					}
					
						
			/*
			 * 
			 * calculate environment similarity
			 * 
			 */ 
						
					if(featuretype.equals("Environments")) {
						
						
					// get the left and right environments of the SM-Score-maximizing binding domains
						
						strtok = new StringTokenizer(domains.get(maxPair[0]));
						strtok.nextToken();
						strtok.nextToken();
						
						
					// start and end indices of domain 1	
						
						int start1 = Integer.parseInt(strtok.nextToken());
						int end1   = Integer.parseInt(strtok.nextToken());
						
						
						strtok = new StringTokenizer(other_domains.get(j).get(maxPair[1]));
						strtok.nextToken();
						strtok.nextToken();
					
						
					// start and end indices of domain 2	
						
						int start2 = Integer.parseInt(strtok.nextToken());
						int end2   = Integer.parseInt(strtok.nextToken());
						
						int left_range  = Math.min(Integer.parseInt(params[0]), Math.min(start1, start2));
						
						int right_range = Math.min(Integer.parseInt(params[0]), Math.min(seq1.length() - end1 , seq2.length() - end2));
						
						
					// left environment of domain 1	
						
						int envstart = start1 - left_range;
						int envend   = start1;
						
						String lenv1 = seq1.substring(envstart, envend);
						
						
					// right environment of domain 1	
						
						envstart = end1;
						envend   = end1 + right_range;
						
						String renv1 = seq1.substring(envstart, envend);
						
						
					// left environment of domain 2	
						
						envstart = start2 - left_range;
						envend   = start2;
						
						String lenv2 = seq2.substring(envstart, envend);
						
						
					// right environment of domain 2	
						
						envstart = end2;
						envend   = end2 + right_range;
						
						String renv2 = seq2.substring(envstart, envend);
						
						
					// calculate the left / right environment based similarity scores	
						
						double left_score  = aligner.getSMBasedSimilarity(lenv1, lenv2);
						double right_score = aligner.getSMBasedSimilarity(renv1, renv2);
						
						
					// calculate the average environment score
						
						double avg_score = 0.5 * (left_score + right_score);
						
						if(lenv1.length() == 0 || lenv2.length() == 0) {
							
							avg_score = right_score;
							
						}
						
						if(renv1.length() == 0 || renv2.length() == 0) {
							
							avg_score = left_score;
							
						}
							
						
						if(lenv1.length() == 0 && renv1.length() == 0 || lenv2.length() == 0 && renv2.length() == 0) {
							
							avg_score = 0.0;
							
						}
						
						
						
					// print score to output file
						
						bw.write(tfname + " vs. " + other_names.get(j) + " :\t" + avg_score + "\n");
					
					}
					
					
				/*
				 * 
				 * calculate secondary structure similarity
				 * 
				*/ 	
					
					
				if(featuretype.equals("SecondaryStructure")) {
							
		
			// look-up the secondary structures
						
					String ss1 = null;
					String ss2 = null;
					
				
				// store the respective protein sequences in seq1 & seq2
					
					strtok = new StringTokenizer(dom1);
					strtok.nextToken();
					
					if(strtok.nextToken().equals("1")) {
						
						if(secstruct1 == null) {
							
							secstruct1 = helper.getSecondaryStructure(sequence1).replace("E", "A").replace("H", "T");
							
							ss1 = secstruct1;
							
						}
						
						else {
						
							ss1 = secstruct1;
						
						}
						
					}
					else {
						
						if(secstruct2 == null) {
							
							secstruct2 = helper.getSecondaryStructure(sequence2).replace("E", "A").replace("H", "T");
							
							ss1 = secstruct2;
							
						}
						
						else {
						
							ss1 = secstruct2;
						
						}
					
					}
					
					
					strtok = new StringTokenizer(dom2);
					strtok.nextToken();
					
					if(strtok.nextToken().equals("1")) {
						ss2 = other_secstructs1.get(j).replace("E", "A").replace("H", "T");
					}
					else {
						ss2 = other_secstructs2.get(j).replace("E", "A").replace("H", "T");
					}
					
	
				// calculate the secondary structure similarity score	
							
					double ss_score = aligner_dna.getSMBasedSecondaryStructureSimilarity(ss1, ss2);
							
							
				// print score to output file
							
					bw.write(tfname + " vs. " + other_names.get(j) + " :\t" + ss_score + "\n");
							
						
						
				}
					
			
			}
			
			bw.flush();
			bw.close();

			
		}
		catch(IOException ioe) {
			System.out.println(ioe.getMessage());
			System.out.println("IOException occurred while comparing sequences.");
		}
		catch(BioException bioe) {
			System.out.println(bioe.getMessage());
			System.out.println("BioException occurred while comparing sequences.");
		}
		
		
	}
	
	
	public void calculateAllSequenceFeatures(String featuretype, String[] params, String matrix, String outfile, ArrayList<String> relevant_pairs) {
		
		
		if(other_names == null) {	
			System.out.println("Fatal Error. Unable to calculate domain features. Global variable \"other_names\" was not initialized. Aborting.");
			System.exit(0);
		}
		
		if(other_domains == null) {	
			System.out.println("Fatal Error. Unable to calculate domain features. Global variable \"other_domains\" was not initialized. Aborting.");
			System.exit(0);
		}
		
		if(other_sequences1 == null) {	
			System.out.println("Fatal Error. Unable to calculate domain features. Global variable \"other_sequences1\" was not initialized. Aborting.");
			System.exit(0);
		}
		
		if(other_sequences2 == null) {	
			System.out.println("Fatal Error. Unable to calculate domain features. Global variable \"other_sequences2\" was not initialized. Aborting.");
			System.exit(0);
		}
		
		
		BufferedWriter bw = null;		
		
		
		try {
				
			
			bw = new BufferedWriter(new FileWriter(new File(outfile)));
			
			if(matrix != null) {
				
				aligner = new SequenceAligner(matrix, "NW");
			
			}
			
			if(featuretype.equals("SecondaryStructure")) {
				
				aligner_dna = new SequenceAligner("DNA", "SS-matrix.dat", "NW");
				
			}
			
			
			double maxScore = 0.0;
			double score = 0.0;	

			
			int[] maxPair = new int[] {0,0} ;
			
			StringTokenizer strtok = null;
			
			
			/*
			 * 
			 * compare domains of all pairs of TFs in training set
			 * 
			 */ 
			
	
			if(!silent) System.out.println("Comparing query");
				
			for(int i=0; i<other_domains.size()-1; i++) {
				
				for(int j=i+1; j<other_domains.size(); j++) {
					
					
					if (! relevant_pairs.contains(i + "_" + j)) continue;
					
					if(!silent) System.out.println("  with " + other_names.get(j) + ".");
						
					maxScore = Double.NEGATIVE_INFINITY; 
						
						
					// determine score-maximizing domain pair	
						
					for(int k=0; k<other_domains.get(i).size(); k++) {
							
						for(int l=0; l<other_domains.get(j).size(); l++) {
								
								
							score = aligner.getSMBasedSimilarity  (new StringTokenizer(other_domains.get(i).get(k)).nextToken(), new StringTokenizer(other_domains.get(j).get(l)).nextToken());
								
								
							if(score > maxScore)  {
									
								maxScore = score;
									
								maxPair[0] = k;
								maxPair[1] = l;
								
								}
								
							}
						
						}
					
					
					// get the SM-Score-maximizing domains	
					
						String dom1 = other_domains.get(i).get(maxPair[0]);
						String dom2 = other_domains.get(j).get(maxPair[1]);
					
						String seq1 = null;
						String seq2 = null;
					
					// store the respective protein sequences in seq1 & seq2
					
						strtok = new StringTokenizer(dom1);
						strtok.nextToken();
					
						if(strtok.nextToken().equals("1")) {
							seq1 = other_sequences1.get(i);
						}
						else {
							seq1 = other_sequences2.get(i);
						}
					
						strtok = new StringTokenizer(dom2);
						strtok.nextToken();
					
						if(strtok.nextToken().equals("1")) {
							seq2 = other_sequences1.get(j);
						}
						else {
							seq2 = other_sequences2.get(j);
						}
					
						
					/*
			 	* 
			 	* calculate environment similarity
			 	* 
			 	*/ 
						
						if(featuretype.equals("Environments")) {
						
						
						// get the left and right environments of the SM-Score-maximizing binding domains
						
							strtok = new StringTokenizer(other_domains.get(i).get(maxPair[0]));
							strtok.nextToken();
							strtok.nextToken();
						
						
						// start and end indices of domain 1	
						
							int start1 = Integer.parseInt(strtok.nextToken());
							int end1   = Integer.parseInt(strtok.nextToken());
						
						
							strtok = new StringTokenizer(other_domains.get(j).get(maxPair[1]));
							strtok.nextToken();
							strtok.nextToken();
					
						
						// start and end indices of domain 2	
						
							int start2 = Integer.parseInt(strtok.nextToken());
							int end2   = Integer.parseInt(strtok.nextToken());
						
							int left_range  = Math.min(Integer.parseInt(params[0]), Math.min(start1, start2));
						
							int right_range = Math.min(Integer.parseInt(params[0]), Math.min(seq1.length() - end1 , seq2.length() - end2));
						
						
						// left environment of domain 1	
						
							int envstart = start1 - left_range;
							int envend   = start1;
						
							String lenv1 = seq1.substring(envstart, envend);
						
						
						// right environment of domain 1	
						
							envstart = end1;
							envend   = end1 + right_range;
						
							String renv1 = seq1.substring(envstart, envend);
						
						
						// left environment of domain 2	
						
							envstart = start2 - left_range;
							envend   = start2;
						
							String lenv2 = seq2.substring(envstart, envend);
						
						
						// right environment of domain 2	
						
							envstart = end2;
							envend   = end2 + right_range;
						
							String renv2 = seq2.substring(envstart, envend);
						
						
						// calculate the left / right environment based similarity scores	
						
							double left_score  = aligner.getSMBasedSimilarity(lenv1, lenv2);
							double right_score = aligner.getSMBasedSimilarity(renv1, renv2);
						
						
						// calculate the average environment score
						
							double avg_score = 0.5 * (left_score + right_score);
						
							if(lenv1.length() == 0 || lenv2.length() == 0) {
							
								avg_score = right_score;
							
							}
						
							if(renv1.length() == 0 || renv2.length() == 0) {
							
								avg_score = left_score;
							
							}
							
						
							if(lenv1.length() == 0 && renv1.length() == 0 || lenv2.length() == 0 && renv2.length() == 0) {
							
								avg_score = 0.0;
							
							}
						
						
						
						// print score to output file
						
							bw.write(other_names.get(i) + " vs. " + other_names.get(j) + " :\t" + avg_score + "\n");
					
						}
					
					
				   /*
				 	* 
				 	* calculate secondary structure similarity
				 	* 
				 	*/ 	
					
					
					if(featuretype.equals("SecondaryStructure")) {
							
		
					// look-up the secondary structures
						
						String ss1 = null;
						String ss2 = null;
					
				
					// store the respective protein sequences in seq1 & seq2
					
						strtok = new StringTokenizer(dom1);
						strtok.nextToken();
					
						if(strtok.nextToken().equals("1")) {
							ss1 = other_secstructs1.get(i).replace("E", "A").replace("H", "T");
						}
						else {
							ss1 = other_secstructs2.get(i).replace("E", "A").replace("H", "T");
						}
					
					
						strtok = new StringTokenizer(dom2);
						strtok.nextToken();
					
						if(strtok.nextToken().equals("1")) {
							ss2 = other_secstructs1.get(j).replace("E", "A").replace("H", "T");
						}
						else {
							ss2 = other_secstructs2.get(j).replace("E", "A").replace("H", "T");
						}
					
	
					// calculate the secondary structure similarity score	
							
						double ss_score = aligner_dna.getSMBasedSecondaryStructureSimilarity(ss1, ss2);
							
							
					// print score to output file
							
						bw.write(other_names.get(i) + " vs. " + other_names.get(j) + " :\t" + ss_score + "\n");
							
						
						
					}
					
			
				}
			}
			bw.flush();
			bw.close();

			
		}
		catch(IOException ioe) {
			System.out.println(ioe.getMessage());
			System.out.println("IOException occurred while comparing sequences.");
		}
		catch(BioException bioe) {
			System.out.println(bioe.getMessage());
			System.out.println("BioException occurred while comparing sequences.");
		}

	}
}

