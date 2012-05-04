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
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.StringTokenizer;

import org.biojava.bio.BioException;

import de.zbit.util.progressbar.AbstractProgressBar;


public class DomainFeatureCalculator {
	
	
	SequenceAligner aligner = null;
	
	public boolean silent = false;
	public boolean predicted_domains = false;
	public String basedir = null;
	
	ArrayList<ArrayList<String>> other_domains = null;
	
	ArrayList<String> other_names = null;
	
	/**
	 * 
	 */
	AbstractProgressBar progress = null;
	
	public void setProgressBar(AbstractProgressBar progress) {
	  this.progress = progress;
	}
	
	
	/*
	 * 
	 * initializes global variables "other_names" and "other_domains"
	 * 
	 */
	
	public ArrayList<ArrayList<String>> get_other_domains() {
		return other_domains;
	}
	
	public void parseRelevantDomains(boolean[] irrelevantPairs, String class_id, String train_dir) {
		
		
		other_domains = new ArrayList<ArrayList<String>>();
		
		other_names = new ArrayList<String>();
		
		
		String line = null;
			
		int entry_counter = 0;
		 
		BufferedReader br = null;
		
		 
		try {
			 
			 br = new BufferedReader(new FileReader(new File(train_dir + "trainingset_" + class_id + ".rawdata")));
			 
			/*
			 * 
			 * parse all relevant domains of the training set
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
						
				other_names.add(line.substring(4));
						
				if(!silent) System.out.println("Parsing domains of " + line.substring(4) + ".");
							
				ArrayList<String> domaingroup = new ArrayList<String>();
						
				line = br.readLine(); 		// XX
						
						
			// parse domains of this tf entry
						
				while((line = br.readLine()).startsWith("DO  ")) {
								
					domaingroup.add(new StringTokenizer(line.substring(4)).nextToken());
							
				}
						
				other_domains.add(domaingroup);
						
						
				while(!line.startsWith("//")) { line = br.readLine(); }
						
						
				line = br.readLine(); 		// XX
							
						
			}	// proceed to next tf entry
					
					
					
			if(!silent) System.out.println("\ndone.\n");
					
			 
		}
		catch(IOException ioe) {
			System.out.println(ioe.getMessage());
			System.out.println("IOException occurred while parsing relevant domains.");
		}
		 
		 
	}
	
	
	public void calculateDomainFeatureFile(String tfname, ArrayList<String> domains1, String featuretype, String[] params, String matrix, String outfile) {
		
	  if (progress!=null) {
	    progress.DisplayBar();
	  }
		
		if(other_names == null) {	
			System.out.println("Fatal Error. Unable to calculate domain features. Global variable \"other_names\" was not initialized. Aborting.");
			System.exit(0);
		}
		
		if(other_domains == null) {	
			System.out.println("Fatal Error. Unable to calculate domain features. Global variable \"other_domains\" was not initialized. Aborting.");
			System.exit(0);
		}
		
		
		ArrayList<String> domains = new ArrayList<String>();
		
		for(int i=0; i<domains1.size(); i++) {
			
			domains.add(new StringTokenizer(domains1.get(i)).nextToken());
			
		}
		
		ArrayList<String> all_domains = new ArrayList<String>();
		if (predicted_domains) {
			for(int i=0; i<domains.size(); i++) {
				all_domains.add(domains.get(i));
			}
		}
		
		
		for(int i=0; i<domains1.size(); i++) {
			
			domains.add(new StringTokenizer(domains1.get(i)).nextToken());
			
		}
		
		BufferedWriter bw = null;
		
		try {			
			
			bw = new BufferedWriter(new FileWriter(new File(outfile)));
			
			
			if(matrix != null) {
				
				aligner = new SequenceAligner(matrix, "NW");
			
			}
			
			
			double maxScore = 0.0;
			double score = 0.0;	
			
			
		/*
		 * 
		 * compare domains of input tf vs. all relevant domains of training set
		 * 
		 */ 
				
			if(!silent) System.out.println("Comparing query" );
				
			for(int j=0; j<other_domains.size(); j++) {
					
				if(!silent) System.out.println("  with " + other_names.get(j) + ".");
				
				maxScore = Double.NEGATIVE_INFINITY; 
				
				if (predicted_domains) {
					domains = new ArrayList<String>();
					domains.add(all_domains.get(j));
				}
					
				for(int k=0; k<domains.size(); k++) {
						
					for(int l=0; l<other_domains.get(j).size(); l++) {
						
						if (featuretype.equals("SequenceIdentity"))  	score = aligner.getSequenceSimilarity (domains.get(k), other_domains.get(j).get(l) );
						if (featuretype.equals("SMBasedIdentity"))   	score = aligner.getSMBasedIdentity    (domains.get(k), other_domains.get(j).get(l), Double.parseDouble(params[0]) );
						if (featuretype.equals("SMBasedSimilarity")) 	score = aligner.getSMBasedSimilarity  (domains.get(k), other_domains.get(j).get(l) );
						if (featuretype.equals("LocalAlignmentKernel"))	score = aligner.getNormalizedLAKernel (domains.get(k), other_domains.get(j).get(l), params[0]);
						if (featuretype.equals("MismatchKernel"))		score = aligner.getNormalizedMMKernel (domains.get(k), other_domains.get(j).get(l), Integer.parseInt(params[0]), Integer.parseInt(params[1]), basedir);
						
						if(score > maxScore) maxScore = score;
					}
				}
				bw.write(tfname + " vs. " + other_names.get(j) + " :\t" + maxScore + "\n");
			}
				
			
			bw.flush();
			bw.close();
		
		}
		catch(IOException ioe) {
			System.out.println(ioe.getMessage());
			System.out.println("IOException occurred while comparing domains.");
		}
		catch(BioException bioe) {
			System.out.println(bioe.getMessage());
			System.out.println("BioException occurred while comparing domains.");
		}
		
	}
	
	public void calculateAllDomainFeatures(String featuretype, String[] params, String matrix, String outfile, ArrayList<String> relevant_pairs) {
		
		
		if(other_names == null) {	
			System.out.println("Fatal Error. Unable to calculate domain features. Global variable \"other_names\" was not initialized. Aborting.");
			System.exit(0);
		}
		
		if(other_domains == null) {	
			System.out.println("Fatal Error. Unable to calculate domain features. Global variable \"other_domains\" was not initialized. Aborting.");
			System.exit(0);
		}
		
		BufferedWriter bw = null;
		
		
		try {
			
			
			bw = new BufferedWriter(new FileWriter(new File(outfile)));
			
			
			if(matrix != null) {
				
				aligner = new SequenceAligner(matrix, "NW");
			
			}
			
			
			double maxScore = 0.0;
			double score = 0.0;	
			
			
		/*
		 * 
		 * compare domains of all pairs of TFs in training set
		 * 
		 */ 
				
			if(!silent) System.out.println("Comparing query" );
			
			for(int i=0; i<other_domains.size()-1; i++) { 
				
				for(int j=i+1; j<other_domains.size(); j++) {
					
					if (! relevant_pairs.contains(i + "_" + j)) continue;
					
					if(!silent) System.out.println("  with " + other_names.get(j) + ".");
				
					maxScore = Double.NEGATIVE_INFINITY; 
					
					for(int k=0; k<other_domains.get(i).size(); k++) {
						
						for(int l=0; l<other_domains.get(j).size(); l++) {
						
							if (featuretype.equals("SequenceIdentity"))  	score = aligner.getSequenceSimilarity (other_domains.get(i).get(k), other_domains.get(j).get(l) );
							if (featuretype.equals("SMBasedIdentity"))   	score = aligner.getSMBasedIdentity    (other_domains.get(i).get(k), other_domains.get(j).get(l), Double.parseDouble(params[0]) );
							if (featuretype.equals("SMBasedSimilarity")) 	score = aligner.getSMBasedSimilarity  (other_domains.get(i).get(k), other_domains.get(j).get(l) );
							if (featuretype.equals("LocalAlignmentKernel"))	score = aligner.getNormalizedLAKernel (other_domains.get(i).get(k), other_domains.get(j).get(l), params[0]);
							if (featuretype.equals("MismatchKernel"))		score = aligner.getNormalizedMMKernel (other_domains.get(i).get(k), other_domains.get(j).get(l), Integer.parseInt(params[0]), Integer.parseInt(params[1]), basedir );
						
						
						if(score > maxScore) maxScore = score;
						
						
						}
					}
					
					bw.write(other_names.get(i) + " vs. " + other_names.get(j) + " :\t" + maxScore + "\n");
				
				}
			}
			
			bw.flush();
			bw.close();
		
		}
		catch(IOException ioe) {
			System.out.println(ioe.getMessage());
			System.out.println("IOException occurred while comparing domains.");
		}
		catch(BioException bioe) {
			System.out.println(bioe.getMessage());
			System.out.println("BioException occurred while comparing domains.");
		}
		
	}
	
	
	public int getTrainingSetSize(String class_id, String train_dir) {
        
		String line = null;
        int entry_counter = 0;
        BufferedReader br = null;
              	
        try {
        	br = new BufferedReader(new FileReader(new File(train_dir + "trainingset_" + class_id + ".rawdata")));
        } catch (FileNotFoundException e) {
        	System.out.println(e.getMessage());
        	System.out.println("Error occured before counting entries in training set");
        }
        try {
        	while((line = br.readLine()) != null) {
        		if(line.startsWith("NA  ")) {
        			entry_counter++;
        		}
        	}
      
        } catch (IOException e) {
        	System.out.println(e.getMessage());
        	System.out.println("Error occured while counting entries in training set");
        }
        return entry_counter;
	}   
}

