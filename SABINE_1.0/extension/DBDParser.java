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
import java.util.StringTokenizer;

public class DBDParser {

	ArrayList<String> gene_IDs = new ArrayList<String>();
	ArrayList<String> prot_IDs = new ArrayList<String>();
	ArrayList<String> sequences = new ArrayList<String>();
	ArrayList<ArrayList<String>> domains = new ArrayList<ArrayList<String>>();
	ArrayList<ArrayList<String>> hmm_IDs = new ArrayList<ArrayList<String>>();
	ArrayList<String[]> uniprot_IDs = new ArrayList<String[]>();
	
	String species_ID = "";
	String species_name = "";
	
	public void parseProtIDs(String infile) {
		
		String line, curr_prot_ID;
		StringTokenizer strtok;
		
		try {
			 
			 BufferedReader br = new BufferedReader(new FileReader(new File(infile)));
			 
			 br.readLine();    // skip first line
			 
			 while ((line = br.readLine()) != null) {
				 
				 strtok = new StringTokenizer(line);
				 strtok.nextToken();
				 curr_prot_ID = strtok.nextToken().trim();
				 
				 
				 if (! prot_IDs.contains(curr_prot_ID)) prot_IDs.add(curr_prot_ID);
			 }
			 br.close();
		}
		catch(IOException ioe) {
			System.out.println(ioe.getMessage());
			System.out.println("IOException occurred while parsing Ensembl protein IDs.");
		}
	}
	
	public void parseSequences(String infile) {
		
		StringTokenizer strtok;
		String line, curr_prot_ID, curr_seq;
		int curr_idx;
		
		// initialize sequences and gene IDs
		for (int i=0; i<prot_IDs.size(); i++) {
			sequences.add("");
			gene_IDs.add("");
		}
		
		try {
			 
			 BufferedReader br = new BufferedReader(new FileReader(new File(infile)));
			 line = br.readLine();
			 
			 while(line != null) {
						
				if (! line.startsWith(">")) {
					System.out.println("Parse Error. \">\" expected at the beginning of the line.");
					System.exit(0);
				}
				 
				 strtok = new StringTokenizer(line.substring(1));   		// skip ">"
				 curr_prot_ID = strtok.nextToken().trim();					// parse protein ID
				 
				 if (! prot_IDs.contains(curr_prot_ID)) {					// no domain available -> skip sequence
					 
					 while (!((line = br.readLine()) == null) && ! line.startsWith(">"));
				 }
				 else { 
					 
					 curr_idx = prot_IDs.indexOf(curr_prot_ID);
					 
					 if (strtok.nextToken().startsWith("pep:")) {			// "pep:"							 
						 strtok.nextToken();								// "chromosome:"
					 
						 // parse Ensembl gene ID
						 System.out.println(strtok.nextToken().trim().split(":")[1]);
						 System.exit(0);
						 gene_IDs.set(curr_idx, strtok.nextToken().trim().split(":")[1]);
					 }
					 
					 // parse sequence
					 curr_seq = "";
					 while (!((line = br.readLine()) == null) && ! line.startsWith(">")) {
						 curr_seq += line; 
					 }
					 sequences.set(curr_idx, curr_seq);
				 }
			 }
		}
		catch(IOException ioe) {
			System.out.println(ioe.getMessage());
			System.out.println("IOException occurred while parsing sequences and DBDs.");
		}
	}

	
	public void parseDomains(String infile) {
		
		if (sequences.size() == 0) {
			System.out.println("Fatal Error. Unable to parse domains. Global variable \"gene_IDs\" was not initialized. Aborting.");
			System.exit(0);
		}
		
		StringTokenizer strtok;
		String line, curr_hmm_ID, curr_seq_ID, curr_domain;
		ArrayList<String> curr_domains, curr_hmm_IDs;
		
		// initialize domains and HMM-IDs
		for (int i=0; i<sequences.size(); i++) {
			hmm_IDs.add(new ArrayList<String>());
			domains.add(new ArrayList<String>());
		}

		int curr_idx;
		
		try {
			 
			 BufferedReader br = new BufferedReader(new FileReader(new File(infile)));
			 line = br.readLine();
			 //int line_counter = 1;
			 String[] splitted_domains;
			 
			 while ((line = br.readLine()) != null) {
				 
				 //if (++line_counter % 1000 == 0)
				 //    System.out.println("  Parsing line: " + line_counter);
				 
				 strtok = new StringTokenizer(line);
				 
				 curr_hmm_ID = strtok.nextToken();    // read HMM-ID
				 curr_seq_ID = strtok.nextToken();    // read Ensembl Protein-ID
				 curr_domain = strtok.nextToken(); 	  // read domain assignment
				 
				 if (prot_IDs.contains(curr_seq_ID)) {
					 curr_idx = prot_IDs.indexOf(curr_seq_ID);

					 curr_domains = domains.get(curr_idx);					// save domain
					 curr_hmm_IDs = hmm_IDs.get(curr_idx);					// save HMM-ID
					 
					 if (curr_domain.matches("\\d*\\p{Punct}\\d*")) {       // only one domain ?
						 curr_domains.add(curr_domain.replace("-", "\t"));
					 	 curr_hmm_IDs.add(curr_hmm_ID);
				 	 }	 
					 else {
						 splitted_domains = curr_domain.split(",");
						 
						 for (int s=0; s<splitted_domains.length; s++){
							 curr_domains.add(splitted_domains[s].replace("-", "\t"));
							 curr_hmm_IDs.add(curr_hmm_ID);
						 }
					 }
					 domains.set(curr_idx, curr_domains);
					 hmm_IDs.set(curr_idx, curr_hmm_IDs);
				 }									  		
				 
			 }
			 
		}
		catch(IOException ioe) {
			System.out.println(ioe.getMessage());
			System.out.println("IOException occurred while parsing domain assignments.");
		}
	}
	
	
	
	
	public static void main(String[] args) {
		
	}

}


