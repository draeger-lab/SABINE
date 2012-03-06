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
import java.util.Iterator;
import java.util.StringTokenizer;

public class UniprotParser {
	
	ArrayList<String> uniprot_IDs = new ArrayList<String>();
	ArrayList<String> sequences = new ArrayList<String>();
	ArrayList<ArrayList<String>> domains = new ArrayList<ArrayList<String>>();
	ArrayList<ArrayList<String>> domain_classes = new ArrayList<ArrayList<String>>();
	
	String species_name = null;
	
	public void parseSequences(String infile) {
		
		String curr_seq, curr_ID, curr_species, line;
		StringTokenizer strtok;
		
		try {
			 
			 BufferedReader br = new BufferedReader(new FileReader(new File(infile)));
			 line = br.readLine();
			 
			 while (line != null) {
				 
				 if (! line.startsWith(">")) {
					 System.out.println("Parse error. \">\" expected at the beginning of the line.");
				 }
				 
				 // parse uniprot ID and species
				 strtok = new StringTokenizer(line, "|");
				 strtok.nextToken();
				 
				 curr_ID = strtok.nextToken().trim();
				 curr_species = line.split(" OS=")[1].split(" ")[0] + " " + line.split(" OS=")[1].split(" ")[1];
				 
				 // skip factors of other species
				 if (! curr_species.toUpperCase().equals(species_name.toUpperCase())) {
					 while((line=br.readLine()) != null && !line.startsWith(">"));
					 continue;
				 }
				 
				 // parse sequence
				 curr_seq = "";
				 while(! (line = br.readLine()).startsWith(">")) {
					curr_seq += line.trim();
				 }
				 sequences.add(curr_seq);
				 uniprot_IDs.add(curr_ID);
			 }
		}
		catch(IOException ioe) {
			System.out.println(ioe.getMessage());
			System.out.println("IOException occurred while parsing sequences.");
		}
		
		
	}
	
	public void parseDomains(String infile) {
		
		ArrayList<String> curr_domains, curr_domain_classes;
		String line;
		boolean[] dbd_found = new boolean[sequences.size()];
		int entry_counter = 0;
		StringTokenizer strtok;
		String curr_ID;
		
		try{
			BufferedReader br = new BufferedReader(new FileReader(new File(infile)));
			
			line = br.readLine();           // skip first line 
			line = br.readLine();
			
			while (line != null && entry_counter < uniprot_IDs.size()) {
				
				if (! line.startsWith("##sequence-region")) { 
					System.out.println("Parse Error. \"##sequence-region\" expected at the beginning of the line.");
				}
				
				strtok = new StringTokenizer(line);
				strtok.nextToken();
				curr_ID = strtok.nextToken().trim();
				    
				// skip entrys corresponding to other organisms
				if (! curr_ID.equals(uniprot_IDs.get(entry_counter))) {
					while (! (line = br.readLine()).startsWith("##sequence-region"));
					continue;
				} 

				curr_domains = new ArrayList<String>();
				curr_domain_classes = new ArrayList<String>();
				
				while (! (line = br.readLine()).startsWith("##sequence-region")) {
					
					if (line.contains("DNA binding")) {
						
						// parse domain 
						curr_domains.add(line.split("\t")[3].trim() + "\t" + line.split("\t")[4].trim());
						
						// parse domain class (if available)
						if (line.contains("Note=")) {
							curr_domain_classes.add(line.split("Note=")[1].trim());
						}
						else {
							curr_domain_classes.add("NA");
						}
						
						dbd_found[entry_counter] = true;
					}
				}
				domains.add(curr_domains);
				domain_classes.add(curr_domain_classes);
				entry_counter++;
			}
			
			int dbd_counter = 0;
			for (int i=0; i<dbd_found.length; i++) {
				if (dbd_found[i]) dbd_counter++;
			}
		}
		
		catch(IOException ioe) {
			System.out.println(ioe.getMessage());
			System.out.println("IOException occurred while parsing DBDs.");
		}
	}
	
	/*
	 *  filter factors with annotated DBD and DBD class
	 */
	
	public void filterDomains() {
		
		boolean[] relevant_entries = new boolean[uniprot_IDs.size()];
		String curr_entry;
		
		// find relevant entries
		for (int i=0; i<uniprot_IDs.size(); i++) {
			if (! domains.get(i).isEmpty()) {
				
				curr_entry = "NA";
				Iterator<String> iter = domain_classes.get(i).iterator();
				while(iter.hasNext() && (curr_entry = iter.next()).equals("NA"));  // skip "NA" entries
				
				if (! curr_entry.equals("NA")) relevant_entries[i] = true;
			}
		}
		
		// remove irrelevant entries
		for (int i=uniprot_IDs.size()-1; i>=0; i--) {
			if (! relevant_entries[i]) {
				uniprot_IDs.remove(i);
				sequences.remove(i);
				domains.remove(i);
				domain_classes.remove(i);
			}
		}
	}
	
	public static void main(String[] args) {
		
		// set paths
		String indir = "/home/jei/SABINE/data/uniprot/";

		
		UniprotParser uniprot_parser = new UniprotParser();
		uniprot_parser.species_name = "Homo sapiens";
		
		uniprot_parser.parseSequences(indir + "uniprot-keyword-805-AND-keyword-DNA-binding-238.fasta");
		uniprot_parser.parseDomains(indir + "uniprot-keyword-805-AND-keyword-DNA-binding-238.gff");
		uniprot_parser.filterDomains();
	}
}

