package util.parsers.uniprobe;
/*
 * ===============================================
 * (C) Florian Topf, University of Tuebingen, 2010
 * ===============================================
 */

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.StringTokenizer;

public class UniProtParser {
	
	ArrayList<String> uniprot_IDs = new ArrayList<String>();
	ArrayList<String> sequences = new ArrayList<String>();
	ArrayList<ArrayList<String>> domains = new ArrayList<ArrayList<String>>();
	
	public void parseSequences(String infile) {
		
		String curr_seq, curr_ID, line;
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
				 
				 // parse sequence
				 curr_seq = "";
				 
				 // while condition changed due to null-pointer exception...
				 while ((line=br.readLine()) != null && !line.startsWith(">")) {
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
		
		ArrayList<String> curr_domains;
		String line;
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
				    
				// skip entries corresponding to other organisms
				if (! curr_ID.equals(uniprot_IDs.get(entry_counter))) {
					while (! (line = br.readLine()).startsWith("##sequence-region"));
					continue;
				} 

				curr_domains = new ArrayList<String>();
				
				// while condition adapted to prevent null-pointer exception...
				while ((line=br.readLine()) != null && !line.startsWith("##sequence-region")) {
					
					if (line.contains("DNA binding")) {
						
						// parse domain 
						curr_domains.add(line.split("\t")[3].trim() + "\t" + line.split("\t")[4].trim());
						
					}
				}
				domains.add(curr_domains);
				entry_counter++;
			}
		}
		
		catch(IOException ioe) {
			System.out.println(ioe.getMessage());
			System.out.println("IOException occurred while parsing DBDs.");
		}
	}
	
	
	
	/*
	 *  filter factors with empty domain
	 */
	public void filterDomains() {
		
		boolean[] relevant_entries = new boolean[uniprot_IDs.size()];
		
		// find relevant entries
		for (int i=0; i<uniprot_IDs.size(); i++) {
			if (! domains.get(i).isEmpty()) {

				relevant_entries[i] = true;
				
			}
		}
		
		// remove irrelevant entries
		for (int i=uniprot_IDs.size()-1; i>=0; i--) {
			if (! relevant_entries[i]) {
				uniprot_IDs.remove(i);
				sequences.remove(i);
				domains.remove(i);
			}
		}
	}
	
	// Startup handler
	public void parse(String fasta, String gff) {
		
		//fasta
		parseSequences(fasta);
		//gff
		parseDomains(gff);
		// filter relevant entries
		filterDomains();
		System.out.println(uniprot_IDs.size()+" sequences with DNA-binding parsed.");
	}
	
	
	// just testing...
	public static void main(String[] args) {
		
		UniProtParser uniprot_parser = new UniProtParser();
		
		// fasta
		uniprot_parser.parseSequences(args[0]);
		
		// gff
		uniprot_parser.parseDomains(args[1]);
		
		uniprot_parser.filterDomains();

		
	}

	public ArrayList<String> getUniprot_IDs() {
		return uniprot_IDs;
	}


	public ArrayList<String> getSequences() {
		return sequences;
	}


	public ArrayList<ArrayList<String>> getDomains() {
		return domains;
	}
}

