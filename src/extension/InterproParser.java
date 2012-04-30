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
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.StringTokenizer;

public class InterproParser {

	ArrayList<String> uniprot_IDs = new ArrayList<String>();
	ArrayList<String> sequences = new ArrayList<String>();
	ArrayList<ArrayList<String>> hmm_IDs = new ArrayList<ArrayList<String>>();
	ArrayList<ArrayList<String>> domains = new ArrayList<ArrayList<String>>();
	
	ArrayList<ArrayList<String>> pfam_IDs = new ArrayList<ArrayList<String>> ();
	ArrayList<ArrayList<String>> prosite_IDs = new ArrayList<ArrayList<String>> ();
	ArrayList<ArrayList<String>> smarts_IDs = new ArrayList<ArrayList<String>> ();
	boolean[] dna_binding;
	
	
	public InterproParser() {}
	
	public InterproParser(ArrayList<String> IDs, ArrayList<String> seqs) {
		uniprot_IDs = IDs;
		sequences = seqs;
	}
	
	
	public ArrayList<String> readPfamIDs(String infile) {
		
		String line, curr_ID;
		StringTokenizer strtok;
		
		ArrayList<String> pfam_ids = new ArrayList<String>();
		
		try {
			 
			BufferedReader br = new BufferedReader(new FileReader(new File(infile)));
			
			while ((line = br.readLine()) != null) {
				
				strtok = new StringTokenizer(line, "\t");
				
				strtok.nextToken();   
				strtok.nextToken();
				curr_ID = strtok.nextToken().trim();
				
				if (curr_ID.matches("PF\\d{5}")) {   // PFAM identifier ?
					
					pfam_ids.add(curr_ID);
				}
			}
		}
		
		catch(IOException ioe) {
			System.out.println(ioe.getMessage());
			System.out.println("IOException occurred while reading PFAM IDs.");
		}
		
		return pfam_ids;
	}

	
	public void parseMapping(String infile) {
		
		String line, curr_pfam_ID, curr_smarts_ID, curr_prosite_ID;
		
		ArrayList<String> curr_pfam_IDs = new ArrayList<String>();
		ArrayList<String> curr_smarts_IDs = new ArrayList<String>();
		ArrayList<String> curr_prosite_IDs = new ArrayList<String>();
		
		try {
			BufferedReader br = new BufferedReader(new FileReader(new File(infile)));
		
			while ((line = br.readLine()) != null) {
				
				if (line.contains(" db=\"PFAM\" ") && line.trim().startsWith("<db_xref ")) {
					
					curr_pfam_IDs = new ArrayList<String>();
					curr_smarts_IDs = new ArrayList<String>();
					curr_prosite_IDs = new ArrayList<String>();
					
					curr_pfam_ID = line.trim().split(" ")[3].substring(7,14);
					curr_pfam_IDs.add(curr_pfam_ID);
					
					while (! (line = br.readLine()).contains("</member_list>")) {
						
						if (line.contains(" db=\"PFAM\" ")){
							curr_pfam_ID = line.trim().split(" ")[3].substring(7,14);
							curr_pfam_IDs.add(curr_pfam_ID);
						}
						else if (line.contains(" db=\"PROSITE\" ")) {
							curr_prosite_ID = line.trim().split(" ")[3].substring(7,14);
							curr_prosite_IDs.add(curr_prosite_ID);
						}
						else if (line.contains(" db=\"SMART\" ")) {
							curr_smarts_ID = line.trim().split(" ")[3].substring(7,14);
							curr_smarts_IDs.add(curr_smarts_ID);
						}
					}
					pfam_IDs.add(curr_pfam_IDs);
					prosite_IDs.add(curr_prosite_IDs);
					smarts_IDs.add(curr_smarts_IDs);
				}
			}
		}
			
		catch(IOException ioe) {
			System.out.println(ioe.getMessage());
			System.out.println("IOException occurred while reading PFAM IDs.");
		}
	}
	
	public void filterBindingDomains(String infile) {
		
		dna_binding = new boolean[pfam_IDs.size()];
		
		ArrayList<String> known_DBD_IDs = readPfamIDs(infile);

		/*
		 *  identify HMM-IDs of DNA binding domains
		 */
		
			
		for (int i=0; i<pfam_IDs.size(); i++) {
			for (int j=0; j<pfam_IDs.get(i).size(); j++) {
				
				if (known_DBD_IDs.contains(pfam_IDs.get(i).get(j)))
					dna_binding[i] = true;
			}	
		}
	
		/*
		 *  filter DNA binding domains
		 */
		
		for (int i=dna_binding.length-1; i>=0; i--) {
			
			if (!dna_binding[i]) {
				
				pfam_IDs.remove(i);
				prosite_IDs.remove(i);
				smarts_IDs.remove(i);
			}
		}
	}
	
	
	
	public void parseDomains(String infile, String outfile) {
		
		String line, curr_uniprot_ID, curr_pfam_ID, domain_start, domain_end;
		StringTokenizer strtok;
		ArrayList<String> curr_domains, curr_pfam_IDs;
		int line_counter = 0;
		int entry_counter = 0;
		//String blank = " ";
		
		// initialize HMM IDs and domains
		for (int i=0; i<uniprot_IDs.size(); i++) {
			hmm_IDs.add(new ArrayList<String>());
			domains.add(new ArrayList<String>());
		}
		
		try {
			 
			BufferedReader br = new BufferedReader(new FileReader(new File(infile)));
			BufferedWriter bw = new BufferedWriter(new FileWriter(new File(outfile)));
			line = br.readLine();
			int idx;
			
			while (line != null) {
			
				// go to next protein ID
				while ((line = br.readLine()) != null && ! line.contains("<protein id=")) line_counter++;
				line_counter++;
				
				if (line == null) return;
				
				// parse UniProt ID
				strtok = new StringTokenizer(line, "\"");
				strtok.nextToken();
				curr_uniprot_ID = strtok.nextToken().trim();
			
				// relevant protein ?
				if (! uniprot_IDs.contains(curr_uniprot_ID)) {    	// irrelevant -> skip entry
					continue;
				}
				else {									  			// relevant -> parse domain assignment
					
					curr_domains = new ArrayList<String>();
					curr_pfam_IDs = new ArrayList<String>();
					
					while (true) {

						// search PFAM ID
						while (! (line = br.readLine()).matches(".*PF\\d{5}.*") &&
								! line.contains("</protein>")) line_counter++;
						line_counter++;
					
						if (line.contains("</protein>")) break;   	// no PFAM ID found -> skip entry
					
						// parse PFAM ID
						strtok = new StringTokenizer(line, "\"");
						strtok.nextToken();
						curr_pfam_ID = strtok.nextToken().trim();
					
						if (! line.contains("dbname=\"PFAM\"")) {
							System.out.println("Parse Error.");
							System.exit(0);
						}
					
						// search domain assignment 
						while (! (line=br.readLine()).contains("<lcn start=") &&
								! line.contains("</protein>")) line_counter++;
						line_counter++;
					
						if (line.contains("</protein>")) break;   	// no domain found -> skip entry
					
						// parse domain assignment
						strtok = new StringTokenizer(line, "\"");
						strtok.nextToken();  						
						domain_start = strtok.nextToken().trim();
						strtok.nextToken();  						
						domain_end = strtok.nextToken().trim();
						
						curr_pfam_IDs.add(curr_pfam_ID);
						curr_domains.add(domain_start + "\t" + domain_end);
					}
					idx = uniprot_IDs.indexOf(curr_uniprot_ID);
					
					hmm_IDs.set(idx, curr_pfam_IDs);
					domains.set(idx, curr_domains);
					
					entry_counter++;
					/*
					if (entry_counter % 1000 == 0){
						blank = "";
					}
					if (entry_counter % 100 == 0 || entry_counter == uniprot_IDs.size()){
						System.out.println("    " + blank + entry_counter + " proteins parsed" );
					}
					*/
					
					bw.write("Uniprot ID " + entry_counter + ": " + 
							 uniprot_IDs.get(idx) + " (line: " + line_counter + ")\n");
					for (int i=0; i<curr_domains.size(); i++) { 
						bw.write("  DBD " + (i+1) + ": " + curr_pfam_IDs.get(i) + 
								 "\t" + curr_domains.get(i) + "\n");
					}
					bw.write("\n");
					
				} 
			}
			br.close();
			bw.flush();
			bw.close();
		}
		catch(IOException ioe) {
			System.out.println(ioe.getMessage());
			System.out.println("IOException occurred while parsing DBDs.");
		}
	}
	
	public static void main(String[] args) {
		
		/*
		// set paths
		String indir_uniprot = "/home/jei/SABINE/data/uniprot/";
		String indir_interpro = "/home/jei/SABINE/data/interpro/";
	
		UniprotParser uniprot_parser = new UniprotParser();
		uniprot_parser.parseSequences(indir_uniprot + "uniprot-keyword-805-AND-keyword-DNA-binding-238.fasta");
		
		InterproParser interpro_parser = new InterproParser(uniprot_parser.uniprot_IDs, 
															uniprot_parser.sequences);
		
		interpro_parser.parseDomains(indir_interpro + "match_complete.xml", "/home/jei/SABINE/tmp/parsed_domains.txt");
		*/
		
		InterproParser interpro_parser = new InterproParser();
		
		interpro_parser.parseMapping("home/jei/SABINE/data/interpro/interpro.xml");
		interpro_parser.filterBindingDomains("/home/jei/SABINE/mappings/TFclassMapping.txt");
	
	}

}

