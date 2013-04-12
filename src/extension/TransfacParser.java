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


import help.FormatConverter;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.StringTokenizer;

public class TransfacParser {
	
	ArrayList<String> tf_names = new ArrayList<String>();
	ArrayList<String> species = new ArrayList<String>();
	ArrayList<String> crossrefs = new ArrayList<String>();
	ArrayList<String> classes = new ArrayList<String>();
	ArrayList<String> sequences1 = new ArrayList<String>();
	ArrayList<String> sequences2 = new ArrayList<String>();
	ArrayList<ArrayList<String>> hmm_IDs = new ArrayList<ArrayList<String>>();
	ArrayList<ArrayList<String>> domains = new ArrayList<ArrayList<String>>();
	ArrayList<ArrayList<String>> pfm_names = new ArrayList<ArrayList<String>>();
	ArrayList<ArrayList<String[]>> pfms = new ArrayList<ArrayList<String[]>>();

	ArrayList<ArrayList<String>> pfam_IDs = new ArrayList<ArrayList<String>>();
	ArrayList<ArrayList<String>> smarts_IDs = new ArrayList<ArrayList<String>>();
	ArrayList<ArrayList<String>> prosite_IDs = new ArrayList<ArrayList<String>>();
	
	ArrayList<ArrayList<String>> pfam_domains = new ArrayList<ArrayList<String>>();
	ArrayList<ArrayList<String>> smarts_domains = new ArrayList<ArrayList<String>>();
	ArrayList<ArrayList<String>> prosite_domains = new ArrayList<ArrayList<String>>();
	
	// lists to store information from factor.dat flatfile
	ArrayList<String> factor_ids = new ArrayList<String>();
	ArrayList<String> factor_names = new ArrayList<String>();
	ArrayList<String> factor_synonyms = new ArrayList<String>();
	
	public boolean silent = false;
	
	
	public ArrayList<String> get_tf_names() {
		return tf_names;
	}
	
	public ArrayList<ArrayList<String>> get_pfm_names() {
		return pfm_names;
	}
	
	public ArrayList<ArrayList<String[]>> get_pfms() {
		return pfms;
	}
	
	public ArrayList<String> get_classes() {
		return classes;
	}
	
	public void parseFactors(String infile) {
		
		String line, curr_name, curr_species, curr_ref, curr_class, seq1, seq2;
		StringTokenizer strtok;
		ArrayList<String> curr_domains, curr_pfm_names;
		ArrayList<String[]> curr_pfms;
		String[] splitted_species, curr_pfm;
		
		try {
			 
			 BufferedReader br = new BufferedReader(new FileReader(new File(infile)));
			 
			 while((line = br.readLine()) != null && line.length() > 0) {
				
				if (! line.startsWith("NA")) {
					System.out.println("Parse Error. \"NA\" expected at the beginning of the line.\nLine: " + line);
					System.exit(0);
				}
				 
				strtok = new StringTokenizer(line.substring(2));
				
				curr_name = strtok.nextToken().trim();                       // parse TRANSFAC name
				
				br.readLine();                                               // XX
				
				curr_species = (line = br.readLine()).substring(4).trim();   // parse species
				
				splitted_species = curr_species.split(",");
				
				if(splitted_species.length == 1) {
					curr_species = splitted_species[0];
				}
				else if(splitted_species.length == 2) {
					curr_species = splitted_species[1].trim();
				}
						
				br.readLine();											     // XX
				
				// parse UniProt ID
				if ((line = br.readLine()).startsWith("RF")) {
					curr_ref = line.substring(4).trim();   	 
				
					br.readLine();											     // XX
					line = br.readLine();										 // CL, S1 or S2
					
					if (! line.startsWith("CL") && ! line.startsWith("S")) {
						System.out.println("Parse Error. CL, S1 or S2 expected.");
						System.exit(0);
					}
				}
				else {
					curr_ref = "NA";
				}
				
				// parse superclass
				if (line.startsWith("CL")) {
					curr_class = line.substring(4).trim(); 
					
					br.readLine();												 // XX
					line = br.readLine();									     // S1 or S2
				}
				else {
					curr_class = "NA";
				}
				                                               
				
				// parse sequence(s)
				
				seq1 = null;
				seq2 = null;
				
				if (line.startsWith("S1  ")) {								 // parse first sequence
					seq1 = line.substring(4);
					while ( (line = br.readLine()).startsWith("S1  ") ) {
						seq1 += line.substring(4);
					}
					line = br.readLine();                                    // XX
				}
				
				if (line.startsWith("S2  ")) {								 // parse second sequence
					seq2 = line.substring(4);
					while ( (line = br.readLine()).startsWith("S2  ") ) {
						seq2 += line.substring(4);		    
					}
					line = br.readLine();								     // XX
				}
				
				// parse domains			
				curr_domains = new ArrayList<String>();
				
				while (line.startsWith("FT")) {
					curr_domains.add(line.substring(4).trim());
					line = br.readLine();
				}
				if (curr_domains.size() > 0) {
					line = br.readLine();
				}
				
				// parse PFMs
				curr_pfm_names = new ArrayList<String>();
				curr_pfms = new ArrayList<String[]>();
				
				while(line.startsWith("MN")) {
					
					strtok = new StringTokenizer(line.substring(4));	   // parse matrix name
					curr_pfm_names.add(strtok.nextToken().trim());
					
					line = br.readLine();										   // XX
					
					curr_pfm = new String[4];							   // parse matrix
					curr_pfm[0] = br.readLine().substring(4).trim();
					curr_pfm[1] = br.readLine().substring(4).trim();
					curr_pfm[2] = br.readLine().substring(4).trim();
					curr_pfm[3] = br.readLine().substring(4).trim();
					curr_pfms.add(curr_pfm);
					
					line = br.readLine();                                         // XX
					line = br.readLine();                                         // MN
				}
				
				if (line.startsWith("//")) {
					
					tf_names.add(curr_name);
					species.add(curr_species);
					crossrefs.add(curr_ref);
					classes.add(curr_class);
					sequences1.add(seq1);
					sequences2.add(seq2);
					domains.add(curr_domains);
					pfm_names.add(curr_pfm_names);
					pfms.add(curr_pfms);
					
					br.readLine();									  	  // XX
				}	
			 }
		}
		catch(IOException ioe) {
			System.out.println(ioe.getMessage());
			System.out.println("IOException occurred while parsing transcription factors .");
		}
	}
	
	
	/*
	 *  returns number of TFs and number of TFs with PFM
	 */
	
	public void getNumTFsWithPFM(String infile) {
		
		
		
		try {
			BufferedReader br = new BufferedReader(new FileReader(new File(infile)));
			
			String line;
			int num_TFs = 0;
			int num_TFs_with_PFM = 0;
			
			while ((line = br.readLine()) != null) {
				
				if (line.startsWith("AC  ")) {
					num_TFs++;
						
					while (!(line = br.readLine()).startsWith("//")) {
						if (line.startsWith("MX  ")) num_TFs_with_PFM++; 
					}
				}
			}
			br.close();
			
			System.out.println(num_TFs_with_PFM + " / " + num_TFs + 
							   " (" + num_TFs_with_PFM*100/num_TFs + "%) TFs" +
							   " are attributed to a PFM.");
		}
		catch (IOException e) {
			e.printStackTrace();
		}
		
		
	}
	
	
	/*
	 *  filters parsed factors with given superclass (e.g. 0, 1, 2, 3, 4)
	 */
	
	
	public void parse_IDs_and_names(String infile) {
		
		String line;
		String curr_id, curr_name, curr_synonyms;
		StringTokenizer strtok;
		
		try {
			 
			 BufferedReader br = new BufferedReader(new FileReader(new File(infile)));
			 curr_id = curr_name = curr_synonyms = "";
			 ArrayList<String> curr_alt_ids = new ArrayList<String>();
			 
			 while((line = br.readLine()) != null && ! line.startsWith("//"));
			 
			 while((line = br.readLine()) != null && line.length() > 0) {
				 
				 if (line.startsWith("AC  ")) {
					 curr_id = line.substring(2).trim();
				 }
				 if (line.startsWith("AS  ")) {
					 strtok = new StringTokenizer(line.substring(2));
					 
					 while (strtok.hasMoreTokens()) {
						 curr_alt_ids.add(strtok.nextToken().trim().substring(0,6));
					 }
				 }
				 if (line.startsWith("FA  ")) {
					 curr_name = line.substring(2).trim();
				 }
				 if (line.startsWith("SY  ")) {
					 curr_synonyms = line.substring(2).trim();
				 }
				 if (line.startsWith("//")) {
					 
					 if (!curr_id.isEmpty() && !curr_name.isEmpty()) {
						 factor_ids.add(curr_id);
						 factor_names.add(curr_name);
						 factor_synonyms.add(curr_synonyms);
						 
						 for (int i=0; i<curr_alt_ids.size(); i++) {
							 
							 curr_id = curr_alt_ids.get(i);
							 
							 factor_ids.add(curr_id);
							 factor_names.add(curr_name);
							 factor_synonyms.add(curr_synonyms);
						 }
					 }
					 else {
						 System.out.println("Parse Error.");
						 System.out.println("  ID:       " + curr_id);
						 System.out.println("  Name:     " + curr_name);
						 System.out.println("  Synonyms: " + curr_synonyms);
						 System.exit(0);
					 }
					 curr_id = curr_name = curr_synonyms = "";
					 curr_alt_ids = new ArrayList<String>();
				 }
			 }
			 if (!silent) System.out.println(factor_ids.size() + " factor names parsed.");
		}
		catch(IOException ioe) {
			System.out.println(ioe.getMessage());
			System.out.println("IOException occurred while parsing transcription factors .");
		}
		
		 // manual mapping
		 factor_ids.add("T05986");
		 factor_names.add("ALC");
		 factor_synonyms.add("ALCATRAZ");
		 
		 factor_ids.add("T06748");
		 factor_names.add("KNAT6");
		 factor_synonyms.add("At1g23380");
		 
		 factor_ids.add("T07191");
		 factor_names.add("WRKY33");
		 factor_synonyms.add("At2g38470");
		 
		 factor_ids.add("T07870");
		 factor_names.add("At5g39610");
		 factor_synonyms.add("");
		 
	}
	
	public int[] map_ids_to_names() {
		
		int[] idx_map = new int[tf_names.size()];
		
		for (int i=0; i<tf_names.size(); i++) {
			idx_map[i] = factor_ids.indexOf(tf_names.get(i));
			
			if (idx_map[i] == -1) {
				System.out.println("Mapping error.");
				System.out.println("  Factor: " + tf_names.get(i));
				System.exit(0);
			}
		}
		return idx_map;
	}
	
	
	public void write_IDs_and_names(String outfile, int[] idx_map) {
		try {
			 
			 BufferedWriter bw = new BufferedWriter(new FileWriter(new File(outfile)));
			 
			 // write header
			 bw.write("TRANSFAC ID\tGene name\tSynonyms\n");
			 for (int i=0; i<tf_names.size(); i++) {
				 
				 bw.write(tf_names.get(i) + "\t" + 
						  factor_names.get(idx_map[i]) + "\t" + 
						  factor_synonyms.get(idx_map[i]) + "\n");
			 }
			 bw.flush();
			 bw.close();
		}
		catch(IOException ioe) {
			System.out.println(ioe.getMessage());
			System.out.println("IOException occurred while parsing transcription factors .");
		}
	}
	
	public void filterFactorsBySuperclass(String class_id) {
		
		String[] curr_class;
		String[] all_classes = new String[classes.size()];
		FormatConverter class_identifier = new FormatConverter();
		
		for (int i=0; i < classes.size(); i++) {
			if (classes.get(i).equals("NA"))
				curr_class = new String[] {"0"};
			else
				curr_class = class_identifier.getTransfacClass(classes.get(i)).split("\\p{Punct}");
			
			all_classes[i] = curr_class[0];
			
		}
		
		for (int i=tf_names.size()-1; i>=0; i--) {
			
			if (! class_id.equals(all_classes[i])) {
				
				tf_names.remove(i);
				species.remove(i);
				crossrefs.remove(i);
				classes.remove(i);
				sequences1.remove(i);
				sequences2.remove(i);
				domains.remove(i);
				pfm_names.remove(i);
				pfms.remove(i);
			}
		}
		
		if (!silent)
			System.out.println("  " + tf_names.size() + " factors from TRANSFAC superclass " + class_id + " found.\n");
	}
	
	public void filterFactorsByOrganism(String organism) {
		
		for (int i=species.size()-1; i>=0; i--) {
			
			if (! species.get(i).equals(organism)) {
				
				tf_names.remove(i);
				species.remove(i);
				crossrefs.remove(i);
				classes.remove(i);
				sequences1.remove(i);
				sequences2.remove(i);
				domains.remove(i);
				pfm_names.remove(i);
				pfms.remove(i);
			}
		}
		
		if (!silent)
			System.out.println("  " + tf_names.size() + " factors from " + organism + " found.\n");
	}
	
	public void filterByBlacklist(String blacklistfile) {
		
		BufferedReader br = null;
		
		ArrayList<String> blacklist = new ArrayList<String>();

		try {
			
			String line = null;
			
			StringTokenizer strtok = null;
			
			br = new BufferedReader(new FileReader(new File(blacklistfile)));
			
			while((line = br.readLine()) != null) {
				
				strtok = new StringTokenizer(line);
				
				blacklist.add(strtok.nextToken());
				
			}
			
			boolean[] irrelevant_tfs = new boolean[tf_names.size()];
			
			for (int i=0; i<tf_names.size(); i++) {
				
				if (blacklist.contains(tf_names.get(i))) {
					irrelevant_tfs[i] = true;
				}
			}
			
			for (int i=irrelevant_tfs.length-1; i>=0; i--) {
				
				if (irrelevant_tfs[i]) {
					
					tf_names.remove(i);
					species.remove(i);
					crossrefs.remove(i);
					classes.remove(i);
					sequences1.remove(i);
					sequences2.remove(i);
					domains.remove(i);
					pfm_names.remove(i);
					pfms.remove(i);
				}
			}
		}
		catch(IOException ioe) {
			System.out.println(ioe.getMessage());
			System.out.println("IOException occurred while filtering transcription factors by black list.");
		}
		
		System.out.println("Number of TFs: " + tf_names.size());
	}
	
	public void removeTFs(boolean[] removeFlags) {
		
		for (int i=tf_names.size()-1; i>=0; i--) {
			
			if (removeFlags[i]) {
				tf_names.remove(i);
				species.remove(i);
				crossrefs.remove(i);
				classes.remove(i);
				sequences1.remove(i);
				sequences2.remove(i);
				domains.remove(i);
				pfm_names.remove(i);
				pfms.remove(i);
			}
		}
	}
		
	public void getNumFactorsPerClass() {
		
		FormatConverter converter = new FormatConverter();
		int curr_class;
		int[] class_cnt = new int[] {0, 0, 0, 0, 0};
		
		for (int i=0; i<classes.size(); i++) {
			
			if (classes.get(i).equals("NA")) {
				curr_class = 0;
			}
			else {
				curr_class = Integer.parseInt(converter.getTransfacClass(classes.get(i)).substring(0,1));
			}
			//System.out.println(tf_names.get(i) + ":  " + curr_class);
			class_cnt[curr_class]++;
		}
		System.out.println("Basic Domain Factors:      " + class_cnt[1]);
		System.out.println("Zinc Finger Factors:       " + class_cnt[2]);
		System.out.println("Helix-Turn-Helix Factors:  " + class_cnt[3]);
		System.out.println("Beta Scaffold Factors:     " + class_cnt[4]);
		System.out.println("Other Factors:             " + class_cnt[0]);
		
		System.out.println("Total:                     " + (class_cnt[0] + class_cnt[1] + class_cnt[2] + class_cnt[3] + class_cnt[4]));
	}
	
	public void renameClasses() {
		
		FormatConverter converter = new FormatConverter();
		String curr_class;
		int class_idx;
		
		final String[] superclasses = new String[] { "other", "basic domain", "zinc finger", "helix-turn-helix", "beta scaffold" };
		
		ArrayList<String> class_list = new ArrayList<String>();
		for (int i=0; i<superclasses.length; i++) {
			class_list.add(superclasses[i]);
		}
		
		for (int i=0; i<classes.size(); i++) {
			
			if (classes.get(i).equals("NA")) {
				curr_class = "other";
			}
			else if (class_list.contains(classes.get(i))) {
				curr_class = classes.get(i);
			}
			else {
				class_idx = Integer.parseInt(converter.getTransfacClass(classes.get(i)).substring(0,1));
				curr_class = superclasses[class_idx];
			}
			classes.set(i, curr_class);
			
			System.out.println(new StringTokenizer(tf_names.get(i)).nextToken() + "  " + curr_class);
		}
	}
	
	public void writeTrainingSetFile(String outfile) {
		
		String REGEX1 = "PF\\d{5}";
		String REGEX2 = "[OPQ]\\d[0-9A-Z]{3}\\d(.\\d)?";
		String REGEX3 = "At\\dg\\d{5}";
		
		try {
			BufferedWriter bw = new BufferedWriter(new FileWriter(new File(outfile)));
			StringTokenizer strtok;
			
			String domain_id, domain_seq, domain_entry;
			int domain_start, domain_stop;
			
			String sequence1, sequence2;
			
			for (int i=0; i<tf_names.size(); i++) {
			
				bw.write("NA  " + new StringTokenizer(tf_names.get(i)).nextToken() + "\n");
				bw.write("XX\n");
				
				sequence1 = sequences1.get(i);
				sequence2 = sequences2.get(i);
				
				for (int j=0; j<domains.get(i).size(); j++) {
					
					// parse current domain
					domain_entry = "";
					strtok = new StringTokenizer(domains.get(i).get(j));
					
					domain_id = strtok.nextToken();
					domain_start = Integer.parseInt(strtok.nextToken().replace("<", ""));
					domain_stop =  Integer.parseInt(strtok.nextToken().replace(">", ""));				
					
					// domain entry comes from Pfam ?	
						
						if(domain_id.matches(REGEX1)) {
							if(domain_stop > sequence1.length()) {
								System.out.println("Problems with domain of " + tf_names.get(i) +"\nLength of S1: " + sequence1.length() + "\nEnd-Index of Domain: " + domain_stop);
							}	
							else {
								
								domain_seq = sequence1.substring(domain_start-1,domain_stop);
								
								domain_entry = domain_seq + " 1 " + " " + (domain_start-1) + " " + domain_stop;
							
							}
						}
						
					// domain entry comes from Entrez !
						
						else {
					
					// sequences from Entrez and Transfac were different --> take S2	

							if (sequence2 != null) {
								if(domain_stop > sequence2.length()) {
									System.out.println("Problems with domain of " + tf_names.get(i) +"\nLength of S2: " + sequence2.length() + "\nEnd-Index of Domain: " + domain_stop);
								}	
								else {
									domain_seq = sequence2.substring(domain_start-1,domain_stop);
									domain_entry = domain_seq + " 2 " + " " + (domain_start-1) + " " + domain_stop;
								}
							}
					
					// sequences from Entrez and Transfac were identical --> take S1	
							
							else {
								if(domain_stop > sequence1.length()) {
									System.out.println("Problems with domain of " + tf_names.get(i));
								}	
								else {
									domain_seq = sequence1.substring(domain_start-1,domain_stop);
									domain_entry = domain_seq + " 1 " + (domain_start-1) + " " + domain_stop;
								}
							}
						}
						
						// write domain
						bw.write("DO  " + domain_entry + "\n");
				}
				bw.write("XX\n");
				
				/*
				*  write sequences
				*/
				
				int SEQLINELENGTH = 60;
					
				if(sequence1 != null) {
							
					for(int j=0; j<(sequence1.length()/SEQLINELENGTH); j++) {
						bw.write("S1  "); 
						bw.write(sequence1.toUpperCase(), j*SEQLINELENGTH, SEQLINELENGTH);
						bw.write("\n");
					}
					if(sequence1.length()-(sequence1.length()/SEQLINELENGTH)*SEQLINELENGTH > 0) {
						bw.write("S1  "); 
						bw.write(sequence1.toUpperCase(), (sequence1.length()/SEQLINELENGTH)*SEQLINELENGTH, sequence1.length()-(sequence1.length()/SEQLINELENGTH)*SEQLINELENGTH);
						bw.write("\n");
					}
							
					bw.write("XX\n");
							
				}
	
				if(sequence2 != null) {
							
					for(int j=0; j<(sequence2.length()/SEQLINELENGTH); j++) {
						bw.write("S2  "); 
						bw.write(sequence2.toUpperCase(), j*SEQLINELENGTH, SEQLINELENGTH);
						bw.write("\n");
					}
					if(sequence2.length()-(sequence2.length()/SEQLINELENGTH)*SEQLINELENGTH > 0) {
						bw.write("S2  "); 
						bw.write(sequence2.toUpperCase(), (sequence2.length()/SEQLINELENGTH)*SEQLINELENGTH, sequence2.length()-(sequence2.length()/SEQLINELENGTH)*SEQLINELENGTH);
						bw.write("\n");
					}
							
					bw.write("XX\n");
				}
					
				/*
				 *  parse matrices
				 */
					
				for (int j=0; j<pfms.get(i).size(); j++) {
					bw.write("MA  " + pfms.get(i).get(j)[0] + "\n");
					bw.write("MA  " + pfms.get(i).get(j)[1] + "\n");
					bw.write("MA  " + pfms.get(i).get(j)[2] + "\n");
					bw.write("MA  " + pfms.get(i).get(j)[3] + "\n");
					bw.write("XX\n");
				}
				bw.write("//\n");
				bw.write("XX\n");
			}
			bw.flush();
			bw.close();
		}
		catch(IOException ioe) {
			System.out.println(ioe.getMessage());
			System.out.println("IOException occurred while writing factors to training set file.");
		}
	}
	
	public void parseSequencesAndDomains(String infile) {
		parseSequencesAndDomains(infile, null);
	}
	
	public void parseSequencesAndDomains(String infile, String species_name) {
		
		String line;
		
		String curr_acc, curr_spec, curr_class, curr_seq, curr_ref, curr_dom, curr_ID, id_type;
		curr_acc = curr_spec = curr_seq = curr_dom = curr_ID = id_type = null;
		curr_ref = curr_class = "NA";
		
		String[] split;
		StringTokenizer strtok;
		
		int entry_counter, line_counter;
		entry_counter = line_counter = 0;
		
		ArrayList<String> curr_pfm_IDs = 			new ArrayList<String>();
		
		ArrayList<String> curr_pfam_domains = 		new ArrayList<String>();
		ArrayList<String> curr_smarts_domains = 	new ArrayList<String>();
		ArrayList<String> curr_prosite_domains = 	new ArrayList<String>();
		ArrayList<String> curr_pfam_IDs = 			new ArrayList<String>();
		ArrayList<String> curr_smarts_IDs = 		new ArrayList<String>();
		ArrayList<String> curr_prosite_IDs = 		new ArrayList<String>();
		
		try {
			
			/*
			 *  count factors
			 */
			
			BufferedReader br = new BufferedReader(new FileReader(new File(infile)));
			
			while ((line=br.readLine()) != null) {
				if (line.startsWith("AC  ")) entry_counter++;
			}
			
				
			boolean[] acc_parsed = 		new boolean[entry_counter];
			boolean[] species_parsed = 	new boolean[entry_counter];
			boolean[] class_parsed = 	new boolean[entry_counter];
			boolean[] seq_parsed = 		new boolean[entry_counter];
			boolean[] ref_parsed = 		new boolean[entry_counter];
			boolean[] dom_parsed = 		new boolean[entry_counter];
			boolean[] pfm_parsed = 		new boolean[entry_counter];
			
			
			/*
			 *  parse factors
			 */
			
			entry_counter = 0;

			br = new BufferedReader(new FileReader(new File(infile)));
			
			br.readLine();  // skip first three lines
			br.readLine();
			br.readLine();
			line_counter += 3;
			
			while ((line=br.readLine()) != null) {
				line_counter++;
				

				// parse accession number
				if (line.startsWith("AC")) {
					strtok = new StringTokenizer(line.substring(4));
					curr_acc = strtok.nextToken().trim();
					acc_parsed[entry_counter] = true;
				}
				
				// parse species
				if (line.startsWith("OS")) {
					curr_spec = line.substring(4).trim();
					if ((split = curr_spec.split(",")).length == 2) 
						curr_spec = split[1].trim();
					if (!line.contains("/")) species_parsed[entry_counter] = true;
					
				}
				
				
				// parse TRANSFAC class
				if (line.startsWith("CL")) {
					curr_class = line.substring(4).trim();
					class_parsed[entry_counter] = true;
				}
				
				// parse sequence
				if (line.startsWith("SQ")) {
					
					curr_seq = line.substring(4).trim();
				
					while ((line = br.readLine()).startsWith("SQ")) {
						curr_seq += line.substring(4).trim();
						line_counter++;
					}
					line_counter++;
					seq_parsed[entry_counter] = true;
				}
				
				// parse reference to UniProt
				if (line.startsWith("SC")) {
					strtok = new StringTokenizer(line.substring(4).replace("#", " "));
					while (strtok.hasMoreTokens()) {
						if (strtok.nextToken().trim().equals("Swiss-Prot")) {
							if (strtok.hasMoreTokens()) {
								curr_ref = strtok.nextToken().trim();
								ref_parsed[entry_counter] = true;
							}
							break;
						}
					}
				}
				
				if (!ref_parsed[entry_counter] && line.startsWith("DR  SWISSPROT")) {
					
					curr_ref = line.replace("DR  SWISSPROT: ", "").substring(0,6);
					ref_parsed[entry_counter] = true;
				}
				
				// parse domains
				if (line.startsWith("FT")) {
					
					curr_pfam_domains = new ArrayList<String>();
					curr_smarts_domains = new ArrayList<String>();
					curr_prosite_domains = new ArrayList<String>();
					
					curr_pfam_IDs = new ArrayList<String>();
					curr_smarts_IDs = new ArrayList<String>();
					curr_prosite_IDs = new ArrayList<String>();
					
					while (line.startsWith("FT")) {
						id_type = "";
						if (line.substring(4).matches(".*PF\\d{5}\\p{Punct}.*")) id_type = "PF";
						if (line.substring(4).matches(".*SM\\d{5}\\p{Punct}.*")) id_type = "SM";
						if (line.substring(4).matches(".*PS\\d{5}\\p{Punct}.*")) id_type = "PS";
						
						if (id_type.length() == 2) {
							strtok = new StringTokenizer(line.substring(4));
							
							curr_dom = strtok.nextToken().trim() + "\t" + strtok.nextToken().trim();
							curr_ID = strtok.nextToken().trim();
							curr_ID = curr_ID.substring(0, curr_ID.length()-1);
							dom_parsed[entry_counter] = true;
							
							if (id_type.equals("PF")) {
								curr_pfam_domains.add(curr_dom);
								curr_pfam_IDs.add(curr_ID);
							}
							else if (id_type.equals("SM")) {
								curr_smarts_domains.add(curr_dom);
								curr_smarts_IDs.add(curr_ID);
							}
							else if (id_type.equals("PS")) {
								curr_prosite_domains.add(curr_dom);
								curr_prosite_IDs.add(curr_ID);
							}
						}
						line = br.readLine();
						line_counter++;
					}
				}
					
				if (line.startsWith("MX")) {
					pfm_parsed[entry_counter] = true;
					curr_pfm_IDs = new ArrayList<String>();
					while (line.startsWith("MX")) {
						curr_pfm_IDs.add(line.substring(4, 10));
						line = br.readLine();
					}
				}
				
				// check parsed information
				if (line.startsWith("//")) {
					
					if (species_name != null && !curr_spec.equals(species_name)) {
						continue;
					}
					
					if (acc_parsed[entry_counter] 
					    && species_parsed[entry_counter]
					    && seq_parsed[entry_counter]) {
						
						tf_names.add(curr_acc);
						species.add(curr_spec);
						crossrefs.add(curr_ref);
						classes.add(curr_class);
						sequences1.add(curr_seq);
						sequences2.add(null);
						
						// add domain (if available)
						if (dom_parsed[entry_counter]) {
							pfam_IDs.add(curr_pfam_IDs);
							smarts_IDs.add(curr_smarts_IDs);
							prosite_IDs.add(curr_prosite_IDs);
							pfam_domains.add(curr_pfam_domains);
							smarts_domains.add(curr_smarts_domains);
							prosite_domains.add(curr_smarts_domains);
							
						} else {
							pfam_IDs.add(null);
							smarts_IDs.add(null);
							prosite_IDs.add(null);
							pfam_domains.add(null);
							smarts_domains.add(null);
							prosite_domains.add(null);
						}

						// add PFM (if available)
						if (pfm_parsed[entry_counter]) {
							pfm_names.add(curr_pfm_IDs);
						} else {
							pfm_names.add(null);
						}
					}
					entry_counter++;
					curr_ref = "NA";
					curr_class = "NA";
				}
			}
			br.close();
			
			/*
			 *  print statistics
			 */
			
			if (!silent) {
			
				int acc_counter, spec_counter, class_counter, seq_counter, ref_counter, dom_counter, pfm_counter;
				acc_counter = spec_counter = class_counter = seq_counter = ref_counter = dom_counter = pfm_counter = 0;
			
				for(int i=0; i<entry_counter; i++) if (acc_parsed[i]) acc_counter++;
				for(int i=0; i<entry_counter; i++) if (species_parsed[i]) spec_counter++;
				for(int i=0; i<entry_counter; i++) if (class_parsed[i]) class_counter++;
				for(int i=0; i<entry_counter; i++) if (seq_parsed[i]) seq_counter++;
				for(int i=0; i<entry_counter; i++) if (ref_parsed[i]) ref_counter++;
				for(int i=0; i<entry_counter; i++) if (dom_parsed[i]) dom_counter++;
				for(int i=0; i<entry_counter; i++) if (pfm_parsed[i]) pfm_counter++;
			
				System.out.println("Number of transcription factors: " + entry_counter);
				System.out.println("Number of parsed factors:        " + tf_names.size() + "\n");
			
				System.out.println("Accession numbers:  " + acc_counter);
				System.out.println("Species:            " + spec_counter);
				System.out.println("TRANSFAC class:     " + class_counter);
				System.out.println("Sequences:          " + seq_counter);
				System.out.println("UniProt References: " + ref_counter);
				System.out.println("Domains:            " + dom_counter);
				System.out.println("Matrices:           " + pfm_counter);
				
			}
		}
		catch(IOException ioe) {
			System.out.println(ioe.getMessage());
			System.out.println("IOException occurred while parsing transcription factors from FACTOR table.");
		}
	}
		

	
	
	
	public void filterDomains(String interpro_infile, String mapping_infile) {

		/*
		 *  obtain HMM-IDs of DNA binding domains
		 */
		
		InterproParser interpro_parser = new InterproParser();

		interpro_parser.parseMapping(interpro_infile);
		
		interpro_parser.filterBindingDomains(mapping_infile);
		
		ArrayList<String> pfam_DBDs = new ArrayList<String>(); 
		for (int i=0; i<interpro_parser.pfam_IDs.size(); i++) {
			for (int j=0; j<interpro_parser.pfam_IDs.get(i).size(); j++) {
				pfam_DBDs.add(interpro_parser.pfam_IDs.get(i).get(j));
			}
		}
		
		ArrayList<String> prosite_DBDs = new ArrayList<String>(); 
		for (int i=0; i<interpro_parser.prosite_IDs.size(); i++) {
			for (int j=0; j<interpro_parser.prosite_IDs.get(i).size(); j++) {
				prosite_DBDs.add(interpro_parser.prosite_IDs.get(i).get(j));
			}
		}
		
		ArrayList<String> smarts_DBDs = new ArrayList<String>(); 
		for (int i=0; i<interpro_parser.smarts_IDs.size(); i++) {
			for (int j=0; j<interpro_parser.smarts_IDs.get(i).size(); j++) {
				smarts_DBDs.add(interpro_parser.smarts_IDs.get(i).get(j));
			}
		}
		
		/*
		 *  filter DNA binding domains
		 */
		
		// initialize HMM-IDs and domains
		for (int i=0; i<tf_names.size(); i++) {
			hmm_IDs.add(new ArrayList<String>());
			domains.add(new ArrayList<String>());
		}
		
		ArrayList<String> curr_domains, curr_IDs;
		
		for (int i=0; i<pfam_IDs.size(); i++) {
			
			curr_domains = domains.get(i);
			curr_IDs = hmm_IDs.get(i);
			
			for (int j=0; j<pfam_IDs.get(i).size(); j++) {
				
				if (pfam_DBDs.contains(pfam_IDs.get(i).get(j))) {
					
					curr_IDs.add(pfam_IDs.get(i).get(j));
					curr_domains.add(pfam_domains.get(i).get(j));
				}
			}
			domains.set(i, curr_domains);
			hmm_IDs.set(i, curr_IDs);
		}
		
		for (int i=0; i<prosite_IDs.size(); i++) {
			
			curr_domains = domains.get(i);
			curr_IDs = hmm_IDs.get(i);
			
			for (int j=0; j<prosite_IDs.get(i).size(); j++) {
				
				if (prosite_DBDs.contains(prosite_IDs.get(i).get(j))) {
					
					curr_IDs.add(prosite_IDs.get(i).get(j));
					curr_domains.add(prosite_domains.get(i).get(j));
				}
			}
			domains.set(i, curr_domains);
			hmm_IDs.set(i, curr_IDs);
		}
		
		for (int i=0; i<smarts_IDs.size(); i++) {
			
			curr_domains = domains.get(i);
			curr_IDs = hmm_IDs.get(i);
			
			for (int j=0; j<smarts_IDs.get(i).size(); j++) {
				
				if (smarts_DBDs.contains(smarts_IDs.get(i).get(j))) {
					
					curr_IDs.add(smarts_IDs.get(i).get(j));
					curr_domains.add(smarts_domains.get(i).get(j));
				}
			}
			domains.set(i, curr_domains);
			hmm_IDs.set(i, curr_IDs);
		}
	}
	
	public int[] countDomains() {
		
		int pfam_counter, prosite_counter, smarts_counter;
		pfam_counter = prosite_counter = smarts_counter = 0;
		
		for (int i=0; i<pfam_IDs.size(); i++) {
			for (int j=0; j<pfam_IDs.get(i).size(); j++) {
				pfam_counter++;
			}
		}
	
		for (int i=0; i<prosite_IDs.size(); i++) {
			for (int j=0; j<prosite_IDs.get(i).size(); j++) {
				prosite_counter++;
			}
		}
		
		for (int i=0; i<smarts_IDs.size(); i++) {
			for (int j=0; j<smarts_IDs.get(i).size(); j++) {
				smarts_counter++;;
			}
		}
		
		return new int[] {pfam_counter , prosite_counter , smarts_counter};
	}
	
	public void writeFactorsToFile(String outfile) {
		writeFactorsToFile(outfile, -1);
	}

	public void writeFactorsToFile(String outfile, int class_id) {
		
		int SEQLINELENGTH = 60;
		String curr_seq;
		FormatConverter converter = new FormatConverter();
	
		try {
			BufferedWriter bw = new BufferedWriter(new FileWriter(new File(outfile)));
	
			for (int i=0; i<tf_names.size(); i++) {
	
				// only write TFs of certain superclass
				if (class_id != -1) {
					if (classes.get(i).equals("NA")) {
						continue;
					}
					int currClass = Integer.parseInt(converter.getTransfacClass(classes.get(i)).substring(0,1));
					if (currClass != class_id) {
						continue;
					}
				}
				
				bw.write("NA  " + tf_names.get(i) + "\n" +
						"XX\n" + 
						"SP  " + species.get(i) + "\n" + 
						"XX\n" + 
						"RF  " + crossrefs.get(i) + "\n" + 
						"XX\n" + 
						"CL  " + classes.get(i) + "\n" + 
						"XX\n");
	
				curr_seq = sequences1.get(i);
	
				if (curr_seq != null) { 		
					for(int j=0; j<(curr_seq.length()/SEQLINELENGTH); j++) {
	
						bw.write("S1  "); 
						bw.write(curr_seq.toUpperCase(), j*SEQLINELENGTH, SEQLINELENGTH);
						bw.write("\n");	
					}
	
					if(curr_seq.length()-(curr_seq.length()/SEQLINELENGTH)*SEQLINELENGTH > 0) {						
	
						bw.write("S1  "); 
						bw.write(curr_seq.toUpperCase(), (curr_seq.length()/SEQLINELENGTH)*SEQLINELENGTH, curr_seq.length()-(curr_seq.length()/SEQLINELENGTH)*SEQLINELENGTH);
						bw.write("\n");	
					}
					bw.write("XX\n");
				}
	
				curr_seq = sequences2.get(i);
	
				if (curr_seq != null) { 
					for(int j=0; j<(curr_seq.length()/SEQLINELENGTH); j++) {
	
						bw.write("S2  "); 
						bw.write(curr_seq.toUpperCase(), j*SEQLINELENGTH, SEQLINELENGTH);
						bw.write("\n");	
					}
	
					if(curr_seq.length()-(curr_seq.length()/SEQLINELENGTH)*SEQLINELENGTH > 0) {
	
						bw.write("S2  "); 
						bw.write(curr_seq.toUpperCase(), (curr_seq.length()/SEQLINELENGTH)*SEQLINELENGTH, curr_seq.length()-(curr_seq.length()/SEQLINELENGTH)*SEQLINELENGTH);
						bw.write("\n");	
					}
					bw.write("XX\n");
				}
	
				// write domains 
				if (domains.size() > 0) {
					for (int j=0; j<domains.get(i).size(); j++) {
						bw.write("FT  " + domains.get(i).get(j) + "\n");
					}
					bw.write("XX\n");
				}
	
	
				// write PFMs
				ArrayList<String> curr_pfm_names = pfm_names.get(i);
				if (curr_pfm_names != null) {
					for (int j=0; j<curr_pfm_names.size(); j++) {
						bw.write("MN  " + pfm_names.get(i).get(j) + "\n");
						bw.write("XX\n");
	
						bw.write("MA  " + pfms.get(i).get(j)[0] + "\n");
						bw.write("MA  " + pfms.get(i).get(j)[1] + "\n");
						bw.write("MA  " + pfms.get(i).get(j)[2] + "\n");
						bw.write("MA  " + pfms.get(i).get(j)[3] + "\n");
						bw.write("XX\n");
					}
				}
				bw.write("//\n" +
						"XX\n");
			}
	
			bw.flush();
			bw.close();
		}
	
		catch(IOException ioe) {
			System.out.println(ioe.getMessage());
			System.out.println("IOException occurred while parsing transcription factors .");
		}
	}
	
	public void writeLabelFile(String outfile) {
		
		try {
			BufferedWriter bw = new BufferedWriter(new FileWriter(new File(outfile)));
			PFMFormatConverter converter = new PFMFormatConverter();
			ArrayList<String> curr_matrix = new ArrayList<String>();
			
			for (int i=0; i<tf_names.size(); i++) {
				
				
				bw.write("NA  " + tf_names.get(i) + "\n" +
						 "XX\n" + 
						 "SP  " + species.get(i) + "\n" + 
						 "XX\n" + 
						 "RF  " + crossrefs.get(i) + "\n" + 
						 "XX\n" + 
						 "CL  " + classes.get(i) + "\n" + 
						 "XX\n");
				
				
				// write PFMs
				for (int j=0; j<pfms.get(i).size(); j++) {
					bw.write("MN  " + pfm_names.get(i).get(j) + "\n");
					bw.write("XX\n");
					
					curr_matrix = converter.convertTransfacToSTAMP(pfms.get(i).get(j));
					
					for (int k=0; k<curr_matrix.size(); k++) {
						bw.write("MA  " + curr_matrix.get(k) + "\n");
					}
					bw.write("XX\n");
				}
				bw.write("//\n" +
						 "XX\n");
			}
			
			bw.flush();
			bw.close();
		}
		
		catch(IOException ioe) {
			System.out.println(ioe.getMessage());
			System.out.println("IOException occurred while parsing transcription factors .");
		}
	}
	
	public static void main(String[] args) {
		
		//TransfacParser tf_parser = new TransfacParser();
		//tf_parser.parseFactors("data/transfac/transpall_interleaved_classes.out");
		//tf_parser.writeFactorsToFile("data/transfac/my_file");
		
		//tf_parser.parseFactors("trainingsets_eucaryotes/transpall_interleaved_classes.out");
		//tf_parser.filterByBlacklist("/home/eichner/Desktop/black_list.out");
		//tf_parser.writeFactorsToFile("/home/eichner/Desktop/basisdatensatz1.txt");
		
		//tf_parser.parseFactors("/home/eichner/Desktop/uniprot_transfac_predictions.txt");
		//tf_parser.getNumFactorsPerClass();
		//tf_parser.renameClasses();
		//tf_parser.writeFactorsToFile("/home/eichner/Desktop/uniprot_transfac_predictions_new.txt");
		
		
		//tf_parser.parseSequencesAndDomains("/data/dat0/transfac/2008.3/dat/factor.dat", "Homo sapiens");
		//tf_parser.filterDomains("/home/jei/SABINE/data/interpro/interpro.xml",
		//						"/home/jei/SABINE/mappings/TFclassMapping.txt");
		
		/*
		tf_parser.silent = false;
		tf_parser.parseFactors("/rahome/eichner/projects/sabine/supplement/basisdatensatz2.txt");
		tf_parser.filterFactorsByOrganism("Arabidopsis thaliana");
		tf_parser.parse_IDs_and_names("/rahome/eichner/data/transfac_professional/cgi-bin/data/factor.dat");
		int[] idx_map = tf_parser.map_ids_to_names();
		tf_parser.write_IDs_and_names("/rahome/eichner/Desktop/A_thaliana_TFs.txt", idx_map);
		*/
		
		// tf_parser.getNumTFsWithPFM("/rahome/eichner/data/transfac_2010.1/dat/factor.dat");
		
		/*
		TransfacParser newDataParser = new TransfacParser();
		newDataParser.parseSequencesAndDomains("/rahome/eichner/data/biobase/transfac_2011.4/dat/factor.dat");
		newDataParser.parseMatrices("/rahome/eichner/data/biobase/transfac_2011.4/dat/matrix.dat");
		newDataParser.writeFactorsToFile("/rahome/eichner/projects/tfpredict/data/tf_pred/sabine_files/transfac_2011.4_flatfile.txt");
		*/
		
		TransfacParser trainSetParser = new TransfacParser();
		trainSetParser.parseFactors("/rahome/eichner/projects/sabine/data/trainingsets/latest/trainingset_private.txt");
		trainSetParser.writeLabelFile("/rahome/eichner/projects/sabine/data/trainingsets/latest/labelfile_private.txt");
	}
}

