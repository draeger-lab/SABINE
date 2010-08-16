package util.parsers.tools.IPRscan;
/*
 * Parser for Input-File ("IP" extended SABINE-FORMAT)
 * Adapted from SABINEInputFileWriter.java
 */

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;

public class InputFileParser {

	ArrayList<String> symbols = new ArrayList<String>();
	ArrayList<String> species = new ArrayList<String>();
	ArrayList<ArrayList<String>> uniprot_ids = new ArrayList<ArrayList<String>>();
	ArrayList<Boolean> flags_arr = new ArrayList<Boolean>();
	ArrayList<String> transfac = new ArrayList<String>();
	
	ArrayList<ArrayList<ArrayList<String>>> matrices = new ArrayList<ArrayList<ArrayList<String>>>();
	ArrayList<ArrayList<String>> sequences = new ArrayList<ArrayList<String>>();
	ArrayList<ArrayList<String>> MN = new ArrayList<ArrayList<String>>();
	ArrayList<ArrayList<String>> bindings = new ArrayList<ArrayList<String>>();

	
	
	
	public void parse(String infile) {
		
		System.out.println("Parsing: "+infile);
		
		try {
			
			BufferedReader br = new BufferedReader(new FileReader(new File(infile)));
			
			String line;
			String curr_name = "";
			String curr_species = "";
			ArrayList<String> curr_ref = null;
			ArrayList<String> curr_seq = null;
			ArrayList<String> curr_domains = null;
			ArrayList<String> curr_pfm_names = null;
			ArrayList<ArrayList<String>> curr_pfms = null;
			
			
			while((line = br.readLine()) != null && (line.length() > 0)) {
				// NA
				if (! line.startsWith("NA")) {
					System.out.println("Parse Error. \"NA\" expected at the beginning of the line.\nLine: " + line);
					System.exit(0);
				} else {
					curr_name = line.substring(4).trim();
				}
				// XX
				br.readLine();
				// SP
				curr_species = (line = br.readLine()).substring(4).trim();
				// XX			
				br.readLine();
				// RF
				curr_ref = new ArrayList<String>();
				while((line = br.readLine()).startsWith("RF")) {
					curr_ref.add(line.substring(4).trim());
				}
				// extension IP
				boolean IP = false;
				if ((line = br.readLine()).startsWith("IP")) {
					IP = true;
					// XX
					br.readLine();
					line = br.readLine();
				}
				// CL
				String curr_class = "";
				if (IP == false && line.startsWith("CL")) {
					curr_class = line.substring(4).trim(); 
					// XX
					br.readLine();
					line = br.readLine();
				}
				// SX
				curr_seq = new ArrayList<String>();
				while (line.substring(0, 2).matches("S[0-9]")) {
					String seq_num = line.substring(1, 2);
					String seq = line.substring(4).trim();
					while ( (line = br.readLine()).startsWith("S"+seq_num) ) {
						seq += line.substring(4).trim();
					}
					curr_seq.add(seq);
					// XX
					line = br.readLine();
				}
				// FT		
				curr_domains = new ArrayList<String>();
				if (IP == false) {
					while (line.startsWith("FT")) {
						curr_domains.add(line.substring(4).trim());
						line = br.readLine();
					}
					line = br.readLine();
				}
				// MN
				curr_pfm_names = new ArrayList<String>();
				// MA
				curr_pfms = new ArrayList<ArrayList<String>>();
				while (line.startsWith("MN")) {
					ArrayList<String> curr_pfm = new ArrayList<String>();
					curr_pfm_names.add(line.substring(4).trim());
					// XX
					br.readLine();
					// MA
					curr_pfm.add(br.readLine().substring(4).trim());
					curr_pfm.add(br.readLine().substring(4).trim());
					curr_pfm.add(br.readLine().substring(4).trim());
					curr_pfm.add(br.readLine().substring(4).trim());
					curr_pfms.add(curr_pfm);

					// XX
					br.readLine();
					line = br.readLine();
				}
				// //
				if (line.startsWith("//")) {
					
					symbols.add(curr_name);
					species.add(curr_species);
					uniprot_ids.add(curr_ref);
					flags_arr.add(IP);
					transfac.add(curr_class);
					sequences.add(curr_seq);
					bindings.add(curr_domains);
					MN.add(curr_pfm_names);
					matrices.add(curr_pfms);
					
					// XX
					br.readLine();
					}
			}
			
		}
		catch(IOException ioe) {
				System.out.println(ioe.getMessage());
				System.out.println("IOException occurred while parsing transcription factors .");
			}
		System.out.println(symbols.size()+" entrie(s) parsed.");
	}

	
	public ArrayList<String> getSymbols() {
		return symbols;
	}

	public ArrayList<String> getSpecies() {
		return species;
	}

	public ArrayList<ArrayList<String>> getUniprot_ids() {
		return uniprot_ids;
	}

	public ArrayList<Boolean> getFlags_arr() {
		return flags_arr;
	}

	public ArrayList<String> getTransfac() {
		return transfac;
	}

	public ArrayList<ArrayList<ArrayList<String>>> getMatrices() {
		return matrices;
	}

	public ArrayList<ArrayList<String>> getSequences() {
		return sequences;
	}

	public ArrayList<ArrayList<String>> getMN() {
		return MN;
	}

	public ArrayList<ArrayList<String>> getBindings() {
		return bindings;
	}
}
