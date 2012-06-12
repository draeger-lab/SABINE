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

package help;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.StringTokenizer;

import main.FBPPredictor;

import org.biojava.bio.BioException;

import help.SequenceAligner;


public class FormatConverter {
	
	public static final String NonTFclassID = "Non-TF";
	
	public String[] convertToInternalFormat(String infile, String outfile) {
		
		
		System.out.println("\nConverting input file to internal format.\n");
		
		
		BufferedReader br = null;
		
		String line = null;
		
		StringTokenizer strtok = null;
		
		
		ArrayList<String> known_class_ids = new ArrayList<String>();
		
		known_class_ids.add("0.");
		known_class_ids.add("1.");
		known_class_ids.add("2.");
		known_class_ids.add("3.");
		known_class_ids.add("4.");
		
		
		/*
		 * 
		 * parse list of possible species
		 * 
		 */
		
		ArrayList<String> known_species = new ArrayList<String>();
		
		try {
			
			br = new BufferedReader(new FileReader(new File("input" + File.separator + "species.dat")));
			
			while((line = br.readLine()) != null) {
				
				known_species.add(line.trim());
				
			}
			
			br.close();
			
		}
		catch(IOException ioe) {
			System.out.println(ioe.getMessage());
			System.out.println("IOException occurred while parsing list of possible species.");
		}
		
		
		
		String name = null;
		
		String class_id = null;
		
		String species = null;
		
		String seq1 = null;
		String seq2 = null;
		
		ArrayList<int[]> domains1 = null;
		ArrayList<int[]> domains2 = null;
		
		/*
		 * 
		 * parse input file
		 * 
		 */
		
		try {
			
			
			br = new BufferedReader(new FileReader(new File(infile)));
			
			
		// parse TF name
			
			name = br.readLine().trim();
			
			if((line = br.readLine()) == null || line.length() > 0) {	
				System.out.println("Parse Error. Blank line expected after TF name. Aborting.");
				System.out.println("Line : " + line);
				System.exit(0);
			}
			
			
		// parse class
			
			class_id = br.readLine().trim();
			
			if(class_id == null || !known_class_ids.contains(class_id)) {
				System.out.println("Parse Error. Unknown TF class. Possible values: 0. , 1. , 2. , 3. , 4. . Aborting.");
				System.out.println("Class : " + line);
				System.exit(0);
			}
			
			if((line = br.readLine()) == null || line.length() > 0) {	
				System.out.println("Parse Error. Blank line expected after TF class. Aborting.");
				System.out.println("Line : " + line);
				System.exit(0);
			}
			
			
		// parse species
			
			species = br.readLine().trim();
			
			if(species == null || !known_species.contains(species.toUpperCase())) {
				System.out.println("Parse Error. Unknown species. List of possible species can be found in \"input/species.dat\". Aborting.");
				System.out.println("Species : " + line);
				System.exit(0);
			}
			
			if((line = br.readLine()) == null || line.length() > 0) {	
				System.out.println("Parse Error. Blank line expected after species. Aborting.");
				System.out.println("Line : " + line);
				System.exit(0);
			}
			
			
		// parse sequence
			
			seq1 = "";
			
			while((line = br.readLine()) != null && line.length() > 0) {
				
				seq1 += line.trim();
				
			}	
			
			if(seq1.length() == 0) {	
				System.out.println("Parse Error. Unable to parse protein sequence. Aborting.");
				System.exit(0);
			}
			
			
		// parse domains of sequence
			
			domains1 = new ArrayList<int[]>();
			
			
			while((line = br.readLine()) != null && line.length() > 0) {
				
				
				
				strtok = new StringTokenizer(line);
				
				int[] startend = new int[2];
				
				
				
				if(strtok.hasMoreTokens()) {
					
					String token = strtok.nextToken();
					
					if(token.matches("\\d*")) {
						
						startend[0] = Integer.parseInt(token);
						
						if(startend[0] < 1) {
							
							System.out.println("Parse Error. Problems while parsing start index of domain. Index is smaller than 1. Aborting.");
							System.out.println("Index : " + token);
							System.exit(0);
							
						}
						
					}
					
					else {
						System.out.println("Parse Error. Problems while parsing start index of domain. Aborting.");
						System.out.println("Index : " + token);
						System.exit(0);
					}
					
				}
				else {
					System.out.println("Parse Error. Problems while parsing start index of domain. Aborting.");
					System.out.println("Line : " + line);
					System.exit(0);
				}
				
				
				
				if(strtok.hasMoreTokens()) {
					
					String token = strtok.nextToken();
					
					if(token.matches("\\d*")) {
						
						startend[1] = Integer.parseInt(token);
						
						if(startend[1] > seq1.length()) {
							
							System.out.println("Parse Error. Problems while parsing end index of domain. Index is greater than " + seq1.length() + ". Aborting.");
							System.out.println("Index : " + token);
							System.exit(0);
							
						}
						
					}
					
					else {
						System.out.println("Parse Error. Problems while parsing end index of domain. Aborting.");
						System.out.println("Index : " + token);
						System.exit(0);
					}
					
				}
				else {
					System.out.println("Parse Error. Problems while parsing end index of domain. Aborting.");
					System.out.println("Line : " + line);
					System.exit(0);
				}
				
				
				domains1.add(startend);
				
			}	
			
			if(domains1.size() == 0) {	
				System.out.println("Parse Error. Unable to parse DNA binding domains. Aborting.");
				System.exit(0);
			}
			
			
		// parse alternative sequence, if available
			
			seq2 = "";
			
			while((line = br.readLine()) != null && line.length() > 0) {
				
				seq2 += line.trim();
				
			}	
			
			
		// parse domains of alternative sequence, if available
			
			domains2 = new ArrayList<int[]>();
			
			
			while((line = br.readLine()) != null && line.length() > 0) {
				
				
				strtok = new StringTokenizer(line);
				
				int[] startend = new int[2];
				
				
				
				if(strtok.hasMoreTokens()) {
					
					String token = strtok.nextToken();
					
					if(token.matches("\\d*")) {
						
						startend[0] = Integer.parseInt(token);
						
						if(startend[0] < 1) {
							
							System.out.println("Parse Error. Problems while parsing start index of alternative domain. Index is smaller than 1. Aborting.");
							System.out.println("Index : " + token);
							System.exit(0);
							
						}
						
					}
					
					else {
						System.out.println("Parse Error. Problems while parsing start index of alternative domain. Aborting.");
						System.out.println("Index : " + token);
						System.exit(0);
					}
					
				}
				else {
					System.out.println("Parse Error. Problems while parsing start index of alternative domain. Aborting.");
					System.out.println("Line : " + line);
					System.exit(0);
				}
				
				
				
				if(strtok.hasMoreTokens()) {
					
					String token = strtok.nextToken();
					
					if(token.matches("\\d*")) {
						
						startend[1] = Integer.parseInt(token);
						
						if(startend[1] > seq2.length()) {
							System.out.println("Parse Error. Problems while parsing end index of alternative domain. Index is greater than " + seq2.length() + ". Aborting.");
							System.out.println("Index : " + token);
							System.exit(0);
						}
						
					}
					
					else {
						System.out.println("Parse Error. Problems while parsing end index of alternative domain. Aborting.");
						System.out.println("Index : " + token);
						System.exit(0);
					}
					
				}
				else {
					System.out.println("Parse Error. Problems while parsing end index of alternative domain. Aborting.");
					System.out.println("Line : " + line);
					System.exit(0);
				}
				
				
				domains2.add(startend);
				
				
			}	
			
			
			br.close();
		
		}
		catch(IOException ioe) {
			System.out.println(ioe.getMessage());
			System.out.println("IOException occurred while transforming input data to internal format.");
		}
		
		
		try {
			
			BufferedWriter bw = new BufferedWriter(new FileWriter(new File(outfile)));
			
			int SEQLINELENGTH = 60;
			
			
		// write name	
			
			bw.write("NA  " + name + "\n");
			
			bw.write("XX\n");
			
		
		// write domain(s)		
			
			for(int i=0; i<domains1.size(); i++) {
				
				bw.write("DO  " + seq1.substring(domains1.get(i)[0] - 1, domains1.get(i)[1]) + " 1 " + " " + (domains1.get(i)[0] - 1) + " " + domains1.get(i)[1] +  "\n");
				
			}
			
			
		// write alternative domain(s), if available
			
			if(seq2.length() > 0) {
			
				for(int i=0; i<domains2.size(); i++) {
					
					bw.write("DO  " + seq2.substring(domains2.get(i)[0] - 1, domains2.get(i)[1]) + " 2 " + " " + (domains2.get(i)[0] - 1) + " " + domains2.get(i)[1] +  "\n");
					
				}
				
				bw.write("XX\n");
			
			}
			
		// write sequence	
			
			for(int i=0; i<(seq1.length()/SEQLINELENGTH); i++) {
				
				bw.write("S1  "); 
				bw.write(seq1.toUpperCase(), i*SEQLINELENGTH, SEQLINELENGTH);
				bw.write("\n");
				
			}
				
			if(seq1.length()-(seq1.length()/SEQLINELENGTH)*SEQLINELENGTH > 0) {
					
				bw.write("S1  "); 
				bw.write(seq1.toUpperCase(), (seq1.length()/SEQLINELENGTH)*SEQLINELENGTH, seq1.length()-(seq1.length()/SEQLINELENGTH)*SEQLINELENGTH);
				bw.write("\n");
				
			}
			
			bw.write("XX\n");
			
			
		// write alternative sequence, if available
		
			if(seq2.length() > 0) {	
			
				for(int i=0; i<(seq2.length()/SEQLINELENGTH); i++) {
					
					bw.write("S2  "); 
					bw.write(seq2.toUpperCase(), i*SEQLINELENGTH, SEQLINELENGTH);
					bw.write("\n");
					
				}
					
				if(seq2.length()-(seq2.length()/SEQLINELENGTH)*SEQLINELENGTH > 0) {
						
					bw.write("S2  "); 
					bw.write(seq2.toUpperCase(), (seq2.length()/SEQLINELENGTH)*SEQLINELENGTH, seq2.length()-(seq2.length()/SEQLINELENGTH)*SEQLINELENGTH);
					bw.write("\n");
					
				}
				
				bw.write("XX\n");
			
			}
			
			bw.write("//\n");
			bw.write("XX\n");
			
			bw.flush();
			bw.close();
			
		}
		catch(IOException ioe) {
			System.out.println(ioe.getMessage());
			System.out.println("IOException occurred while writing transformed input data to new input file.");
		}
		
		
		String[] res = new String[5 + domains1.size() + domains2.size()] ;	
		
		res[0] = name;
		res[1] = "class" + class_id.charAt(0);
		res[2] = species;
		res[3] = seq1; 
		res[4] = seq2;
		
		for(int i=0; i<domains1.size(); i++) {
			
			res[5 + i] = seq1.substring(domains1.get(i)[0] - 1, domains1.get(i)[1]) + " 1 " + " " + (domains1.get(i)[0] - 1) + " " + domains1.get(i)[1];
			
		}
		
		for(int i=0; i<domains2.size(); i++) {
			
			res[5 + domains1.size() + i] = seq2.substring(domains2.get(i)[0] - 1, domains2.get(i)[1]) + " 2 " + " " + (domains2.get(i)[0] - 1) + " " + domains2.get(i)[1];
			
		}
		
		return res;
	}
	
	
	
	public String[] parseInternalFormat(String infile, String species_file) {
		
		
		String REGEX1 = "PF\\d{5}";
		// String REGEX2 = "[OPQ]\\d[0-9A-Z]{3}\\d(.\\d)?";
		// String REGEX3 = "At\\dg\\d{5}";
		
		
		String[] res = null;
		
		BufferedReader br = null;
		
		try {
			
			br = new BufferedReader(new FileReader(new File(infile)));
			
		// parse name
			
			StringTokenizer strtok = new StringTokenizer(br.readLine());
			
			if (!strtok.nextToken().equals("NA")) {
				System.out.println("Parse Error. Invalid format. \"NA\" expected at the beginning of the line. Aborting.");
				System.out.println("Please specify the name of the transcription factor in the given input file.");
				System.exit(0);
			}	
			
			String name = strtok.nextToken().trim();
			br.readLine();					// XX
			
			
		// parse species and check it
			
			String spec = br.readLine();	
			if (!spec.startsWith("SP")) {	
				System.out.println("Parse Error. Invalid format. \"SP\" expected at the beginning of the line. Aborting.");
				System.out.println("Please specify the species in the given input file.");
				System.exit(0);
			}	
			
			spec = spec.substring(2, spec.length()).trim();
			
			String[] split = spec.split(",");
				
			if(split.length == 1) {
				spec = split[0];
			} else if(split.length == 2) {
				spec = split[1].substring(0,split[1].length()).trim();
			} else {	
				System.out.println("Error while parsing species of " + name + "\nSpecies: " + spec);
			}

			
			ArrayList<String> known_species = new ArrayList<String>();
			
			try {
				
				BufferedReader br_spec = new BufferedReader(new FileReader(new File(species_file)));
				
				StringTokenizer spec_tok = new StringTokenizer(br_spec.readLine());
				
				while (spec_tok.hasMoreTokens()) {
					known_species.add(spec_tok.nextToken(",").trim());
				}	
				
				br_spec.close();
				
			}
			catch(IOException ioe) {
				System.out.println(ioe.getMessage());
				System.out.println("IOException occurred while parsing list of possible species.");
			}
				
			if(spec == null || !known_species.contains(spec.toUpperCase())) {
				System.out.println("Parse Error. Unknown species: \"" + spec + "\".\nList of possible species can be found in \"" + species_file + "\". Aborting.");
				System.exit(0);
			}
			
			
			br.readLine();			// XX
			
			String line = br.readLine();
			
			
		// skip cross-reference, if available	
			
			if(line.startsWith("RF")) {
				
				line = br.readLine();
				line = br.readLine();
				
			}
			
			
			String class_id = null;
			
		// parse class, if available	
			
			if(line.startsWith("CL")) {
				
				// handle Non-TF proteins contained in TFpredict output file
				if (line.startsWith("CL  Unknown") || line.startsWith("CL  Non-TF")) {
					class_id = NonTFclassID;
					return(new String[] { name, class_id, spec });
				}
				
				class_id = "class" + getTransfacClass(line.substring(2).trim()).substring(0,1);
				line = br.readLine();
				line = br.readLine(); 
			}

		// parse sequence(s)
			
			String sequence1 = null;
			String sequence2 = null;
			
			if(line.startsWith("S1")) {
				sequence1 = line.substring(2,line.length()).trim();
				while ((line = br.readLine()) != null && line.startsWith("S1")) {
					sequence1 += line.substring(2,line.length()).trim();
				}
				line = br.readLine();
			}
			
			
			
			if(line != null && line.startsWith("S2")) {
				sequence2 = line.substring(2,line.length()).trim();
				while ((line = br.readLine()) != null && line.startsWith("S2") ) {
					sequence2 += line.substring(2,line.length()).trim();
				}
				line = br.readLine();
			}
			
			
			
			ArrayList<String> pfam_domains = new ArrayList<String>();
			ArrayList<String> other_domains = new ArrayList<String>();
			
			String domain_id = null;
			
			String domain = null;
			
			int start = 0;
			int end = 0;
			
			
		// parse domain(s)	
			
			while (line != null && line.startsWith("FT") ) {
				
				strtok = new StringTokenizer(line.substring(2,line.length()).trim());
				
				domain_id = strtok.nextToken();
				
				start = Integer.parseInt(strtok.nextToken().replace("<", ""));
				end   = Integer.parseInt(strtok.nextToken().replace(">", ""));
				
			
			// domain entry comes from Pfam ?	
				
				if(domain_id.matches(REGEX1)) {
					if(end > sequence1.length()) {
						System.out.println("Problems with domain of " + name +"\nLength of S1: " + sequence1.length() + "\nEnd-Index of Domain: " + end);
					}	
					else {
						
						domain = sequence1.substring(start-1,end);
						
						pfam_domains.add(domain + " 1 " + " " + (start-1) + " " + end);
					
					}
				}
				
			// domain entry comes from Entrez !
				
				else {

			
			// sequences from Entrez and Transfac were different --> take S2	

					if (sequence2 != null) {
						if(end > sequence2.length()) {
							System.out.println("Problems with domain of " + name +"\nLength of S2: " + sequence2.length() + "\nEnd-Index of Domain: " + end);
						}	
						else {
							domain = sequence2.substring(start-1,end);
							other_domains.add(domain + " 2 " + " " + (start-1) + " " + end);
						}
					}
			
			// sequences from Entrez and Transfac were identical --> take S1	
					
					else {
						if(end > sequence1.length()) {
							System.out.println("Problems with domain of " + name);
						}	
						else {
							domain = sequence1.substring(start-1,end);
							other_domains.add(domain + " 1 " + (start-1) + " " + end);
						}
					}
				}
				
				line = br.readLine();
			}
			
			
			br.close();
			
			
			
			ArrayList<String> domains = compareDomains(other_domains, pfam_domains);
			
			
			res = new String[5 + domains.size()];	
			
			res[0] = name;
			res[1] = class_id;
			res[2] = spec;
			res[3] = sequence1; 
			res[4] = sequence2;
			
			for(int i=0; i<pfam_domains.size(); i++) {
				
				res[5 + i] = domains.get(i);
				
			}
			
			return res;
			
			
			
		}
		catch(IOException ioe) {
			System.out.println(ioe.getMessage());
			System.out.println("IOException occurred while parsing input file.");
		}
		
		
		
		return res;
		
	}
	 
	
	public String getTransfacClass(String class_id) {
		
		String CLASS_FORMAT =    "[0-4]\\p{Punct}[0-9][0-2]?\\p{Punct}";
		String CLASS_ID_FORMAT = "C00\\d\\d\\p{Punct}";
		
	// Method 1: check whether last token is of the form X.Y.Z. ...
		
		StringTokenizer strtok_class = new StringTokenizer(class_id);
		
		String class_token = null;
		
		while (strtok_class.hasMoreTokens()) {
			
			class_token = strtok_class.nextToken();
			
		}
		
		if((class_token.length() > 3 && class_token.substring(0,4).matches(CLASS_FORMAT))
		   || (class_token.length() > 4 && class_token.substring(0,5).matches(CLASS_FORMAT))) {
			
			return class_token;
			
		}
		
		
	//	Method 2: check whether class is obtainable by its class_id ...	
		
		if(class_id.length() >= 6 && class_id.substring(0,6).matches(CLASS_ID_FORMAT)) {
			
			class_id = class_id.substring(0,5);
			
		}
	
	
		BufferedReader br = null;
		
		StringTokenizer strtok = null;
		
		String line = null;
		
		String res = null;
		
		boolean found = false;
		
		try {
			
			br = new BufferedReader(new FileReader(new File(FBPPredictor.classMappingFile)));
			
			
			while((line = br.readLine()) != null) {
				
				strtok = new StringTokenizer(line);
				
				res = strtok.nextToken();
				
				if(class_id.equals(strtok.nextToken())) {
				
					found = true;
					break;
					
				}	
				
			}
			
			br.close();
			
		}
		catch (IOException e) {		
			System.out.println("IOException caught while fetching class."); 
		}
		
		if(!found) {
			System.out.println("No classification found for \"" + class_id + "\". Aborting.");
			System.exit(0);
		}
		
		return res;
		
	}
	
	
	private ArrayList<String> compareDomains(ArrayList<String> other_domains, ArrayList<String> pfam_domains) {
		
		int length_threshold = 15;
		double similarity_threshold = 0.9;
		
		String seq1 = "";
		String seq2 = "";
		double score = 0;
		double similarity = 0;
		
		StringTokenizer strtok = null;
		
		boolean match = false;
		
		try {
			
			SequenceAligner globalAligner  = new SequenceAligner(1, -1, "NW");
			
			for(int i=0; i<other_domains.size(); i++) {
				
				strtok = new StringTokenizer(other_domains.get(i));
				
				seq1 = strtok.nextToken();
				
			
				
				for(int j=0; j<pfam_domains.size(); j++) {
					
					strtok = new StringTokenizer(pfam_domains.get(j));
					
					seq2 = strtok.nextToken();
					
					score = 0.0;
					similarity = 0.0;
					
					if(Math.max(seq1.length(), seq2.length()) - Math.min(seq1.length(), seq2.length()) <= length_threshold) {
						
						score = globalAligner.getAlignmentScore(seq1, seq2);
						
						similarity = score / Math.max(seq1.length(), seq2.length());
						
						
						if(similarity >= similarity_threshold) {
						
							match = true;
							break;
							
						}	
					}
				}
				
			// domain was very similar to another domain? --> ignore this domain.	
				
				if(match) {
					match = false;
				}
				
			// domain was unique? --> add domain to non-redundant domain set (pfam-domains)	
				else {
					pfam_domains.add(other_domains.get(i));
				}
					
			}	
		}
		catch (BioException e) {System.out.println("Exception caught while comparing domains.");} 
		catch (IOException e) {	System.out.println("IOException caught."); } 
		catch (Exception e) { System.out.println("Exception caught.");}
		
		
		return pfam_domains;
	}
	
}

