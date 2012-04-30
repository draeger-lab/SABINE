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
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.StringTokenizer;

import org.biojava.bio.BioException;


public class DomainMatrixSequenceExtractor {
	
	boolean considerClasses = false;
	String superclass = "1.";
	
	boolean considerQualities = false;
	int quality_threshold = 4;
	
	boolean control_quality_annotations = false;
	
	
	public void extractDomainsAndMatrices(String infile, String outfile, String species) {
		
		BufferedReader br = null;
		BufferedWriter bw = null;
		
		String domains_known = "";
		
		String name = null;
		String spec = null;
		String domain = null;
		String line = null;
		String sequence1 = null;
		String sequence2 = null;
		String domain_id = null;
		
		ArrayList<String> pfam_domains = null;
		ArrayList<String> other_domains = null;
		ArrayList<String> matrices = null;
		ArrayList<String> domains = null;
		
		String REGEX1 = "PF\\d{5}";
		String REGEX2 = "[OPQ]\\d[0-9A-Z]{3}\\d(.\\d)?";
		String REGEX3 = "At\\dg\\d{5}";
		
		String CLASS_FORMAT =    "[0-4]\\p{Punct}[0-9][0-2]?\\p{Punct}";
		String CLASS_ID_FORMAT = "C00\\d\\d\\p{Punct}";
		
		int start = 0;
		int end = 0;
		
		String[] split = null;
		
		StringTokenizer strtok = null;
		StringTokenizer strtok_class = null;
		
		boolean obtain = false;
		
		try {

			br = new BufferedReader(new FileReader(new File(infile)));
			bw  = new BufferedWriter(new FileWriter(new File(outfile)));
			

			
			while((line = br.readLine()) != null) {
			
			
				
				obtain = false;	
				
			// parse name
				
				if(!line.startsWith("NA  ")) {
					System.out.println("Parse Error. NA expected.\nLine = " + line);
					System.exit(0);
				}
				else {
					name = new StringTokenizer(line.substring(4,line.length())).nextToken();		// SwissProt ID
				}
				
				line = br.readLine();
				
			// parse domain-knowledge-flag
				
				line = br.readLine();
				domains_known = line;
				
				line = br.readLine();
				spec = br.readLine();
			
			// parse species	
				
				spec = spec.substring(4, spec.length());
				
				split = spec.split(",");
				if(split.length == 1) {
					spec = split[0];
				}
				else {
					if(split.length == 2) {
						spec = split[1].substring(1,split[1].length());
					}
					else{
						System.out.println("Error while parsing species of " + name + "\nSpecies: " + spec);
					}
				}
			
			
			
				//	species OK? --> store domains and matrices in output file 	
				
					if(spec.matches(species)) {
						
						obtain = true;
					
					}
		
				// if classes are considered --> entry must only be obtained if class info is available	
					
					if(considerClasses) {
							
						obtain = false;
						
					}
						
				// goto sequence
					
					while (!(line = br.readLine()).startsWith("S1  ") && !line.startsWith("S2  ")) {
							
						if(considerClasses && line.startsWith("CL  ")) {
								
							
					// parse class		
							
							String factor_class = null;
							
							
						// Method 1: check whether last token is of the form X.Y.Z. ...

							strtok_class = new StringTokenizer(line.substring(4));
							
							String class_token = null;
							
							while (strtok_class.hasMoreTokens()) {
								class_token = strtok_class.nextToken();
							}
							
							if(class_token.length() > 3 && class_token.substring(0,4).matches(CLASS_FORMAT)) {
								
								factor_class = class_token;
								
							}
							
							
						//	Method 2: check whether class is obtainable by its class_id ...	
							
							if(factor_class == null && line.substring(4,10).matches(CLASS_ID_FORMAT)) {
								
								factor_class = getTransfacClass(line.substring(4,9));
								
							}
							
							if(factor_class == null) {
								
								System.out.println("No class found for " + name + ". Aborting.");
								System.exit(0);
								
							}
							
						// obtain domains and matrices of this factor if its class matches the TRANSFAC-superclass
							
							if(factor_class.startsWith(superclass)) {
								
								obtain = true;
								
							}
							
							if(obtain) System.out.println("\t OBTAINED.");
							else       System.out.println("\t REJECTED.");							
						}
							
							
					}
				
			// obtain this entry
					
					if(obtain) {
					
					System.out.println(name);
					
				
				// parse sequence	
					
					sequence1 = null;
					sequence2 = null;
					
					if(line.startsWith("S1  ")) {
						sequence1 = line.substring(4,line.length());
						while ( (line = br.readLine()).startsWith("S1  ") ) {
							sequence1 += line.substring(4,line.length());
						}
						line = br.readLine();
					}
					
					if(line.startsWith("S2  ")) {
						sequence2 = line.substring(4,line.length());
						while ( (line = br.readLine()).startsWith("S2  ") ) {
							sequence2 += line.substring(4,line.length());
						}
						line = br.readLine();
					}
					
					pfam_domains = new ArrayList<String>();	
					other_domains = new ArrayList<String>();
					
					
				// generate domain sequences
					
					while ( line.startsWith("FT  ") ) {
						
						strtok = new StringTokenizer(line.substring(4,line.length()));
						
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
							if(!domain_id.matches(REGEX2)) {
								
								if(!domain_id.equals("AAH36092") && !domain_id.matches(REGEX3) && !spec.equals("Drosophila melanogaster")) {
									System.out.println("Fatal Error. Domain ID " + domain_id + " does not match regular expressions.");
									System.exit(1);
								}	
							}
					
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
					
		//			System.out.println("#Pfam  Domains: " + pfam_domains.size());
		//			System.out.println("#Other Domains: " + other_domains.size());
					
					matrices = new ArrayList<String>();	
					
				// parse matrices
					
					while(!(line = br.readLine()).startsWith("//")) {
						
						if(!line.startsWith("MN  ")) {
							System.out.println("Parse Error. MN expected.\nLine = " + line);
							System.exit(0);
						}
						
				// check quality of this matrix
						
						if(!considerQualities || checkQuality(line.substring(4))) {
						
							line = br.readLine();	// XX
						
							matrices.add ( 	  br.readLine() + "\n"
											+ br.readLine() + "\n"
											+ br.readLine() + "\n"
											+ br.readLine() + "\n" );
						
						}
						
						else {
							
							line = br.readLine();	// XX
							
							br.readLine();
							br.readLine();
							br.readLine();
							br.readLine();
							
						}
						
						line = br.readLine();

					}
					
					
					
					domains = compareDomains(other_domains, pfam_domains);
					
					int SEQLINELENGTH = 60;
					
		/*
		 * 
		 * write entry to output file
		 * 
		 */			
					
					if(domains.size() > 0 && matrices.size() > 0) {
						bw.write("NA  " + name + "\n");
						bw.write("XX\n");
						
						bw.write(domains_known + "\n");
						bw.write("XX\n");
					
						for(int i=0; i<domains.size(); i++) {
							bw.write("DO  " + domains.get(i) + "\n");
						}
						
						bw.write("XX\n");
						
						
						if(sequence1 != null) {
							
							for(int i=0; i<(sequence1.length()/SEQLINELENGTH); i++) {
								bw.write("S1  "); 
								bw.write(sequence1.toUpperCase(), i*SEQLINELENGTH, SEQLINELENGTH);
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
							
							for(int i=0; i<(sequence2.length()/SEQLINELENGTH); i++) {
								bw.write("S2  "); 
								bw.write(sequence2.toUpperCase(), i*SEQLINELENGTH, SEQLINELENGTH);
								bw.write("\n");
							}
							if(sequence2.length()-(sequence2.length()/SEQLINELENGTH)*SEQLINELENGTH > 0) {
								bw.write("S2  "); 
								bw.write(sequence2.toUpperCase(), (sequence2.length()/SEQLINELENGTH)*SEQLINELENGTH, sequence2.length()-(sequence2.length()/SEQLINELENGTH)*SEQLINELENGTH);
								bw.write("\n");
							}
							
							bw.write("XX\n");
							
						}
						
						
						for(int i=0; i<matrices.size(); i++) {
							bw.write(matrices.get(i));
							bw.write("XX\n");
						}
						
						
						bw.write("//\n");
						bw.write("XX\n");
					}
					
					line = br.readLine();
				
				}
					
		//  obtain == false --> proceed to next entry	
				
				else {
					
					while(!line.startsWith("//"))line = br.readLine();
					line = br.readLine();
				
				}
				
			}
			
			bw.flush();
			bw.close();
			
			
		}
		catch(FileNotFoundException fnfe) {
			System.out.println(fnfe.getMessage());
			System.out.println("File not found.");
		}
		catch(IOException ioe) {
			System.out.println(ioe.getMessage());
			System.out.println("IO Exception occured.");
		}
		
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
						
						
		//				System.out.println("Comparing Domains ");
		//				System.out.println(seq1);
		//				System.out.println(seq2);
		//				System.out.println("Similarity: " + similarity + "\n");
						
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
	
	public String getTransfacClass(String class_id) {
		
		BufferedReader br = null;
		
		StringTokenizer strtok = null;
		
		String line = null;
		
		String res = null;
		
		boolean found = false;
		
		try {
			
			br = new BufferedReader(new FileReader(new File("interleaved" + File.separator + "Classes")));
			
			
			while((line = br.readLine()) != null) {
				
				strtok = new StringTokenizer(line);
				
				res = strtok.nextToken();
				
				if(class_id.equals(strtok.nextToken())) {
				
					found = true;
					break;
					
				}	
				
			}
			
		}
		catch (IOException e) {		
			System.out.println("IOException caught while fetching class."); 
		}
		
		if(!found) {
			System.out.println("No classification found for " + class_id + ". Aborted.");
			System.exit(0);
		}
		
		return res;
		
	}
	
	/*
	 * 
	 * gets the MN-line of a matrix entry without the "MN  "-part
	 * returns true iff the matrix quality exceeds the specified threshold
	 * 
	 */
	
	public boolean checkQuality(String matrix_name) {
		
		
		String TRANSFAC_REGEX    = "M\\d{5}";
		String JASPAR_REGEX      = "MA\\d{4}";
		String YEASTRACT_REGEX   = "MY\\d{4}";
		
		String QUALITY_REGEX = "Q\\d";
		
		boolean res = false;
		
		
		StringTokenizer strtok = new StringTokenizer(matrix_name);
		
		
		String matrix_identifier = strtok.nextToken();
		
		
	// matrix comes from transfac ?
		if(matrix_identifier.matches(TRANSFAC_REGEX)) {
			
			StringTokenizer strtok2 = new StringTokenizer(strtok.nextToken(), "_");
			strtok2.nextToken();
			String quality = strtok2.nextToken();
			
			if(control_quality_annotations) {
				System.out.println(quality);
			}
			
			if(		quality.matches(QUALITY_REGEX) 
					&& 
					Integer.parseInt("" + quality.charAt(1)) <= quality_threshold
			) {
				
				res = true;
				
			}
			
		}
		
	// matrix comes from jaspar ?
		else {
				
			if(matrix_identifier.matches(JASPAR_REGEX)) {
					
				res = true;
				
			}
			
	// matrix comes from yeastract ?
			else {
				
				if(matrix_identifier.matches(YEASTRACT_REGEX)) {
						
					res = false;
					
				}
				
	// matrix comes from another database ?			
				else {
					
					System.out.println("Source of " + matrix_identifier + " could not be identified.");
				
				}
			}
		}	
		
		return res;
	}
	
	

	public static void main(String[] args) {
		
		DomainMatrixSequenceExtractor extractor = new DomainMatrixSequenceExtractor();
		
		extractor.considerClasses = true;
		extractor.superclass = "4.";
		
		extractor.considerQualities = false;
		extractor.quality_threshold = 5;
		
		extractor.control_quality_annotations = false;
		
		extractor.extractDomainsAndMatrices("trainingsets_new/transpall_interleaved_classes.out", "trainingsets_new/sequences_matrices_class4.out", ".*");
		
	}
	
}


