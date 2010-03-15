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

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

public class DBDEnsemblParser {
	
	
	public void writeSabineInfile() {
		
		DBDParser dbd_parser =  new DBDParser();
		
		// set species name and ID
		dbd_parser.species_name = "Streptomyces coelicolor";
		dbd_parser.species_ID = "sf";
	
		// parse datasets and save variables to disk
		String ensembl_dir = "../data/ensembl/";
		String seq_filename = ensembl_dir + dbd_parser.species_ID + ".fa";
		String dom_filename = ensembl_dir + dbd_parser.species_ID + ".tf.ass";
		
		System.out.println("  Parsing Ensembl sequences.\n    (file: " + seq_filename + ")\n");
		dbd_parser.parseProtIDs(dom_filename);

		dbd_parser.parseSequences(seq_filename);
		System.out.println("  " + dbd_parser.sequences.size() + " sequences parsed.\n");
	
		
		System.out.println("  Parsing DBD domain assignments.\n    (file: " + dom_filename + ")\n");
		dbd_parser.parseDomains(dom_filename);
		
		// transfer fields to DatabaseParser object
		DatabaseParser data_parser = new DatabaseParser(); 
		data_parser.transferDBDFields(dbd_parser);
		
		int old_num_proteins = data_parser.domains.size();
		data_parser.filterDomains();
		System.out.println("  Domain assignment found for " + data_parser.domains.size() + " / "
						   + old_num_proteins + " proteins.\n");
		
		// write input file for Standalone predictor
		int[] all_factors = new int[data_parser.sequences.size()];
		for (int i=0; i<all_factors.length; i++) all_factors[i] = i;
		
		String outfile_all_factors = "input/" + dbd_parser.species_ID + "_all_factors.input";
		System.out.println("  Writing input file for SABINE which contains all factors.\n" +
						   "    (file: " + outfile_all_factors + ")\n");
		 
		writeInputFile(outfile_all_factors, data_parser, all_factors);
	}
	
	
	
	public void writeInputFile(String outfile, DatabaseParser parser, int[] relevant_factors) {
		
		String curr_seq;
		int SEQLINELENGTH = 60;
		int i;
		
		try {
			BufferedWriter bw = new BufferedWriter(new FileWriter(new File(outfile)));
			
			for (int x=0; x<relevant_factors.length; x++) {
				
				i = relevant_factors[x];
				
				/*
				 *    write name and species
				 */
				
				bw.write("NA  " + parser.prot_IDs.get(i) + "\n");
				bw.write("XX\n" +   
						 "SP  Streptomyces Coelicolor\n" + 
						 "XX\n"); 
				
				
				/*
				 *    write TRANSFAC superclass
				 */
				
				bw.write("CL  0.0.0.0.0.\n" + 
						 "XX\n");
				
				/*
				 *    write sequence 
				 */
				
				curr_seq = parser.sequences.get(i);
				
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
				
				/*
				 *    write domains 
				 */
				
				for (int j=0; j<parser.domains.get(i).size(); j++) {
					bw.write("FT  " + parser.hmm_IDs.get(i).get(j) + "\t" + parser.domains.get(i).get(j) + "\n");
				}
				
				bw.write("XX\n" + 
						 "//\n" +
						 "XX\n");	
			}
			
			bw.flush();
			bw.close();
		}
		
		catch(IOException ioe) {
			System.out.println(ioe.getMessage());
			System.out.println("IOException occurred while writing input file for SABINE.");
		}
	}
	
	
	
	
	public static void main(String[] args) {
		
		DBDEnsemblParser parser = new DBDEnsemblParser();
		parser.writeSabineInfile();
	}
}

