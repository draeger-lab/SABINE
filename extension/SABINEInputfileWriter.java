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

import help.FileFilter;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.StringTokenizer;

public class SABINEInputfileWriter {
	
	ArrayList<String> tf_names = new ArrayList<String>();
	ArrayList<String> uniprot_IDs = new ArrayList<String>();
	ArrayList<String> species = new ArrayList<String>();
	ArrayList<String[]> pwms = new ArrayList<String[]>();
	ArrayList<String> sequences = new ArrayList<String>();
	ArrayList<ArrayList<String>> domains = new ArrayList<ArrayList<String>>();
	
	ArrayList<String> unique_tf_names = new ArrayList<String>();
	ArrayList<String> unique_species = new ArrayList<String>();
	ArrayList<IDMappingParser> id_parser = new ArrayList<IDMappingParser>();

	ArrayList<UniProtParser4Procaryotes> sequence_parser = new ArrayList<UniProtParser4Procaryotes>();
	
	
	public void readInfile(String infile) {
		
		String line;
		String[] curr_matrix;
		
		try {
			 
			 BufferedReader bw = new BufferedReader(new FileReader(new File(infile)));
			 
			 while ((line = bw.readLine()) != null) {

				 if (line.startsWith("NA")) {
					 tf_names.add(line.substring(3).trim());
				 }
				 else {
					 System.out.println("Parse Error. \"NA\" expected at beginning of the line.");
					 System.exit(0);
				 }
				 bw.readLine();
				 if ((line = bw.readLine()).startsWith("SP")) {
					 species.add(line.substring(3).trim());
				 }
				 bw.readLine();
				 if ((line = bw.readLine()).startsWith("MA")) {
					curr_matrix = new String[4];
					curr_matrix[0] = line.substring(3).trim();
					curr_matrix[1] = bw.readLine().substring(3).trim();
					curr_matrix[2] = bw.readLine().substring(3).trim();
					curr_matrix[3] = bw.readLine().substring(3).trim();
					pwms.add(curr_matrix);
				 }
				 bw.readLine();		// XX
				 bw.readLine();		// //
				 bw.readLine();		// XX
			 }
			 bw.close();
		}
		catch(IOException ioe) {
			System.out.println(ioe.getMessage());
			System.out.println("IOException occurred while writing species mapping.");
		}
	}
	
	public void prepareIDMapping(String mapdir) {
		
		StringTokenizer strtok;
		String curr_species;
		
		// generate unique species list
		for (int i=0; i<species.size(); i++) {
			
			strtok = new StringTokenizer(species.get(i));
			curr_species = strtok.nextToken() + " " + strtok.nextToken();
			
			if (! unique_species.contains(curr_species)) {
				unique_species.add(curr_species);
			}
		}
		
		// compute mapping (gene name -> UniProt ID)
		IDMappingParser curr_id_parser;
		String curr_infile;
		String species_part2;
		
		for (int i=0; i<unique_species.size(); i++) { 
			
			strtok = new StringTokenizer(unique_species.get(i)); 
			curr_infile = mapdir + strtok.nextToken();
			species_part2 = strtok.nextToken();
			curr_infile += species_part2.substring(0,1).toUpperCase() + species_part2.substring(1) + "_MappingFile.txt";
			curr_id_parser = new IDMappingParser();
			curr_id_parser.parseIDMapping(curr_infile);
			id_parser.add(curr_id_parser);
		}
	}
		
	public void addUniProtIDs() {
		
		int name_idx, spec_idx;
		StringTokenizer strtok;
		String curr_species;
		String curr_name;
		
		for (int i=0; i<tf_names.size(); i++) {
			
			strtok = new StringTokenizer(species.get(i));
			curr_species = strtok.nextToken() + " " + strtok.nextToken();
			spec_idx = unique_species.indexOf(curr_species);
			
			curr_name = new StringTokenizer(tf_names.get(i)).nextToken().toUpperCase();
			name_idx = id_parser.get(spec_idx).gene_names.indexOf(curr_name);
			
			if (name_idx == -1) {
				uniprot_IDs.add("NA");
				// System.out.println("No UniProt ID found for " + tf_names.get(i) + " of " + curr_species + ".");
			}
			else {
				uniprot_IDs.add(id_parser.get(spec_idx).uniprot_IDs.get(name_idx));
			}
		}
	}	
		
	public void writeInfile(String outfile) {
		
		try {
			 BufferedWriter bw = new BufferedWriter(new FileWriter(new File(outfile)));
			 
			 for (int i=0; i<tf_names.size(); i++) {
			 
				 bw.write("NA  " + tf_names.get(i) + "\n");
				 bw.write("XX  \n");
				 bw.write("SP  " + species.get(i) + "\n");
				 bw.write("XX  \n");
				 bw.write("RF  " + uniprot_IDs.get(i) + "\n");
				 bw.write("XX  \n");
				 bw.write("MA  " + pwms.get(i)[0] + "\n");
				 bw.write("MA  " + pwms.get(i)[1] + "\n");
				 bw.write("MA  " + pwms.get(i)[2] + "\n");
				 bw.write("MA  " + pwms.get(i)[3] + "\n");
				 bw.write("XX  \n");
				 bw.write("//  \n");
				 bw.write("XX  \n");
			 }
			 bw.flush();
			 bw.close();
		}
		catch(IOException ioe) {
			System.out.println(ioe.getMessage());
			System.out.println("IOException occurred while writing species mapping.");
		}
	}
	
	public void mapFactors2UniProt(String indir, String mapdir, String outdir) {
		
		/*
		 *  list pwm files in given path
		 */		
		
		FileFilter filter = new FileFilter();
		filter.setFormat(".*.pwms");
		filter.setDirectory(indir);
		String[] files = filter.listFilesFullPath();
		
		String infile;
		for (int i=0; i<files.length; i++) {
			infile = indir + files[i];
			
			readInfile(infile);
			prepareIDMapping(mapdir);
			addUniProtIDs();
			writeInfile(outdir + files[i]);
			
			tf_names.clear();
			uniprot_IDs.clear();
			species.clear();
			pwms.clear();
			unique_species.clear();
			id_parser.clear();
		}
	}
	
	public int addInfile(String infile, String database) {
		
		String line, curr_name;
		String[] curr_matrix;
		int tf_counter = 0;
		
		try {
			 
			 BufferedReader bw = new BufferedReader(new FileReader(new File(infile)));
			 
			 while ((line = bw.readLine()) != null) {

				 if (line.startsWith("NA")) {
					 curr_name = line.substring(3).trim();
					 
					 if (!unique_tf_names.contains(new StringTokenizer(curr_name.toUpperCase()).nextToken())) {
						 unique_tf_names.add(new StringTokenizer(curr_name.toUpperCase()).nextToken());
						 tf_names.add(curr_name + "  [" + database + "]");
						 tf_counter++;
					 }
					 else {
						 while ((line = bw.readLine()) != null && !line.startsWith("//"));
						 bw.readLine();
						 continue;
					 }
				 }
				 else {
					 System.out.println("Parse Error. \"NA\" expected at beginning of the line.");
					 System.exit(0);
				 }
				 bw.readLine();
				 if ((line = bw.readLine()).startsWith("SP")) {
					 species.add(line.substring(3).trim());
				 }
				 bw.readLine();
				 if ((line = bw.readLine()).startsWith("RF")) {
					 uniprot_IDs.add(line.substring(3).trim());
				 }
				 bw.readLine();
				 if ((line = bw.readLine()).startsWith("MA")) {
					curr_matrix = new String[4];
					curr_matrix[0] = line.substring(3).trim();
					curr_matrix[1] = bw.readLine().substring(3).trim();
					curr_matrix[2] = bw.readLine().substring(3).trim();
					curr_matrix[3] = bw.readLine().substring(3).trim();
					pwms.add(curr_matrix);
				 }
				 bw.readLine();		// XX
				 bw.readLine();		// //
				 bw.readLine();		// XX
			 }
			 bw.close();
		}
		catch(IOException ioe) {
			System.out.println(ioe.getMessage());
			System.out.println("IOException occurred while writing species mapping.");
		}
		return tf_counter;
	}
	
	public void mergeFactorsByName(ArrayList<String> infiles, String outfile) {
		
		//clear class variables
		tf_names.clear();
		uniprot_IDs.clear();
		species.clear();
		pwms.clear();
		unique_species.clear();
		id_parser.clear();
		
		// read input files
		int sum_counter = 0;
		int tf_counter = 0;
		String database;
		
		for (int i=0; i<infiles.size(); i++) {
			
			database = infiles.get(i).substring(infiles.get(i).lastIndexOf("/")+1,infiles.get(i).lastIndexOf("."));
			tf_counter = addInfile(infiles.get(i), database);
			sum_counter += tf_counter;
			System.out.println(tf_counter + " PWMs found in file: " + infiles.get(i).substring(infiles.get(i).lastIndexOf("/")+1, infiles.get(i).length()));
		}
		System.out.println("\n" + sum_counter + " PWMs found.");
		writeInfile(outfile);
	}
	
	public void filterByUniprotID(String infile, String outfile) {
		
		boolean[] id_found = new boolean[tf_names.size()];
		for (int i=0; i<tf_names.size(); i++) {
			if (! uniprot_IDs.get(i).equals("NA")) {
				id_found[i] = true;
			}
		}
		
		int num_factors = tf_names.size();
		for (int i=tf_names.size()-1; i>0; i--) {
			if (! id_found[i]) {
				tf_names.remove(i);
				uniprot_IDs.remove(i);
				species.remove(i);
				pwms.remove(i);
			}
		}
		System.out.println("UniProt IDs found for " + tf_names.size() + "/" + num_factors + " factors.");
		writeInfile(outfile);
	}
	
	public void init_sequence_parser(String indir) {
		
		FileFilter filter = new FileFilter();
		filter.setFormat(".*.dat");
		filter.setDirectory(indir);
		String[] files = filter.listFilesFullPath();
		
		StringTokenizer strtok;
		String species_part1, species_part2;
		UniProtParser4Procaryotes uniprot_parser;
		
		String curr_species;
		
		for (int i=0; i<species.size(); i++) {
			
			strtok = new StringTokenizer(species.get(i));
			curr_species = strtok.nextToken() + " " + strtok.nextToken();
			
			if (! unique_species.contains(curr_species)) {
				unique_species.add(curr_species);
			}
		}
		
		for (int i=0; i<unique_species.size(); i++) {
			for (int j=0; j<files.length; j++) {
			
				strtok = new StringTokenizer(unique_species.get(i)); 
				species_part1 = strtok.nextToken();
				species_part2 = strtok.nextToken();
				species_part2 = species_part2.substring(0,1).toUpperCase() + species_part2.substring(1);
				
				if (files[j].startsWith(species_part1 + species_part2)) {
					uniprot_parser = new UniProtParser4Procaryotes();
					uniprot_parser.parseSequencesAndDomains(indir + files[j]);
					sequence_parser.add(uniprot_parser);
					// System.out.println("Processing: " + species_part1 + " " + species_part2);
					break;
				}
			}
		}
	}
	
	public void add_sequences() {
		
		int spec_idx = 0;
		int all_id_idx;
		int id_idx;
		
		String curr_seq;
		ArrayList<String> curr_domains;
		
		String curr_species;
		StringTokenizer strtok;
		
		for (int i=0; i<uniprot_IDs.size(); i++) {
			
			strtok = new StringTokenizer(species.get(i));
			curr_species = strtok.nextToken() + " " + strtok.nextToken();
			
			spec_idx = unique_species.indexOf(curr_species);
			curr_seq = "NA";
			curr_domains = new ArrayList<String>();
			
			all_id_idx = sequence_parser.get(spec_idx).all_uniprot_IDs.indexOf(uniprot_IDs.get(i));
			if (all_id_idx != -1) {
				id_idx = sequence_parser.get(spec_idx).id_map.get(all_id_idx)[0].intValue();
				
				curr_seq = sequence_parser.get(spec_idx).sequences.get(id_idx);
				curr_domains = sequence_parser.get(spec_idx).domains.get(id_idx);
			}
			sequences.add(curr_seq);
			domains.add(curr_domains);
		}
	}
	
	public void writeOutfile(String outfile) {
		
		int SEQLINELENGTH = 60;
		StringTokenizer strtok;
		
		try {
			BufferedWriter bw = new BufferedWriter(new FileWriter(new File(outfile)));
			
			int num_factors = tf_names.size();
			int num_domains = 0;
			
			for (int i=0; i<tf_names.size(); i++) {
				
				// skip entry if domain annotation or pwm is missing
				if (domains.get(i).isEmpty() || pwms.get(i)[0].isEmpty()) {
					continue;
				}
				num_domains++;
				
				strtok = new StringTokenizer(species.get(i));
				
				bw.write("NA  " + new StringTokenizer(tf_names.get(i)).nextToken() + "\n");
				bw.write("XX  \n");
				bw.write("SP  " + strtok.nextToken() + " " + strtok.nextToken() + "\n");
				bw.write("XX  \n");
				bw.write("RF  " + uniprot_IDs.get(i) + "\n");
				bw.write("XX  \n");
				bw.write("CL  6.0.0.0.0. \n");
				bw.write("XX  \n");
				
				// write sequence
				for (int j=0; j<(sequences.get(i).length()/SEQLINELENGTH); j++) {
					
					bw.write("S1  "); 
					bw.write(sequences.get(i).toUpperCase(), j*SEQLINELENGTH, SEQLINELENGTH);
					bw.write("\n");
				}
				if(sequences.get(i).length()-(sequences.get(i).length()/SEQLINELENGTH)*SEQLINELENGTH > 0) {
					
					bw.write("S1  "); 
					bw.write(sequences.get(i).toUpperCase(), (sequences.get(i).length()/SEQLINELENGTH)*SEQLINELENGTH, sequences.get(i).length()-(sequences.get(i).length()/SEQLINELENGTH)*SEQLINELENGTH);
					bw.write("\n");
				}
				bw.write("XX  \n");
				
				// write domains
				for (int j=0; j<domains.get(i).size(); j++) {
					bw.write("FT  " + uniprot_IDs.get(i) + "    " + domains.get(i).get(j) + "\n");
				}
				bw.write("XX  \n");

				// write matrix
				bw.write("MN  " + tf_names.get(i) + "\n");
				bw.write("XX\n");
				for (int j=0; j<4; j++) {
					bw.write("MA  " + pwms.get(i)[j] + "\n");
				}
				bw.write("XX  \n");
				bw.write("//  \n");
				bw.write("XX  \n");
			}
			bw.flush();
			bw.close();
			
			System.out.println("Domains found for " + num_domains + "/" + num_factors + " factors.");
		}
		catch(IOException ioe) {
			System.out.println(ioe.getMessage());
			System.out.println("IOException occurred while writing output file");
		}
	}
	
	public void checkUniprotIDs() {
		
		ArrayList<String> unique_uniprot_IDs = new ArrayList<String>();
		int idx;
		
		for (int i=0; i<uniprot_IDs.size(); i++) {
			
			if (! unique_uniprot_IDs.contains(uniprot_IDs.get(i))) {
				unique_uniprot_IDs.add(uniprot_IDs.get(i));
			}
			else {
				idx = unique_uniprot_IDs.indexOf(uniprot_IDs.get(i));
				System.out.println("Factor " + tf_names.get(idx) + " and factor " + tf_names.get(i) + " have identical IDs.");
				System.exit(0);
			}
		}
		System.out.println("Number of entrys:      " + uniprot_IDs.size());
		System.out.println("Number of UniProt IDs: " + unique_uniprot_IDs.size());
	}
	
	public static void main(String[] args) {
		
		String path = "/afs/cs/ra/share/eichner/pwm_procaryotes/";
		//String path = "S:/share/eichner/pwm_procaryotes/";
		
		String indir = path + "pwms/";
		String mapdir = path + "MappingFiles/";
		String outdir = path + "mapped_pwms/";
		String seqdir = path + "Uniprot_Annotations/";
		
		SABINEInputfileWriter infile_writer = new SABINEInputfileWriter();
		infile_writer.mapFactors2UniProt(indir, mapdir, outdir);
		
		ArrayList<String> infiles = new ArrayList<String>();
		infiles.add(outdir + "regulondb.pwms");
		infiles.add(outdir + "prodoric.pwms");    
		infiles.add(outdir + "dp_interact.pwms");  
		infiles.add(outdir + "swiss_regulon.pwms");
		String outfile = outdir + "all_databases.pwms";
		
		infile_writer.mergeFactorsByName(infiles, outfile);
		infile_writer.filterByUniprotID(outfile, outfile);
		infile_writer.init_sequence_parser(seqdir);
		infile_writer.add_sequences();
		infile_writer.writeOutfile(path + "training_factors_procaryotes.txt");
		infile_writer.checkUniprotIDs();
		
		for (int i=0; i<infile_writer.unique_species.size(); i++) {
			System.out.println(infile_writer.unique_species.get(i));
		}
		
		/*
		String indir = "/afs/cs/ra/share/eichner/pwm_procaryotes/pwms/";
		String infile = "dp_interact.pwms";
		
		SABINEInputfileWriter infile_writer = new SABINEInputfileWriter();
		infile_writer.readInfile(indir + infile);
		
		String mapdir = "/afs/cs/ra/share/eichner/pwm_procaryotes/MappingFiles/";
		
		infile_writer.prepareIDMapping(mapdir);
		infile_writer.addUniProtIDs();
		
		infile_writer.writeInfile(indir + outfile);
		*/
	}
}

