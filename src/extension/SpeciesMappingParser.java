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


public class SpeciesMappingParser {

	ArrayList<String> species_IDs = new ArrayList<String>();
	ArrayList<String> species_names = new ArrayList<String>();
	
	public void parseSpeciesMapping(String infile) {
		
		StringTokenizer strtok;
		String line;
		BufferedReader br;
		
		try {
			 
			 br = new BufferedReader(new FileReader(new File(infile)));
			 
			 while((line = br.readLine()) != null) {
						
				strtok = new StringTokenizer(line);
				 
				species_IDs.add(strtok.nextToken());					 	// parse species ID
				species_names.add(strtok.nextToken().toUpperCase() + " " + 
								  strtok.nextToken().toUpperCase()); 		// parse species name
			 }
			
			 
		}
		catch(IOException ioe) {
			System.out.println(ioe.getMessage());
			System.out.println("IOException occurred while parsing species mapping.");
		}
	}
	
	public void filterRelevantSpecies(String infile) {
		
		if ((species_IDs.size() == 0) || (species_names.size() == 0)) {
			System.out.println("Fatal Error. Unable to filter species mapping. Global variable \"species_names\" or \"species_IDs\" was not initialized. Aborting.");
			System.exit(0);
		}
		
		ArrayList<String> species = null;
			
		try {
			 
			 BufferedReader br = new BufferedReader(new FileReader(new File(infile)));
			 String line = br.readLine();
				
			 StringTokenizer strtok = new StringTokenizer(line, ",");
			 species = new ArrayList<String>();
				
			 while(strtok.hasMoreTokens()) {
				 
				 species.add(strtok.nextToken().trim());
					
			 }		
		}
		catch(IOException ioe) {
			System.out.println(ioe.getMessage());
			System.out.println("IOException occurred while filtering species mapping.");
		}
		 
		for (int i=species_names.size()-1; i>=0; i--) {
			if (! species.contains(species_names.get(i))) {
				species_names.remove(i);
				species_IDs.remove(i);
			}
		}
	}
	
	
	public void write_mapping(String outfile) {
		
		if ((species_IDs.size() == 0) || (species_names.size() == 0)) {
			System.out.println("Fatal Error. Unable to filter species mapping. Global variable \"species_names\" or \"species_IDs\" was not initialized. Aborting.");
			System.exit(0);
		}
		
		try {
			 
			 BufferedWriter bw = new BufferedWriter(new FileWriter(new File(outfile)));
			 
			 for (int i=0; i<species_names.size(); i++) {
				 bw.write(species_IDs.get(i) + "\t" + species_names.get(i) + "\n");
			 }
			 bw.flush();
			 bw.close();
			 
		}
		catch(IOException ioe) {
			System.out.println(ioe.getMessage());
			System.out.println("IOException occurred while writing species mapping.");
		}	
	}
	
	private void writeSpeciesFile(String infile, String outfile) {
		
		String line;
		ArrayList<Integer> species_counter = new ArrayList<Integer>();
		int idx;
		
		try {
			 
			 BufferedReader br = new BufferedReader(new FileReader(new File(infile)));
			 
			 while ((line = br.readLine()) != null) {
				 
				 if (line.startsWith("SP  ") && !species_names.contains(line.substring(4).trim())) {
					 species_names.add(line.substring(4).trim());
					 species_counter.add(1);
				 }
				 else if (line.startsWith("SP  ")) {
					 idx = species_names.indexOf(line.substring(4).trim());
					 species_counter.set(idx, species_counter.get(idx) + 1);
				 }
			 }
			 br.close();
			 
			 BufferedWriter bw = new BufferedWriter(new FileWriter(new File(outfile)));
			 int sum = 0;
			 for (int i=0; i<species_names.size(); i++) {
				 bw.write(species_counter.get(i) + "\t" + species_names.get(i) + "\n");
				 sum += species_counter.get(i);
			 }
			 bw.flush();
			 bw.close();
			 
			 System.out.println(sum + " entries found !");
		}
		catch(IOException ioe) {
			System.out.println(ioe.getMessage());
			System.out.println("IOException occurred while writing species mapping.");
		}	
	}
	
	
	public static void main(String[] args) {
		
		SpeciesMappingParser species_parser = new SpeciesMappingParser();
		
		String indir = "/afs/cs/ra/share/eichner/pwm_procaryotes/pwms/";
		
		species_parser.writeSpeciesFile(indir + "prodoric.pwms", indir + "prodoric.species");
		
		/*
		String infile = "/home/jei/SABINE/mappings/SpeciesMapping.old";
		species_parser.parseSpeciesMapping(infile);
		
		String species_file = "/home/jei/SABINE/trainingsets/new_phylogenetic_distances.out";
		species_parser.filterRelevantSpecies(species_file);
		
		String outfile = "/home/jei/SABINE/mappings/SpeciesMapping.txt";
		species_parser.write_mapping(outfile);
		*/
	}

}


