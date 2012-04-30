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

package core;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.StringTokenizer;

public class SpeciesFeatureCalculator {
	
	
	ArrayList<String> species = null;
	
	double[][] distances = null;
	
	
	public void calculatePhylogeneticDistances(String tfname, String inputspecies, String class_id,  boolean[] irrelevantPairs, String lookupfile, String outfile, String train_dir) {
		
		BufferedReader br = null;
		
		ArrayList<String> ordered_species = new ArrayList<String>();
		 
		ArrayList<String> ordered_names = new ArrayList<String>();
		
		
		try {
			 
			String line = null;
			
			
			/*
			 * 
			 * parse species of tfs of the training set
			 * 
			 */
			
			br = new BufferedReader(new FileReader(new File(train_dir + "trainingset_" + class_id + ".species")));
			
			int line_count = 0;
			
			while((line = br.readLine()) != null) {
				
				if(!irrelevantPairs[line_count]) {
					
					ordered_species.add(line.trim());
					
				}
				
				line_count++;
				
			}
			
			br.close();	
			
			
			/*
			 * 
			 * parse names of tfs of the training set
			 * 
			 */
			
			br = new BufferedReader(new FileReader(new File(train_dir + "trainingset_" + class_id + ".tfnames")));
			
			line_count = 0;
			
			while((line = br.readLine()) != null) {
				
				if(!irrelevantPairs[line_count]) {
					
					ordered_names.add(line.trim());
					
				}
				
				line_count++;
				
			}
			
			br.close();	
			
			
			/*
			 * 
			 * init phylogenetic distance matrix for look-up
			 *
			 */
			
			br = new BufferedReader(new FileReader(new File(lookupfile)));
			
			line = br.readLine();
		
			StringTokenizer strtok = new StringTokenizer(line, ",");
			
			species = new ArrayList<String>();
			
			while(strtok.hasMoreTokens()) {
				
				species.add(strtok.nextToken().trim());
				
			}
			
			distances = new double[species.size()][species.size()];
			
			for(int i=0; i<distances.length; i++) {
				
				strtok = new StringTokenizer(br.readLine());
				
				for(int j=0; j<distances.length; j++) {
					
					distances[i][j] = Double.parseDouble(strtok.nextToken());
					
				}
				
			}
			
			br.close();
			
			
			BufferedWriter bw = new BufferedWriter(new FileWriter(new File(outfile)));
			
			for(int i=0; i<ordered_names.size(); i++) {
				
				bw.write(tfname + " vs. " + ordered_names.get(i) + " : " + getDistance(inputspecies, ordered_species.get(i)) + "\n");
				
			}
			
			bw.flush();
			bw.close();
			 
		}
		catch(IOException ioe) {
			System.out.println(ioe.getMessage());
			System.out.println("IOException occurred while looking up phylogenetic distances.");
		}
		
	}
	
	
public void calculateAllPhylogeneticDistances(String class_id,  boolean[] irrelevantPairs, String lookupfile, String outfile, ArrayList<String> relevant_pairs, String train_dir) {
		
		BufferedReader br = null;
		
		ArrayList<String> ordered_species = new ArrayList<String>();
		 
		ArrayList<String> ordered_names = new ArrayList<String>();
		
		
		try {
			 
			String line = null;
			
			
			/*
			 * 
			 * parse species of tfs of the training set
			 * 
			 */
			
			br = new BufferedReader(new FileReader(new File(train_dir + "trainingset_" + class_id + ".species")));
			
			int line_count = 0;
			
			while((line = br.readLine()) != null) {
				
				if(!irrelevantPairs[line_count]) {
					
					ordered_species.add(line.trim());
					
				}
				
				line_count++;
				
			}
			
			br.close();	
			
			
			/*
			 * 
			 * parse names of tfs of the training set
			 * 
			 */
			
			br = new BufferedReader(new FileReader(new File(train_dir + "trainingset_" + class_id + ".tfnames")));
			
			line_count = 0;
			
			while((line = br.readLine()) != null) {
				
				if(!irrelevantPairs[line_count]) {
					
					ordered_names.add(line.trim());
					
				}
				
				line_count++;
				
			}
			
			br.close();	
			
			
			/*
			 * 
			 * init phylogenetic distance matrix for look-up
			 *
			 */
			
			br = new BufferedReader(new FileReader(new File(lookupfile)));
			
			line = br.readLine();
		
			StringTokenizer strtok = new StringTokenizer(line, ",");
			
			species = new ArrayList<String>();
			
			while(strtok.hasMoreTokens()) {
				
				species.add(strtok.nextToken().trim());
				
			}
			
			distances = new double[species.size()][species.size()];
			
			for(int i=0; i<distances.length; i++) {
				
				strtok = new StringTokenizer(br.readLine());
				
				for(int j=0; j<distances.length; j++) {
					
					distances[i][j] = Double.parseDouble(strtok.nextToken());
					
				}
				
			}
			br.close();
			
			
			BufferedWriter bw = new BufferedWriter(new FileWriter(new File(outfile)));
			
			for(int i=0; i<ordered_names.size()-1; i++) {
				
				for(int j=i+1; j<ordered_names.size(); j++) {
					
					if (! relevant_pairs.contains(i + "_" + j)) continue;
					
					bw.write(ordered_names.get(i) + " vs. " + ordered_names.get(j) + " : " + getDistance(ordered_species.get(i), ordered_species.get(j)) + "\n");
				
				}
			}
			
			bw.flush();
			bw.close();
			 
		}
		catch(IOException ioe) {
			System.out.println(ioe.getMessage());
			System.out.println("IOException occurred while looking up phylogenetic distances.");
		}
		
	}
	
	/*
	 * 
	 * returns the distance for a given pair of species
	 * 
	 */
	public double getDistance(String species1, String species2) {
		
		if(species == null || distances == null) {
			
			System.out.println("List of species and distance matrix must be initialized before using this function. Aborting.");
			System.exit(0);
			
		}
		
		species1 = species1.replace("(", "").replace(")", "").toUpperCase();
		species2 = species2.replace("(", "").replace(")", "").toUpperCase();
		
		int index1 = species.indexOf(species1);
		int index2 = species.indexOf(species2);
		
		// System.out.println(species1 + " : " + index1);
		// System.out.println(species2 + " : " + index2 + "\n");
		
		return distances[index1][index2];
		
	}
	
}

