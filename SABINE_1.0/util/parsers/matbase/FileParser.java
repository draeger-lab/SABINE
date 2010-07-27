package util.parsers.matbase;
/*
 * ===============================================
 * (C) Florian Topf, University of Tuebingen, 2010
 * ===============================================
 */

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;

/*
 * File Parser (inherits OrganismListParser)
 */

public class FileParser {
	
	Filter filter = new Filter();
	
	// parse organism_list.txt as lowercase
	public ArrayList<String> parseOrgFile(String infile) {
		
		System.out.println("Parsing: "+infile);
		
		ArrayList<String> organisms = new ArrayList<String>();
		
		String line = null, curr_organism = null;
		
		try {
			 
			 BufferedReader br = new BufferedReader(new FileReader(infile));
			 line = br.readLine();
			 
			 while (line != null) {
				 				  
				 // get organisms
				 curr_organism = line.trim();
				 organisms.add(filter.filtername(curr_organism.toLowerCase()));
				 
				 line = br.readLine();
			}
			 			 
			 br.close();
		
		}
		
		catch(IOException ioe) {
			System.out.println(ioe.getMessage());
			System.out.println("IOException occurred while parsing...");
		}
		
		return organisms;
		
	}

	
	// parse organism_map.txt file
	public StringMap parseSpecFile(String infile) {
		
		System.out.println("Parsing: " +infile);
		
		StringMap hm = new StringMap(true);
			
			String line = null;
		
			try {
				 
				 BufferedReader br = new BufferedReader(new FileReader(infile));
				 line = br.readLine();
				 
				 while (line != null) {
					 
					 String[] stray = line.split("\t");
					 
					 //System.out.println(stray[0] + "\t" + stray[1]);					 
					 
					 hm.put(stray[0].trim(), stray[1].trim());
					 
					 line = br.readLine();
					 
				 }
				 
				 br.close();
				 
			} 						
		
			catch(IOException ioe) {
				System.out.println(ioe.getMessage());
				System.out.println("IOException occurred while parsing...");
			}
			
	return hm;
			
	}
	
	
	public HashMap<String, String> annoteHM(String infile, HashMap<String, String> transHM) {
		
		int count = 0;
		try {
			
			BufferedReader br = new BufferedReader(new FileReader(infile));
			String line = br.readLine();
		
			 while (line != null) {
				 
				 line = line.trim();
				 
				 if (line.isEmpty() != true) {
					 
					 String[] tmp = line.split("\t");
					 
					 transHM.put(tmp[1].trim().toUpperCase(), tmp[0].trim());
				     
				     count++;
				 }
				 line = br.readLine();
			 }
			
		 
		} catch (IOException ioe) {
			ioe.printStackTrace();
		}
		
		System.out.println(count+" manual annotations added to transfac-mapping.");
		return transHM;
		
	}
	
	
}
