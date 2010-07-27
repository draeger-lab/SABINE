package util.parsers.uniprobe;
/*
 * ===============================================
 * (C) Florian Topf, University of Tuebingen, 2010
 * ===============================================
 */

/*
 * Parser for Organism List
 */

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;

public class OrganismListParser {
	
	static ArrayList<String> Organism = new ArrayList<String>();

	public void parseFile(File infile) {
		
		String line = null, curr_organism = null;
		
		try {
			 
			 BufferedReader br = new BufferedReader(new FileReader(infile));
			 line = br.readLine();
			 
			 while (line != null) {
				 				  
				 // get organisms
				 curr_organism = line.trim();
				 Organism.add(curr_organism);
				 
				 line = br.readLine();
			}
			 			 
			 br.close();
		
		}
		
		catch(IOException ioe) {
			System.out.println(ioe.getMessage());
			System.out.println("IOException occurred while parsing...");
		}
		
	}
	
	public static void main(String[] args) throws Exception {
		
		// organism_list.txt
		File infile = new File(args[0]);
		
		OrganismListParser organism_parser = new OrganismListParser();
		
		System.out.println("Parsing: "+args[0]);
		
		organism_parser.parseFile(infile);
		
		System.out.println("SUCCESS!!");

    	
	}

}
