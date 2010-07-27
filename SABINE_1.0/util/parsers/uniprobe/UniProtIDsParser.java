package util.parsers.uniprobe;
/*
 * ===============================================
 * (C) Florian Topf, University of Tuebingen, 2010
 * ===============================================
 */

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.StringTokenizer;

public class UniProtIDsParser {

	// Hash Map for secondary -> primary lookup
	HashMap<String,String> hm = new HashMap<String,String>();
	
	public HashMap<String, String> parseFile(File infile) {
		
		String line = null;
	
		try {
			 
			 BufferedReader br = new BufferedReader(new FileReader(infile));
			 line = br.readLine();
			 			 
			 StringTokenizer strtok;
			 
			 while (line != null) {
				 
				 strtok = new StringTokenizer(line);
				 
				 if (strtok.countTokens() == 2) {
					 
					 //System.out.println(line);
					 
					 if (line.contains(":")) {
						
					 } else {
						 
						 hm.put(strtok.nextToken().toUpperCase().trim(), strtok.nextToken().toUpperCase().trim());
						 
					 }
					 
				 } else {
					 
					 line = br.readLine();
					 strtok = new StringTokenizer(line);
					 
				 }
				 
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
	
	public static void main(String[] args) throws Exception {
		
		File infile = new File(args[0]);
		
		UniProtIDsParser uniprot_parser = new UniProtIDsParser();
		
		System.out.println("Parsing: "+args[0]);
		
		uniprot_parser.parseFile(infile);
		
		System.out.println("SUCCESS!!");
		
	}
}
