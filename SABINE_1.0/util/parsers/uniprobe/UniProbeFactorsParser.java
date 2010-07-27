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
import java.util.ArrayList;
import java.util.regex.Pattern;
import java.util.Arrays;

public class UniProbeFactorsParser {
	
	static ArrayList<String> Pnames = new ArrayList<String>();
	static ArrayList<String> UniProtIDs = new ArrayList<String>();
	// 2D-ArrayList for Domains
	static ArrayList<ArrayList<String>> Domainholder = new ArrayList<ArrayList<String>>();
	static ArrayList<String> Species = new ArrayList<String>();

	public void parseFile(File infile) {
		
		String line = null, curr_Prot = null, curr_ID = null, curr_species = null;
		
		final String STRPATTERN = "\t";
		Pattern p = Pattern.compile(STRPATTERN);
		
		// Sub-pattern for Domains
		final String SUBPATTERN = ",";
		Pattern s = Pattern.compile(SUBPATTERN);
						
		try {
			 
			 BufferedReader br = new BufferedReader(new FileReader(infile));
			 line = br.readLine(); // skip first line
			 line = br.readLine();
			 
			 while (line != null) {
				 				  
				 String[] strData = p.split(line);
				 
				 // get protein names
				 curr_Prot = strData[0].trim();
				 //System.out.println(curr_Prot);
				 Pnames.add(curr_Prot);
				 
				 // get UniProtID
				 curr_ID = strData[1].trim();
				 //System.out.println(curr_ID);
				 UniProtIDs.add(curr_ID);
			 
				 // get Domains
				 String[] strSub = s.split(strData[2].trim(), 0);
				 // trimming
				 for (int i = 0 ; i < strSub.length; i++) {
					 strSub[i] = strSub[i].trim();
					 // only grab the first string, not the extension (HTH APSES-type)
					 if (strSub[i].contains(" ")) {
						 String[] subsub = strSub[i].split(" ");
						 strSub[i] = subsub[0];
					 }

				 }
				 ArrayList<String> arrSub = new ArrayList<String>(Arrays.asList(strSub));
				 Domainholder.add(arrSub);
				 
				 // get Species
				 curr_species = strData[3].trim();
				 //System.out.println(curr_species);
				 Species.add(curr_species);
				 
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
		
		// uniprobe_factors.txt
		File infile = new File(args[0]);
		
		UniProbeFactorsParser uniprobe_parser = new UniProbeFactorsParser();
		
		System.out.println("Parsing: "+args[0]);
		
		uniprobe_parser.parseFile(infile);
		
		System.out.println("SUCCESS!!");
		
	}
}
