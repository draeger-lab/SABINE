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
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class TransfacClassificationParser {

	// Hash Map for name -> number lookup
	HashMap<String,String> hm = new HashMap<String,String>();
	
	public void parseFile(File infile) {
		
		String line = null;
	
		try {
			 
			 BufferedReader br = new BufferedReader(new FileReader(infile));
			 line = br.readLine();
			 
			 final String STRPATTERN = " ";
			 Pattern p = Pattern.compile(STRPATTERN);
			 
			 Pattern re = Pattern.compile("^[0-4][.][0-9]+");
			 			 
			 while (line != null) {
				 
				 String key = null;
				 String val = null;
				 
				 Matcher ma = re.matcher(line.trim());
				 
				 if (ma.find()) {
					 
					 String[] strData = p.split(line);
					 
					 for (int i = 0 ; i < strData.length; i++) {
						 
						 strData[i] = strData[i].trim();
						 Matcher sub = re.matcher(strData[i]);
						 
						 if (sub.find()) {
	 
							 if ((strData.length == i+1)  != true ) {
								 
								 Matcher subsub = re.matcher(strData[i+1]);
								 if (subsub.find() != true ) {
									
									if (strData[i].startsWith("0") || strData[i].startsWith("1") || strData[i].startsWith("2") || strData[i].startsWith("3") || strData[i].startsWith("4")) {
										
										if (strData[i+1].contains("Class")) {
											
											if (strData[strData.length-1].trim().contains("factor") || strData[strData.length-1].trim().contains("domain")) {
												
												//System.out.println(strData[i].trim()+" "+strData[strData.length-2].trim());
												val = strData[i].trim();
												key = strData[strData.length-2].trim();
												
											} else {
												
												//System.out.println(strData[i].trim()+" "+strData[strData.length-1].trim());
												val = strData[i].trim();
												key = removeBrackets(strData[strData.length-1].trim());
												
											}

										} else if (strData[i+1].contains("Family") || strData[i+1].contains("Subfamily")) {
											
											//System.out.println(strData[i].trim()+" "+strData[i+2].trim());
											val = strData[i].trim();
											key = removeLike(removeBrackets(strData[i+2].trim()));
											
										} else {
											
											//System.out.println(strData[i].trim()+" "+strData[i+1].trim());
											val = strData[i].trim();
											key = strData[i+1].trim() ;
											
										}
										
									}
									
								 } 
							 } 	
						 }
					 }
				 }
				
				if (key != null && val != null) {
					System.out.println(val+" "+key);
					hm.put(key, val);
				}

				line = br.readLine();
				 
			 	}
			 			 
			 br.close();
		
		}
		
		catch(IOException ioe) {
			System.out.println(ioe.getMessage());
			System.out.println("IOException occurred while parsing...");
		}

	}
	
	private String removeBrackets(String in) {
		String out = null;
		
		out = in.replaceAll("[()]", "");
		out = out.replaceAll("[.]", "");
		
		return out;
	}
	
	private String removeLike(String in) {
		String out = null;
		
		out = in.replaceAll("-like", "");
				
		return out;
	}
	
	
	public HashMap<String, String> gethm() {
		return hm;
	}
	
	public static void main(String[] args) throws Exception {
		
		File infile = new File(args[0]);
		
		TransfacClassificationParser transfac_parser = new TransfacClassificationParser();
		
		System.out.println("Parsing: "+args[0]);
		
		transfac_parser.parseFile(infile);
		
	}
}
