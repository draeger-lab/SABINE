package util.parsers.matbase;
/*
 * ===============================================
 * (C) Florian Topf, University of Tuebingen, 2010
 * ===============================================
 */

import java.io.BufferedReader;
import java.io.ByteArrayInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.StringTokenizer;

/*
 * Converter UniProtApiResults -> Java
 */

public class UniProtResultRefiner {
	
	
	// process IDMapping to HashMap (dirty string bashing)
	public HashMap<String, ArrayList<String>> processIDs(String result) {
		
		String[] lines = result.split("\n");
		
		ArrayList<String> from = new ArrayList<String>();
		ArrayList<String> to = new ArrayList<String>();
		
		// convert result to 2 arraylists
		for (int i = 1; i < lines.length; i++) {
			
			String[] curr = lines[i].split("\t");
			
			from.add(curr[0]);
			to.add(curr[1]);
		}
		
		HashMap<String, ArrayList<String>> idHM = new HashMap<String, ArrayList<String>>();
		
		// assemble hashmap-entries
		for (int i = 0; i < from.size(); i++) {
			
			ArrayList<String> values = new ArrayList<String>();
			
			String key = from.get(i);
			
			if (idHM.containsKey(key)) {
				// add entry to arraylist from hashmap
				values = idHM.get(key);
				values.add(to.get(i));
				idHM.put(key, values);
			} else {
				// add new entry and put it on hashmap
				values.add(to.get(i));
				idHM.put(key, values);
			}
			
		}
		
		return idHM;
		
	}
	
	
	// process Sequences to HashMap
	public HashMap<String, String> processSeqs(String result) {
		
		String curr_seq, curr_ID, line;
		StringTokenizer strtok;
		
		HashMap<String,String> seqHM = new HashMap<String,String>();
		
		try {
			// convert string to bufferedreader
			InputStream is = new ByteArrayInputStream(result.getBytes("UTF-8"));
			InputStreamReader isr = new InputStreamReader(is);
			BufferedReader br = new BufferedReader(isr);
			
			line = br.readLine();
			
			while (line != null) {
				
				if (! line.startsWith(">")) {
					System.out.println("Parse error. \">\" expected at the beginning of the line.");
				}
			
				// parse uniprot ID
				strtok = new StringTokenizer(line, "|");
				strtok.nextToken();
				
				curr_ID = strtok.nextToken().trim();
			
				// parse sequence
				curr_seq = "";
				
				// while condition changed due to null-pointer exception...
				while ((line=br.readLine()) != null && !line.startsWith(">")) {
					curr_seq += line.trim();
				}
				// put it on!
				seqHM.put(curr_ID, curr_seq);
			}
		}	 
		catch(IOException ioe) {
			System.out.println(ioe.getMessage());
			System.out.println("IOException occurred while parsing sequences.");
		}
				
		return seqHM;
		
	}
	
	
	//  process Bindings to HashMap
	public HashMap<String, ArrayList<String>> processBinds(String result) {
		
		HashMap<String,ArrayList<String>> bindsHM = new HashMap<String,ArrayList<String>>();
		
		ArrayList<String> curr_domains;
		String line;
		StringTokenizer strtok;
		String curr_ID;
		
		try {
			// convert string to bufferedreader
			InputStream is = new ByteArrayInputStream(result.getBytes("UTF-8"));
			InputStreamReader isr = new InputStreamReader(is);
			BufferedReader br = new BufferedReader(isr);
			
			line = br.readLine();		// skip first line 
			line = br.readLine();
			
			while (line != null) {
				
				if (! line.startsWith("##sequence-region")) { 
					System.out.println("Parse Error. \"##sequence-region\" expected at the beginning of the line.");
				}
				
				strtok = new StringTokenizer(line);
				strtok.nextToken();
				curr_ID = strtok.nextToken().trim();

				curr_domains = new ArrayList<String>();
				
				// while condition adapted to prevent null-pointer exception...
				while ((line=br.readLine()) != null && !line.startsWith("##sequence-region")) {
					
					if (line.contains("DNA binding")) {
						
						// parse domain 
						curr_domains.add(line.split("\t")[3].trim() + "\t" + line.split("\t")[4].trim());
						
					}
				}
				bindsHM.put(curr_ID, curr_domains);
				// debugging
				//System.out.println(curr_ID +" "+ curr_domains);
			}

		}
		
		catch(IOException ioe) {
			System.out.println(ioe.getMessage());
			System.out.println("IOException occurred while parsing DBDs.");
		}
		
		return bindsHM;
		
	}
	
}
