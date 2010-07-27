package util.parsers.matbase;
/*
 * ===============================================
 * (C) Florian Topf, University of Tuebingen, 2010
 * ===============================================
 */

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.StringTokenizer;

/*
 * OutputFileWriter to conform to SABINE-InputFileFormat
 */

public class OutputFileWriter {
	
	// This function was taken and adapted from SABINEInputFileWriter.java
	public void writeOutfile(ArrayList<String> symbols, ArrayList<ArrayList<String>> uniprot_ids, ArrayList<String> species, ArrayList<ArrayList<ArrayList<String>>> matrices, ArrayList<ArrayList<String>> sequences, ArrayList<ArrayList<String>> bindings, ArrayList<String> transfac, boolean[] flags, String outfile) {
		
		// statistics...
		int count = symbols.size();
		
		int SEQLINELENGTH = 60;
		StringTokenizer strtok;
		
		try {
			BufferedWriter bw = new BufferedWriter(new FileWriter(new File(outfile)));
			
			//for (int i = 0; i < 5; i++) {
			for (int i=0; i<symbols.size(); i++) {
				
				// skip missing entries
				//if (uniprot_ids.get(i).isEmpty() || matrices.get(i).isEmpty() || sequences.get(i).isEmpty() || bindings.get(i).isEmpty() || transfac.get(i).isEmpty()) {
				//	count--;
				//	continue;
				//}
				
				strtok = new StringTokenizer(species.get(i));
				
				bw.write("NA  " + new StringTokenizer(symbols.get(i)).nextToken() + "\n");
				bw.write("XX  \n");
				bw.write("SP  " + strtok.nextToken() + " " + strtok.nextToken() + "\n");
				bw.write("XX  \n");
				
				// write multiple IDs
				Iterator<String> itu = uniprot_ids.get(i).iterator();
			    while (itu.hasNext()) {
			    	String element = itu.next();
			    	bw.write("RF  " + element + "\n");
			    }

			    bw.write("XX  \n");
			    
				if (!flags[i]) {
					bw.write("CL  "+ convert2full(transfac.get(i)) +"\n");
					bw.write("XX  \n");
				}
				
				// paste IP for flagged entries
				else {
					bw.write("IP  "+"\n");
					bw.write("XX  \n");
				}
				
				// write sequence
				Iterator<String> its = sequences.get(i).iterator();
				int num = 1;
			    while (its.hasNext()) {
			    	String element = its.next();
			    	
			    	for (int j=0; j<(element.length()/SEQLINELENGTH); j++) {
			    		bw.write("S" + num + "  ");
			    		bw.write(element.toUpperCase(), j*SEQLINELENGTH, SEQLINELENGTH);
			    		bw.write("\n");
			    	}
			    	if (element.length()-(element.length()/SEQLINELENGTH)*SEQLINELENGTH > 0) {
			    		bw.write("S" + num + "  ");
			    		bw.write(element.toUpperCase(), (element.length()/SEQLINELENGTH)*SEQLINELENGTH, element.length()-(element.length()/SEQLINELENGTH)*SEQLINELENGTH);
						bw.write("\n");
			    	}
			    	bw.write("XX  \n");
			    	num++;
			    }
				
			    
				// write DNAbinding
				String tab = "\t";
				
			    if (!flags[i]) {
					
					for (int j=0; j<bindings.get(i).size(); j++) {
						
						String curr = bindings.get(i).get(j);
						String[] arr = curr.split(";");
						
						bw.write("FT  " + arr[0].trim() + "    " + arr[1].trim().replaceAll(tab, "   ") + "\n");
					}
					
					bw.write("XX  \n");
			    	
			    }

			    
				// write matrices
				ArrayList<ArrayList<String>> currmat = matrices.get(i);
				
				
				for (int j = 0; j < currmat.size(); j++) {
					
					bw.write("MN  " + "MatBase" + "\n");
					bw.write("XX\n");
					
					ArrayList<String> tmpmat = currmat.get(j);
					
					for (int k=0; k < 4; k++) {
						bw.write("MA  " + tmpmat.get(k).trim().replaceAll(tab, "   ") + "\n");
					}
					bw.write("XX  \n");
					
				}
				
				bw.write("//  \n");
				bw.write("XX  \n");
			}
			bw.flush();
			bw.close();
		}
		catch(IOException ioe) {
			System.out.println(ioe.getMessage());
			System.out.println("IOException occurred while writing output file");
		}
		
		System.out.println(count + " entrie(s) written.");
		
	}


	// convert tranfac-classification to SABINE-InputFileFormat
	private String convert2full(String input) {
		String output = input+".";
		
		while (output.length() <= 9) {
			output = output.concat("0.");
		}
		
		return output;
	}
	
}
