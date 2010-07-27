package util.parsers.uniprobe;
/*
 * ===============================================
 * (C) Florian Topf, University of Tuebingen, 2010
 * ===============================================
 */

/*
 * OutputFileWriter to conform to SABINE-InputFileFormat
 */

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.StringTokenizer;

public class OutputFileWriter {
	
	// This function was taken and adapted from SABINEInputFileWriter.java
	public void writeOutfile(ArrayList<String> pnames, ArrayList<String> uniProtIDs, ArrayList<ArrayList<String>> domainholder, ArrayList<String> species, ArrayList<ArrayList<String>> pwms, ArrayList<String> transfac, ArrayList<String> sequences, ArrayList<ArrayList<String>> DNAbinding, ArrayList<String> MNs, String outfile) {
		
		int SEQLINELENGTH = 60;
		StringTokenizer strtok;
		
		try {
			BufferedWriter bw = new BufferedWriter(new FileWriter(new File(outfile)));
			
			for (int i=0; i<pnames.size(); i++) {
				
				// skip entry if domain annotation or pwm is missing
				if (domainholder.get(i).isEmpty() || pwms.get(i).isEmpty()) {
					continue;
				}
				
				strtok = new StringTokenizer(species.get(i));
				
				bw.write("NA  " + new StringTokenizer(pnames.get(i)).nextToken() + "\n");
				bw.write("XX  \n");
				bw.write("SP  " + strtok.nextToken() + " " + strtok.nextToken() + "\n");
				bw.write("XX  \n");
				bw.write("RF  " + uniProtIDs.get(i) + "\n");
				bw.write("XX  \n");
				bw.write("CL  "+ convert2full(transfac.get(i)) +"\n");
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
				
				// write DNAbinding
				
				String tab = "\t";
				
				for (int j=0; j<DNAbinding.get(i).size(); j++) {
					bw.write("FT  " + uniProtIDs.get(i) + "    " + DNAbinding.get(i).get(j).replaceAll(tab, "   ") + "\n");
				}
				
				bw.write("XX  \n");

				// write pwm-matrix
				bw.write("MN  " + filterMN(MNs.get(i)) + "\n");
				bw.write("XX\n");
				for (int j=0; j<4; j++) {
					bw.write("MA  " + pwms.get(i).get(j).substring(3).replaceAll(tab, "   ") + "\n");
				}
				
				bw.write("XX  \n");
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
	}


	// convert tranfac-classification to SABINE-InputFileFormat
	private String convert2full(String input) {
		String output = input+".";
		
		while (output.length() <= 9) {
			output = output.concat("0.");
		}
		
		return output;
	}
	
	
	// filters MN out of FileEntryString
	private String filterMN(String input) {
		String[] tmp = input.split("/");
		
		String output = tmp[1] + "_" + tmp[2].replaceAll(".pwm", "").replaceAll(".txt", "");
		
		return output;
	}


}
