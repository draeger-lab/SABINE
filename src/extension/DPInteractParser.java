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
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.ArrayList;

public class DPInteractParser {
	
	ArrayList<String> tf_names = new ArrayList<String>();
	ArrayList<ArrayList<String>> sequences = new ArrayList<ArrayList<String>>();
	ArrayList<String[]> pwms = new ArrayList<String[]>();
 	
	public void parseNamesAndSequences(String infile) {
		
		String line;
		ArrayList<String> curr_seqs;
		
		try {
			 
			 BufferedReader br = new BufferedReader(new FileReader(new File(infile)));
			 
			 while((line = br.readLine()) != null) {
					
				while (! line.startsWith("*")) {
						line = br.readLine();
				}
				 
				if (! line.startsWith("*")) {
					System.out.println("Parse Error. \"'*\" expected at the beginning of the line.\nLine: " + line);
					System.exit(0);
				}
				
				line = br.readLine();
				tf_names.add(line.substring(0, line.indexOf(".dat")));
				
				curr_seqs = new ArrayList<String>();
				br.readLine();        // ***
				br.readLine();        //
				while ((line = br.readLine()) != null && line.length() > 0) {      // > ... (fasta header)
					line = br.readLine();										   // binding site
					curr_seqs.add(line.trim().toUpperCase());					   
				}
				sequences.add(curr_seqs);
			 }
			 br.close();
			 System.out.println("Number of matrices: " + tf_names.size());
		}
		catch(IOException ioe) {
			System.out.println(ioe.getMessage());
			System.out.println("IOException occurred while parsing matrices from DPInteract.");
		}
	}
	
	public void compute_pwms() {
		
		String[] curr_pwm;
		String curr_line;
		int num_seqs = 0;
		double[][] base_counter = null;
		
		DecimalFormat fmt = new DecimalFormat();
		DecimalFormatSymbols symbs = fmt.getDecimalFormatSymbols();
        symbs.setDecimalSeparator('.');
        fmt.setDecimalFormatSymbols(symbs);
        fmt.setMaximumFractionDigits(4);
        fmt.setMinimumFractionDigits(4);
		
		for (int i=0; i<sequences.size(); i++) {
			
			base_counter = new double[4][sequences.get(i).get(0).length()];
			num_seqs = sequences.get(i).size();
			
			for (int j=0; j<sequences.get(i).size(); j++) {
				for (int k=0; k<sequences.get(i).get(j).length(); k++) {
					if (sequences.get(i).get(j).charAt(k) == 'A') {
						base_counter[0][k]++;
					}
					else if (sequences.get(i).get(j).charAt(k) == 'C') {
						base_counter[1][k]++;
					}
					else if (sequences.get(i).get(j).charAt(k) == 'G') {
						base_counter[2][k]++;
					}
					else if (sequences.get(i).get(j).charAt(k) == 'T') {
						base_counter[3][k]++;
					}
					else {
						System.out.println("Parse Error. A,C,G or T expected.\n");
						System.exit(0);
					}
				}
				
			}
			curr_pwm = new String[4];
			
			for (int j=0; j<4; j++) {
				curr_line = "";
				
				for (int k=0; k<base_counter[0].length; k++) {
					curr_line += fmt.format(base_counter[j][k] / num_seqs) + "   ";
				}
				curr_pwm[j] = curr_line;
			}
			pwms.add(curr_pwm);
		}
	}
	
	public void write_outfile(String outfile) {
		
		try {
			
			 BufferedWriter bw = new BufferedWriter(new FileWriter(new File(outfile)));
		
			 for (int i=0; i<tf_names.size(); i++) {
				 
				 bw.write("NA  " + tf_names.get(i) + "\n");
				 bw.write("XX  \n");
				 bw.write("SP  Escherichia coli\n");
				 bw.write("XX  \n");
				 bw.write("MA  " + pwms.get(i)[0] + "\n");
				 bw.write("MA  " + pwms.get(i)[1] + "\n");
				 bw.write("MA  " + pwms.get(i)[2] + "\n");
				 bw.write("MA  " + pwms.get(i)[3] + "\n");
				 bw.write("XX  \n");
				 bw.write("//  \n");
				 bw.write("XX  \n");
			 }
			 bw.flush();
			 bw.close();
		}
		
		
		catch(IOException ioe) {
			System.out.println(ioe.getMessage());
			System.out.println("IOException occurred while writing matrices from DPInteract to output file.");
		}
	}
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		
		String infile = "/afs/cs/ra/share/eichner/pwm_procaryotes/raw_data/all_matrices.dat";
		String outfile = "/afs/cs/ra/share/eichner/pwm_procaryotes/pwms/dp_interact.pwms";
		
		DPInteractParser dp_parser = new DPInteractParser();
		
		dp_parser.parseNamesAndSequences(infile);
		dp_parser.compute_pwms();
		dp_parser.write_outfile(outfile);		
	}
}

