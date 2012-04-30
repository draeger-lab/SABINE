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
import java.util.ArrayList;

import help.FileFilter;

public class RegTransBaseParser {
	
	ArrayList<String> tf_names = new ArrayList<String>();
	ArrayList<String> organisms = new ArrayList<String>();
	ArrayList<ArrayList<String>> sequences = new ArrayList<ArrayList<String>>();
	ArrayList<String[]> pwms = new ArrayList<String[]>();
	
	private void parseMSA(String infile) {
		
		// parse factor name and organism from filename
		int pos = infile.substring(infile.lastIndexOf('/')).indexOf('_') + infile.lastIndexOf('/');
		
		tf_names.add(infile.substring(infile.lastIndexOf('/')+1, pos));
		organisms.add(infile.substring(pos+1, infile.lastIndexOf(".fasta")).replace('_', ' '));
		
		//System.out.println("Factor Name: " + tf_names.get(tf_names.size()-1));
		//System.out.println("Organism:    " + organisms.get(tf_names.size()-1));
		
		try {
			BufferedReader br = new BufferedReader(new FileReader(new File(infile)));
			ArrayList<String> curr_seqs = new ArrayList<String>();
			
			while (br.readLine() != null) {					// fasta header
				curr_seqs.add(br.readLine().trim().toUpperCase());		// binding site
				//System.out.println(curr_seqs.get(curr_seqs.size()-1));
			}
			sequences.add(curr_seqs);
			br.close();
		}
		catch(IOException ioe) {
			System.out.println(ioe.getMessage());
			System.out.println("IOException occurred while parsing matrices from RegTransBase.");
		}
	}
	
	private void parseAllMSAs(String data_path) {
		
		/*
		 *  list fasta files in given path
		 */		
		
		FileFilter filter = new FileFilter();
		filter.setFormat(".*.fasta");
		filter.setDirectory(data_path);
		String[] files = filter.listFilesFullPath();
		
		System.out.println(files.length + " FASTA files found.");
		
		String infile;
		
		/*
		 *  parse sequences of all MSAs
		 */		
		
		for (int i=0; i<files.length; i++) {
			infile = data_path + files[i];
			parseMSA(infile);
		}
		
		/*
		 *  compute PFMs
		 */
		
		DPInteractParser dp_parser = new DPInteractParser();
		dp_parser.tf_names = tf_names;
		dp_parser.sequences = sequences;
		dp_parser.compute_pwms();
		pwms = dp_parser.pwms;
		
	}
	
	public void write_outfile(String outfile) {
		
		try {
			
			 BufferedWriter bw = new BufferedWriter(new FileWriter(new File(outfile)));
		
			 for (int i=0; i<tf_names.size(); i++) {
				 
				 bw.write("NA  " + tf_names.get(i) + "\n");
				 bw.write("XX  \n");
				 bw.write("SP  " + organisms.get(i) + "\n");
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
			System.out.println("IOException occurred while writing matrices from RegTransBase to output file.");
		}
	}

	public static void main(String[] args) {
		
		RegTransBaseParser rtb_parser = new RegTransBaseParser();
		
		String data_path = "/afs/cs/ra/share/eichner/pwm_procaryotes/raw_data/regtransbase_alignments_download/";
		String outfile = "/afs/cs/ra/share/eichner/pwm_procaryotes/pwms/regtransbase.pwms";
		
		rtb_parser.parseAllMSAs(data_path);
		rtb_parser.write_outfile(outfile);
	}
}

