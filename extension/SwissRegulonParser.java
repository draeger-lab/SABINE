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
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.StringTokenizer;

public class SwissRegulonParser {

	ArrayList<String> tf_names = new ArrayList<String>();
	ArrayList<ArrayList<int[]>> matrices = new ArrayList<ArrayList<int[]>>();
	ArrayList<String[]> pwms = new ArrayList<String[]>();
	
	
	private void parseMatrices(String infile) {
		
		String line;
		StringTokenizer strtok;
		ArrayList<int[]> curr_matrix = new ArrayList<int[]>();
		int[] base_counter;
		
		try {
			BufferedReader br = new BufferedReader(new FileReader(new File(infile)));

			while ((line = br.readLine()) != null) {
				
				if (line.startsWith("NA")) {
					tf_names.add(line.substring(3).trim().substring(0,line.substring(3).trim().indexOf('_')));
				}
				if (line.startsWith("PO")) {
					
					curr_matrix = new ArrayList<int[]>();
					while (! (line = br.readLine()).startsWith("//")) {
						strtok = new StringTokenizer(line.substring(3));
						base_counter = new int[4];
						
						base_counter[0] = Integer.parseInt(strtok.nextToken());
						base_counter[1] = Integer.parseInt(strtok.nextToken());
						base_counter[2] = Integer.parseInt(strtok.nextToken());
						base_counter[3] = Integer.parseInt(strtok.nextToken());
						
						curr_matrix.add(base_counter);
					}
					matrices.add(curr_matrix);
				}
			}
			br.close();
		}
		catch(IOException ioe) {
			System.out.println(ioe.getMessage());
			System.out.println("IOException occurred while parsing matrices from PRODORIC.");
		}
	}	
	
	private void computePWMs() {
		
		ProdoricParser prod_parser = new ProdoricParser();
		prod_parser.tf_names = tf_names;
		prod_parser.matrices = matrices;
		prod_parser.computePWMs();
		pwms = prod_parser.pwms;
	}
	
	private void writeOutfile(String outfile) {
		
		DPInteractParser dp_parser = new DPInteractParser();
		dp_parser.tf_names = tf_names;
		dp_parser.pwms = pwms;
		dp_parser.write_outfile(outfile);
	}

	public static void main(String[] args) {
		
		SwissRegulonParser sr_parser = new SwissRegulonParser();
		
		String infile = "/afs/cs/ra/share/eichner/pwm_procaryotes/raw_data/Swiss_Regulon_e.coli_PWMs.txt";
		String outfile = "/afs/cs/ra/share/eichner/pwm_procaryotes/pwms/swiss_regulon.pwms";
		
		sr_parser.parseMatrices(infile);
		sr_parser.computePWMs();
		sr_parser.writeOutfile(outfile);
	}

}

