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
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.ArrayList;
import java.util.StringTokenizer;

public class RegulonDBParser {
	
	ArrayList<String> tf_names = new ArrayList<String>();
	ArrayList<int[][]> matrices = new ArrayList<int[][]>();
	ArrayList<String[]> pwms = new ArrayList<String[]>();
	
	private void parseMatrices(String infile) {
		
		String line;
		StringTokenizer strtok;
		int[][] curr_matrix = null;
		int j;
		
		try {
			BufferedReader br = new BufferedReader(new FileReader(new File(infile)));

			while ((line = br.readLine()) != null) {
				
				if (line.startsWith("Transcription Factor Name")) {
					tf_names.add(line.substring(26).trim());
				}
				if (line.startsWith("Matrix")) {

					for (int i=0; i<4; i++) {
						strtok = new StringTokenizer(br.readLine().substring(2));
						
						if (i == 0) curr_matrix = new int[4][strtok.countTokens()];
						j = 0;
						while (strtok.hasMoreTokens()) {
							curr_matrix[i][j] = Integer.parseInt(strtok.nextToken());
							j++;
						}
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
		
		int[][] curr_matrix;
		String[] curr_pwm;
		double sum;
		
		DecimalFormat fmt = new DecimalFormat();
		DecimalFormatSymbols symbs = fmt.getDecimalFormatSymbols();
        symbs.setDecimalSeparator('.');
        fmt.setDecimalFormatSymbols(symbs);
        fmt.setMaximumFractionDigits(4);
        fmt.setMinimumFractionDigits(4);
        
		for (int i=0; i<matrices.size(); i++) {	
			curr_matrix = matrices.get(i);
			curr_pwm = new String[] {"","","",""};
			
			for(int j=0; j<curr_matrix[0].length; j++) {
				
				sum = curr_matrix[0][j] + curr_matrix[1][j] + curr_matrix[2][j] + curr_matrix[3][j];
				
				curr_pwm[0] += fmt.format(curr_matrix[0][j]/sum) + "   ";
				curr_pwm[1] += fmt.format(curr_matrix[1][j]/sum) + "   ";
				curr_pwm[2] += fmt.format(curr_matrix[2][j]/sum) + "   ";
				curr_pwm[3] += fmt.format(curr_matrix[3][j]/sum) + "   ";
			}
			pwms.add(curr_pwm);
		}
	}
	
	private void writeOutfile(String outfile) {
		
		DPInteractParser dp_parser = new DPInteractParser();
		dp_parser.tf_names = tf_names;
		dp_parser.pwms = pwms;
		dp_parser.write_outfile(outfile);
	}
	
	public static void main(String[] args) {
		
		String infile = "/afs/cs/ra/share/eichner/pwm_procaryotes/raw_data/Matrix_AlignmentSet.txt";
		String outfile = "/afs/cs/ra/share/eichner/pwm_procaryotes/pwms/regulondb.pwms";
		
		RegulonDBParser rdb_parser = new RegulonDBParser();
		rdb_parser.parseMatrices(infile);
		rdb_parser.computePWMs();
		rdb_parser.writeOutfile(outfile);
	}
}

