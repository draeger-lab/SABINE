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

public class ProdoricParser {
	
	ArrayList<String> tf_names = new ArrayList<String>();
	ArrayList<String> unique_names = new ArrayList<String>();
	ArrayList<String> organisms = new ArrayList<String>();
	ArrayList<ArrayList<int[]>> matrices = new ArrayList<ArrayList<int[]>>();
	ArrayList<String[]> pwms = new ArrayList<String[]>();
	ArrayList<Integer> num_sites = new ArrayList<Integer>();
	
	private void parseMatrices(String infile) {
		
		String line, curr_name;
		StringTokenizer strtok;
		ArrayList<int[]> curr_matrix = new ArrayList<int[]>();
		int[] base_counter;
		int site_counter = 0;
		
		boolean name_parsed = false;
		boolean species_parsed = false;
		boolean matrix_parsed = false;
		
		try {
			BufferedReader br = new BufferedReader(new FileReader(new File(infile)));

			
			while ((line = br.readLine()) != null) {
				
				if (line.startsWith("NA")) {
					tf_names.add(line.substring(3).trim());
					name_parsed = true;
					
					curr_name = new StringTokenizer(line.substring(3)).nextToken();
					if (! unique_names.contains(curr_name)) {
						unique_names.add(curr_name);
					}
				}
				if (line.startsWith("BF")) {
					organisms.add(line.substring(line.indexOf("Species:")+9));
					species_parsed = true;
				}
				if (line.startsWith("P0")) {
					
					curr_matrix = new ArrayList<int[]>();
					while (! (line = br.readLine()).startsWith("XX")) {
						strtok = new StringTokenizer(line.substring(3));
						base_counter = new int[4];
						
						base_counter[0] = Integer.parseInt(strtok.nextToken());
						base_counter[1] = Integer.parseInt(strtok.nextToken());
						base_counter[2] = Integer.parseInt(strtok.nextToken());
						base_counter[3] = Integer.parseInt(strtok.nextToken());
						
						curr_matrix.add(base_counter);
					}
					matrices.add(curr_matrix);
					matrix_parsed = true;
				}
				
				// count number of binding sites
				if (line.startsWith("BS")) {
					site_counter++;
				}
				num_sites.add(site_counter);
				site_counter = 0;
				
				if (line.startsWith("//")) {
					if (name_parsed && species_parsed && matrix_parsed) {
						name_parsed = species_parsed = matrix_parsed = false;
					}
					else {
						System.out.println("Parse Error. Error occurred while parsing matrices from PRODORIC.\n");
						System.out.println("Name:   " + tf_names.get(tf_names.size()-1) + "  (" + name_parsed + ")");
						System.exit(0);
					}
				}
			}
			br.close();
			//System.out.println(unique_names.size() + " factor names found.");
		}
		catch(IOException ioe) {
			System.out.println(ioe.getMessage());
			System.out.println("IOException occurred while parsing matrices from PRODORIC.");
		}
	}	
	
	public void computePWMs() {
		
		DecimalFormat fmt = new DecimalFormat();
		DecimalFormatSymbols symbs = fmt.getDecimalFormatSymbols();
        symbs.setDecimalSeparator('.');
        fmt.setDecimalFormatSymbols(symbs);
        fmt.setMaximumFractionDigits(4);
        fmt.setMinimumFractionDigits(4);
        
		String[] curr_pwm;
		double sum;
			
		for (int i=0; i<matrices.size(); i++) {
			curr_pwm = new String[] {"", "", "", ""};
			
			for (int j=0; j<matrices.get(i).size(); j++) {
				sum = matrices.get(i).get(j)[0] + matrices.get(i).get(j)[1] + matrices.get(i).get(j)[2] + matrices.get(i).get(j)[3];
				curr_pwm[0] += fmt.format(matrices.get(i).get(j)[0] / sum) + "   ";
				curr_pwm[1] += fmt.format(matrices.get(i).get(j)[1] / sum) + "   ";
				curr_pwm[2] += fmt.format(matrices.get(i).get(j)[2] / sum) + "   ";
				curr_pwm[3] += fmt.format(matrices.get(i).get(j)[3] / sum) + "   ";
			}
			pwms.add(curr_pwm);
		}
	}
	
	private void writeOutfile(String outfile) {
		
		RegTransBaseParser rtb_parser = new RegTransBaseParser();
		rtb_parser.tf_names = tf_names;
		rtb_parser.organisms = organisms;
		rtb_parser.pwms = pwms;
		rtb_parser.write_outfile(outfile);
	}
	
	
	private void filterMultipleEntrys() {
		
		String curr_name, next_name;
		int curr_max, max_idx;
		boolean[] relevant_entrys = new boolean[tf_names.size()];
		
		int i = 0;
		while (i < tf_names.size()-1) {
			
			curr_name = new StringTokenizer(tf_names.get(i)).nextToken();
			next_name = new StringTokenizer(tf_names.get(i+1)).nextToken();
			
			if (curr_name.equals(next_name)) {
				
				curr_max = num_sites.get(i);
				max_idx = i;
				if (num_sites.get(++i) > curr_max) {
					curr_max = num_sites.get(i);
					max_idx = i;
				}
				
				next_name = new StringTokenizer(tf_names.get(i+1)).nextToken();
				while (next_name.equals(curr_name) && i<tf_names.size()-2) {
				
					if (num_sites.get(++i) > curr_max) {
						curr_max = num_sites.get(i);
						max_idx = i;
					}
					next_name = new StringTokenizer(tf_names.get(i+1)).nextToken();
				}
				relevant_entrys[max_idx] = true;
			}
			else {
				relevant_entrys[i] = true;
			}
			i++;
		}
		
		int tf_counter = 0;
		for (int j=tf_names.size()-1; j >= 0; j--) {
			
			if (! relevant_entrys[j]) {
				
				tf_names.remove(j);
				organisms.remove(j);
				matrices.remove(j);
				pwms.remove(j);
				num_sites.remove(j);
			}
			else {
				tf_counter++;
			}
		}
		// System.out.println(tf_counter + " factors filtered.");
	}
	
	
	public static void main(String[] args) {
		
		//String path = "/afs/cs/ra/";
		String path = "S:/";
		
		String infile = path + "share/eichner/pwm_procaryotes/raw_data/PRODORIC8.9.txt";
		String outfile = path + "share/eichner/pwm_procaryotes/pwms/prodoric2.pwms";
		
		ProdoricParser prod_parser = new ProdoricParser();
		prod_parser.parseMatrices(infile);
		prod_parser.computePWMs();
		prod_parser.filterMultipleEntrys();
		prod_parser.writeOutfile(outfile);
	}
}

