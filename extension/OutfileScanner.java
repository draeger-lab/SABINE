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

public class OutfileScanner {
	
	ArrayList<String> transfac_factor_ids = new ArrayList<String>();
	ArrayList<String> public_factor_ids = new ArrayList<String>();
	ArrayList<String> professional_factor_ids = new ArrayList<String>();
	
	ArrayList<String> transfac_pwm_ids = new ArrayList<String>();
	ArrayList<String> public_pwm_ids = new ArrayList<String>();
	ArrayList<String> professional_pwm_ids = new ArrayList<String>();
	
	/*
	 *  generates list of factor or matrix IDs and returns IDs as ArrayList
	 */
	
	public void generateIDList(String infile, String outfile) {
		try {
			
			BufferedReader br = new BufferedReader(new FileReader(new File(infile)));
			BufferedWriter bw = new BufferedWriter(new FileWriter(new File(outfile)));
			
			String line;
			int idx=1;
			while ((line = br.readLine()) != null) {
				if (line.startsWith(idx + ":")) {
					
					bw.write(line.substring(line.indexOf(':')+1).trim() + "\n");
					idx++;
				}
			}
			br.close();
			bw.flush();
			bw.close();
			System.out.println((idx-1) + " TRANSFAC IDs parsed.");
		}
		catch(IOException ioe) {
			System.out.println(ioe.getMessage());
			System.out.println("IOException occurred while generating TRANSFAC ID lists.");
		}
	}
	 
	public void getProfessionalIDs(String training_factor_file, String public_factorID_file, String public_matrixID_file) {
		
		TransfacParser tf_parser = new TransfacParser();
		tf_parser.parseFactors(training_factor_file);
		
		// parse all TRANSFAC factor names
		for (int i=0; i<tf_parser.tf_names.size(); i++) {
			if (!transfac_factor_ids.contains(tf_parser.tf_names.get(i)) && tf_parser.tf_names.get(i).matches("T\\d{5}")) {
				transfac_factor_ids.add(tf_parser.tf_names.get(i));
			}
		}
		
		// parse all TRANSFAC matrix names
		for (int i=0; i<tf_parser.pfm_names.size(); i++) { 
			for (int j=0; j<tf_parser.pfm_names.get(i).size(); j++) { 
				
				if (!transfac_pwm_ids.contains(tf_parser.pfm_names.get(i).get(j)) && tf_parser.pfm_names.get(i).get(j).matches("M\\d{5}")) {
					transfac_pwm_ids.add(tf_parser.pfm_names.get(i).get(j));
				}
			}
		}
		
		try {
			// parse all TRANSFAC Public factor names
			BufferedReader br = new BufferedReader(new FileReader(new File(public_factorID_file)));
			
			String line;
			while ((line = br.readLine()) != null) {
					public_factor_ids.add(line.trim());
			}
			br.close();
			
			// parse all TRANSFAC Public matrix names
			br = new BufferedReader(new FileReader(new File(public_matrixID_file)));
			while ((line = br.readLine()) != null) {
					public_pwm_ids.add(line.trim());
			}
			br.close();
			
			// get TRANSFAC Professional Factor IDs
			for (int i=0; i<transfac_factor_ids.size(); i++) {
				
				if (! public_factor_ids.contains(transfac_factor_ids.get(i))) {
					professional_factor_ids.add(transfac_factor_ids.get(i));
				}
			}
			
			// get TRANSFAC Professional Matrix IDs
			for (int i=0; i<transfac_pwm_ids.size(); i++) {
				
				if (! public_pwm_ids.contains(transfac_pwm_ids.get(i))) {
					professional_pwm_ids.add(transfac_pwm_ids.get(i));
				}
			}
		}
		catch(IOException ioe) {
			System.out.println(ioe.getMessage());
			System.out.println("IOException occurred while calling SABINE.");
		}
		
		// print statistics
		System.out.println(transfac_factor_ids.size() + " TRANSFAC Factor IDs parsed.");
		System.out.println(public_factor_ids.size() + " TRANSFAC Public Factor IDs parsed.");
		System.out.println(professional_factor_ids.size() + " TRANSFAC Professional Factor IDs parsed.");
		
		System.out.println(transfac_pwm_ids.size() + " TRANSFAC Matrix IDs parsed.");
		System.out.println(public_pwm_ids.size() + " TRANSFAC Public Matrix IDs parsed.");
		System.out.println(professional_pwm_ids.size() + " TRANSFAC Professional Matrix IDs parsed.");
	}
	
	
	/*
	 *  remove factor IDs and matrix IDs from TRANSFAC Professional   
	 */
	
	public void recomputeTrainingFile(String training_infile, String training_outfile) {
		
		TransfacParser tf_parser = new TransfacParser();
		tf_parser.parseFactors(training_infile);
		
		int del_factor_cnt = 0;
		int del_pfm_cnt = 0;
		
		int num_factors = tf_parser.tf_names.size();
		
		for (int i=tf_parser.tf_names.size()-1; i>=0; i--) {
			
			// remove factor IDs from TRANSFAC Professional
			if (professional_factor_ids.contains(tf_parser.tf_names.get(i))) {
				
				del_factor_cnt++;
				
				tf_parser.tf_names.remove(i);
				tf_parser.species.remove(i);
				tf_parser.crossrefs.remove(i);
				tf_parser.classes.remove(i);
				tf_parser.sequences1.remove(i);
				tf_parser.sequences2.remove(i);
				tf_parser.domains.remove(i);
				tf_parser.pfm_names.remove(i);
				tf_parser.pfms.remove(i);
				
				continue;
			}
			
			// remove matrix IDs from TRANSFAC Professional
			for (int j=tf_parser.get_pfm_names().get(i).size()-1; j>=0; j--) {
				
				if (professional_pwm_ids.contains(tf_parser.get_pfm_names().get(i).get(j))) {
					tf_parser.pfm_names.get(i).remove(j);
					tf_parser.pfms.get(i).remove(j);
				}
			}
			
			// remove factor, if no matrix is left
			if (tf_parser.pfm_names.get(i).isEmpty()) {
				
				System.out.println(del_pfm_cnt + ": " + tf_parser.tf_names.get(i));
				
				del_pfm_cnt++;
				
				tf_parser.tf_names.remove(i);
				tf_parser.species.remove(i);
				tf_parser.crossrefs.remove(i);
				tf_parser.classes.remove(i);
				tf_parser.sequences1.remove(i);
				tf_parser.sequences2.remove(i);
				tf_parser.domains.remove(i);
				tf_parser.pfm_names.remove(i);
				tf_parser.pfms.remove(i);
			}
		}
		
		// print statistics
		System.out.println("\n" + tf_parser.tf_names.size() + " / " + num_factors +  " are publicly available.");
		System.out.println(del_factor_cnt + " factors with IDs from TRANSFAC Professional were deleted.");
		System.out.println(del_pfm_cnt + " factors with matrices from TRANSFAC Professional were deleted.");
		
		// write training file
		tf_parser.writeFactorsToFile(training_outfile);
	}
	
	
	public static void main(String[] args) {
		
		
		
		/*
		String factor_table_file = "/home/eichner/projects/sabine/supplement/transfac_public_factortable.txt";
		String matrix_table_file = "/home/eichner/projects/sabine/supplement/transfac_public_matrixtable.txt";
		*/
		String factor_table_list = "/home/eichner/projects/sabine/supplement/transfac_public_factor_IDs.txt";
		String matrix_table_list = "/home/eichner/projects/sabine/supplement/transfac_public_matrix_IDs.txt";
		
		String training_infile = "/home/eichner/projects/sabine/supplement/basisdatensatz1_ohne_regulonDB.txt";
		String training_outfile = "/home/eichner/projects/sabine/supplement/basisdatensatz1_public.txt";
		
		
		OutfileScanner scanner = new OutfileScanner();
		/*
		scanner.generateIDList(factor_table_file, factor_table_list);
		scanner.generateIDList(matrix_table_file, matrix_table_list);
		*/
		
		scanner.getProfessionalIDs(training_infile, factor_table_list, matrix_table_list);
		scanner.recomputeTrainingFile(training_infile, training_outfile);
	}
}

