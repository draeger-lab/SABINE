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
import java.util.StringTokenizer;

public class PWMModelFileParser {
	
	
	ArrayList<String> header = new ArrayList<String>();
	ArrayList<String> factor_names = new ArrayList<String>();
	ArrayList<String> species = new ArrayList<String>();
	ArrayList<ArrayList<String>> pwm_names = new ArrayList<ArrayList<String>>();
	ArrayList<ArrayList<String[]>> pwms = new ArrayList<ArrayList<String[]>>();
	
	public void parsePWMModels(String infile) {
		
		try {
			BufferedReader br = new BufferedReader(new FileReader(new File(infile)));
			String line;
			ArrayList<String> curr_pwm_names = new ArrayList<String>();
			ArrayList<String[]> curr_pwms = new ArrayList<String[]>();
			String[] curr_pwm = new String[4];
			
			// parse header
			while ((line = br.readLine()) != null && line.startsWith("#")) {
				header.add(line.trim());
			}
			
			// parse entries
			while ((line = br.readLine()) != null) {
				
				if (line.startsWith("NA")) {
					factor_names.add(line.substring(4).trim());
				}
				else {
					System.out.println(line);
					System.out.println("Parse Error. \"NA\" expected.");
					System.exit(0);
				}
				line = br.readLine();
				
				if (line.startsWith("SP")) {
					species.add(line.substring(4).trim());
				}
				else {
					System.out.println("Parse Error. \"SP\" expected.");
					System.exit(0);
				} 
				line = br.readLine();
				
				
				curr_pwm_names = new ArrayList<String>();
				curr_pwms = new ArrayList<String[]>();
				
				while (line.startsWith("MN") || line.startsWith("MA")) {
					
					// skip matrix with missing ID
					if (line.startsWith("MA")) {
						
						curr_pwm_names.add("");
						
						curr_pwm = new String[4];
						curr_pwm[0] = line.substring(4).trim();
						curr_pwm[1] = br.readLine().substring(4).trim();
						curr_pwm[2] = br.readLine().substring(4).trim();
						curr_pwm[3] = br.readLine().substring(4).trim();
						
						line = br.readLine();
					}
					else {
						curr_pwm_names.add(line.substring(4).trim());
						line = br.readLine();
						
						curr_pwm = new String[4];
						if (line.startsWith("MA")) {
							curr_pwm[0] = line.substring(4).trim();
							curr_pwm[1] = br.readLine().substring(4).trim();
							curr_pwm[2] = br.readLine().substring(4).trim();
							curr_pwm[3] = br.readLine().substring(4).trim();
							
							line = br.readLine();
						}
					}
					curr_pwms.add(curr_pwm);
					
					if (line == null) break;
					
				}
 				if (line == null || line.startsWith("XX")) {
					pwm_names.add(curr_pwm_names);
					pwms.add(curr_pwms);
				}
				else {
					System.out.println("Parse Error. \"XX\" expected.");
					System.exit(0);
				}

			}
			br.close();
		}
		catch(IOException ioe) {
			System.out.println(ioe.getMessage());
			System.out.println("IOException occurred while reading PWM model file.");
		}
	}
	
	/*
	public void addPublicPWMs(String training_factor_file, String transfac_public_file) {
		
		OutfileScanner scanner = new OutfileScanner();
		scanner.getProfessionalIDs(training_factor_file, transfac_public_file);
		
		TransfacParser tf_parser = new TransfacParser();
		tf_parser.parseFactors(training_factor_file);
		
		ArrayList<String> all_pfm_names = new ArrayList<String>();
		ArrayList<String[]> all_pfms = new ArrayList<String[]>();
		
		for (int i=0; i<tf_parser.pfm_names.size(); i++) {
			for (int j=0; j<tf_parser.pfm_names.get(i).size(); j++) {
				
				if (scanner.public_pwm_ids.contains(tf_parser.pfm_names.get(i).get(j)) &&
				    ! all_pfm_names.contains(tf_parser.pfm_names.get(i).get(j))) {
					
					all_pfm_names.add(tf_parser.pfm_names.get(i).get(j));
					all_pfms.add(tf_parser.pfms.get(i).get(j));
				}
			}
		}
		
		System.out.println(all_pfm_names.size() + " TRANSFAC Public Matrices found.");
		
		ArrayList<String[]> curr_pwms = new ArrayList<String[]>();
		*/
		/*
		for (int i=0; i<all_pfm_names.size(); i++) {
			System.out.println(all_pfm_names.get(i));
		}
		
		for (int i=0; i<pwm_names.size(); i++) {
			for (int j=0; j<pwm_names.get(i).size(); j++) {
				System.out.println(new StringTokenizer(pwm_names.get(i).get(j)).nextToken());
			}
		}
		*/
		/*
		for (int i=0; i<pwm_names.size(); i++) {
		
			curr_pwms = pwms.get(i);
			int idx;
			
			for (int j=0; j<curr_pwms.size(); j++) {
				
				
				
				if (! pwm_names.get(i).get(j).isEmpty() &&
					(idx = all_pfm_names.indexOf(new StringTokenizer(pwm_names.get(i).get(j)).nextToken())) != -1 &&
					curr_pwms.get(j)[0] == null)	{
					
					curr_pwms.set(j, all_pfms.get(idx));
				}
			}
			pwms.set(i, curr_pwms);
			curr_pwms = new ArrayList<String[]>();
		}
	}
	*/
	public void writePWMModelFile(String outfile) {
		
		try {
			BufferedWriter bw = new BufferedWriter(new FileWriter(new File(outfile)));
			
			// write header
			for (int i=0; i<header.size(); i++) {
				bw.write(header.get(i) + "\n");
			}
			bw.write("\n");
			
			// write  entries
			for (int i=0; i<factor_names.size(); i++) {
				
				bw.write("NA  " + factor_names.get(i) + "\n");
				bw.write("SP  " + species.get(i) + "\n");
				
				for (int j=0; j<pwm_names.get(i).size(); j++) {
					
					if (! pwm_names.get(i).get(j).isEmpty()) {
						bw.write("MN  " + pwm_names.get(i).get(j) + "\n");
					}
					
					if (pwms.get(i).get(j)[0] != null) {
						bw.write("MA  " + pwms.get(i).get(j)[0] + "\n");
						bw.write("MA  " + pwms.get(i).get(j)[1] + "\n");
						bw.write("MA  " + pwms.get(i).get(j)[2] + "\n");
						bw.write("MA  " + pwms.get(i).get(j)[3] + "\n");
					}
				}
				bw.write("XX\n");
			}
			bw.flush();
			bw.close();
		}
		catch(IOException ioe) {
			System.out.println(ioe.getMessage());
			System.out.println("IOException occurred while writing PWM model file.");
		}
	}
	
	
	public void filterPWMModels(String dataset_file) {
		
		TransfacParser tf_parser = new TransfacParser();
		tf_parser.parseFactors(dataset_file);
		
		boolean[] relevant_factors = new boolean[factor_names.size()];
		
		for (int i=0; i<factor_names.size(); i++) {
			
			if (tf_parser.tf_names.contains(new StringTokenizer(factor_names.get(i)).nextToken())) {
				
				relevant_factors[i] = true;
			}
		}
		
		for (int i=factor_names.size()-1; i>=0; i--) {
			
			if (! relevant_factors[i]) {
				factor_names.remove(i);
				species.remove(i);
				pwm_names.remove(i);
				pwms.remove(i);
			}
		}
		System.out.println(factor_names.size() + " factors remaining.");
	}
	
	public static void main(String[] args) {
		
		PWMModelFileParser parser = new PWMModelFileParser();
		/* 
		
		String model_file = "/home/eichner/Desktop/supplement/DNA-bindingSpecificityModels.txt";
		String outfile = "/home/eichner/Desktop/supplement/pwm_models.txt";
		
		String training_factor_file = "trainingsets_eucaryotes/transpall_interleaved_classes.out";
		String transfac_public_file = "/home/eichner/Desktop/supplement/transfac_public.txt";
		
		parser.parsePWMModels(model_file);
		parser.addPublicPWMs(training_factor_file, transfac_public_file);
		parser.writePWMModelFile(outfile);
		*/
		
		String model_file = "/home/eichner/Desktop/supplement/pwm_models.txt";
		String dataset_file = "/home/eichner/Desktop/supplement/dataset1.txt";
		String outfile = "/home/eichner/Desktop/supplement/pwm_models_new.txt";
		
		
		parser.parsePWMModels(model_file);
		parser.filterPWMModels(dataset_file);
		parser.writePWMModelFile(outfile);
		
	}
}

