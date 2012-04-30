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

import help.FormatConverter;

import java.io.BufferedWriter;
import java.util.ArrayList;

/*
 *  parses SABINE *.input and *.output files
 */

public class SABINEInputOutputParser implements SABINEParser {

	public ArrayList<String> tf_names = new ArrayList<String>();
	public ArrayList<String> classes = new ArrayList<String>();
	ArrayList<String> sequences1 = new ArrayList<String>();
	ArrayList<String> sequences2 = new ArrayList<String>();
	ArrayList<ArrayList<String>> domains = new ArrayList<ArrayList<String>>();
	  
	ArrayList<String> species = new ArrayList<String>();
	public ArrayList<ArrayList<String>> matrices = new ArrayList<ArrayList<String>>();
	
	
	public ArrayList<ArrayList<String>> get_domains() {
		return domains;
	}

	public ArrayList<ArrayList<String>> get_matrices() {
		return matrices;
	}

	public ArrayList<String> get_sequences1() {
		return sequences1;
	}

	public ArrayList<String> get_sequences2() {
		return sequences2;
	}

	public ArrayList<String> get_species() {
		return species;
	}

	public ArrayList<String> get_tf_names() {
		return tf_names;
	}
	
	
	public void parseInputFile(String inputfile) {
		
		TransfacParser transfac_parser = new TransfacParser();
		transfac_parser.parseFactors(inputfile);
		tf_names = transfac_parser.tf_names;
		sequences1 = transfac_parser.sequences1;
		sequences2 = transfac_parser.sequences2;
		domains = transfac_parser.domains;
		species = transfac_parser.species;
		classes = transfac_parser.classes;
	}
	
	public void parseLabelFile(String labelfile) {
		
		PredictionEvaluator evaluator = new PredictionEvaluator();
		evaluator.parse_labels(labelfile);
		matrices = evaluator.annotated_matrices;
	}
	
	
	public void writeInputFile(BufferedWriter bw, int[] relevant_factors, int class_id) {
		
		/*
		 *  use method implemented in class SABINETrainingsetParser
		 */
		
		SABINETrainingsetParser trainingset_parser = new SABINETrainingsetParser();
		trainingset_parser.tf_names = tf_names;
		trainingset_parser.sequences1 = sequences1;
		trainingset_parser.sequences2 = sequences2;
		trainingset_parser.domains = domains;
		trainingset_parser.species = species;
		trainingset_parser.matrices = matrices;
		
		trainingset_parser.writeInputFile(bw, relevant_factors, class_id);
	}
	
	
	public void parseAll(int class_id, String inputfile, String labelfile) {
		
		parseInputFile(inputfile);
		if (labelfile != null) parseLabelFile(labelfile); 
		
		FormatConverter converter = new FormatConverter();
		
		boolean[] right_class = new boolean[classes.size()];
		int curr_class;
		
		for (int i=0; i<classes.size(); i++) {
			
			if (classes.get(i).equals("NA")) continue;
			
			curr_class = Integer.parseInt(converter.getTransfacClass(classes.get(i)).substring(0,1));
			
			if (curr_class == class_id) right_class[i] = true;
		}
		
		for (int i=classes.size()-1; i>=0; i--) {
			
			if (! right_class[i]) {
				tf_names.remove(i);
				sequences1.remove(i);
				sequences2.remove(i);
				domains.remove(i);
				species.remove(i);
				classes.remove(i);
				if (matrices.size() > 0) matrices.remove(i); 
			}
		}
	}
	
	
	public static void main(String[] args) {
		
		SABINEInputOutputParser sabine_parser = new SABINEInputOutputParser();
		
		String inputfile = "testdir/hs_testset.input";
		String labelfile = "testdir/hs_testset.labels";
		
		sabine_parser.parseInputFile(inputfile);
		sabine_parser.parseLabelFile(labelfile);
	}
}

