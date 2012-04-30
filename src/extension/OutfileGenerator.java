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

import java.util.ArrayList;
import java.util.StringTokenizer;

public class OutfileGenerator {
	
	public void rewritePredictionOutfile(String infile, String outfile) {
		
		TransfacParser tf_parser = new TransfacParser();
		tf_parser.parseFactors(infile);
		tf_parser.writeFactorsToFile(outfile);
	}
	
	public void addBestMatches(String infile, String predfile, String outfile) {
		
		TransfacParser tf_parser = new TransfacParser();
		tf_parser.parseFactors(infile);
		
		PredictionEvaluator evaluator = new PredictionEvaluator(); 
		evaluator.parse_predictions(predfile);
		
		int idx;
		String curr_factor;
		
		for (int i=0; i<tf_parser.tf_names.size(); i++) {
			
			curr_factor = new StringTokenizer(tf_parser.tf_names.get(i)).nextToken().substring(1);
			idx = evaluator.tf_names.indexOf(curr_factor);
			

			
			ArrayList<String> curr_name = new ArrayList<String>();
			
			if (evaluator.best_matches.get(idx).get(0).getLabel() >= 0.95) {
				
				curr_name.add(evaluator.best_matches.get(idx).get(0).getName());
				System.out.println(curr_factor + " -> " + evaluator.best_matches.get(idx).get(0).getName());
				tf_parser.pfm_names.set(i, curr_name);
				curr_name = new ArrayList<String>();
			}
		}
		
		tf_parser.writeFactorsToFile(outfile);
	}
	
	public static void main(String[] args) {
		

		
		//OutfileGenerator generator = new OutfileGenerator();
		
		//String infile = "/home/eichner/Desktop/supplement/SupplementaryFile3_predictedTFs.txt";	
		//String outfile = "/home/eichner/Desktop/supplement/old_predictions.txt";
		//generator.rewritePredictionOutfile(infile, outfile);
		
		/*
		String infile = "/home/eichner/Desktop/supplement/Rattus_Uniprot_NewPredictedTFs.txt";
		String predfile = "/home/eichner/Desktop/supplement/rn_all_factors.out";
		String outfile = "/home/eichner/Desktop/supplement/new_predictions_rat.txt";
		generator.addBestMatches(infile, predfile, outfile);
		*/

	}
}

