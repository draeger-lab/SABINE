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

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;

public class OutfileWriter {

	
	/*
	 *  removes redundant entries in SABINE output file
	 */
	
	public void filterLabelfile(String infile, String outfile) {

		TransfacParser transfac_parser = new TransfacParser();
			
		transfac_parser.parseFactors(infile);
		transfac_parser.writeLabelFile(outfile);
	}
	
	
	public void filterOutputfile(String infile, String outfile) {
		
		try {
			
			/*
			 *  obtain unique list of factor names
			 */
			
			PredictionEvaluator evaluator = new PredictionEvaluator();
			evaluator.parse_predictions(infile);
			
			ArrayList<String> unique_names = new ArrayList<String>();
			boolean[] irrelevant_elements = new boolean[evaluator.tf_names.size()];
			
			for (int i=0; i<evaluator.tf_names.size(); i++) {
				
				if (unique_names.contains(evaluator.tf_names.get(i))) {
					irrelevant_elements[i] = true;
				}
				else {
					unique_names.add(evaluator.tf_names.get(i));
				}
			}
			
			/*
			 *  filter factors
			 */
			
			for (int i=evaluator.tf_names.size()-1; i>=0; i--) {
				
				if (irrelevant_elements[i]) {
					
					evaluator.tf_names.remove(i);
					evaluator.best_matches.remove(i);
					evaluator.predicted_matrices.remove(i);
				}
			}
			
			/*
			 *  write filtered outfile
			 */
			
			BufferedWriter bw = new BufferedWriter(new FileWriter(new File(outfile)));
			
			for (int i=0; i<evaluator.tf_names.size(); i++) {
				
				// write name
				bw.write("NA  " + evaluator.tf_names.get(i) + "\nXX\n");
				
				// write best matches
				if (evaluator.best_matches.get(i).isEmpty())
					bw.write("BM  none\n");
				
				for (int j=0; j<evaluator.best_matches.get(i).size(); j++) {
					bw.write("BM  " + evaluator.best_matches.get(i).get(j).getName() + "  " + 
									+ evaluator.best_matches.get(i).get(j).getLabel() + "\n");
				}
				bw.write("XX\n");
				
				// write matrix
				if (evaluator.predicted_matrices.get(i).isEmpty())
					bw.write("MA  none\n");
				
				for (int j=0; j<evaluator.predicted_matrices.get(i).size(); j++) {	
					bw.write("MA  " + evaluator.predicted_matrices.get(i).get(j) + "\n");
				}
				bw.write("XX\n//\nXX\n");
			}
			bw.flush();
			bw.close();
		}
		catch (IOException e) {
			System.out.println(e.getMessage());
			System.out.println("IOException occurred while writing output file for SABINE.");
		}
	}
	
	
	public static void main(String[] args) {
		
		OutfileWriter writer = new OutfileWriter();
		
		String infile = "../backup/res_07_12/ensembl/rat/rn_testset.labels";
		String outfile = infile;
		writer.filterLabelfile(infile, outfile);
	}
}

