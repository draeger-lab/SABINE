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

package model;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.ArrayList;

public class SingleFeatureEvaluator {
	
	String base_dir = "";
	String train_dir = "trainingsets/";
	int class_id = 0; 
	
	public void generateSingleFeatureLibSVMFiles(String infile, String outdir) {
		
		try {
			
			/*
			 *  parse feature names
			 */
			
			ArrayList<String> features = new ArrayList<String>();
			String line;
			
			String featurenamesfile = train_dir + "class" + class_id + ".30featurenames";
			
			BufferedReader br = new BufferedReader(new FileReader(new File(featurenamesfile)));
			
			while ((line = br.readLine()) != null) {
				
				features.add(line.trim());
			}
			
			/*
			 *  generate single feature model files
			 */
			
			FeatureSelector feat_sel = new FeatureSelector();
			String outfile;
			
			for (int i=0; i<features.size(); i++) {
				
				outfile = "trainingset_" + features.get(i) + ".lp.class" + class_id + ".att";	
				feat_sel.selectFeatures(new int[] {i+1}, infile, outdir + outfile);
			}
			
			
			/*
			 *  perform cross-validation for each single-feature model
			 */

			
			ModelGenerator model_gen = new ModelGenerator("class" + class_id);
			model_gen.base_dir = base_dir;
			model_gen.train_dir = train_dir;
			model_gen.silent = true;
			
			String pwd = System.getProperty("user.dir");
			String rel_libsvm_dir = outdir;
			if (outdir.startsWith("/") && outdir.startsWith(pwd)) {
				rel_libsvm_dir = outdir.substring(pwd.length());
				//System.out.println(libsvm_dir);
			}
			
			double[] MSE = new double[features.size()];
			double[] AAE = new double[features.size()];
			
			DecimalFormat fmt = new DecimalFormat();
			DecimalFormatSymbols symbs = fmt.getDecimalFormatSymbols();
	        symbs.setDecimalSeparator('.');
	        fmt.setDecimalFormatSymbols(symbs);
	        fmt.setMaximumFractionDigits(2);
	        fmt.setMinimumFractionDigits(2);
			
			for (int i=0; i<features.size(); i++) {
				
				
				System.out.println("Performing CV on training set: trainingset_" + features.get(i) + ".lp.class" + class_id + ".att");
				model_gen.runCrossvalidation(rel_libsvm_dir + "trainingset_" + features.get(i) + ".lp.class" + class_id + ".att");
				
				
				double[] params = model_gen.parseCVResults(outdir + "trainingset_" + features.get(i) + ".lp.class" + class_id + ".att.dat",
												 		   outdir + "trainingset_" + features.get(i) + ".lp.class" + class_id + ".par");
				
				MSE[i] = params[5];
				AAE[i] = params[3];
				
				System.out.println("  MSE: " + fmt.format(MSE[i]));
				System.out.println("  AAE: " + fmt.format(AAE[i]) + "\n");
			}

			/*
			 *  write MSE and AAE for all single feature models to text file
			 */
			
			BufferedWriter bw = new BufferedWriter(new FileWriter(new File(outdir + "result.out")));
				
			for (int i=0; i<MSE.length; i++) {
					
				bw.write(MSE[i] + "\t" + AAE[i] + "\n");
			}
			bw.flush();
			bw.close();
			
		}
		catch(IOException ioe) {
			System.out.println(ioe.getMessage());
			System.out.println("IOException occurred while selecting PCA-based features.");
		}
	}
	
	public static void main(String args[]) {
		
		SingleFeatureEvaluator evaluator = new SingleFeatureEvaluator();
		evaluator.base_dir = "svr_models/my_basedir/";
		evaluator.train_dir = "trainingsets_procaryotes/";
		
		String infile = evaluator.base_dir + "libsvm_files/trainingset.lp.class0.att";
		String outdir = evaluator.base_dir + "libsvm_files/";
		
		evaluator.generateSingleFeatureLibSVMFiles(infile, outdir);
	}
}

