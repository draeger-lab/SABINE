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

package help;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;

import main.FBPPredictor;


public class FeatureFileGenerator {
	
	public boolean silent = false;
	
	public void generateOverallFeatureFile(String featuredir, String featurenamesfile, String outfile) {
		
		LibSVMFileGenerator generator = new LibSVMFileGenerator();
		
		ArrayList<String> features = new ArrayList<String>();
		
		try {
			
			BufferedReader br = new BufferedReader(new FileReader(new File(featurenamesfile)));
			
			String line = null;
			
			int count = 0;
			
			while((line = br.readLine()) != null) {
				
				features.add(line.trim());
				
				if (! silent) System.out.println("  [" + (++count) + "]\t" + line.trim());
				
			}
			
			if (! silent) System.out.println();
			
			br.close();
		}
		catch (IOException ioe) {
			System.out.println(ioe.getMessage());
			System.out.println("IOException occurred while parsing feature names.");
		}
		
		
		for(int i=0; i<features.size(); i++) {
			
			if(i == 0) {
			
				if (! silent) System.out.print("  Generating feature file for " + features.get(i) + "...");
			
				generator.generateUnlabeledLibSVMFile(
						
						featuredir + "/domain_scores_" + features.get(i) + ".out", 
						
						outfile
						
				);

				if (! silent) System.out.println("done.\n");
				
			}
			
			else {
				if (! silent) System.out.print("  Adding feature " + features.get(i) + "...");
				
				generator.addFeature(outfile, featuredir + "/domain_scores_" + features.get(i) + ".out");
				
				if (! silent) System.out.println("done.");

			}
			
		}
	}
		
		
		public void generateLibSVMTrainingSet(String featuredir, String traindir, String class_id, String labelfile, String outfile) {
			
			ArrayList<String> features = new ArrayList<String>();

			try {
				
				/*
				 *  write feature names to file in training directory
				 */
				
				FileFilter filter = new FileFilter();
				filter.setFormat(".*.out");
				filter.setDirectory(featuredir);
				String[] files = filter.listFilesFullPath();
				String curr_feature;
				int count = 0;
				
				for (int i=0; i<files.length; i++) {
					if (files[i].startsWith("domain_scores_")) {
						curr_feature = files[i].substring(14 ,files[i].lastIndexOf(".out"));
						features.add(curr_feature);
					}
				}
				
				Collections.sort(features);
				
				String featurenamesfile = traindir + class_id + FBPPredictor.featureNamesFileSuffix;
				BufferedWriter bw = new BufferedWriter(new FileWriter(new File(featurenamesfile)));
				
				for (int i=0; i<features.size(); i++) {
					bw.write(features.get(i) + "\n");
					if (! silent) System.out.println("  [" + (++count) + "]\t" + features.get(i));
				}
				bw.flush();
				bw.close();
			}
			catch (IOException ioe) {
				System.out.println(ioe.getMessage());
				System.out.println("IOException occurred while parsing feature names.");
			}
			
			
			LibSVMFileGenerator generator = new LibSVMFileGenerator();
			for(int i=0; i<features.size(); i++) {
				
				if(i == 0) {
				
					if (! silent) System.out.print("  Generating feature file for " + features.get(i) + "...");
				
					generator.generateLabeledLibSVMFile(
							
							featuredir + "/domain_scores_" + features.get(i) + ".out", 
							labelfile,
							outfile
							
					);
					
					if (! silent) System.out.println("done.\n");
					
				}
				
				else {
					
					if (! silent) System.out.print("  Adding feature " + features.get(i) + "...");
					
					generator.addFeature(outfile, featuredir + "/domain_scores_" + features.get(i) + ".out");
					
					if (! silent) System.out.println("done.");
				
				}
			}
		}
	}
	
	
	

