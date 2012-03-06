package libsvm;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.PriorityQueue;
import java.util.StringTokenizer;


import optimization.MatlignOptimizer;
import optimization.MoStaOptimizer;
import optimization.Optimizer;


public class SVMPredictor {
	
	public boolean dynamic_threshold = false;
	double sim_thres = 0.0;
	int max_num = 0;
	public boolean silent = false; 
	public boolean gui_output_mode = false;
	public String basedir = "";
	
	public double high_conf_bmt = 0.0;
	public double medium_conf_bmt = 0.0;
	public double low_conf_bmt = 0.0;
	
	/*
	 * global variables that are initialized by extractBestHits()
	 */
	
	
	public void predictLabels(String modelfile, String testset, String outfile) {
		
	
		try {
			
			if (gui_output_mode) {
				System.out.print("done.\nPredicting PFM similarities for candidate factors...");
			}
			
			if (! silent) { 
				System.out.println("\n\nPredicting PFM similarities for all relevant TF pairs.\n");
			
				System.out.println("  Model File       : " + modelfile.substring(modelfile.lastIndexOf("/") + 1));
				System.out.println("  Input File       : " + testset.substring(testset.lastIndexOf("/") + 1));
				System.out.println("  Predicted Labels : " + outfile.substring(outfile.lastIndexOf("/") + 1) + "\n");
			}
				
			svm_predict.main( new String[] { testset, modelfile, outfile } );
			
			

		}
		catch(IOException ioe) {
			System.out.println(ioe.getMessage());
		}
		
	}
	
	
	
	public ArrayList<LabeledTF> extractBestHits(String tfnamepairsfile, String predictedlabels, double threshold, int maxNumPairs) {
		
		ArrayList<LabeledTF> best_matches = new ArrayList<LabeledTF>();		
		
		
		BufferedReader br1 = null;
		BufferedReader br2 = null;
		
		if (! silent) System.out.println("\nExtracting best matches.");
		if (gui_output_mode) System.out.print("done.\n\nExtracting best matches:");
		
		try {
			
			br1 = new BufferedReader(new FileReader(new File(tfnamepairsfile)));
			br2 = new BufferedReader(new FileReader(new File(predictedlabels)));
			
			StringTokenizer strtok = null;
			
			String line1 = null;
			String line2 = null;
			
			String train_tf = null;
			
			double predicted_label = 0.0;
			
	
			PriorityQueue<LabeledTF> matches = new PriorityQueue<LabeledTF>();
			
			while((line1 = br1.readLine()) != null) {
				
				line2 = br2.readLine();
				
				
				strtok = new StringTokenizer(line1);
				
				strtok.nextToken();
				
				strtok.nextToken(); 			// vs.
				
				train_tf = strtok.nextToken();
				
				
				predicted_label = Double.parseDouble(line2);
				
				
				
				matches.add(new LabeledTF(train_tf, predicted_label));
			
			}
			
			
			if (! silent) System.out.println("\n  Number of candidate TFs : " + matches.size());
			
			if (gui_output_mode) System.out.println("\n  Number of candidate TFs : " + matches.size());
			
			int size = matches.size();
			
			boolean first = true;
			
			for(int i=0; i<Math.min(size, maxNumPairs); i++) {
							
				LabeledTF tf = matches.poll();
				
				if (first && dynamic_threshold) {
					
					if (tf.label > high_conf_bmt)
						threshold = high_conf_bmt;
					else if (tf.label > medium_conf_bmt)
						threshold = medium_conf_bmt;
					else 
						threshold = low_conf_bmt;
					
					first = false;
				}
				
				if(tf.label > threshold) {
									
					best_matches.add(tf);
												
				}		
			}
							
			br1.close();
			br2.close();
			
		}
		catch(IOException ioe) {
			System.out.println(ioe.getMessage());
			System.out.println("IOException occurred while calculating best matches.");
		}

		
		if (! silent) System.out.println("  Number of best matches  : " + best_matches.size());
		if (gui_output_mode) System.out.println("  Number of best matches  : " + best_matches.size());
		
	 	return best_matches;
	 	
	}
	
	
	public ArrayList<String> generateFBP(ArrayList<LabeledTF> tfid_list, String fbplookupfile, Optimizer optimizer, double outlier_tolerance) {
		
		ArrayList<String> res = null;
		
		ArrayList<ArrayList<String>> fbps = new ArrayList<ArrayList<String>>();
		
		
		
		if (! silent) System.out.println("\nGenerating Meta-FBP for " + tfid_list.size() + " FBP(s).\n");
		
		
		
		if (! silent) System.out.println("---PARSE---");		
		for(int i=0; i<tfid_list.size(); i++) {
			fbps.add(getFBP(tfid_list.get(i), fbplookupfile));
			
		}
		if (! silent) System.out.println("");
		
		
		
	// filter outlier fbps	
		
		if (! silent) System.out.println("\n  Filtering outlier FBPs. (tolerance = " + (1 + outlier_tolerance) + " * mean distance)\n" );
		
		PWMOutlierFilter filter = new PWMOutlierFilter();
		
		filter.tolerance = outlier_tolerance;
		filter.silent = silent;
		
		fbps = filter.filterOutliers(fbps, optimizer);

		
		if (! silent) System.out.println("  " + fbps.size() + " FBP(s) remaining.\n");
		
	// generate overall fbp	
		
		if (fbps.size() > 1) {
			
			STAMPRunner runner = new STAMPRunner();
			runner.basedir = basedir;
			res = runner.runSTAMP(fbps);
			
		}
		
		else {
			
			res = fbps.get(0);
			
		}
		
		return res;
		
	}
	
	
	
	
	public ArrayList<String> getFBP(LabeledTF tfid, String fbplookupfile) {
		
		ArrayList<String> res = new ArrayList<String>();
		
		BufferedReader br = null;
		if (! silent) System.out.println("SC  FBPs " + tfid.name + " (score = " + tfid.label + ")");
		
		try {
			
			br = new BufferedReader(new FileReader(new File(fbplookupfile)));
			
			String line = null;
			
			
		// goto respective factor	
			
			while((line = br.readLine()) != null && !line.startsWith("DE\t" + tfid.name)) ;
			
			
		// parse its FBP	
			
			while(!(line = br.readLine()).startsWith("XX")) {
				
				res.add(line);
				
			}
			
			br.close();
			
		}
		catch(IOException ioe) {
			System.out.println(ioe.getMessage());
			System.out.println("IOException occurred while looking up FBP for " + tfid.name + ".");
		}
		return res;
		
	}
	
	
	public void performPFMTransfer(ArrayList<LabeledTF> best_matches, String fbplookupfile, Optimizer optimizer, double outlier_tolerance, String outfile) {
			
		if (optimizer instanceof MatlignOptimizer) {
			((MatlignOptimizer) optimizer).basedir = basedir;
		}
		else if (optimizer instanceof MoStaOptimizer) {
			((MoStaOptimizer) optimizer).basedir = basedir;
		}
		
		if(best_matches.size() == 0)  {
					
			if (! silent) System.out.println("\nNo PFM-Transfer possible.");
			if (gui_output_mode) System.out.println("\nNo PFM transfer possible.");
			
			try {
				BufferedWriter bw = new BufferedWriter(new FileWriter(new File(outfile)));
		
				bw.write("BM  none\nXX\nMA  none\nXX");
				
				bw.flush();
				bw.close();
			}
			catch(IOException ioe) {
				System.out.println(ioe.getMessage());
				System.out.println("IOException occurred while writing output file.");
			}
			
			return;
					
		}
		
		/*
		 *  merge at most the PFMs of 5 Factors
		 */
		ArrayList<LabeledTF> merged_factors = new ArrayList<LabeledTF>();
		int max_number_merged_matrices = 5; 
		
		for(int i=0; i<Math.min(best_matches.size(), max_number_merged_matrices); i++) {
			merged_factors.add(best_matches.get(i));
		}
		
		ArrayList<String> predicted_fbp = generateFBP(merged_factors, fbplookupfile, optimizer, outlier_tolerance);
		
		if (! silent) {
			System.out.println("\n-------------");
			System.out.println("Predicted FBP");
			System.out.println("-------------\n");
		
			for(int i=0; i<predicted_fbp.size(); i++) {
			
				System.out.println("MA  " + predicted_fbp.get(i));
			
			}
		
			System.out.println("\ndone.");
		}
		
		// write results to Outfile
		try {
			BufferedWriter bw = new BufferedWriter(new FileWriter(new File(outfile)));
			
			// write best matches and prediction scores to outfile
			for (int i=0; i<best_matches.size(); i++) {
				bw.write("BM  " + best_matches.get(i).name + "  " + best_matches.get(i).label + "\n");
			}
			bw.write("XX\n");
			
			// write PFM to outfile
			for (int i=0; i<predicted_fbp.size(); i++) {
				bw.write("MA  " + predicted_fbp.get(i) + "\n");
			}
			bw.write("XX\n");
			bw.flush();
			bw.close();
		}
		catch(IOException ioe) {
			System.out.println(ioe.getMessage());
			System.out.println("IOException occurred while writing output file.");
		}
		
	}
	
	
	
}
