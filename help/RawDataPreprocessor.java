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
import java.util.PriorityQueue;
import java.util.StringTokenizer;

import core.DomainFeatureCalculator;
import extension.PFMFormatConverter;
import extension.TransfacParser;

/*
 * 
 * reduces, sorts and generates data in folder "trainingset" so that
 * 
 *   - only TFs of the training sets are included in these files
 *   
 *   - all TFs and FBPs are ordered alphabetically
 * 
 * 	 - all secondary structures of the tfs are precalculated
 * 
 * 	 - all species of the tfs are looked up in advance
 * 
 */
public class RawDataPreprocessor {
	
	public String basedir = "./PSIPRED/";
	

	public void filterAndSortDataset(String infile, String trainingsetfactorfile, String outfile) {
		
		
		System.out.println("\nRemoving all TFs that are not part of the training set.\n");
		
		
		BufferedReader br = null;
		
		try {
			
			/*
			 * 
			 * parse ordered list of tf names
			 * 
			 */
			
			br = new BufferedReader(new FileReader(new File(trainingsetfactorfile)));
		
			String line = null;
			
			StringTokenizer strtok = null;
			
			ArrayList<String> orderednames = new ArrayList<String>();
			
			
			while((line = br.readLine()) != null) {
					
				strtok = new StringTokenizer(line);
					
				orderednames.add(strtok.nextToken());
				
			}
			
			br.close();
			
			
			
			System.out.println("  Leaving entries of " + orderednames.size() + " TFs.\n");
			
			
			/*
			 * 
			 * parse unordered list of tfs
			 * 
			 */
			
			int num_cut = 0;
			int num_retained = 0;
			
			br = new BufferedReader(new FileReader(new File(infile)));
			
			String name = null;
			
			boolean cut = false;
			
			ArrayList<ArrayList<String>> tfs = new ArrayList<ArrayList<String>>();
			
			ArrayList<String> tf_names = new ArrayList<String>();
			
			
			
			while((line = br.readLine()) != null) {
				
				if(!line.startsWith("NA  ")) {
					System.out.println("Parse Error. NA expected. Aborting.");
					System.exit(0);
				}
				
				name = new StringTokenizer(line.substring(4)).nextToken();
				
				
				cut = !orderednames.contains(name); 
				
				
				if(cut) num_cut++;
				if(!cut) num_retained++;
				
				
				ArrayList<String> entry = new ArrayList<String>();
				
				
				if(!cut) entry.add(line + "\n");						// name
				
				
				while(!(line = br.readLine()).startsWith("//")) {
					
					if(!cut) entry.add(line + "\n");					// data
						
				}
				
				if(!cut) entry.add(line + "\n");						// //
				
				line = br.readLine();
				
				if(!cut) entry.add(line + "\n");						// XX
				
				
				if(!cut) tf_names.add(name);
				if(!cut) tfs.add(entry);
				
			}
			
			br.close();
			
			
			/*
			 * 
			 * write ordered data to output file
			 * 
			 */
			
			BufferedWriter bw = new BufferedWriter(new FileWriter(new File(outfile)));
			
			for(int i=0; i<orderednames.size(); i++) {
				
				int index = tf_names.indexOf(orderednames.get(i));
				
				ArrayList<String> entry = tfs.get(index);
				
				for(int j=0; j<entry.size(); j++) {
					
					bw.write(entry.get(j));
					
				}
				
			}
			
			bw.flush();
			bw.close();
		
			
			System.out.println("  Number of TFs : " + (num_cut + num_retained));
			System.out.println("  Remaining TFs : " + num_retained);
			System.out.println("  Deleted   TFs : " + num_cut);
			
		
		}
		catch(IOException ioe) {
			System.out.println(ioe.getMessage());
			System.out.println("IOException occurred while cutting dataset.");
		}
	
	}
	
	public void sortDataset(String infile, String namesfile, String outfile) {
		
		BufferedReader br = null;
		
		try {
			
			/*
			 * 
			 * parse ordered list of tf names
			 * 
			 */
			
			br = new BufferedReader(new FileReader(new File(infile)));
		
			String line = null;
			
			StringTokenizer strtok = null;
			
			PriorityQueue<String> orderednames = new PriorityQueue<String>();
			
			
			while((line = br.readLine()) != null) {
				
				if (line.startsWith("NA  ")) {
					strtok = new StringTokenizer(line.substring(4).trim());
					orderednames.add(strtok.nextToken());
				}	
			}
			
			br.close();
			
			
			/*
			 * 
			 * parse unordered list of tfs
			 * 
			 */
			
			
			br = new BufferedReader(new FileReader(new File(infile)));
			
			String name = null;
			
			ArrayList<ArrayList<String>> tfs = new ArrayList<ArrayList<String>>();
			
			ArrayList<String> tf_names = new ArrayList<String>();
			

			while((line = br.readLine()) != null) {
				
				if(!line.startsWith("NA  ")) {
					System.out.println("Parse Error. NA expected. Aborting.");
					System.exit(0);
				}
				
				name = new StringTokenizer(line.substring(4)).nextToken();
				
				ArrayList<String> entry = new ArrayList<String>();
				
				
				entry.add(line + "\n");						// name
				
				
				while(!(line = br.readLine()).startsWith("//")) {
					
					entry.add(line + "\n");					// data
						
				}
				
				entry.add(line + "\n");						// //
				
				line = br.readLine();
				
				entry.add(line + "\n");						// XX
				
				
				tf_names.add(name);
				tfs.add(entry);
				
			}
			
			br.close();
			
			
			/*
			 * 
			 * write ordered data to output files
			 * 
			 */
			
			BufferedWriter bw = new BufferedWriter(new FileWriter(new File(outfile)));
			BufferedWriter bw_names = new BufferedWriter(new FileWriter(new File(namesfile)));
			int num_factors = orderednames.size();
			
			for(int i=0; i<num_factors; i++) {
				
				bw_names.write(orderednames.peek() + "\n");
				
				int index = tf_names.indexOf(orderednames.poll());
				
				ArrayList<String> entry = tfs.get(index);
				
				for(int j=0; j<entry.size(); j++) {
					
					bw.write(entry.get(j));
					
				}
			}
			bw.flush();
			bw.close();
			bw_names.flush();
			bw_names.close();
		}
		catch(IOException ioe) {
			System.out.println(ioe.getMessage());
			System.out.println("IOException occurred while cutting dataset.");
		}
	
	}
	
	
public void sortFBPs(String infile, String namesfile, String outfile) {
		

		BufferedReader br = null;
		
		try {
			
			/*
			 * 
			 * parse ordered list of tf names
			 * 
			 */
			
			br = new BufferedReader(new FileReader(new File(namesfile)));
		
			String line = null;
			
			StringTokenizer strtok = null;
			
			ArrayList<String> orderednames = new ArrayList<String>();
			
			
			while((line = br.readLine()) != null) {
				
				strtok = new StringTokenizer(line);
					
				orderednames.add(strtok.nextToken());
			}
			br.close();
			
			/*
			 * 
			 * parse unordered list of tfs
			 * 
			 */
			

			
			br = new BufferedReader(new FileReader(new File(infile)));
			
			String name = null;
			
			ArrayList<ArrayList<String>> fbps = new ArrayList<ArrayList<String>>();
			
			ArrayList<String> fbp_names = new ArrayList<String>();
			
			
			
			while((line = br.readLine()) != null) {
				
				if(!line.startsWith("DE\t")) {
					System.out.println("Parse Error. NA expected. Aborting.");
					System.exit(0);
				}
				
				name = new StringTokenizer(line.substring(3)).nextToken();


				ArrayList<String> entry = new ArrayList<String>();
				
				
				entry.add(line + "\n");						// name
				
				
				while(!(line = br.readLine()).startsWith("XX")) {
					
					 entry.add(line + "\n");					// data
						
				}
				
				entry.add(line + "\n");						// XX
				
				
				fbp_names.add(name);
				fbps.add(entry);
				
			}
			
			br.close();
			
			
			/*
			 * 
			 * write ordered data to output file
			 * 
			 */
			
			BufferedWriter bw = new BufferedWriter(new FileWriter(new File(outfile)));
			
			for(int i=0; i<orderednames.size(); i++) {
				
				int index = fbp_names.indexOf(orderednames.get(i));
				
				ArrayList<String> entry = fbps.get(index);
				
				for(int j=0; j<entry.size(); j++) {
					
					bw.write(entry.get(j));
					
				}
				
			}
			
			bw.flush();
			bw.close();
		}
		catch(IOException ioe) {
			System.out.println(ioe.getMessage());
			System.out.println("IOException occurred while cutting dataset.");
		}
	
	}
	
	
	public void precalculateSecondaryStructures(String infile, String outfile, String namesfile) {
		
		
		int SEQLINELENGTH = 60;
		
		System.out.println("  Predicting secondary structures.");

		BufferedReader br = null;
		BufferedWriter bw = null;
		
		try {
			
			br = new BufferedReader(new FileReader(new File(infile)));
			bw = new BufferedWriter(new FileWriter(new File(outfile)));
		
			String line = null;
			
			String name = null;
			
			int count_tfs = 0;
			
			int num_factors = getNumFactors(namesfile);
			
			while((line = br.readLine()) != null) {
				
				if(!line.startsWith("NA  ")) {
					System.out.println("Parse Error. NA expected. Aborting.");
					System.exit(0);
				}
				
				name = new StringTokenizer(line.substring(4)).nextToken();
				
				System.out.print("    Processing factor " + ++count_tfs + " / " + num_factors  +  "\t(" + name + ")         \r");
				//System.out.print("  [" + (++count_tfs) + "]\tPrecalculating secondary structure of " + name + "...");
				
				
				bw.write(line + "\n");					// name
				
				line = br.readLine(); 		
				
				bw.write(line + "\n");					// XX
				
			// parse domains
				
				while((line = br.readLine()).startsWith("DO  ")) {
						
					bw.write(line + "\n");				// DNA-binding domain(s)
					
				}
				
				bw.write(line + "\n");					// XX
				
				
				line = br.readLine();
				
				
			// parse first sequence	
				
				String seq1 = "";
				
				if(line.startsWith("S1")) {
					
					while(line.startsWith("S1")) {
						
						seq1 += line.substring(4);
						
						line = br.readLine();
						
					}
					
					String secstr1 = getSecondaryStructure(seq1);
					
					
				// write secondary structure string of this sequence to output file
					
					for(int i=0; i<(secstr1.length()/SEQLINELENGTH); i++) {
						
						bw.write("S1  "); 
						bw.write(secstr1.toUpperCase(), i*SEQLINELENGTH, SEQLINELENGTH);
						bw.write("\n");
					
					}
					
					if(secstr1.length()-(secstr1.length()/SEQLINELENGTH)*SEQLINELENGTH > 0) {
						
						bw.write("S1  "); 
						bw.write(secstr1.toUpperCase(), (secstr1.length()/SEQLINELENGTH)*SEQLINELENGTH, secstr1.length()-(secstr1.length()/SEQLINELENGTH)*SEQLINELENGTH);
						bw.write("\n");
					
					}
					
					bw.write(line + "\n");				// XX
					
					line = br.readLine();
					
				}
				
				
			// parse second sequence
				
				String seq2 = "";
					
				if(line.startsWith("S2")) {

					while(line.startsWith("S2")) {
						
						seq2 += line.substring(4);
						
						line = br.readLine();
						
					}
					
					
					String secstr2 = getSecondaryStructure(seq2);
					
					
				// write secondary structure string of this sequence to output file
						
					for(int i=0; i<(secstr2.length()/SEQLINELENGTH); i++) {
							
						bw.write("S2  "); 
						bw.write(secstr2.toUpperCase(), i*SEQLINELENGTH, SEQLINELENGTH);
						bw.write("\n");
						
					}
						
					if(secstr2.length()-(secstr2.length()/SEQLINELENGTH)*SEQLINELENGTH > 0) {
							
						bw.write("S2  "); 
						bw.write(secstr2.toUpperCase(), (secstr2.length()/SEQLINELENGTH)*SEQLINELENGTH, secstr2.length()-(secstr2.length()/SEQLINELENGTH)*SEQLINELENGTH);
						bw.write("\n");
						
					}
					
					
					bw.write(line + "\n");				// XX
					
					line = br.readLine();
					
				}
				
				
				while(!line.startsWith("//")) {
					
					bw.write(line + "\n");				// PFM(s)
					
					line = br.readLine();
				
				}
				
				bw.write(line + "\n");					// //
				
				
				line = br.readLine(); 
				
				bw.write(line + "\n");					// XX
				
				
				//System.out.println("done.");
				
				bw.flush();
				
			}
			
			br.close();
			bw.flush();
			bw.close();
		
		}
		catch(IOException ioe) {
			System.out.println(ioe.getMessage());
			System.out.println("IOException occurred while cutting dataset.");
		}
		System.out.println("\n");
	}
	
	
	/*
	 * 
	 *	predicts the secondary structure of a given sequence
	 *
	 *	sequence-string should contain ~60 chars per line
	 * 
	 */
	
	public String getSecondaryStructure(String seq) {
		
		BufferedReader br = null;
		BufferedWriter bw = null;
		
		String res = "";
		String line = null;
		String infile = "seq.fasta";
		
		// set temporary base directory and input file name
		String psidir = basedir;
		if (! psidir.equals("./PSIPRED/")) {
			psidir += "psipred/";
			
			infile = RandomStringGenerator.randomString(10) + ".fasta";
		}
		
		try {
			
			/*
			 * 
			 * generate input file
			 * 
			 */
			File input_file = new File(psidir + infile);
			bw = new BufferedWriter(new FileWriter(input_file));
			
			bw.write(">seq\n");
			
			bw.write(seq);
			bw.flush();
			bw.close();
			
			/*
			 * 
			 * generate shell script
			 * 
			 */
			
			File exec_file = new File(psidir + "launchPSIPRED");
			bw = new BufferedWriter(new FileWriter(exec_file));
			
			bw.write("cd ./PSIPRED\n");
			if (psidir.startsWith("/")) {
				bw.write("./runpsipred_single " + psidir + infile + "\n");
			}
			else {
				bw.write("./runpsipred_single ../" + psidir + infile + "\n");
			}
			bw.flush();
			bw.close();
			exec_file.setExecutable(true);
			
			/*
			 * 
			 * run external SS-prediction program (PSIPRED)
			 * 
			 */
			
			String cmdString = psidir + "launchPSIPRED";
		
			Process proc = Runtime.getRuntime().exec(cmdString);
	
			proc.waitFor();
			
			proc.destroy();
			
			/*
			 * 
			 * extract PSIPRED results from output file
			 * 
			 */
			
			String outfilename = infile.substring(0, infile.indexOf(".")) + ".horiz";
			String outfilename_ss = infile.substring(0, infile.indexOf(".")) + ".ss";
			String outfilename_ss2 = infile.substring(0, infile.indexOf(".")) + ".ss2";
			
			File outfile = new File("./PSIPRED/" + outfilename);
			File outfile_ss = new File("./PSIPRED/" + outfilename_ss);
			File outfile_ss2 = new File("./PSIPRED/" + outfilename_ss2);
			
			br = new BufferedReader(new FileReader(outfile));
			
			br.readLine();				// # PSIPRED HFORMAT (PSIPRED V2.6 by David Jones)
			br.readLine();				// 
			
			while((line = br.readLine()) != null) {			// Conf.
				
				line = br.readLine();						// Pred.
								
				res += line.substring(6);
				
				line = br.readLine();						// AA
				line = br.readLine();						// Pos.
				line = br.readLine();						//	
				line = br.readLine();						//
				
			}
			
			
			br.close();
			input_file.delete();
			exec_file.delete();
			outfile.delete();
			outfile_ss.delete();
			outfile_ss2.delete();
		
		}
		catch(IOException ioe) {
			System.out.println(ioe.getMessage());
			System.out.println("IOException occurred while calculating secondary structure.");
		}	
		catch (InterruptedException e) {
			System.out.println(e.getMessage());
			System.out.println("InterruptedException while executing PSIPRED!");
		}
		
		return res;
	}
	
	
	public void lookUpAllSpecies(String infile, String lookupfile, String outfile) {
		
		BufferedReader br = null;
		BufferedWriter bw = null;
		
		// System.out.println("  Looking up all species of the factors in the training set.\n");
		
		try {
			
			
			br = new BufferedReader(new FileReader(new File(infile)));
			bw = new BufferedWriter(new FileWriter(new File(outfile)));
		
			String line = null;
			
			StringTokenizer strtok = null;
			
			String tfid = null;
			
			
			while((line = br.readLine()) != null) {
				
				strtok = new StringTokenizer(line);
					
				tfid = strtok.nextToken();
				bw.write(lookUpSingleSpecies(tfid, lookupfile) + "\n");
			}
			
			br.close();
			bw.flush();
			bw.close();
			
		}	
		catch(IOException ioe) {
			System.out.println(ioe.getMessage());
			System.out.println("IOException occurred while looking up species.");
		}
		
		
	}
	
	public String lookUpSingleSpecies(String tfid, String lookupfile) {
		
		BufferedReader br = null;
		
		String spec = null;
		
		try {
			
			br = new BufferedReader(new FileReader(new File(lookupfile)));
			
			String line = null;
			
			
		// goto tf entry	
			
			while((line = br.readLine()) != null && !line.startsWith("NA  " + tfid));
			
			line = br.readLine();	// XX
			
			spec = br.readLine();
			
			
		// parse species	
				
			spec = spec.substring(4, spec.length());
				
			String[] split = spec.split(",");
				
			if(split.length == 1) {
			
				spec = split[0];
			
			}
			
			else {
					
				if(split.length == 2) {
					
					spec = split[1].substring(1,split[1].length());
					
				}
					
				else {
						
					System.out.println("Error while parsing species of " + tfid + "\nSpecies: " + spec);
				}
					
			}

			
			br.close();
			
		}
		catch(IOException ioe) {
			System.out.println(ioe.getMessage());
			System.out.println("IOException occurred while looking up species of " + tfid + ".");
		}
		
		return spec;
		
	}
	
	
	
	public void precalculateSVMPairwiseHelpFiles(String infile, String tfnamesfile, String outfile) {
		
		BufferedReader br = null;
		BufferedWriter bw = null;
		
		
		ArrayList<String> names = new ArrayList<String>();
		
		
		try {
			
			
		// read tf names	
			
			br = new BufferedReader(new FileReader(new File(tfnamesfile)));
			
			String line = null;
			
			while((line = br.readLine()) != null) {
				
				names.add(line.trim());
				
			}
			
			br.close();
			
		
		// init scores-matrix
			
			double[][] scores = new double[names.size()][names.size()];
			
			for(int i=0; i<scores.length; i++) {
				
				for(int j=0; j<scores.length; j++) {
					
					if(i == j) 	{ scores[i][j] = 1.0; }
					else		{ scores[i][j] = 0.0; } 
				
				}
				
			}
			
		
		// fill scores matrix	
			
			br = new BufferedReader(new FileReader(new File(infile)));
			
			StringTokenizer strtok = null;
			
			String tf1 = null;
			String tf2 = null;
			
			int index1 = 0;
			int index2 = 0;
			
			String score = null;
			
			while((line = br.readLine()) != null) {
				
				strtok = new StringTokenizer(line);
				
				tf1 = strtok.nextToken();
				
				strtok.nextToken();			// vs.
				
				tf2 = strtok.nextToken();
				
				strtok.nextToken();			// :
				
				score = strtok.nextToken();
				
				
				if((index1 = names.indexOf(tf1)) != -1 && (index2 = names.indexOf(tf2)) != -1) {
					
					scores[index1][index2] = Double.parseDouble(score);
					scores[index2][index1] = Double.parseDouble(score);
					
				}
				
				
			}
			
			br.close();
			
			
		// write results to output file
			
			bw = new BufferedWriter(new FileWriter(new File(outfile)));
			
			
			for(int i=0; i<scores.length; i++) {
				
				for(int j=0; j<scores.length; j++) {
					
					bw.write(names.get(i) + " vs. " + names.get(j) + " : " + scores[i][j] + "\n");
					
				}
				
			}
			
			bw.flush();
			bw.close();
			
		}
		catch(IOException ioe) {
			System.out.println(ioe.getMessage());
			System.out.println("IOException occurred while precalculating SVMPairwise Helpfiles.");
		}
	}
	
	
	
	public void adjustTrainingSet(String infile, String outfile) {
		
		LibSVMFeatureScaler scaler = new LibSVMFeatureScaler();
		
		scaler.adjustSingleFeatureFile(infile, outfile);
		
	}
	
	public int getNumFactors(String infile) {
        
        int entry_counter = 0;
        BufferedReader br = null;
              	
        try {
        	br = new BufferedReader(new FileReader(new File(infile)));
        	
        	while(br.readLine() != null) {
        			entry_counter++;
        	}
      
        } catch (IOException e) {
        	System.out.println(e.getMessage());
        	System.out.println("Error occured while counting entries in training set");
        }
        return entry_counter;
	}   
	
	public void computeDomainScores(String class_id, String train_dir, String matrix, String namesfile, String outfile) {
		
		System.out.println("  Calculating substitution matrix based alignment scores (matrix file: " + matrix + ").");
		
		DomainFeatureCalculator domaincalculator = new DomainFeatureCalculator();
		domaincalculator.silent = true;
		
		int num_entries = getNumFactors(namesfile);
	    boolean[] irrelevantPairs = new boolean[num_entries]; 
		  
	    domaincalculator.parseRelevantDomains(irrelevantPairs, class_id, train_dir);
	    
	    ArrayList<String> relevant_pairs = new ArrayList<String>();
	    for (int i=0; i<num_entries-1; i++) {
	    	for (int j=i+1; j<num_entries; j++) { 
	    		
	    		relevant_pairs.add(i + "_" + j);
	    	}
	    }
	    
	    domaincalculator.calculateAllDomainFeatures("SMBasedSimilarity", null, matrix , outfile, relevant_pairs);
	}
	
	public void extractFBPs(String class_id, String infile, String outfile) {
		
		// parse factors
		TransfacParser tf_parser = new TransfacParser();
		tf_parser.parseFactors(infile);
		tf_parser.filterFactorsBySuperclass(class_id.substring(5));
		
		// convert and merge PWMs
		PFMFormatConverter pfm_converter = new PFMFormatConverter();
		pfm_converter.basedir = basedir;
		
		ArrayList<ArrayList<String>> stamp_pfms = new ArrayList<ArrayList<String>>();
		ArrayList<ArrayList<String>> merged_pfms = new ArrayList<ArrayList<String>>();
		ArrayList<String> merged_pfm = new ArrayList<String>();
		ArrayList<String> factor_names = new ArrayList<String>();
		//int matrix_counter = 0;
		
		System.out.println("  Converting PFMs to STAMP format.");
		for (int i=0; i<tf_parser.get_tf_names().size(); i++) {
			 
			// System.out.print("    Converting matrix \t" + ++matrix_counter + " / " + factor_counter + "\r");
			
			// more than one PFM ?
			if (tf_parser.get_pfms().get(i).size() > 1) {
				
				stamp_pfms = pfm_converter.convertAllTransfacToSTAMP(tf_parser.get_pfms().get(i));
				merged_pfm = pfm_converter.mergePFMs(stamp_pfms);
			}
			else {
				merged_pfm = pfm_converter.convertTransfacToSTAMP(tf_parser.get_pfms().get(i).get(0));
			}
			merged_pfms.add(merged_pfm);
			factor_names.add(tf_parser.get_tf_names().get(i));
		}
		
		// write FBP file
		try {
			BufferedWriter bw = new BufferedWriter(new FileWriter(new File(outfile)));
			
			System.out.println("\n  Merging PFMs.\n");
			for(int i=0; i<merged_pfms.size(); i++) {
					
				bw.write("DE\t" + factor_names.get(i) + "\n");

				for(int j=0; j<merged_pfms.get(i).size(); j++) {
						
					bw.write(merged_pfms.get(i).get(j) + "\n");
						
				}	
				bw.write("XX\n");
			}
			bw.flush();
			bw.close();
		}
		
		
		catch(IOException ioe) {
			System.out.println(ioe.getMessage());
			System.out.println("IOException occurred while writing FBPs to training set file.");
		}
	}
	
	
	public void prepareTrainingDirectory(String train_dir, String class_id) {
		
		// copy missing files to training directory
		if (! (new File(train_dir + "Classes")).exists()) {
			FileCopier.copy("trainingsets/Classes", train_dir + "Classes");
		}
		if (! (new File(train_dir + "new_phylogenetic_distances.out")).exists()) {
			FileCopier.copy("trainingsets/new_phylogenetic_distances.out", train_dir + "new_phylogenetic_distances.out");
		}
		
		// create directory for temporary files
		basedir = train_dir + "temp_" + class_id + "/";
		String psipred_dir = basedir + "psipred/";
		String stamp_dir = basedir + "stamp/";
		String mosta_dir = basedir + "mosta/";
		
		if ((! new File(basedir).exists() && ! new File(basedir).mkdir()) |
			(! new File(mosta_dir).exists() && ! new File(mosta_dir).mkdir()) |
			(! new File(psipred_dir).exists() && ! new File(psipred_dir).mkdir()) |
			(! new File(stamp_dir).exists() && ! new File(stamp_dir).mkdir())) {
			
			System.out.println("\nInvalid training directory. Aborting.");
			System.out.println("Training directory: " + train_dir + "\n");
			System.exit(0);
		}
	}
	
	
	public void extractSequencesAndMatrices(String class_id, String infile, String outfile) {
		
		// parse and filter factors
		TransfacParser tf_parser = new TransfacParser();
		tf_parser.silent = false;
		
		System.out.println("  Parsing factors from input file.");
		tf_parser.parseFactors(infile);
		tf_parser.filterFactorsBySuperclass(class_id.substring(5));
		
		// write factors to output file
		tf_parser.writeTrainingSetFile(outfile);
	}
	

	public void runPreprocessor(String input_file, String class_id, String train_dir) {
		
		RawDataPreprocessor preprocessor = new RawDataPreprocessor();
		
		/*
		 * missing files:
		 * 		- sequences_matrices_class0.out			
		 * 		- FBPs_class0.out						done !
		 * 		- domain_scores_BLOSUM_62_class0.out 	done !
		 * 		- domain_scores_PAM_80_class0.out		done !
		 * 		- class0.30featurenames
		 * 		- mosta_30features_class0.out
		 */
		
		preprocessor.prepareTrainingDirectory(train_dir, class_id);
		
		
		preprocessor.extractSequencesAndMatrices(class_id, input_file, train_dir + "sequences_matrices_" + class_id + ".out");

		preprocessor.sortDataset(train_dir + "sequences_matrices_" + class_id + ".out", train_dir + "/trainingset_" + class_id + ".tfnames", train_dir + "/trainingset_" + class_id + ".rawdata");
		
		preprocessor.extractFBPs(class_id, input_file, train_dir + "FBPs_" + class_id + ".out");
		
		preprocessor.sortFBPs(train_dir + "FBPs_" + class_id + ".out", train_dir + "trainingset_" + class_id + ".tfnames", train_dir + "trainingset_" + class_id + ".fbps");
		
		preprocessor.precalculateSecondaryStructures(train_dir + "trainingset_" + class_id + ".rawdata", train_dir + "trainingset_" + class_id + ".secstruct", train_dir + "/trainingset_" + class_id + ".tfnames");
		
		preprocessor.lookUpAllSpecies(train_dir + "/trainingset_" + class_id + ".tfnames", input_file, train_dir + "trainingset_" + class_id + ".species");
	
		preprocessor.computeDomainScores(class_id, train_dir, "BLOSUM_62.dat", train_dir + "trainingset_" + class_id + ".tfnames", train_dir + "domain_scores_BLOSUM_62_" + class_id + ".out");
		preprocessor.computeDomainScores(class_id, train_dir, "PAM_080.dat", train_dir + "trainingset_" + class_id + ".tfnames", train_dir + "domain_scores_PAM_080_" + class_id + ".out");
		
		preprocessor.precalculateSVMPairwiseHelpFiles(train_dir + "domain_scores_BLOSUM_62_" + class_id + ".out", train_dir + "trainingset_" + class_id + ".tfnames", train_dir + "trainingset_" + class_id + ".blo62");
		preprocessor.precalculateSVMPairwiseHelpFiles(train_dir + "domain_scores_PAM_080_" + class_id + ".out", train_dir + "trainingset_" + class_id + ".tfnames", train_dir + "trainingset_" + class_id + ".pam80");
		
		System.out.println();
		// preprocessor.adjustTrainingSet(train_dir + "mosta_30features_" + class_id + ".out", train_dir + "feature_scaling_orientation_" + class_id + ".out");
		
	}
	/*
	 *  default options 
	 */
	
	public static void main(String[] args) {
		
		String inputfile = null;
		String traindir = "trainingsets/";
		int superclass = 1;
		
		if (args.length < 1) usage();
		
		if (args.length == 1 && (args[0].equals("-help") || args[0].equals("--help"))) {
			usage();
		}
		
		inputfile = args[0];
	
		for(int i=1; i<args.length-1; i+=2) {
		
			if(args[i].equals("-d")) { traindir 		 = args[i+1]; 						continue; }
			if(args[i].equals("-c")) { superclass 		 = Integer.parseInt(args[i+1]);		continue; }
		
			if( !args[i].equals("-d") && !args[i].equals("-c")) {	
				
				System.out.println("\n  Invalid argument: " + args[i]);
				usage();
			}
		}	
		
		if (! traindir.endsWith("/")) {
			traindir += "/";
		}
		
		System.out.println();
		System.out.println("  ----------------------------------------------------------------------------");
		System.out.println("  Training Set Generator for SABINE - StandAlone BINding specificity Estimator");
		System.out.println("  ----------------------------------------------------------------------------");
		System.out.println("\n");
		System.out.println("  Superclass:  " + superclass + "\n");
	
		RawDataPreprocessor processor = new RawDataPreprocessor();
		processor.runPreprocessor(inputfile, "class" + superclass, traindir);
	}
	
	private static void usage() {
			
		System.out.println();
		System.out.println("  ----------------------------------------------------------------------------");
		System.out.println("  Training Set Generator for SABINE - StandAlone BINding specificity Estimator");
		System.out.println("  ----------------------------------------------------------------------------");
		System.out.println("\n");
		System.out.println("  Usage   : traintool <input_file> [OPTIONS]\n");
		System.out.println("  OPTIONS : -d <training_set_dir>      (directory to save training set)");
		System.out.println("            -c <transfac_superclass>   (e.g. 0, 1, 2, 3, 4 )                    default = 1\n\n");

			
		System.exit(0);
			
	}
}

