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

package main;

import java.awt.image.BufferedImage;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.StringTokenizer;

import libsvm.LabeledTF;

import org.biojava.bio.gui.DistributionLogo;
import org.biojava.bio.seq.ProteinTools;
import org.biojava.bio.symbol.IllegalSymbolException;

import extension.DirectoryRemover;
import extension.UniProtClient;
import gui.SequenceLogo;
import help.FileCopier;

public class SABINE_Galaxy {
   
	boolean silent = true;
	boolean batchMode = false;
	
	private static final String[] superclassNames = new String[] {"Other", "Basic domain", "Zinc finger", "Helix-turn helix", "Beta scaffold", "Auto-detect"};
	private static final int unknownProteinClass = -1;
	
	int tf_num_limit = 10;
	
	String organism;
	int superclass;
	String sequence;
	ArrayList<String> domains = new ArrayList<String>();
	
	String basedir;
	String inputfile;
	String outputfile;
	String sabineInputFile;
	String sabineOutputFile;
	String trainingset = FBPPredictor.public_trainingset;
	boolean dyn_bmt = true;
	double bmt;
	int mnb; 
	double oft;
	
	// objects to save results
	ArrayList<String> inputTFnames = new ArrayList<String>();
	ArrayList<Boolean> predictionPossible = new ArrayList<Boolean>();
	ArrayList<Integer> numCandidates = new ArrayList<Integer>();
	ArrayList<ArrayList<LabeledTF>> bestMatches = new ArrayList<ArrayList<LabeledTF>>();
	ArrayList<ArrayList<Double[]>> predPFM = new ArrayList<ArrayList<Double[]>>();
	ArrayList<DistributionLogo[]> logoArray = new ArrayList<DistributionLogo[]>();
	ArrayList<BufferedImage> seqLogo = new ArrayList<BufferedImage>();
	
	public void parseInput(String[] args) {
		
		basedir = args[0];
		sabineInputFile = basedir + "infile.tmp";
		int i = 1;
		String curr_dom;
		
		while (i < args.length) {
			
			if (!batchMode && args[i].equals("-o")) {
				organism = args[++i].toUpperCase();
				while (! args[++i].equals("-p")) {
					organism += " " + args[i].toUpperCase();
				}
			}
			
			if (!batchMode && args[i].equals("-p")) {
				sequence = "";
				while (! args[++i].equals("-s")) {
					sequence += args[i].replaceAll("X", "").toUpperCase();
				}
			}
			
			if (!batchMode && args[i].equals("-u")) {
				String uniprot_id = args[++i].toUpperCase();
				UniProtClient uniprot_client = new UniProtClient();
				String fasta_seq = uniprot_client.getUniProtSequence(uniprot_id.toUpperCase(), true);
				
				organism = fasta_seq.substring(fasta_seq.indexOf("OS=")+3, fasta_seq.indexOf("GN=")-1);
				sequence = fasta_seq.replaceFirst(">.*\\n", "").replaceAll("\\n", "");
				i++;
			}
			
			if (!batchMode && args[i].equals("-s")) {
				superclass = Integer.parseInt(args[++i]);
				i++;
			}
			
			if (batchMode && args[i].equals("-i")) {
				inputfile = args[++i];
				i++;
			}
			
			if (batchMode && args[i].equals("--sabine-output-file")) {
				sabineOutputFile = args[++i];
				i++;
			}
			
			if (args[i].equals("-f")) {
				outputfile = args[++i];
				i++;
			}
			
			if (!batchMode && args[i].equals("-d")) {
				curr_dom = args[i+1] + "  " + args[i+2];
				domains.add(curr_dom);
				i = i + 3;
			}
			
			if (args[i].equals("-b")) {
				dyn_bmt = false;
				bmt = Double.parseDouble(args[++i]);
				i++;
			}
			
			if (args[i].equals("-m")) {
				mnb = Integer.parseInt(args[++i]);
				i++;
			}
			
			if (args[i].equals("-t")) {
				oft = Double.parseDouble(args[++i]);
				i++;
			}
			
			if (i < args.length && args[i].equals("--biobase-data")) {
				trainingset = FBPPredictor.biobase_trainingset;
				tf_num_limit = Integer.MAX_VALUE;
				i++;
			}
		
			if (i < args.length && args[i].equals("--batch-mode")) {
				break;
			}
		}
		
		if (! silent) {
			if (batchMode) {
				System.out.println("Input File:         " + inputfile);
				System.out.println("HTML report:        " + outputfile);
				System.out.println("SABINE output file: " + sabineOutputFile);
				
			} else {
			
				System.out.println("Organism:   " + organism);
				System.out.println("Superclass: " + superclassNames[superclass]);
				System.out.println("Sequence:   " + sequence);
				if (! domains.isEmpty()) System.out.println("Domains:    " + domains.get(0));
				for (int j=1; j<domains.size(); j++) {
					System.out.println("            " + domains.get(j));
				}
			}
			if (!dyn_bmt) {
				System.out.println("BMT:        " + bmt);
			} else {
				System.out.println("BMT:        dynamic");
			}
			System.out.println("MNB:        " + mnb);
			System.out.println("OFT:        " + oft + "\n");
		}
	}
	
	public void checkInput() {
		
		// check organism
		if (organism.equals("?")) { 
			System.out.println("Error. No organism selected.");
			System.exit(0);
		}
		
		// check protein sequence
		try {
			ProteinTools.createProtein(sequence);
			
		} catch (IllegalSymbolException e) {
			System.out.println("Error. Illegal symbols in protein sequence.");
			System.exit(0);
		}
		
		// check domains
		int domStart, domEnd;
		StringTokenizer strtok;
		if (domains.isEmpty()) {
			System.out.println("Error. No DNA-binding domain specified.");
			System.exit(0);
		}
		else {
			
			for (int i=0; i<domains.size(); i++) {
				strtok = new StringTokenizer(domains.get(i));
				domStart = Integer.parseInt(strtok.nextToken());
				domEnd = Integer.parseInt(strtok.nextToken());
				
				if ((domEnd < domStart) || domStart < 0 || domEnd < 0 || domEnd > sequence.length()) {
					System.out.println("Error. Invalid DNA-binding domain.");
					System.exit(0);
				}
			}
		}
		
		// check parameters
		if (bmt < 0 || bmt > 1) {
			System.out.println("Error. Invalid Value for Best Match Threshold. Choose value between 0 and 1.");
			System.exit(0);
		}
		if (mnb < 1 || mnb > 1000) {
			System.out.println("Error. Invalid Value for Max. number of best matches. Choose value between 1 and 1000.");
			System.exit(0);
		}
		if (oft < 0 || oft > 1) {
			System.out.println("Error. Invalid Value for Best Match Threshold. Choose value between 0 and 1.");
			System.exit(0);
		}
	} 
	
	public int countTFs(String infile) {
		
		BufferedReader br;
		String line;
		int tf_counter = 0;
		
		try {
			br = new BufferedReader(new FileReader(new File(infile)));
	
			while ((line = br.readLine()) != null) {
				if (line.startsWith("NA")) 
					tf_counter++;
			}
			br.close();
			 
		} catch (IOException e) {
			e.printStackTrace();
		}
		return tf_counter;
	}
	
	public void createBasedir() {
		
		
		if (basedir == null) {

			File basedir_path = new File(basedir);
			
			if (!basedir_path.exists() && !basedir_path.mkdir()) {
				System.out.println("\nInvalid base directory. Aborting.");
				System.out.println("Base directory: " + basedir + "\n");
				System.exit(0);
			}
		}
		
		if (!silent) {
			System.out.println("Basedir:    " + basedir);
		}
	}
	
	public void deleteBasedir() {
		DirectoryRemover.removeDirectory(basedir);
	}
	
	public void writeInputFile(String infile) {

		BufferedWriter bw;
		try {
			bw = new BufferedWriter(new FileWriter(new File(infile)));

			bw.write("NA  QueryTF\n");
			bw.write("XX\n");
			bw.write("SP  " + organism + "\n");
			bw.write("XX\n");
			if (superclass < 5) {
				bw.write("CL  " + superclass + ".0.0.0.0.\n");
				bw.write("XX\n");
			} 
				
			// write sequence
			int SEQLINELENGTH = 60;

			for (int i = 0; i < (sequence.length() / SEQLINELENGTH); i++) {
				bw.write("S1  ");
				bw.write(sequence.toUpperCase(), i * SEQLINELENGTH, SEQLINELENGTH);
				bw.write("\n");
			}

			if (sequence.length() - (sequence.length() / SEQLINELENGTH) * SEQLINELENGTH > 0) {
				bw.write("S1  ");
				bw.write(sequence.toUpperCase(), (sequence.length() / SEQLINELENGTH)
						* SEQLINELENGTH, sequence.length()
						- (sequence.length() / SEQLINELENGTH) * SEQLINELENGTH);
				bw.write("\n");
			}
			bw.write("XX\n");

			// write domains
			for (int i = 0; i < domains.size(); i++) {
				bw.write("FT  POS  " + domains.get(i) + "\n");
			}
			bw.write("XX\n");
			bw.write("//\n");
			bw.write("XX\n");

			bw.flush();
			bw.close();

		} catch (IOException e) {
			e.printStackTrace();
		}
		
		if (! silent) {
			System.out.println("Input file: " + infile);
		}
	}
	
	
	public void runSabine() {
		
		// create directories for temporary files
		SABINE_Caller sabine_caller = new SABINE_Caller();
		sabine_caller.createTempDirectories(basedir);

		// run SABINE on generated input file
		if (batchMode) {
			// set parameters
			sabine_caller.dynamic_threshold = dyn_bmt;
			if (!dyn_bmt) sabine_caller.best_match_threshold = bmt;
			sabine_caller.max_number_of_best_matches = mnb;
			sabine_caller.outlier_filter_threshold = oft;
			sabine_caller.silent = silent;
			
			sabine_caller.launch_SABINE(sabineInputFile, basedir + "outfile.tmp", "n", basedir, trainingset, null);
		}
		else { 
			// set parameters
			FBPPredictor predictor = new FBPPredictor();
			predictor.dynamic_threshold = dyn_bmt;
			if (!dyn_bmt) sabine_caller.best_match_threshold = bmt;
			predictor.max_number_of_best_matches = mnb;
			predictor.outlier_filter_threshold = oft;
			predictor.silent = silent;
			
			predictor.predictFBP(sabineInputFile, basedir, trainingset, null);
		}
	}
	
	
	public void readOutfile(String infile) {

		try {
			BufferedReader br = new BufferedReader(new FileReader(new File(infile)));
			
			StringTokenizer strtok;
			String line, name;
			double A, C, G, T, pos, score;
			
			while ((line=br.readLine()) != null) {
				
				if (line.startsWith("//") || line.startsWith("XX")) {
					continue;
				}
				
				// read name of Input TF
				if (batchMode) {
					if (line.startsWith("NA")) {
						inputTFnames.add(line.replaceFirst("NA  ", "").trim());
						line = br.readLine();	// XX
						line = br.readLine();	// BM
						
					} else {
						
						System.out.println("Parse Error. \"NA\" expected.");
						System.out.println("Line: " + line);
						System.exit(0);
					}
				}
				
				// read Best Matches
				ArrayList<LabeledTF> currBestMatches = new ArrayList<LabeledTF>();
				
				boolean first = true;
				while (line.startsWith("BM")) {
					
					if (line.startsWith("BM  none")) {
						if (line.startsWith("BM  none (Unknown")) {
							numCandidates.add(unknownProteinClass);
							predictionPossible.add(false);
							
						} else {
							numCandidates.add(Integer.parseInt(line.substring(10, line.indexOf(" ", 10))));
							predictionPossible.add(false);
						}
						
					} else {
						if (first) {
							numCandidates.add(0);
							predictionPossible.add(true);
							first = false;
						}

						// add Best Match to list for current TF 
						strtok = new StringTokenizer(line.substring(4));
						name = strtok.nextToken();
						score = Double.parseDouble(strtok.nextToken());
						currBestMatches.add(new LabeledTF(name, score));
					}
					line = br.readLine();
				}
				bestMatches.add(currBestMatches);
				
				line = br.readLine();	// XX
				
				// read PFM
				ArrayList<Double[]> currPFM = null;
				ArrayList<DistributionLogo> currLogoData = null;
				
				if (!line.startsWith("MA  none")) {
					currPFM = new ArrayList<Double[]>();
					currLogoData = new ArrayList<DistributionLogo>();
					
					while (line.startsWith("MA")) {
						
						strtok = new StringTokenizer(line.substring(4));
						
						pos = Integer.parseInt(strtok.nextToken()); 
						A = Double.parseDouble(strtok.nextToken()) / 100;
						C = Double.parseDouble(strtok.nextToken()) / 100;
						G = Double.parseDouble(strtok.nextToken()) / 100;
						T = Double.parseDouble(strtok.nextToken()) / 100;
						
						currPFM.add(new Double[] {pos, A, C, G, T});
						currLogoData.add(SequenceLogo.getLogo(new double[] {A, C, G, T}));
		
						line = br.readLine();
					}
				}
				predPFM.add(currPFM);
				
				DistributionLogo[] currLogoArray = null;
				BufferedImage currLogo = null;
				if (currLogoData != null) {
					currLogoArray = currLogoData.toArray(new DistributionLogo[] {});
					currLogo = SequenceLogo.drawSequenceLogo(currLogoArray, 200);
				} 
				logoArray.add(currLogoArray);
				seqLogo.add(currLogo);
				
				
				
			}

		} catch (IOException e) {
			e.printStackTrace();
		} catch (Exception e1) {
			e1.printStackTrace();
		}
	}
	
	// write HTML header
	public static void writeHTMLheader(BufferedWriter bw) throws IOException {
		
		bw.write("<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.01 Transitional//EN\"\n");
		bw.write("     \"http://www.w3.org/TR/html4/loose.dtd\">\n");
		bw.write("<html>\n");
		bw.write("<head>\n");
		bw.write("<title>SABINE Result</title>\n");
		bw.write("<style type=\"text/css\">\n");
		bw.write("  h1 { font-size: 150%;color: #002780; }\n");
		bw.write("  h2 { font-size: 135%;color: #002780; }\n");
		bw.write("  table { width: 500px; background-color: #E6E8FA; border: 1px solid black; padding: 3px; vertical-align: middle;}\n");
		bw.write("  tr.secRow { background-color:#FFFFFF; margin-bottom: 50px; vertical-align: middle;}\n");
		bw.write("  th { padding-bottom: 8px; padding-top: 8px; text-align: center}\n");
		bw.write("  td { padding-bottom: 8px; padding-top: 8px; text-align: center}\n");
		bw.write("</style>\n");
		bw.write("</head>\n");
		bw.write("<body style=\"padding-left: 30px\">\n");
	}

	
	public void writeHTMLOutput(String result_file) {

		DecimalFormat fmt2 = new DecimalFormat();
		fmt2.setMaximumFractionDigits(2);
		fmt2.setMinimumFractionDigits(2);
		
		DecimalFormat fmt4 = new DecimalFormat();
		fmt4.setMaximumFractionDigits(4);
		fmt4.setMinimumFractionDigits(4);
		
		// write best matches and PFM to HTML tables
		try {
			BufferedWriter bw = new BufferedWriter(new FileWriter(new File(outputfile)));
		
			writeHTMLheader(bw);
			
			int numInputTFs = 1;
			if (batchMode) {
				numInputTFs = inputTFnames.size();
			}
			
			for (int tfIndex=0; tfIndex<numInputTFs; tfIndex++) {
				
				if (tfIndex > 0) {
					bw.write("<br><hr>\n\n");
				}
			
				if (batchMode) {
					bw.write("<h1><span style=\"color:#000000\">Results report: </span>" + inputTFnames.get(tfIndex) + "</h1>\n");
				}
				
				// write SUCCESS output 
				if (predictionPossible.get(tfIndex)) {
				
					if (!batchMode) {
						System.out.println("PFM transfer was possible.");
					}
					
					// write best matches
					bw.write("<h2>Best Matches</h2>\n");
					bw.write("<table>\n");
					bw.write("  <tr><th> Transcription Factor </th><th> PFM similarity </th></tr>\n");
					
					double lowest_score = 0;
					for (int i=0; i<bestMatches.get(tfIndex).size(); i++) {
						if ((i%2) == 1) { 
							bw.write("  <tr>");
						}
						else {
							bw.write("  <tr class=\"secRow\">");
						}
						bw.write("<td> " + bestMatches.get(tfIndex).get(i).getName() + " </td>");
						bw.write("<td> " + fmt4.format(bestMatches.get(tfIndex).get(i).getLabel()) + " </td></tr>\n");
						
						lowest_score = bestMatches.get(tfIndex).get(i).getLabel();
					}
					bw.write("</table>\n\n");
					bw.write("<br><br>\n\n");
					
					// write PFM
					bw.write("<h2>Predicted PFM</h2>\n");
					bw.write("<table>\n");
					bw.write("  <tr><th> Pos. </th><th> A </th><th> C </th><th> G </th><th> T </th></tr>\n");
					
					for (int i=0; i<predPFM.get(tfIndex).size(); i++) {
						if ((i%2) == 1) { 
							bw.write("<tr>");
						}
						else {
							bw.write("<tr class=\"secRow\">");
						}
						bw.write("<td> " +  Math.round(predPFM.get(tfIndex).get(i)[0]) + " </td>");
						bw.write("<td> " +  fmt2.format(predPFM.get(tfIndex).get(i)[1]) + " </td>");
						bw.write("<td> " +  fmt2.format(predPFM.get(tfIndex).get(i)[2]) + " </td>");
						bw.write("<td> " +  fmt2.format(predPFM.get(tfIndex).get(i)[3]) + " </td>");
						bw.write("<td> " +  fmt2.format(predPFM.get(tfIndex).get(i)[4]) + " </td></tr>\n");
					}
					bw.write("</table>\n");
					bw.write("<br>\n\n");
					
					// get threshold values for dynamic BMT
					double[] thresholds = FBPPredictor.getThresholdValues();
					double high_conf_bmt = thresholds[0];
					double medium_conf_bmt = thresholds[1];
					
					// append confidence sign
					String conf_img_path = "../../../static/images/static_code/sabine_static/";
					
					if (lowest_score > high_conf_bmt) {
						conf_img_path = conf_img_path + "high_conf.png";
					} else if (lowest_score > medium_conf_bmt) {
						conf_img_path = conf_img_path + "medium_conf.png";
					} else {
						conf_img_path = conf_img_path + "low_conf.png";
					}
					bw.write("<img src=\"" + conf_img_path + "\" height=\"50px\">\n\n");
					bw.write("<br><br>\n\n");
					
					// append sequence logo
					String png_file;
					if (batchMode) {
						png_file = basedir + "seq_logo_" + inputTFnames.get(tfIndex) + ".png";
					} else {
						png_file = basedir + "seq_logo.png";
					}
					javax.imageio.ImageIO.write(SequenceLogo.drawSequenceLogo(logoArray.get(tfIndex), 100), "png", new File(png_file));
					png_file = png_file.substring(png_file.lastIndexOf('/')+1, png_file.length());
					
					bw.write("<h2>Sequence logo</h2>\n");
					bw.write("<img src=\"" + png_file + "\" height=\"100px\">\n\n");		
					
				// write FAILURE output 
				} else {
					
					if (! batchMode) {
						System.out.println("PFM transfer was not possible.");
					}
					
					bw.write("<h2>Best Matches</h2>\n");
					double thres = bmt;
					if (dyn_bmt) {
						bmt = FBPPredictor.low_conf_bmt;
					}
					if (numCandidates.get(tfIndex) == unknownProteinClass) {
						bw.write("<h3>Input protein is not labeled as Transcription Factor.</h3>\n" +
								"The input protein was either predicted as a Non-TF or could not be classified using the" +
								" tool TFpredict. As SABINE requires a transcription factor as input, no prediction was made for this protein.");
						
					} else if (numCandidates.get(tfIndex) == 0) {
						bw.write("<h3>No candidates for PFM transfer were found.</h3>\nPlease note, that SABINE " +
								 "can only predict a PFM for a given transcription factor if at least one " +
								 "candidate factor is found in the training set which has a normalized domain sequence similarity score " +
								 "of at least 0.3 with respect to the BLOSUM62 substitution matrix. As no candidate factor with sufficient domain " +
								 "sequence similarity was found, SABINE could not predict a PFM.");
					
					} else {
						bw.write("<h3>No Best Matches were found.</h3>\nSABINE compared the input transcription factor to all transcription factors with known DNA motif in the training set. " +
								 numCandidates.get(tfIndex) + " candidate factors with sufficiently high domain sequence similarity to the input factor were found. " +
								 "Unfortunately, none of these candidate factors was predicted to have a PFM similarity greater than the best match threshold (" + thres + "). " +
								 "Please note, that this threshold can be adjusted in the options of SABINE.");
					}
				}
			}
			
			// close HTML file
			bw.write("</body>\n");
			bw.write("</html>\n");
			
			bw.flush();
			bw.close();

		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	public void writeFailureOutput(String result_file) {
		
		try {
			
			BufferedWriter bw = new BufferedWriter(new FileWriter(new File(result_file)));
			writeHTMLheader(bw);
		

			// close HTML file
			bw.write("</body>\n");
			bw.write("</html>\n");
			
			bw.flush();
			bw.close();

		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	public boolean checkBatchModeLimit() {
		int num_factors = countTFs(inputfile);
		
		if (num_factors > tf_num_limit) {
			System.out.println("Number of transcription factors in input file must not exceed " + tf_num_limit + ".");
			
			// generate error message for SABINE output file 
			try {
				BufferedWriter bw = new BufferedWriter(new FileWriter(new File(outputfile)));
				
				bw.write("Limit for number of transcription factors in input file exceeded\n");
				bw.write("================================================================\n");
				bw.write("Limit for number of transcription factors : " + tf_num_limit + "\n");
				bw.write("Number of factors in given input file     : " + num_factors);
				
				bw.flush();
				bw.close();
				
			} catch (IOException e) {
				e.printStackTrace();
			}
			
			// TODO: generate error message for HTML output file 
			
			return(true);
			
		} else {
			return(false);
		}
	}
	
	public static void main(String[] args) {
		
		SABINE_Galaxy starter = new SABINE_Galaxy();
		
		if (args[args.length-1].equals("--batch-mode")) {
			starter.batchMode = true;
		}
		
		// parse input passed from Galaxy webinterface
		starter.parseInput(args);
		if (! starter.batchMode) {
			starter.checkInput();
			starter.writeInputFile(starter.sabineInputFile);
		}
		else {
			
			// quit if too many Input TFs are given
			boolean tooManyInputTFs = starter.checkBatchModeLimit();
			if (tooManyInputTFs) {
				return;
			}
			
			// move input file to temporary working directory
			FileCopier.copy(starter.inputfile, starter.sabineInputFile);
		}
		starter.runSabine();
		
		String sabineResultFile;
		if (starter.batchMode) {
			sabineResultFile = starter.basedir + "outfile.tmp";
			
			// move output file from base directory to destination given from Galaxy
			FileCopier.copy(sabineResultFile, starter.sabineOutputFile);
		}
		else {
			sabineResultFile = starter.basedir + "prediction.out";
		}
		
		starter.readOutfile(sabineResultFile);
		starter.writeHTMLOutput(starter.outputfile);
	}
}

