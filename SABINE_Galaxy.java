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
import gui.SequenceLogo;
import help.FileCopier;

public class SABINE_Galaxy {
   
	boolean silent = true;
	boolean batchMode = false;
	
	final String public_trainingset = "trainingsets_public/";
	final String biobase_trainingset = "trainingsets_biobase/";
	
	int tf_num_limit = 10;
	
	String organism;
	int superclass;
	String sequence;
	ArrayList<String> domains = new ArrayList<String>();
	
	String basedir;
	String inputfile;
	String outputfile;
	String trainingset = public_trainingset;
	boolean dyn_bmt = true;
	double bmt;
	int mnb; 
	double oft;
	
	DistributionLogo[] logoArray;
	BufferedImage seqLogo;
	
	ArrayList<LabeledTF> bestMatches= new ArrayList<LabeledTF>();
	ArrayList<Double[]> predPFM = new ArrayList<Double[]>(); 
	
	public void parseInput(String[] args) {
		
		basedir = args[0];
		int i = 1;
		String curr_dom;
		
		while (i < args.length) {
			
			if (!batchMode && args[i].equals("-o")) {
				organism = args[++i].toUpperCase();
				while (! args[++i].equals("-s")) {
					organism += " " + args[i].toUpperCase();
				}
			}
			
			if (!batchMode && args[i].equals("-s")) {
				superclass = Integer.parseInt(args[++i]);
				i++;
			}
			
			if (!batchMode && args[i].equals("-p")) {
				sequence = "";
				while (! args[++i].equals("-f")) {
					sequence += args[i].replaceAll("X", "").toUpperCase();
				}
			}
			
			if (batchMode && args[i].equals("-i")) {
				inputfile = args[++i];
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
				trainingset = biobase_trainingset;
				tf_num_limit = Integer.MAX_VALUE;
				i++;
			}
		
			if (i < args.length && args[i].equals("--batch-mode")) {
				break;
			}
		}
		
		if (! silent) {
			System.out.println("Organism:   " + organism);
			System.out.println("Superclass: " + superclass);
			System.out.println("Sequence:   " + sequence);
			if (! domains.isEmpty()) System.out.println("Domains:    " + domains.get(0));
			for (int j=1; j<domains.size(); j++) System.out.println("            " + domains.get(j));
			System.out.println("BMT:        " + bmt);
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
			sabine_caller.silent = true;
			
			sabine_caller.launch_SABINE(basedir + "infile.tmp", basedir + "outfile.tmp", "n", basedir, trainingset, null);
		}
		else { 
			// set parameters
			FBPPredictor predictor = new FBPPredictor();
			predictor.best_match_threshold = bmt;
			predictor.max_number_of_best_matches = mnb;
			predictor.outlier_filter_threshold = oft;
			predictor.silent = true;
			
			predictor.predictFBP(basedir + "infile.tmp", basedir, trainingset, null);
		}
	}
	
	
	public boolean readOutfile(String infile) {
		
		try {
			BufferedReader br = new BufferedReader(new FileReader(new File(infile)));
			
			StringTokenizer strtok;
			String line, name;
			double A, C, G, T, pos, score;
			ArrayList<DistributionLogo> logoData = new ArrayList<DistributionLogo>();
			
			// skip name
			if (! silent) System.out.println("Output file:" + infile);
			while ((line = br.readLine()).startsWith("BM")) {
				
				if (line.startsWith("BM  none"))
					return false;
				else { 
					strtok = new StringTokenizer(line.substring(4));
					name = strtok.nextToken();
					score = Double.parseDouble(strtok.nextToken());
					bestMatches.add(new LabeledTF(name, score));
				}
			}
			line = br.readLine();	// XX
			
			while (line.startsWith("MA")) {
				
				strtok = new StringTokenizer(line.substring(4));
				
				pos = Integer.parseInt(strtok.nextToken()); 
				A = Double.parseDouble(strtok.nextToken()) / 100;
				C = Double.parseDouble(strtok.nextToken()) / 100;
				G = Double.parseDouble(strtok.nextToken()) / 100;
				T = Double.parseDouble(strtok.nextToken()) / 100;
				
				predPFM.add(new Double[] {pos, A, C, G, T});
				logoData.add(SequenceLogo.getLogo(new double[] {A, C, G, T}));

				line = br.readLine();
			}
			logoArray = logoData.toArray(new DistributionLogo[] {});
			seqLogo = SequenceLogo.drawSequenceLogo(logoArray, 200); 

		} catch (IOException e) {
			e.printStackTrace();
		} catch (Exception e1) {
			e1.printStackTrace();
		}
		return true;
	}
	
	public void writeSuccessOutput(String result_file) {
		
		// copy output file to galaxy database
		System.out.println("PFM transfer was possible.");
		//FileCopier.copy(result_file, outfile);
		
		DecimalFormat fmt2 = new DecimalFormat();
		fmt2.setMaximumFractionDigits(2);
		fmt2.setMinimumFractionDigits(2);
		
		DecimalFormat fmt4 = new DecimalFormat();
		fmt4.setMaximumFractionDigits(4);
		fmt4.setMinimumFractionDigits(4);
		
		// write best matches and PFM to HTML tables
		try {
			BufferedWriter bw = new BufferedWriter(new FileWriter(new File(outputfile)));
			
			// write HTML header
			bw.write("<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.01 Transitional//EN\"\n");
			bw.write("     \"http://www.w3.org/TR/html4/loose.dtd\">\n");
			bw.write("<html>\n");
			bw.write("<head>\n");
			bw.write("<title>SABINE Result</title>\n");
			bw.write("<style type=\"text/css\">\n");
			bw.write("  h1 { font-size: 150%;color: #002780; }\n");
			bw.write("  table { width: 500px; background-color: #E6E8FA; border: 1px solid black; padding: 3px; vertical-align: middle;}\n");
			bw.write("  tr.secRow { background-color:#FFFFFF; margin-bottom: 50px; vertical-align: middle;}\n");
			bw.write("  th { padding-bottom: 8px; padding-top: 8px; text-align: center}\n");
			bw.write("  td { padding-bottom: 8px; padding-top: 8px; text-align: center}\n");
			bw.write("</style>\n");
			bw.write("</head>\n");
			bw.write("<body style=\"padding-left: 30px\">\n");

			
			// write best matches
			bw.write("<h1>Best Matches</h1>\n");
			bw.write("<table>\n");
			bw.write("  <tr><th> Transcription Factor </th><th> PFM similarity </th></tr>\n");
			
			double lowest_score = 0;
			for (int i=0; i<bestMatches.size(); i++) {
				if ((i%2) == 1) { 
					bw.write("  <tr>");
				}
				else {
					bw.write("  <tr class=\"secRow\">");
				}
				bw.write("<td> " + bestMatches.get(i).getName() + " </td>");
				bw.write("<td> " + fmt4.format(bestMatches.get(i).getLabel()) + " </td></tr>\n");
				
				lowest_score = bestMatches.get(i).getLabel();
			}
			bw.write("</table>\n\n");
			bw.write("<br><br>\n\n");
			
			// write PFM
			bw.write("<h1>Predicted PFM</h1>\n");
			bw.write("<table>\n");
			bw.write("  <tr><th> Pos. </th><th> A </th><th> C </th><th> G </th><th> T </th></tr>\n");
			
			for (int i=0; i<predPFM.size(); i++) {
				if ((i%2) == 1) { 
					bw.write("<tr>");
				}
				else {
					bw.write("<tr class=\"secRow\">");
				}
				bw.write("<td> " +  Math.round(predPFM.get(i)[0]) + " </td>");
				bw.write("<td> " +  fmt2.format(predPFM.get(i)[1]) + " </td>");
				bw.write("<td> " +  fmt2.format(predPFM.get(i)[2]) + " </td>");
				bw.write("<td> " +  fmt2.format(predPFM.get(i)[3]) + " </td>");
				bw.write("<td> " +  fmt2.format(predPFM.get(i)[4]) + " </td></tr>\n");
			}
			bw.write("</table>\n");
			bw.write("<br>\n\n");
			
			// get threshold values for dynamic BMT
			double[] thresholds = FBPPredictor.getThresholdValues();
			double high_conf_bmt = thresholds[0];
			double medium_conf_bmt = thresholds[1];
			
			// append confidence sign
			String path2database = "../../../static/images/sabine/";
			
			
			String conf_img_path = "";
			
			if (lowest_score > high_conf_bmt)
				conf_img_path = "sabine_static/high_conf.png";
			else if (lowest_score > medium_conf_bmt) {
				conf_img_path = "sabine_static/medium_conf.png";
			}
			else {
				conf_img_path = "sabine_static/low_conf.png";
			}
			bw.write("<img src=\"" + path2database + conf_img_path + "\" height=\"50px\">\n\n");
			bw.write("<br><br>\n\n");
			
			// append sequence logo
			String png_file = outputfile.substring(0, outputfile.length()-4) + ".png";
			javax.imageio.ImageIO.write(SequenceLogo.drawSequenceLogo(logoArray, 100), "png", new File(png_file));
			png_file = png_file.substring(png_file.lastIndexOf('/')+1, png_file.length());
			

			
			bw.write("<h1>Sequence logo</h1>\n");
			bw.write("<img src=\"" + path2database + png_file + "\" height=\"100px\">\n\n");
			
			// close HTML file
			bw.write("</body>\n");
			bw.write("</html>\n");
			
			bw.flush();
			bw.close();

		} catch (Exception e) {
			e.printStackTrace();
		}
		
		
	}
	
	public void writeFailureOutput() {
		
		try {
			BufferedWriter bw = new BufferedWriter(new FileWriter(new File(outputfile)));
			
			System.out.println("PFM transfer was not possible.");
			bw.write("PFM transfer was not possible.");

			bw.flush();
			bw.close();

		} catch (IOException e) {
			e.printStackTrace();
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
			starter.writeInputFile(starter.basedir + "infile.tmp");
		}
		else {
			
			int num_factors = starter.countTFs(starter.inputfile);
			
			if (num_factors > starter.tf_num_limit) {
				System.out.println("Number of transcription factors in input file must not exceed " + starter.tf_num_limit + ".");
				
				try {
					BufferedWriter bw = new BufferedWriter(new FileWriter(new File(starter.outputfile)));
					
					bw.write("Limit for number of transcription factors in input file exceeded\n");
					bw.write("================================================================\n");
					bw.write("Limit for number of transcription factors : " + starter.tf_num_limit + "\n");
					bw.write("Number of factors in given input file     : " + num_factors);
					
					bw.flush();
					bw.close();
					
				} catch (IOException e) {
					e.printStackTrace();
				}
				return;
			}
			
			// move input file to temporary working directory
			FileCopier.copy(starter.inputfile, starter.basedir + "infile.tmp");
		}
		
		starter.runSabine();
		
		if (starter.batchMode) {
			
			// copy output file to Galaxy database
			FileCopier.copy(starter.basedir + "outfile.tmp", starter.outputfile);
		}
		else {
			boolean pfm_transferred = starter.readOutfile(starter.basedir + "prediction.out");
		
			if (pfm_transferred) {
				starter.writeSuccessOutput(starter.basedir + "prediction.out");
			}
			else {
				starter.writeFailureOutput();
			}
			//starter.deleteBasedir();
		}
	}
}

