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

import help.FeatureFileGenerator;
import help.FormatConverter;
import help.LibSVMFeatureScaler;
import help.ModelFileObtainer;
import help.TFNamePairsFileGenerator;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;

import libsvm.LabeledTF;
import libsvm.SVMPredictor;

import org.biojava.bio.BioException;

import core.DomainFeatureCalculator;
import core.FeatureRepairer;
import core.IrrelevantPairIdentifier;
import core.SVMPairwiseFeatureCalculator;
import core.SequenceAligner;
import core.SequenceFeatureCalculator;
import core.SpeciesFeatureCalculator;
import de.zbit.util.progressbar.AbstractProgressBar;


public class FBPPredictor {
	
	public static final String public_trainingset = "data/trainingsets_public/";
	public static final String private_trainingset = "data/trainingsets_private/";
	public static final String classMappingFile = public_trainingset + "Classes";
	
	public static final String matrix_dir = "data/substitutionMatrices/";
	public static final String model_dir = "data/models/";
	public static final String MMkernelDir = "lib/MismatchKernel/";
	public static final String LAkernelDir = "lib/LAKernel/";
	public static final String PsiPredDir = "lib/PSIPRED/";
	public static final String MoStaDir = "lib/MoSta/";
	public static final String MatlignDir = "lib/Matlign/";
	public static final String StampDir = "lib/STAMP/";
	
	double best_match_threshold = 0.95;
	
	int max_number_of_best_matches = 5;
	
	double outlier_filter_threshold = 0.5;
	
	
	boolean silent = false;
	
	boolean minimal_output_mode = false;
	
	boolean gui_output_mode = false;
	
	boolean dynamic_threshold = true;
	
	public int num_candidates;
	
	public static final double high_conf_bmt = 0.95;
	public static final double medium_conf_bmt = 0.8;
	public static final double low_conf_bmt = 0.5;
	
	/*
	 * parses an input file and converts it to the internally used file format
	 */
	FormatConverter converter = new FormatConverter();
	
	/**
	 * calculates alignment- and kernel-scores of two binding-domains
	 */
	DomainFeatureCalculator domaincalculator = new DomainFeatureCalculator();
	
	
	/**
	 * calculates alignment scores of domain-environments and secondary structures
	 */
	SequenceFeatureCalculator sequencecalculator = new SequenceFeatureCalculator();
	
	
	/**
	 * looks up phylogentic distances between species
	 */
	SpeciesFeatureCalculator speciescalculator = new SpeciesFeatureCalculator();
	
	/**
	 * calculates scores according to the svm-pairwise approach 
	 */
	SVMPairwiseFeatureCalculator svmpairwisecalculator = new SVMPairwiseFeatureCalculator();
	
	
	/**
	 * eliminates numerical instabilities within certain features (MEHP950101, MEHP950103, MIYT790101)
	 */
	FeatureRepairer repairer = new FeatureRepairer();
	
	
	/**
	 * produces an unlabeled feature file in libsvm-format  
	 */
	FeatureFileGenerator libsvmfilegenerator = new FeatureFileGenerator();
	
	/**
	 * scales features in unlabeled feature file with respect to a given training set
	 */
	LibSVMFeatureScaler featurescaler = new LibSVMFeatureScaler();
	
	
	/**
	 * predicts FBP-similarity scores of TF-pairs, extracts best matches and performs the FBP-transfer
	 */
	SVMPredictor predictor = new SVMPredictor();
	
	
	/**
	 * looks up the optimal libsvm-modelfile for the predictions
	 */
	ModelFileObtainer obtainer = new ModelFileObtainer();
	
	
	/*
	 * generates a list of pairs of tf-names for the FBP-transfer
	 */
	TFNamePairsFileGenerator namepairsfilegenerator = new TFNamePairsFileGenerator();
	
	 /**
   * 
   */
  AbstractProgressBar progress = null;
  
  public void setProgressBar(AbstractProgressBar progress) {
    this.progress = progress;
  }
	
	public static double[] getThresholdValues() {
		return new double[] {high_conf_bmt , medium_conf_bmt , low_conf_bmt};
	}
	
	
	private static int getTrainingSetSize(String class_id, String train_dir) {
		
		String line = null;
		
		int entry_counter = 0;
		 
		BufferedReader br = null;
			 
		try {
			br = new BufferedReader(new FileReader(new File(train_dir + "trainingset_" + class_id + ".rawdata")));
		} catch (FileNotFoundException e) {
			System.out.println(e.getMessage());
			System.out.println("Error occured before counting entries in training set");
		}
					
		try {
			while((line = br.readLine()) != null) {
							
				if(line.startsWith("NA  ")) {
						entry_counter++;
				} 
			}
		} catch (IOException e) {
			System.out.println(e.getMessage());
			System.out.println("Error occured while counting entries in training set");
		}
		return entry_counter;
	}	
	
	
	public void calculateFBP(String name, String class_id, String species, String sequence1, String sequence2, ArrayList<String> domains, String base_dir, String train_dir, String model_file) {

		/*
		 *  predict superclass if necessary
		 */
		
		if (class_id == null) {
			class_id = "class0";
			if (!silent) System.out.println("\nSuperclass not found. Predicting superclass based on sequence homology.");
			if (gui_output_mode) System.out.print("\nPredicting superclass based on sequence homology.");
			
			String[] all_class_ids = new String[]{"class0", "class1", "class2", "class3", "class4"};
			String[] all_class_names = new String[]{"Other", "Basic Domain", "Zinc Finger", "Helix-Turn-Helix", "Beta Scaffold"};
			
			try {
				// read available protein sequences for all superclasses
				
				double max_score = Double.NEGATIVE_INFINITY;
				for (int c=0; c<all_class_ids.length; c++) {
					
					SequenceFeatureCalculator seqCalc = new SequenceFeatureCalculator();
					seqCalc.silent = true;
					boolean [] irrelPairs = new boolean[getTrainingSetSize(all_class_ids[c], train_dir)];
					seqCalc.parseRelevantDomainsAndSequences(irrelPairs, all_class_ids[c], train_dir);
					SequenceAligner class_predictor = new SequenceAligner(matrix_dir + "BLOSUM_62.dat", "NW");
					
					HashSet<String> all_seqs = new HashSet<String>();
					all_seqs.addAll(seqCalc.get_other_sequences1());
					all_seqs.addAll(seqCalc.get_other_sequences2());
					all_seqs.remove("NO SEQUENCE.");
					String[] class_seqs = all_seqs.toArray(new String[]{});
	
					for(int i=0; i<class_seqs.length; i++) {	
						
						double blosum_score = class_predictor.getSMBasedSimilarity(sequence1, class_seqs[i]);
						
						if (blosum_score > max_score) {
							max_score = blosum_score;
							class_id = "class" + c;
						}
					}
				}
			} catch(IOException ioe) {
				System.out.println(ioe.getMessage());
				System.out.println("IOException occurred while predicting superclass.");
				
			} catch(BioException bioe) {
				System.out.println(bioe.getMessage());
				System.out.println("BioException occurred while predicting superclass.");
			}
			
			int class_num = Integer.parseInt(class_id.substring(5));
			System.out.println("\n  Predicted class: " + class_num + " (" + all_class_names[class_num] + ")");
		}

		/*
		 * 
		 * calculate 2 features for all pairs
		 * 
		 */
		
		int num_entries = getTrainingSetSize(class_id, train_dir);
		boolean[] irrelevantPairs = new boolean[num_entries];

		domaincalculator.silent = true;
		domaincalculator.basedir = base_dir;
		sequencecalculator.silent = true;
		sequencecalculator.basedir = base_dir;
		domaincalculator.setProgressBar(null);
		
		
		if (! silent) System.out.println("\n  Calculating features for all tf pairs.");
		
		domaincalculator.parseRelevantDomains(irrelevantPairs, class_id, train_dir);

		
		/*
		 *  predict DNA binding domain if no domain assignment is available
		 */
		
		double max_blosum_score = Double.NEGATIVE_INFINITY;
		String[] res, res1, res2;
		double score, score1, score2;
		int start_pos, end_pos, seq_idx;
		start_pos = end_pos = seq_idx = 1;
		
		try {
			if (domains.isEmpty()) {
				if (!silent) System.out.println("\nDNA-binding domain not found. Predicting domain based on sequence homology.");
				SequenceAligner dom_predictor;
				domaincalculator.predicted_domains = true;
				sequencecalculator.predicted_domains = true;
				
				dom_predictor = new SequenceAligner(matrix_dir + "BLOSUM_62.dat", "SW");

				for(int i=0; i<domaincalculator.get_other_domains().size(); i++) {	
					for(int j=0; j<domaincalculator.get_other_domains().get(i).size(); j++) {
						
						res1 = new String[3];
						res2 = new String[3];
						score1 = score2 = Double.NEGATIVE_INFINITY;
						
						String curr_domain_seq = domaincalculator.get_other_domains().get(i).get(j); 
						
						if (sequence1 != null) {
							res1 = dom_predictor.getMatchingRegionAndScore(sequence1, curr_domain_seq);
							score1 = Double.parseDouble(res1[0]);
							seq_idx = 1;
						}
						if (sequence2 != null) {
							res2 = dom_predictor.getMatchingRegionAndScore(sequence2, curr_domain_seq);
							score2 = Double.parseDouble(res2[0]);
							seq_idx = 2;
						}
						
						if (score1 >= score2) res = res1;
						else 				  res = res2; 
						
						if ((score = Double.parseDouble(res[0])) > max_blosum_score) {
							max_blosum_score = score;
							start_pos = Integer.parseInt(res[1]);
							end_pos = Integer.parseInt(res[2]);
						}
					}
					if (seq_idx == 1) {
						domains.add(sequence1.substring(start_pos-1,end_pos) + " " + seq_idx + "  " + start_pos + " " + end_pos);
					}
					else {
						domains.add(sequence2.substring(start_pos-1,end_pos) + " " + seq_idx + "  " + start_pos + " " + end_pos);
					}
				}
			}
		}
		catch(IOException ioe) {
			System.out.println(ioe.getMessage());
			System.out.println("IOException occurred while predicting domains.");
		}
		catch(BioException bioe) {
			System.out.println(bioe.getMessage());
			System.out.println("BioException occurred while predicting domains.");
		}
		
		domaincalculator.calculateDomainFeatureFile(name, domains, "SMBasedSimilarity", null, matrix_dir + "BLOSUM_62.dat", base_dir + "allpairs/domain_scores_BLOSUM_62.out");
		domaincalculator.calculateDomainFeatureFile(name, domains, "SMBasedSimilarity", null, matrix_dir + "PAM_080.dat"  , base_dir + "allpairs/domain_scores_PAM_080.out");
		
		IrrelevantPairIdentifier identifier = new IrrelevantPairIdentifier();
		
		irrelevantPairs = identifier.identifyIrrelevantPairs(base_dir + "allpairs/domain_scores_BLOSUM_62.out", 0.3, irrelevantPairs.length);
		
		// count number of candidate TFs for PFM transfer
		int numRelevantPairs = 0;
		for (boolean isIrrelevant: irrelevantPairs) {
			if (!isIrrelevant) {
				numRelevantPairs++;
			}
		}
		num_candidates = numRelevantPairs;
		
		/*
		 * 
		 * calculate 30 feature values for all relevant tf pairs
		 * 
		 */
		if (! silent) System.out.println("\nCalculating features.");
		if (gui_output_mode) System.out.println("Calculating features:");
		if (! silent) System.out.println("\n  Calculating features for relevant tf pairs.");
		
		domaincalculator.parseRelevantDomains(irrelevantPairs, class_id, train_dir);
		
		sequencecalculator.parseRelevantDomainsAndSequences(irrelevantPairs, class_id, train_dir);
		sequencecalculator.parseRelevantSecondaryStructures(irrelevantPairs, class_id, train_dir);
		
		// Configure the progress bar
    if (progress!=null) {
      progress.setNumberOfTotalCalls(30);
      progress.setEstimateTime(false);
      domaincalculator.setProgressBar(progress);
    }
		
		if (! silent) System.out.println("\n    Calculating substitution matrix based alignment scores.");
		if (gui_output_mode) System.out.print("  Calculating substitution matrix based alignment scores...");
		
		domaincalculator.calculateDomainFeatureFile(name, domains, "SMBasedSimilarity", null, matrix_dir + "BLOSUM_62.dat"   , base_dir + "relevantpairs/domain_scores_BLOSUM_62.out");
		domaincalculator.calculateDomainFeatureFile(name, domains, "SMBasedSimilarity", null, matrix_dir + "PAM_080.dat"     , base_dir + "relevantpairs/domain_scores_PAM_080.out");
		
		domaincalculator.calculateDomainFeatureFile(name, domains, "SMBasedSimilarity", null, matrix_dir + "PAM_010.dat"     , base_dir + "relevantpairs/domain_scores_PAM_010.out");
		domaincalculator.calculateDomainFeatureFile(name, domains, "SMBasedSimilarity", null, matrix_dir + "WEIL970101.dat"  , base_dir + "relevantpairs/domain_scores_WEIL970101.out");
		domaincalculator.calculateDomainFeatureFile(name, domains, "SMBasedSimilarity", null, matrix_dir + "MEHP950101.dat"  , base_dir + "relevantpairs/domain_scores_MEHP950101.out");
		domaincalculator.calculateDomainFeatureFile(name, domains, "SMBasedSimilarity", null, matrix_dir + "MEHP950102.dat"  , base_dir + "relevantpairs/domain_scores_MEHP950102.out");
		domaincalculator.calculateDomainFeatureFile(name, domains, "SMBasedSimilarity", null, matrix_dir + "MEHP950103.dat"  , base_dir + "relevantpairs/domain_scores_MEHP950103.out");
		domaincalculator.calculateDomainFeatureFile(name, domains, "SMBasedSimilarity", null, matrix_dir + "LUTR910102.dat"  , base_dir + "relevantpairs/domain_scores_LUTR910102.out");
		domaincalculator.calculateDomainFeatureFile(name, domains, "SMBasedSimilarity", null, matrix_dir + "NIEK910102.dat"  , base_dir + "relevantpairs/domain_scores_NIEK910102.out");
		domaincalculator.calculateDomainFeatureFile(name, domains, "SMBasedSimilarity", null, matrix_dir + "RISJ880101.dat"  , base_dir + "relevantpairs/domain_scores_RISJ880101.out");
		domaincalculator.calculateDomainFeatureFile(name, domains, "SMBasedSimilarity", null, matrix_dir + "MIYS930101.dat"  , base_dir + "relevantpairs/domain_scores_MIYS930101.out");
		domaincalculator.calculateDomainFeatureFile(name, domains, "SMBasedSimilarity", null, matrix_dir + "MIYT790101.dat"  , base_dir + "relevantpairs/domain_scores_MIYT790101.out");
		
		if (! silent) System.out.println("    Calculating sequence identity based alignment scores.");
		if (gui_output_mode) System.out.print("done.\n  Calculating sequence identity based alignment scores...");
		
		domaincalculator.calculateDomainFeatureFile(name, domains, "SequenceIdentity" , null, matrix_dir + "BLOSUM_62.dat"   , base_dir + "relevantpairs/domain_scores_BLOSUM_62_si.out");
		
		if (! silent) System.out.println("    Calculating sequence similarity based alignment scores.");
		if (gui_output_mode) System.out.print("done.\n  Calculating sequence similarity based alignment scores...");
		
		domaincalculator.calculateDomainFeatureFile(name, domains, "SMBasedIdentity" , new String[] {"1.0"} , matrix_dir + "BLOSUM_62.dat"   , base_dir + "relevantpairs/domain_scores_BLOSUM_62_t=1.0.out");
		domaincalculator.calculateDomainFeatureFile(name, domains, "SMBasedIdentity" , new String[] {"3.0"} , matrix_dir + "BLOSUM_62.dat"   , base_dir + "relevantpairs/domain_scores_BLOSUM_62_t=3.0.out");
		domaincalculator.calculateDomainFeatureFile(name, domains, "SMBasedIdentity" , new String[] {"5.0"} , matrix_dir + "BLOSUM_62.dat"   , base_dir + "relevantpairs/domain_scores_BLOSUM_62_t=5.0.out");
		
		if (! silent) System.out.println("    Calculating local alignment kernel scores.");
		if (gui_output_mode) System.out.print("done.\n  Calculating local alignment kernel scores...");
		
		domaincalculator.calculateDomainFeatureFile(name, domains, "LocalAlignmentKernel", new String[] {"GCBopt.dat"} 		, null, base_dir + "relevantpairs/domain_scores_lak_GCBopt.out");
		domaincalculator.calculateDomainFeatureFile(name, domains, "LocalAlignmentKernel", new String[] {"JTTopt.dat"} 		, null, base_dir + "relevantpairs/domain_scores_lak_JTTopt.out");
		domaincalculator.calculateDomainFeatureFile(name, domains, "LocalAlignmentKernel", new String[] {"BLOSUM_62opt.dat"} 	, null, base_dir + "relevantpairs/domain_scores_lak_BLOSUM_62opt.out");
		domaincalculator.calculateDomainFeatureFile(name, domains, "LocalAlignmentKernel", new String[] {"PAM_250opt.dat"} 	, null, base_dir + "relevantpairs/domain_scores_lak_PAM_250opt.out");
		domaincalculator.calculateDomainFeatureFile(name, domains, "LocalAlignmentKernel", new String[] {"blosum62.dat"} 		, null, base_dir + "relevantpairs/domain_scores_lak_blosum62.out");
		
		if (! silent) System.out.println("    Calculating mismatch kernel scores.");
		if (gui_output_mode) System.out.print("done.\n  Calculating mismatch kernel scores...");
		
		domaincalculator.calculateDomainFeatureFile(name, domains, "MismatchKernel", new String[] {"4", "1"}, null, base_dir + "relevantpairs/domain_scores_mmk_4_1.out");
		domaincalculator.calculateDomainFeatureFile(name, domains, "MismatchKernel", new String[] {"5", "1"}, null, base_dir + "relevantpairs/domain_scores_mmk_5_1.out");
		domaincalculator.calculateDomainFeatureFile(name, domains, "MismatchKernel", new String[] {"6", "1"}, null, base_dir + "relevantpairs/domain_scores_mmk_6_1.out");
		
		if (! silent) System.out.println("    Calculating secondary structure scores.");
		if (gui_output_mode) System.out.print("done.\n  Calculating secondary structure scores...");
		
		sequencecalculator.calculateSequenceFeatureFile(name, domains, sequence1, sequence2, "SecondaryStructure", null, matrix_dir + "BLOSUM_62.dat", base_dir + "relevantpairs/domain_scores_secstr_blo62.out");
		if (progress!=null) progress.DisplayBar(); // 25.
		
		if (! silent) System.out.println("    Calculating DNA-binding domain environment scores.");
		if (gui_output_mode) System.out.print("done.\n  Calculating DNA-binding domain environment scores...");
		
		sequencecalculator.calculateSequenceFeatureFile(name, domains, sequence1, sequence2, "Environments", new String[] {"25"}, matrix_dir + "BLOSUM_62.dat", base_dir + "relevantpairs/domain_scores_env_25_BLOSUM_62.out");
		if (progress!=null) progress.DisplayBar();  // 26.
		sequencecalculator.calculateSequenceFeatureFile(name, domains, sequence1, sequence2, "Environments", new String[] {"50"}, matrix_dir + "BLOSUM_62.dat", base_dir + "relevantpairs/domain_scores_env_50_BLOSUM_62.out");
		if (progress!=null) progress.DisplayBar();  // 27.
		
		if (! silent) System.out.println("    Calculating phylogenetic distance based scores.");
		if (gui_output_mode) System.out.print("done.\n  Calculating phylogenetic distance based scores...");
		
		speciescalculator.calculatePhylogeneticDistances(name, species, class_id, irrelevantPairs, train_dir + "new_phylogenetic_distances.out", base_dir + "relevantpairs/domain_scores_phyl_dist.out", train_dir);
		if (progress!=null) progress.DisplayBar();  // 28.
		
		if (! silent) System.out.println("    Calculating SVM pairwise scores.");
		if (gui_output_mode) System.out.print("done.\n  Calculating SVM pairwise scores...");
		
		svmpairwisecalculator.calculateSVMPairwiseScores(name, base_dir + "relevantpairs/domain_scores_BLOSUM_62.out", base_dir + "allpairs/domain_scores_BLOSUM_62.out", train_dir + "trainingset_" + class_id + ".blo62", base_dir + "relevantpairs/domain_scores_svm_pairwise_BLOSUM_62.out");
		if (progress!=null) progress.DisplayBar();  // 29.
		svmpairwisecalculator.calculateSVMPairwiseScores(name, base_dir + "relevantpairs/domain_scores_PAM_080.out", base_dir + "allpairs/domain_scores_PAM_080.out", train_dir + "trainingset_" + class_id + ".pam80", base_dir + "relevantpairs/domain_scores_svm_pairwise_PAM_080.out");
		if (progress!=null) {
		  progress.DisplayBar();  // 30.
		  progress.finished();
		}
		
		if (! silent) System.out.println("\nScreening features for numerical instabilities.");
		
	// repair feature files	
		repairer.silent = silent;
		repairer.repairAllFeatureFiles(base_dir + "relevantpairs");
		
		
		if (! silent) System.out.println("\nGenerating input file in LIBSVM format.\n");
		if (gui_output_mode) System.out.print("done.\n\nGenerating input file in LIBSVM format...");
		
	// generate libsvm-file with all features	
		libsvmfilegenerator.silent = silent;
		libsvmfilegenerator.generateOverallFeatureFile(base_dir + "relevantpairs", train_dir + class_id + ".30featurenames", base_dir + "libsvmfiles/unlabeled_testset.out");
		
	// scale features in this file	
		
		featurescaler.scaleFeatureFile(train_dir + "feature_scaling_orientation_" + class_id + ".out", base_dir + "libsvmfiles/unlabeled_testset.out", base_dir + "libsvmfiles/scaled_unlabeled_testset.out");
	
		
	// predict PFM-similarities for all relevant TF-pairs	
		predictor.silent = silent;
		predictor.gui_output_mode = gui_output_mode;
		predictor.basedir = base_dir;
		predictor.num_candidates = num_candidates;
		
		predictor.dynamic_threshold = dynamic_threshold;
		predictor.high_conf_bmt = high_conf_bmt;
		predictor.medium_conf_bmt = medium_conf_bmt;
		predictor.low_conf_bmt = low_conf_bmt;

		
		if (model_file == null) {
			model_file = obtainer.obtainModelFile(class_id);
		}
		
		predictor.predictLabels(model_file, base_dir + "libsvmfiles/scaled_unlabeled_testset.out", base_dir + "libsvmfiles/predicted_labels.out");
		
		
	// generate list of all relevant pairs of tfnames 	
		
		namepairsfilegenerator.generateTFNamePairsFile(base_dir + "relevantpairs/domain_scores_BLOSUM_62.out", base_dir + "libsvmfiles/tfnamepairs.out");
		
		
	// identify best matching TFs
		ArrayList<LabeledTF> best_matches = predictor.extractBestHits(base_dir + "libsvmfiles/tfnamepairs.out", base_dir + "libsvmfiles/predicted_labels.out", best_match_threshold, max_number_of_best_matches);
		
		
	// construct and transfer FBP to input tf
		
		predictor.performPFMTransfer(best_matches, train_dir + "FBPs_" + class_id + ".out", obtainer.obtainOptimizer(class_id), outlier_filter_threshold, base_dir + "prediction.out");
		
		
	}
	
	
	
	
	
	public void predictFBP(String inputfile, String base_dir, String train_dir, String model_file) {
		
		/*
		 * 
		 * parse input data
		 * 
		 */
		
		//		String[] data = converter.convertToInternalFormat(inputfile, "internal/input.tf");
		
		
		if (! silent) System.out.println("Parsing input file.\n");
		
		
		String[] data = converter.parseInternalFormat(inputfile, train_dir + "new_phylogenetic_distances.out");
		
		//System.exit(0);
		
		String name = data[0];
		
		String class_id = data[1];		// class0
		
		String species = data[2].trim();
		
		
		if (class_id != null && class_id.equals(FormatConverter.NonTFclassID)) {
			if (! silent) {
				System.out.println("  The protein \"" + name + "\" will be omitted.");
			}
			try {
				BufferedWriter bw = new BufferedWriter(new FileWriter(new File(base_dir + "prediction.out")));
				bw.write("BM  none (Unknown protein class or Non-TF)\nXX\nMA  none\nXX");
				bw.flush();
				bw.close();
			} catch (IOException e) {
				e.printStackTrace();
			}
			return;
		}
		
		if (! silent) {
			System.out.println("  Name    : "  + name);
			if (class_id != null) {
				System.out.println("  Class   : "  + class_id.substring(5));
			}
			System.out.println("  Species : "  + species);
		}
		
		if (minimal_output_mode) {
			System.out.println("    Factor: " + name + 
							   "    Species: " + species); 
		}
		
		String sequence1 = data[3];
		String sequence2 = null;
		

		
		if(data[4] != null && data[4].length() > 0) sequence2 = data[4];
		
		
		// remove gaps in input sequences
		sequence1 = sequence1.replaceAll("\\s", "");
		if (sequence2 != null) sequence2 = sequence2.replaceAll("\\s", "");
		
		ArrayList<String> domains = new ArrayList<String>();
		
		for(int i=5; i<data.length; i++) {
			
			domains.add(data[i]);
			
		}
		
		
		/*
		 * 
		 * calculate features
		 * 
		 */
		
		calculateFBP(name, class_id, species, sequence1, sequence2, domains, base_dir, train_dir, model_file);
		
		
	}
	
	
	public static void main(String[] args) {
		
		FBPPredictor predictor = new FBPPredictor();

		if(args.length == 0) {
			
			SABINE_Main.usage();
			
		}

		/*
		 * 
		 * parse options
		 *
		 */
		String verbose_option = "y";
		String train_dir = "trainingsets_public/";  // directory that contains training factors (= candidates for PWM transfer)
		String base_dir = null;
		String model_file = null;
		
		for(int i=1; i<args.length-1; i+=2) {
			
			if(args[i].equals("-s")) { predictor.best_match_threshold 			= Double.parseDouble (args[i+1]);   continue; }
			if(args[i].equals("-m")) { predictor.max_number_of_best_matches 	= Integer.parseInt   (args[i+1]); 	continue; }
			if(args[i].equals("-o")) { predictor.outlier_filter_threshold	 	= Double.parseDouble (args[i+1]); 	continue; }
			if(args[i].equals("-v")) { verbose_option	 						= args[i+1]; 						continue; }
			if(args[i].equals("-b")) { base_dir		 							= args[i+1]; 						continue; }
			if(args[i].equals("-t")) { train_dir		 						= args[i+1]; 						continue; }
			if(args[i].equals("-c")) { model_file		 						= args[i+1]; 						continue; }
			if(args[i].equals("-d")) { predictor.dynamic_threshold				= Boolean.parseBoolean(args[i+1]); 	continue; }
			
			if( !args[i].equals("-s") && !args[i].equals("-m") && !args[i].equals("-o") && !args[i].equals("-v") && !args[i].equals("-c") && !args[i].equals("-t") && !args[i].equals("-b") && !args[i].equals("-d")) {	
				System.out.println("  Illegal argument: " + args[i] + "\n");
				SABINE_Main.usage();
			}
		}
		
		if (verbose_option.equals("n") || verbose_option.equals("no") || verbose_option.equals("h")) 
			predictor.silent = true; 
		if (verbose_option.equals("h"))
			predictor.minimal_output_mode = true;
		
		if (! predictor.silent) {
			
			System.out.println("");
			System.out.println("  Min. best match similarity  :  " + predictor.best_match_threshold);
			System.out.println("  Max. number best matches    :  " + predictor.max_number_of_best_matches);
			System.out.println("  Max. outlier-FBP deviation  :  " + predictor.outlier_filter_threshold + "\n\n");
		}
		
		predictor.predictFBP(args[0], base_dir, train_dir, model_file);
		
		
	}
}

