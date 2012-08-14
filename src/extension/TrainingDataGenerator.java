package extension;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.StringTokenizer;

import main.SABINE_Main;

/*
 *  All other functions required for generating trainingsets for TFpredict & SABINE can be found 
 *  in the TFpredict project (same class name) 
 */
public class TrainingDataGenerator {
	
	// merges PFMs in given SABINE training set file
	private static void mergePFMsInFlatfile(String flatFile, String outputFile) {
		
		// parse SABINE training set
		TransfacParser tf_parser = new TransfacParser();
		tf_parser.parseFactors(flatFile);
	
		
		// filter outlier PFMs and merge them using STAMP
		PFMFormatConverter pfm_converter = new PFMFormatConverter();
		String base_dir = SABINE_Main.createBaseDir();
		SABINE_Main.createTempDirectories(base_dir);
		pfm_converter.basedir = base_dir;
		pfm_converter.silent = true;
		
		ArrayList<ArrayList<String>> stamp_pfms = new ArrayList<ArrayList<String>>();
		ArrayList<String> merged_pfm = new ArrayList<String>();
		ArrayList<String[]> merged_transfac_pfms = new ArrayList<String[]>();
		//ArrayList<ArrayList<String[]>> merged_pfms = new ArrayList<ArrayList<String[]>>();
		
		for (int i=0; i<tf_parser.tf_names.size(); i++) {
			
			System.out.println("Processing the " + tf_parser.pfms.get(i).size() + 
							   " PFM(s) of TF " + tf_parser.tf_names.get(i) + 
							   "\t(" + (i+1) + " / " + tf_parser.tf_names.size() + ")");
			
			// more than one PFM ?
			if (tf_parser.pfms.get(i).size() > 1) {
				
				// only use first 5 PFMs for merging
				ArrayList<String[]> curr_pfms = new ArrayList<String[]>();
				for (int j=0; j<Math.min(5, tf_parser.pfms.get(i).size()); j++) {
					curr_pfms.add(tf_parser.pfms.get(i).get(j));
				}
				stamp_pfms = pfm_converter.convertAllTransfacToSTAMP(curr_pfms);
				merged_pfm = pfm_converter.mergePFMs(stamp_pfms);
				
				String[] merged_converted_pfm = pfm_converter.convertSTAMPToTransfac(merged_pfm);
				merged_transfac_pfms.add(merged_converted_pfm);	
				//tf_parser.pfms.set(i, merged_transfac_pfm);
				
			} else {
				merged_transfac_pfms.add(tf_parser.pfms.get(i).get(0));
			}
		}
		
		// replace PFMs by merged FBPs
		for (int i=0; i<tf_parser.tf_names.size(); i++) {
			if (tf_parser.pfms.get(i).size() > 1) {
				ArrayList<String[]> currPFMs = new ArrayList<String[]>();
				currPFMs.add(merged_transfac_pfms.get(i));
				tf_parser.pfms.set(i, currPFMs);
			}
		}
		
		// adjust matrix IDs
		tf_parser.pfm_names = pfm_converter.mergePFMnames(tf_parser.pfm_names);
		
		// write SABINE training set file
		tf_parser.writeFactorsToFile(outputFile);
	}
	
	
	// merges PFMs in given SABINE training set file
	private static void truncatePFMnames(String flatFile, String outputFile) {
		
		// parse SABINE training set
		TransfacParser tf_parser = new TransfacParser();
		tf_parser.parseFactors(flatFile);
		
		// truncate PFM names (only include first 5 entries)
		for (int i=0; i<tf_parser.pfm_names.size(); i++) {
		
			int[] pos = getAllIndicesOf(tf_parser.pfm_names.get(i).get(0), ",");
			if (pos.length > 4) {
				int last_pos = pos[4];
				String pfmName = tf_parser.pfm_names.get(i).get(0).substring(0, last_pos);
				tf_parser.pfm_names.get(i).set(0, pfmName);
			} 
		}
		
		// write SABINE training set file
		tf_parser.writeFactorsToFile(outputFile);
	}
	
	private static void filterBySpecies(String flatFile, String supportedSpeciesFile, String outputFile) {
		
		// parse SABINE training set
		TransfacParser tf_parser = new TransfacParser();
		tf_parser.parseFactors(flatFile);
		
		// read list of supported species
		ArrayList<String> supportedSpecies = readFile2List(supportedSpeciesFile, false);
		boolean[] incompatibleTFs = new boolean[tf_parser.species.size()];
		
		for (int i=0; i<tf_parser.species.size(); i++) {
			String currSpecies = formatSpecies(tf_parser.species.get(i));
			
			if (!supportedSpecies.contains(currSpecies)) {
				incompatibleTFs[i] = true;
				System.out.println("TF \"" + tf_parser.tf_names.get(i) + "\" has incompatible species \"" + currSpecies + "\".");
			} else {
				tf_parser.species.set(i, currSpecies);
			}
		}
		
		tf_parser.removeTFs(incompatibleTFs);
		
		// write SABINE training set file
		tf_parser.writeFactorsToFile(outputFile);
	}
	
	private static String formatSpecies(String speciesString) {
		
		StringTokenizer strtok = new StringTokenizer(speciesString.replaceAll("[\\(\\)]", ""));
		StringBuffer formattedSpeciesString = new StringBuffer();
		int tokenCnt = 0;
		while (strtok.hasMoreTokens() && tokenCnt++ < 2) {
			String currToken = strtok.nextToken();
			formattedSpeciesString.append(currToken.substring(0, 1).toUpperCase() + currToken.substring(1).toLowerCase() + " ");
		}
		return(formattedSpeciesString.toString().trim());
	}
	
	public static int[] getAllIndicesOf(String string, String substring) {
		
		ArrayList<Integer> indices = new ArrayList<Integer>();
		int index = string.indexOf(substring);
		while(index >= 0) {
			indices.add(index);
			index = string.indexOf(substring, index+1);
		}
		return Integer2int(indices.toArray(new Integer[]{}));
	}
	
	public static int[] Integer2int(Integer[] intArray) {
		
		int[] res = new int[intArray.length];
		for (int i=0; i<res.length; i++) {
			res[i] = intArray[i].intValue();
		}
		return(res);
	}
	
	public static ArrayList<String> readFile2List(String fileName, boolean upperCase) {
		
		ArrayList<String> fileContent = null;
		
		try {
			fileContent = readStream2List(new FileInputStream(new File(fileName)), upperCase);
		
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
		return fileContent;
	}
	
	public static ArrayList<String> readStream2List(InputStream stream, boolean upperCase) {
		
		ArrayList<String> fileContent = new ArrayList<String>();
		
		try {
			BufferedReader br = new BufferedReader(new InputStreamReader(stream));
			
			String line;
			while ((line = br.readLine()) != null) {
				if (upperCase) {
					fileContent.add(line.trim().toUpperCase());
				} else {
					fileContent.add(line.trim());
				}
			}
			br.close();
			
		} catch (IOException ioe) {
			ioe.printStackTrace();
		}
		
		return(fileContent);
	}
	
	public static void main(String[] args) {
		
		String flatfileDir = "/rahome/eichner/projects/tfpredict/data/tf_pred/sabine_files/latest/";
		String mergedFlatfileDBDs = flatfileDir + "transfac2012.2_matbase8.2_flatfile_with_dbd.txt";
		String mergedFlatfileDBDsFBPs = flatfileDir + "transfac2012.2_matbase8.2_flatfile_with_dbd_and_mergedPFM.txt";
		String mergedFlatfileDBDsFBPsValidSpecies = flatfileDir + "transfac2012.2_matbase8.2_flatfile_with_dbd_and_mergedPFM_validSpecies.txt";
		String supportedSpeciesList = "/rahome/eichner/workspace/SABINE/src/resources/txt/organism_list.txt";
		
		// merge PFMs using STAMP
		//mergePFMsInFlatfile(mergedFlatfileDBDs, mergedFlatfileDBDsFBPs);
		//truncatePFMnames(mergedFlatfileDBDsFBPs, mergedFlatfileDBDsFBPsNew);
		
		filterBySpecies(mergedFlatfileDBDsFBPs, supportedSpeciesList, mergedFlatfileDBDsFBPsValidSpecies);
	}
}
