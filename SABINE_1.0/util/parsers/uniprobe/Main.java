package util.parsers.uniprobe;
/*
 * ===============================================
 * (C) Florian Topf, University of Tuebingen, 2010
 * ===============================================
 */

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.zip.ZipEntry;

/*
 * Calls the subroutines
 * 
 * INPUT:
 * 
 * uniprobe_factors.txt
 * transfac_classification.txt 
 * organism_list.txt 
 * uniprot_prim2sec_id.txt 
 * uniprobe_pwms.zip
 * manual_annotation.txt
 * .fasta
 * .gff
 * 
 */

public class Main {

	public static void main(String[] args) throws IOException {
		
		/*
		// Check input
		for (int i=0; i<args.length; i++) {
			System.out.println(i+" "+args[i]);
		}
		*/
		
		///////
		// Parsing txt files
		// uniprobe_factors.txt
		File in0 = new File(args[0]);
		System.out.println("Parsing: "+args[0]);
		UniProbeFactorsParser uniprobe_parser = new UniProbeFactorsParser();
		uniprobe_parser.parseFile(in0);
		
		// transfac_classification.txt
		File in1 = new File(args[1]);
		System.out.println("Parsing: "+args[1]);
		TransfacClassificationParser transfac_parser = new TransfacClassificationParser();
		transfac_parser.parseFile(in1);
		
		// organism_list.txt
		File in2 = new File(args[2]);
		System.out.println("Parsing: "+args[2]);
		OrganismListParser organism_parser = new OrganismListParser();
		organism_parser.parseFile(in2);

		// uniprot_prim2sec_id.txt
		File in3 = new File(args[3]);
		System.out.println("Parsing: "+args[3]);
		UniProtIDsParser uniprot_parser = new UniProtIDsParser();
		// Create HashMap for UniProtIDs
		HashMap<String,String> hm_uniprot = uniprot_parser.parseFile(in3);
		//
		///////
		
		
		// Clean things up...
		Cleaner cleaner = new Cleaner();
		cleaner.cleanup();
		System.out.println("Cleaned up "+((uniprobe_parser.Pnames.size())-(cleaner.Pnames_cln.size()))+" of "+(uniprobe_parser.Pnames.size())+" entries.");
		

		//////
		// Allocate processed Objects
		ArrayList<String> Pnames = cleaner.getPnames();
		ArrayList<String> UniProtIDs = cleaner.getUniProtIDs();
		ArrayList<ArrayList<String>> Domainholder = cleaner.getDomainholder();
		ArrayList<String> Species = cleaner.getSpecies();
		
		HashMap<String,String> hm_transfac = transfac_parser.gethm();
		//
		//////

		
		// Map UniProtIDs ...
		UniProtIDsMapper uniprot_mapper = new UniProtIDsMapper();
		System.out.println("Mapping UniProtIDs...");
		UniProtIDs = uniprot_mapper.mapIDs(UniProtIDs, hm_uniprot);
		System.out.println("UniProtIDs mapped.");
		
		
		// uniprobe_pwms.zip
		File in4 = new File(args[4]);
		System.out.println("Reading: "+args[4]);
		ZipHelper zip_helper = new ZipHelper();
		// Get entries of zip-file
		ArrayList<ZipEntry> ZipEntries = zip_helper.getZipEntries(in4);
		
		
		// Parse zip-entries
		PWMParser pwm_parser = new PWMParser();
		pwm_parser.parse(ZipEntries, in4);
		// Gets the result
		ArrayList<String> PWMnames = pwm_parser.getPWMnames();
		ArrayList<ArrayList<String>> Pwms = pwm_parser.getPwms();
		// MN (filenames of pwms)
		ArrayList<String> MNs = pwm_parser.getMNs();
				
		
		// Map PWMs
		System.out.println("Mapping PWMs...");
		PWMMapper pwm_mapper = new PWMMapper();
		pwm_mapper.map(Pnames, UniProtIDs, Domainholder, Species, PWMnames, Pwms, MNs);
		
		
		// Refresh Objects
		Pnames = pwm_mapper.getPnames_mapped();
		UniProtIDs = pwm_mapper.getUniProtIDs_mapped();
		Domainholder = pwm_mapper.getDomainholder_mapped();
		Species = pwm_mapper.getSpecies_mapped();
		Pwms = pwm_mapper.getPwms_mapped();
		//MN
		MNs = pwm_mapper.getMNs_mapped();
		
		
		// Add manual annotation entries from file to HM
		File in5 = new File(args[5]);
		System.out.println("Parsing: "+args[5]);
		ManualAnnotation manual_annotation = new ManualAnnotation();
		hm_transfac = manual_annotation.annote(in5, hm_transfac);
		// Map Transfac-Classification
		TransfacClassificationMapper transfac_mapper = new TransfacClassificationMapper();
		transfac_mapper.map(hm_transfac, Pnames, UniProtIDs, Domainholder, Species, Pwms, MNs);

		
		// Refresh Objects
		Pnames = transfac_mapper.getPnames_mapped();
		UniProtIDs = transfac_mapper.getUniProtIDs_mapped();
		Domainholder = transfac_mapper.getDomainholder_mapped();
		Species = transfac_mapper.getSpecies_mapped();
		Pwms = transfac_mapper.getPwms_mapped();
		ArrayList<String> Transfac = transfac_mapper.getTransfac_mapped();
		// MN
		MNs = transfac_mapper.getMNs_mapped();
		
		// Parse UniProt fasta/gff data
		UniProtParser raw_uniprot_parser = new UniProtParser();
		System.out.println("Parsing: "+args[6]);
		System.out.println("Parsing: "+args[7]);
		raw_uniprot_parser.parse(args[6], args[7]);
				
		ArrayList<String> raw_uniprot_IDs = raw_uniprot_parser.getUniprot_IDs();
		ArrayList<String> raw_sequences = raw_uniprot_parser.getSequences();
		ArrayList<ArrayList<String>> raw_DNAbinding = raw_uniprot_parser.getDomains();
		
		
		// Merging of the sequence-data and DNA-bindings with the rest
		Merger data_merger = new Merger();
		System.out.println("Merging data...");
		data_merger.merge(Pnames, UniProtIDs, Domainholder, Species, Pwms, Transfac, raw_uniprot_IDs, raw_sequences, raw_DNAbinding, MNs);
				
		// Refresh Objects
		Pnames = data_merger.getPnames_merged();
		UniProtIDs = data_merger.getUniProtIDs_merged();
		Domainholder = data_merger.getDomainholder_merged();
		Species = data_merger.getSpecies_merged();
		Pwms = data_merger.getPwms_merged();
		Transfac = data_merger.getTransfac_merged();
		ArrayList<String> Sequences = data_merger.getSequences_merged();
		ArrayList<ArrayList<String>> DNAbinding = data_merger.getDNAbinding_merged();
		// MN
		MNs = data_merger.getMNs_merged();
		
		// Writing Output File
		OutputFileWriter output_writer = new OutputFileWriter();
		System.out.println("Writing output file: "+args[8]);
		output_writer.writeOutfile(Pnames, UniProtIDs, Domainholder, Species, Pwms, Transfac, Sequences, DNAbinding, MNs, args[8]);
		
		
		// Finished! 
		System.out.println("Done.");
	
	}

}
