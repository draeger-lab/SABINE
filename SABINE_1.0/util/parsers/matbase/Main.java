package util.parsers.matbase;
/*
 * ===============================================
 * (C) Florian Topf, University of Tuebingen, 2010
 * ===============================================
 */

import java.sql.Connection;
import java.util.ArrayList;
import java.util.HashMap;

/*
 * Main method
 * 
 * -d (dump) needs the following input parameters:
 * url (jdbc:mysql://jdbc:mysql://radb.cs.uni-tuebingen.de/MATBASE_82)
 * user
 * password
 *  
 */

import org.apache.commons.cli.*;

public class Main {
	
	@SuppressWarnings({ "unchecked" })
	public static void main(String[] args) {
		
		///
		// command line stuff...
		/*
		 * -d (dump) requires:
		 * url (jdbc:mysql://radb.cs.uni-tuebingen.de/MATBASE_82)
		 * username
		 * password
		 */
		// create options object
		Options options = new Options();
		
		// add d option
		options.addOption("d", false, "dump data");
		
		CommandLine cmd = null;
		CommandLineParser parser = new PosixParser();
		try {
			cmd = parser.parse(options, args);
		} catch (ParseException e) {
			e.printStackTrace();
		}
		////
		
		
		////
		// initialize all
		ArrayList<String> gene_ids = null;
		ArrayList<String> symbols = null;
		ArrayList<Integer> taxon_ids = null;
		ArrayList<ArrayList<String>> domains = null;
		HashMap<Integer, String> specmap = null;
		ArrayList<ArrayList<ArrayList<String>>> matrices = null;
		////
		
		Filter filter = new Filter();
		
		////
		// dump data
		if(cmd.hasOption("d")) {
			System.out.println("Dumping data...");
			
			
			// fetch command-line: url, username, password
			String sDbUrl = args[1];
			String sUsr = args[2];
			String sPwd = args[3];
			
			// establish database connection
			MySQLAPI mysqlapi = new MySQLAPI();
			Connection cn = mysqlapi.connect(sDbUrl, sUsr, sPwd);
			
			
			// dump data
			MySQLHandler mysqlhandler = new MySQLHandler();
			mysqlhandler.dumpData(cn);
			
			// get objects
			gene_ids = mysqlhandler.getGene_id();
			symbols = mysqlhandler.getSymbol();
			taxon_ids = mysqlhandler.getTaxon_id();
			// get map of taxon_id's
			specmap = mysqlhandler.getSpecmap();

			
			// filter redundant entries
			filter.checkDub(gene_ids, symbols, taxon_ids);
			
			// refresh objects
			gene_ids = filter.getGene_ids_fi();
			symbols = filter.getSymbols_fi();
			taxon_ids = filter.getTaxon_ids_fi();
			
			
			// fetch domains
			domains = mysqlhandler.dumpDomains(cn, gene_ids);
			
			// get matrices
			matrices = mysqlhandler.dumpMatrices(cn, gene_ids);
			
			// refresh objects
			gene_ids = filter.getGene_ids_fi();
			symbols = filter.getSymbols_fi();
			taxon_ids = filter.getTaxon_ids_fi();
			domains = filter.getDomains_fi();
			matrices = filter.getMatrices_fi();

			
			////
			// write objects
			ObjectRW objectrw = new ObjectRW();
					
			objectrw.write(gene_ids,"gene_ids");
			objectrw.write(symbols,"symbols");
			objectrw.write(taxon_ids,"taxon_ids");
			objectrw.write(domains,"domains");
			
			objectrw.write(specmap,"specmap");
			
			objectrw.write(matrices, "matrices");
			////
			
		}
		
		// or read from file
		else {
			System.out.println("Loading data from files...");
			
			// read objects
			ObjectRW objectrw = new ObjectRW();
			
			gene_ids = (ArrayList<String>) objectrw.read("gene_ids");
			symbols = (ArrayList<String>) objectrw.read("symbols");
			taxon_ids = (ArrayList<Integer>) objectrw.read("taxon_ids");
			domains = (ArrayList<ArrayList<String>>) objectrw.read("domains");
			
			specmap = (HashMap<Integer, String>) objectrw.read("specmap");
			
			matrices = (ArrayList<ArrayList<ArrayList<String>>>) objectrw.read("matrices");
			
		}
		////
		
		printstatus(symbols);
		
		// filter entries without corresponding matrix
		filter.checkMat(gene_ids, symbols, taxon_ids, domains, matrices);
		
		printstatus(symbols);
		
		
		// map taxon_ids -> species
		Mapper mapper = new Mapper();
		ArrayList<String> species = mapper.mapSpec(taxon_ids, specmap);
		
		// parse mapping of species names
		FileParser fileparser = new FileParser();
		StringMap matspecmap = fileparser.parseSpecFile("organism_map.txt");
		// map names of matbase-species -> organism_list.txt species
		species = mapper.mapMatSpec(species, matspecmap);
		
		// parse organism_list.txt
		ArrayList<String> organisms = fileparser.parseOrgFile("organism_list.txt");
		// filter unsupported species
		filter.checkOrg(gene_ids, symbols, species, domains, matrices, organisms);
		
		printstatus(symbols);


		// idmapping geneid -> uniprotid
		ArrayList<ArrayList<String>> uniprot_ids = mapper.mapGeneIDs(gene_ids,species);

		//filter entries without uniprot_id (null)
		filter.checkIDs(uniprot_ids, symbols, species, domains, matrices);

		printstatus(symbols);

		
		// map sequences from uniprot
		ArrayList<ArrayList<String>> sequences = mapper.mapSeqs(uniprot_ids);
		filter.checkSeq(uniprot_ids, symbols, species, domains, matrices, sequences);

		printstatus(symbols);

		
		// map binding from uniprot
		ArrayList<ArrayList<String>> bindings = mapper.mapBinds(uniprot_ids);


		// map transfac classification 
		ObjectRW objectrw = new ObjectRW();
		HashMap<String,String> transHM = (HashMap<String, String>) objectrw.read("transmap");
		// add manual entries to HM
		transHM = fileparser.annoteHM("manual_annotation.txt", transHM);
		ArrayList<String> transfac = mapper.mapTrans(domains, transHM);

		
		// flag entries of missing domain/transfac and binding
		boolean[] flags = filter.setflags(uniprot_ids, symbols, species, transfac, matrices, sequences, bindings);
		
		
		//// TODO: InterProScan
		// use flagged "IP" entries only.

		////
		printstatus(symbols);
		

		// write output file
		OutputFileWriter output_writer = new OutputFileWriter();
		System.out.println("Writing output file...");
		output_writer.writeOutfile(symbols, uniprot_ids, species, matrices, sequences, bindings, transfac, flags, "output.txt");
		
		
		// Finished!
		System.out.println("Done.");
		
	}
	
	private static void printstatus(ArrayList<String> symbols){
		System.out.println("Processing " + symbols.size() + " entries...");
	}
	
}
