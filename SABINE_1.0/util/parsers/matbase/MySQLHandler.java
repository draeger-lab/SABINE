package util.parsers.matbase;
/*
 * ===============================================
 * (C) Florian Topf, University of Tuebingen, 2010
 * ===============================================
 */

import java.io.Serializable;
import java.sql.*;
import java.util.ArrayList;
import java.util.HashMap;

public class MySQLHandler implements Serializable{
	
	// serializable id
	private static final long serialVersionUID = 2645541236557306118L;
	
	ArrayList<String> gene_id = new ArrayList<String>();
	ArrayList<String> symbol = new ArrayList<String>();
	ArrayList<Integer> taxon_id = new ArrayList<Integer>();
	
	HashMap<Integer,String> specmap = new HashMap<Integer,String>();
	
	
	// initialize database support
	MySQLAPI mysqlapi = new MySQLAPI();
	
	// initialize database converter
	MySQLResultSetRefiner mysqlresultsetrefiner = new MySQLResultSetRefiner();
	
	
	// get tf data
	public void dumpData(Connection cn) {
		
		// get gene_id, symbol, taxon_id
		ResultSet rs1 = mysqlapi.query(cn, "SELECT transcription_factor.gene_id, transcription_factor.symbol, transcription_factor.taxon_id FROM transcription_factor");
		
		// process it!
		gene_id = mysqlresultsetrefiner.processStr(rs1, "gene_id");
		symbol = mysqlresultsetrefiner.processStr(rs1, "symbol");
		taxon_id = mysqlresultsetrefiner.processInt(rs1, "taxon_id");
				
		// get mapping: taxon_id -> organism (species)
		ResultSet rs2 = mysqlapi.query(cn, "SELECT organisms.taxon_id, organisms.name FROM organisms");
				
		// process it!
		specmap = mysqlresultsetrefiner.processHM(rs2, "taxon_id", "name");
		
	}

	
	// get domain of tf
	public ArrayList<ArrayList<String>> dumpDomains(Connection cn, ArrayList<String> gene_ids) {
		
		ArrayList<ArrayList<String>> domains = new ArrayList<ArrayList<String>>();
		
		System.out.println("Fetching domains...");
		
		for (int i = 0; i < gene_ids.size(); i++) {
			
			ResultSet rs = mysqlapi.query(cn, "SELECT domain.name FROM transcription_factor JOIN tf2family JOIN family JOIN family2domain JOIN domain ON transcription_factor.id = tf2family.tf_id AND tf2family.family_id = family.id AND family.id = family2domain.family_id AND family2domain.domain_id = domain.id WHERE transcription_factor.gene_id = \"" + gene_ids.get(i) + "\"");
			
			domains.add(mysqlresultsetrefiner.processDom(rs));
			
		}
		
		return domains;
	}
	
	
	// get data for a specific tf.symbol and tf.taxon_id
	public ArrayList<ArrayList<ArrayList<String>>> dumpMatrices(Connection cn, ArrayList<String> gene_ids) {
		
		ArrayList<ArrayList<ArrayList<String>>> matrices = new ArrayList<ArrayList<ArrayList<String>>>();
		
		System.out.println("Dumping matrices...");
		
		int count = 0;
		
		for (int i = 0; i < gene_ids.size(); i++) {
			// get id, position, a, c, g, t
			ResultSet rs = mysqlapi.query(cn, "SELECT distributionmatrix.position, distributionmatrix.a, distributionmatrix.c, distributionmatrix.g, distributionmatrix.t FROM transcription_factor JOIN tf2family JOIN family JOIN matrix JOIN distributionmatrix ON transcription_factor.id = tf2family.tf_id AND tf2family.family_id = family.id AND family.id = matrix.family_id AND matrix.id = distributionmatrix.matrix_id WHERE transcription_factor.gene_id = \"" + gene_ids.get(i) + "\"");
			
			// add processed entries ...
			matrices.add(mysqlresultsetrefiner.processMat(rs));
			
			count += matrices.get(i).size();
			
		}
		
		System.out.println("Dumped " + count + " matrices.");
		
		return matrices;
	}
	

	public ArrayList<String> getGene_id() {
		return gene_id;
	}

	public ArrayList<String> getSymbol() {
		return symbol;
	}

	public ArrayList<Integer> getTaxon_id() {
		return taxon_id;
	}

	public HashMap<Integer, String> getSpecmap() {
		return specmap;
	}
	




}
