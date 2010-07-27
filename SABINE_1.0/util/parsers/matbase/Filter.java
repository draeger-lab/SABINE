package util.parsers.matbase;
/*
 * ===============================================
 * (C) Florian Topf, University of Tuebingen, 2010
 * ===============================================
 */

import java.util.ArrayList;
import java.util.Iterator;

/*
 * diverse filter
 */

public class Filter {
	
	private ArrayList<String> gene_ids_fi = null;
	private ArrayList<String> symbols_fi = null;
	private ArrayList<Integer> taxon_ids_fi = null;
	private ArrayList<ArrayList<String>> domains_fi = null;
	private ArrayList<ArrayList<ArrayList<String>>> matrices_fi = null;

	////
	// filters applied during dump
	////
	
	// filter redundant entries
	public void checkDub(ArrayList<String> gene_ids, ArrayList<String> symbols, ArrayList<Integer> taxon_ids) {
		
		ArrayList<String> gene_ids_tmp = new ArrayList<String>();
		ArrayList<String> symbols_tmp = new ArrayList<String>();
		ArrayList<Integer> taxon_ids_tmp = new ArrayList<Integer>();
		
		int count = 0;
		
		for (int i=0; i<gene_ids.size(); i++) {
			
			String curr = gene_ids.get(i);
			
			if (! gene_ids_tmp.contains(curr)) {
				
				gene_ids_tmp.add(curr);
				symbols_tmp.add(symbols.get(i));
				taxon_ids_tmp.add(taxon_ids.get(i));
			}
			else {
				count++;
			}
			
		}
		
		// refresh global
		gene_ids_fi = gene_ids_tmp;
		symbols_fi = symbols_tmp;
		taxon_ids_fi = taxon_ids_tmp;
		
		System.out.println(count + " dublicate(s) removed.");
	}
	
	////
	// filters applied after dump
	////
	
	// filter entries without a corresponding matrix
	public void checkMat(ArrayList<String> gene_ids, ArrayList<String> symbols, ArrayList<Integer> taxon_ids, ArrayList<ArrayList<String>> domains, ArrayList<ArrayList<ArrayList<String>>> matrices) {
		
		boolean[] relevant = new boolean[gene_ids.size()];
		
		// find relevant entries
		for (int i=0; i<gene_ids.size(); i++) {
			
			if (matrices.get(i).get(0).get(0).length() != 0) {
				relevant[i] = true;
			}
			
		}
		
		int count = 0;
		
		// remove irrelevant entries
		for (int i=gene_ids.size()-1; i>=0; i--) {
			if (! relevant[i]) {
				
				count++;
				
				gene_ids.remove(i);
				symbols.remove(i);
				taxon_ids.remove(i);
				domains.remove(i);
				matrices.remove(i);
			}
		}
		
		// refresh global
		gene_ids_fi = gene_ids;
		symbols_fi = symbols;
		taxon_ids_fi = taxon_ids;
		domains_fi = domains;
		matrices_fi = matrices;
		
		System.out.println(count + " entrie(s) without matrix removed.");
	}
		
		
	// filter unsupported species
	public void checkOrg(ArrayList<String> gene_ids, ArrayList<String> symbols, ArrayList<String> species, ArrayList<ArrayList<String>> domains, ArrayList<ArrayList<ArrayList<String>>> matrices, ArrayList<String> organisms) {

		ArrayList<String> unsupported = new ArrayList<String>();
		
		boolean[] relevant = new boolean[species.size()];
		
		// find relevant entries
		for (int i=0; i<species.size(); i++) {
			if (organisms.contains(species.get(i).toLowerCase())) {
				
				relevant[i] = true;
			}
		}
		
		int count = 0;
		
		// remove irrelevant entries
		for (int i=species.size()-1; i>=0; i--) {
			if (! relevant[i]) {
				
				
				// add unsupported to a list
				String tmp = species.get(i);
				if (! unsupported.contains(tmp)) {
					unsupported.add(tmp);
				}
				
				count++;
				
				gene_ids.remove(i);
				symbols.remove(i);
				species.remove(i);
				domains.remove(i);
				matrices.remove(i);
			}
		}
		
		// print found unsupported species
		System.out.println("Unsupported species found:");
		Iterator<String> itr0 = unsupported.iterator();
	    while (itr0.hasNext()) {
	      String element = itr0.next();
	      System.out.println(element);
	    }
		System.out.println(count + " entrie(s) of unsupported species removed.");
	}
	
	
	// filter entries without uniprotid
	public void checkIDs(ArrayList<ArrayList<String>> uniprot_ids, ArrayList<String> symbols, ArrayList<String> species, ArrayList<ArrayList<String>> domains, ArrayList<ArrayList<ArrayList<String>>> matrices) {
		
		boolean[] relevant = new boolean[uniprot_ids.size()];
		
		// find relevant entries
		for (int i=0; i<uniprot_ids.size(); i++) {
	
			if (uniprot_ids.get(i) != null) {
				
				relevant[i] = true;
			}
		}
		
		int count = 0;
		
		// remove irrelevant entries
		for (int i=uniprot_ids.size()-1; i>=0; i--) {
			if (! relevant[i]) {
				
				count++;
				
				uniprot_ids.remove(i);
				symbols.remove(i);
				species.remove(i);
				domains.remove(i);
				matrices.remove(i);
			}
		}
		System.out.println(count + " entrie(s) without UniProtID removed.");
	}
	
	
	// returns only the first two words
	public String filtername(String input) {
		String output = null;
		
		input = input.trim();
		
		if (input.contains(" ")) {
			String tmp[] = input.split(" ");
			output = tmp[0] + " " + tmp[1];
		}
		else {
			output = input;
		}

		return output;
	}
	
	
	public void checkSeq(ArrayList<ArrayList<String>> uniprot_ids, ArrayList<String> symbols, ArrayList<String> species, ArrayList<ArrayList<String>> domains, ArrayList<ArrayList<ArrayList<String>>> matrices, ArrayList<ArrayList<String>> sequences) {
		
		boolean[] relevant = new boolean[uniprot_ids.size()];
		
		// find relevant entries
		for (int i=0; i<uniprot_ids.size(); i++) {
	
			if (! sequences.get(i).isEmpty()) {
				
				relevant[i] = true;
			}
		}
		
		int count = 0;
		
		// remove irrelevant entries
		for (int i=uniprot_ids.size()-1; i>=0; i--) {
			if (! relevant[i]) {
				
				count++;
				
				uniprot_ids.remove(i);
				symbols.remove(i);
				species.remove(i);
				domains.remove(i);
				matrices.remove(i);
				sequences.remove(i);
			}
		}
		System.out.println(count + " entrie(s) without sequence information removed.");
	}
	
	
	// flag missing entries for later processing
	public boolean[] setflags(ArrayList<ArrayList<String>> uniprot_ids, ArrayList<String> symbols, ArrayList<String> species, ArrayList<String> transfac, ArrayList<ArrayList<ArrayList<String>>> matrices, ArrayList<ArrayList<String>> sequences, ArrayList<ArrayList<String>> bindings) {
		
		boolean[] flags = new boolean[uniprot_ids.size()];
		
		int count = 0;
		
		// find missing entries
		for (int i=0; i<uniprot_ids.size(); i++) {
				
			if ((bindings.get(i).isEmpty()) || (transfac.get(i).isEmpty())) {
				flags[i] = true;
				count++;
			}
		}
		
		System.out.println(count+" entrie(s) flagged.");
		return flags;
	}

	
	// depreciated
	public ArrayList<String> getGene_ids_fi() {
		return gene_ids_fi;
	}
	
	public ArrayList<String> getSymbols_fi() {
		return symbols_fi;
	}

	public ArrayList<Integer> getTaxon_ids_fi() {
		return taxon_ids_fi;
	}
	
	public ArrayList<ArrayList<String>> getDomains_fi() {
		return domains_fi;
	}
	
	public ArrayList<ArrayList<ArrayList<String>>> getMatrices_fi() {
		return matrices_fi;
	}
}

