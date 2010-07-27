package util.parsers.matbase;
/*
 * ===============================================
 * (C) Florian Topf, University of Tuebingen, 2010
 * ===============================================
 */

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;

/*
 * map diverse
 */

public class Mapper {
	
	Filter filter = new Filter();
	
	// map taxon_ids to species
	public ArrayList<String> mapSpec(ArrayList<Integer> taxon_ids, HashMap<Integer, String> specmap) {
		
		ArrayList<String> species = new ArrayList<String>();
		
		for (int i = 0; i < taxon_ids.size(); i++) {
			
			species.add(filter.filtername(specmap.get(taxon_ids.get(i))));
			
		}
		
		return species;
		
	}
	
	
	// map matbase-species names to organism_list.txt names
	public ArrayList<String> mapMatSpec(ArrayList<String> species, StringMap matspecmap) {
		
		ArrayList<String> species_conv = new ArrayList<String>();
		
		for (int i = 0; i < species.size(); i++) {
			
			String curr = species.get(i);
			
			if (matspecmap.containsKey(curr)) {
				
				species_conv.add((String)matspecmap.get(curr));
				
			}
			else {
			
				species_conv.add(curr);				
				
			}
			
		}
		
		return species_conv;
		
	}
	
	
	// map sequences to uniprot_ids
	public ArrayList<ArrayList<String>> mapSeqs(ArrayList<ArrayList<String>> uniprot_ids) {
		
		ArrayList<ArrayList<String>> sequences = new ArrayList<ArrayList<String>>();
		
		// build query
		String query = "";
		for (int i = 0; i < uniprot_ids.size(); i++) {
			Iterator<String> itr = uniprot_ids.get(i).iterator();
		    while (itr.hasNext()) {
		    	String element = itr.next();
		    	query = query.concat(" " + element);
		    }
		}
		
		UniProtAPI uniprotapi = new UniProtAPI();
		UniProtResultRefiner uniprotresultrefiner = new UniProtResultRefiner();
		
		// debugging
		//System.out.println(query.trim());
		
		System.out.println("Fetching Sequences from UniProt...");
		HashMap<String,String> seqHM = uniprotresultrefiner.processSeqs(uniprotapi.getSeqs(query.trim()));
		
		
		// assemble sequences
		for (int i = 0; i < uniprot_ids.size(); i++) {
			
			ArrayList<String> subseqs = new ArrayList<String>();
			
			Iterator<String> itr = uniprot_ids.get(i).iterator();
		    while (itr.hasNext()) {
		    	String element = itr.next();
		    	String seq = seqHM.get(element.toUpperCase());
		    	
		    	if (seq != null) {
		    		
			    	// only add new sequences
			    	if (!subseqs.contains(seq)) {
			    		subseqs.add(seq);
			    	} 
		    		
		    	}
		    	
		    }
		    
		    sequences.add(subseqs);
			
		}
		
		return sequences;
	}


	public ArrayList<ArrayList<String>> mapBinds(ArrayList<ArrayList<String>> uniprot_ids) {
		
		ArrayList<ArrayList<String>> binds = new ArrayList<ArrayList<String>>();
		
		// build query
		String query = "";
		for (int i = 0; i < uniprot_ids.size(); i++) {
			Iterator<String> itr = uniprot_ids.get(i).iterator();
		    while (itr.hasNext()) {
		    	String element = itr.next();
		    	query = query.concat(" " + element);
		    }
		}
		
		//System.out.println(query.trim());
		
		UniProtAPI uniprotapi = new UniProtAPI();
		UniProtResultRefiner uniprotresultrefiner = new UniProtResultRefiner();
		
		System.out.println("Fetching DNA Bindings from UniProt...");
		HashMap<String,ArrayList<String>> bindsHM = uniprotresultrefiner.processBinds(uniprotapi.getBinds(query.trim()));
		
		// assemble binds
		for (int i = 0; i < uniprot_ids.size(); i++) {

			// debugging
			//System.out.println(uniprot_ids.get(i));
			
			ArrayList<String> subbinds = new ArrayList<String>();
			
			ArrayList<String> tmp = new ArrayList<String>();
			
			Iterator<String> itr = uniprot_ids.get(i).iterator();
		    while (itr.hasNext()) {
		    	String id = itr.next();
		    	tmp = bindsHM.get(id);
		    	
		    	// debugging
		    	//System.out.println(tmp);
		    	
		    	if (tmp != null) {
		    		
		    		for (int j = 0; j < tmp.size(); j++) {
			    		
			    		// debugging
			    		//System.out.println(tmp.get(j));
			    		
			    		String entry = id+";"+tmp.get(j);
			    		
			    		boolean add = true;
			    		
			    		if (subbinds.isEmpty()) {
			    		}
			    		else {
			    			for (int k = 0; k < subbinds.size(); k++) {
				    			String curr = subbinds.get(k);
				    			String[] arr = curr.split(";");  
				    			if (arr[1] == tmp.get(j)) {
				    				add = false;
				    			}
			    			}
			    		}
			    		if (add) {
			    			subbinds.add(entry);
			    		}
				    }
		    	}
		    }
		    binds.add(subbinds);
		}
		return binds;
	}
	
	
	// map domains to transfac-classification
	public ArrayList<String> mapTrans(ArrayList<ArrayList<String>> domains, HashMap<String,String> transHM) {

		ArrayList<String> transfac = new ArrayList<String>();
		
		for (int i = 0; i < domains.size(); i++) {
			
			ArrayList<String> tmp = domains.get(i);
			
			String max = "";
			String hit = "";
			
			for (int j = 0 ; j < tmp.size(); j++) {
				
				String entry = tmp.get(j).replace(" domain", " ").replace("-domain", " ").replaceAll("[()]", "").toUpperCase().replace("WINGED HELIX", "").replaceAll("-", " ").trim();
				//System.out.println(entry);
				
				if (transHM.containsKey(entry)) {
					hit = transHM.get(entry);
					if (hit.length() > max.length()) {
						max = hit;
					}
				}
			}
			transfac.add(max);
		}
		return transfac;
	}
	
	
	// map gene_ids to uniprot using uniprot
	// map gene_ids to uniprot using entrez
	// map gene_ids to uniprot using synergizer (species needed)
	public ArrayList<ArrayList<String>> mapGeneIDs(ArrayList<String> gene_ids, ArrayList<String> species){
		
		// final result
		ArrayList<ArrayList<String>> uniprot_ids = new ArrayList<ArrayList<String>>();
		
		// gene_id -> uniprot_ids HM
		HashMap<String,ArrayList<String>> idHM = null;
		
		
		////
		// UniProtAPI
		// build query
		String query = "";
		for (int i = 0; i < gene_ids.size(); i++) {
			query = query.concat(" " + gene_ids.get(i));
		}
		
		// retrieve results
		UniProtAPI uniprotapi = new UniProtAPI();
		UniProtResultRefiner uniprotresultrefiner = new UniProtResultRefiner();
		
		System.out.println("Fetching IDs from UniProt...");
		idHM = uniprotresultrefiner.processIDs(uniprotapi.getIDs(query.trim()));
		
	
		////
		// Entrez Gene Wget
		EntrezGeneWget wget = new EntrezGeneWget();
		
		System.out.println("Getting IDMapping from Entrez Gene...");
		// query entrez gene
		for (int i = 0; i < gene_ids.size(); i++) {
			String id = gene_ids.get(i);
			
			ArrayList<String> destArr = new ArrayList<String>();
			
			if (!idHM.containsKey(id) && !id.startsWith("ME:") && !id.startsWith("PL:")) {
				destArr = wget.retrieveIds(id);
			}
			
			if (!destArr.isEmpty()) {
				idHM.put(id, destArr);
			}
		}

		
		////
		// BridgeDBAPI (Synergizer)
		BridgeDBAPI bridgedbapi = new BridgeDBAPI();

		System.out.println("Querying Synergizer for IDMapping...");
		// query synergizer
		for (int i = 0; i < gene_ids.size(); i++) {
			String id = gene_ids.get(i);
			String spec = species.get(i);
			
			ArrayList<String> destArr = new ArrayList<String>();
			
			if (!idHM.containsKey(id) && !id.startsWith("ME:") && !id.startsWith("PL:")) {
				destArr = bridgedbapi.retrieve(id, spec);
			}
			
			if (!destArr.isEmpty()) {
				idHM.put(id, destArr);
			}
		}

		
		////
		// assemble final result
		for (int i = 0; i < gene_ids.size(); i++) {
			
			ArrayList<String> tmp = idHM.get(gene_ids.get(i));
			
			uniprot_ids.add(tmp);
		}
		
		return uniprot_ids;
	}

}
