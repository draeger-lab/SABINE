package util.parsers.tools.IPRscan;
/*
 * ===============================================
 * (C) Florian Topf, University of Tuebingen, 2010
 * ===============================================
 */

import java.util.ArrayList;

/*
 * diverse filter
 */

public class Filter {
	
	public void checkFlag(ArrayList<String> symbols,
			ArrayList<ArrayList<String>> uniprot_ids, ArrayList<String> species,
			ArrayList<ArrayList<ArrayList<String>>> matrices,
			ArrayList<ArrayList<String>> sequences,
			ArrayList<ArrayList<String>> bindings, ArrayList<String> transfac,
			ArrayList<Boolean> flags_arr, ArrayList<ArrayList<String>> MN) {

			boolean[] relevant = new boolean[uniprot_ids.size()];
			
			// find relevant entries
			for (int i=0; i<uniprot_ids.size(); i++) {
		
				if (! flags_arr.get(i)) {
					
					relevant[i] = true;
				}
			}
			
			int count = 0;
			
			// remove irrelevant entries
			for (int i=uniprot_ids.size()-1; i>=0; i--) {
				if (! relevant[i]) {
					
					count++;
					
					symbols.remove(i);
					uniprot_ids.remove(i);
					species.remove(i);
					matrices.remove(i);
					sequences.remove(i);
					bindings.remove(i);
					transfac.remove(i);
					flags_arr.remove(i);
					MN.remove(i);
				}
			}
			System.out.println(count + " entrie(s) with flag removed.");
		}
	
	
	// flag missing entries for later processing
	public ArrayList<Boolean> setflags(ArrayList<ArrayList<String>> bindings, ArrayList<String> transfac, ArrayList<Boolean> flags_arr) {
		
		int count = 0;
		
		// find missing entries
		for (int i=0; i<flags_arr.size(); i++) {
				
			if ((bindings.get(i).isEmpty()) || (transfac.get(i).isEmpty())) {
				flags_arr.set(i, true);
				count++;
			}
			else {
				flags_arr.set(i, false);
			}
		}
		
		System.out.println(count+" entrie(s) flagged.");
		return flags_arr;
	}
	
	
	// pick relevant uniprot_ids / sequences
	public void picker(ArrayList<ArrayList<String>> uniprot_ids, ArrayList<ArrayList<String>> sequences, ArrayList<ArrayList<String>> bindings) {
		
		System.out.println("Picking...");
		
		for (int i = 0; i < uniprot_ids.size(); i++) {
		
			// final result
			ArrayList<String> new_ids = new ArrayList<String>();
			ArrayList<String> new_seqs = new ArrayList<String>();
			
			// id <-> sequence
			ArrayList<String> curr_ids = new ArrayList<String>();
			ArrayList<String> curr_seqs = new ArrayList<String>();
			
			for (int j = 0; j < sequences.get(i).size(); j++) {
				
				String seq_tmp = sequences.get(i).get(j);
				String[] seq_tmp_arr = seq_tmp.split(";");
				String seq_id = seq_tmp_arr[0].trim();
				String seq = seq_tmp_arr[1].trim();
				
				curr_ids.add(seq_id);
				curr_seqs.add(seq);
			}
			
			
			
			// relevant ids (those with corresponding binding)
			ArrayList<String> rel_ids = new ArrayList<String>();
			
			// search bindings for id match
			if (!bindings.get(i).isEmpty()) {
				
				for (int j = 0; j < bindings.get(i).size(); j++) {
					
					String bind_tmp = bindings.get(i).get(j);
					
					String[] bind_tmp_arr = bind_tmp.split("    ");
					String bind_id = bind_tmp_arr[0].trim();
					
					if (curr_ids.contains(bind_id)) {
						// relevant id and sequence found
						if (!rel_ids.contains(bind_id)) {
							rel_ids.add(bind_id);
						}
					} else {
						System.err.println("Dataset inconsistent!");
						System.exit(1);
					}
				}
				for (int j = 0; j < rel_ids.size(); j++) {
					
					// get rel_id in curr_ids
					for (int k = 0; k < curr_ids.size(); k++) {
						
						String rel_id = rel_ids.get(j);
						
						if (curr_ids.get(k).equals(rel_id)) {
							
							new_ids.add(curr_ids.get(k));
							new_seqs.add(curr_seqs.get(k));
						}
					}
				}
				
				// TODO optional
				// only two ids, sequences and bindings
			
				/*
				if (rel_ids.size() == 1) {
					
					String rel_id = rel_ids.get(0);
					
					// get rel_id in curr_ids
					for (int j = 0; j < curr_ids.size(); j++) {
						
						if (curr_ids.get(j).equals(rel_id)) {
							
							new_ids.add(curr_ids.get(j));
							new_seqs.add(curr_seqs.get(j));
							
						}
					}
				}
				else if (rel_ids.size() > 1) {
					
					for (int j = 0; j < rel_ids.size(); j++) {
						
						// get rel_id in curr_ids
						for (int k = 0; k < curr_ids.size(); k++) {
							
							String rel_id = rel_ids.get(j);
							
							if (curr_ids.get(k).equals(rel_id)) {
								
								new_ids.add(curr_ids.get(k));
								new_seqs.add(curr_seqs.get(k));
							}
						}
					}
				}
				else {
					System.err.println("Dataset inconsistent!");
					System.exit(1);
				}
				 */
			}	
			else {
				// no binding information
				// pick longest sequence
				
				String max1 = "";
				String id1 = "";
				
				// more than one sequence
				if (curr_seqs.size() > 1) {
					
					for (int j = 0; j < curr_seqs.size(); j++) {
						
						String seq = curr_seqs.get(j);
						
						if (seq.length() > max1.length()) {
							
							max1 = seq;
							id1 = curr_ids.get(j);
						}
					}
					
					new_ids.add(id1);
					new_seqs.add(max1);
					
				}
				else {
					
					new_ids.add(curr_ids.get(0));
					new_seqs.add(curr_seqs.get(0));
					
				}
				
			}
			
			
			/*
			else {
			// no binding information
			// pick longest sequences
			
			String max1 = "";
			String max2 = null;
			
			String id1 = "";
			String id2 = null;
			
			// more than one sequence
			if (curr_seqs.size() > 1) {
				
				for (int j = 0; j < curr_seqs.size(); j++) {
					
					String seq = curr_seqs.get(j);
					
					if (seq.length() > max1.length()) {
						
						max2 = max1;
						max1 = seq;
						
						id2 = id1;
						id1 = curr_ids.get(j);
					}
				}
				
				new_ids.add(id1);
				new_seqs.add(max1);
				
				if ((!max2.isEmpty()) && (!new_seqs.contains(max2))) {
					
					new_ids.add(id2);
					new_seqs.add(max2);
					
				}
				
			}
			else {
				
				new_ids.add(curr_ids.get(0));
				new_seqs.add(curr_seqs.get(0));
				
			}
			
		}
		*/
			
			uniprot_ids.set(i, new_ids);
			sequences.set(i, new_seqs);
		}
			
	}

}

