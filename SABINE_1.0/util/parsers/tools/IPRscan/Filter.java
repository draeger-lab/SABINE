package util.parsers.tools.IPRscan;
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
			System.out.println(count + " entrie(s) having flag set removed.");
		}
}

