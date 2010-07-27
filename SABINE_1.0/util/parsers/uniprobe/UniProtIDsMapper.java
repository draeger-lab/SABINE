package util.parsers.uniprobe;
/*
 * ===============================================
 * (C) Florian Topf, University of Tuebingen, 2010
 * ===============================================
 */

import java.util.ArrayList;
import java.util.HashMap;

/*
 * Maps UniProtIDs to primary ones
 */
public class UniProtIDsMapper {
	
	ArrayList<String> UniProtIDs_mapped = new ArrayList<String>();
	
	public ArrayList<String> mapIDs(ArrayList<String> ids, HashMap<String, String> hm) {
		
		for (int i = 0; i < ids.size(); i++) {
			
			String tmp = ids.get(i).toUpperCase().trim();
			String pri = tmp;
			
			// already primary
			if (hm.containsValue(tmp) != true) {
				
				// map primary
				if (hm.containsKey(tmp)) {
					
					pri = hm.get(tmp);
					
					System.out.println(tmp+" --> "+pri);
				
				}	
			}
			UniProtIDs_mapped.add(pri);
		}
		return UniProtIDs_mapped;
	}

}
