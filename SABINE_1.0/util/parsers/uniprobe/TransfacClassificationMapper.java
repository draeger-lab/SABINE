package util.parsers.uniprobe;
/*
 * ===============================================
 * (C) Florian Topf, University of Tuebingen, 2010
 * ===============================================
 */

import java.util.ArrayList;
import java.util.HashMap;


public class TransfacClassificationMapper {
	
	ArrayList<String> Pnames_mapped = new ArrayList<String>();
	ArrayList<String> UniProtIDs_mapped = new ArrayList<String>();
	ArrayList<ArrayList<String>> Domainholder_mapped = new ArrayList<ArrayList<String>>();
	ArrayList<String> Species_mapped = new ArrayList<String>();
	ArrayList<ArrayList<String>> Pwms_mapped = new ArrayList<ArrayList<String>>();
	ArrayList<String> Transfac_mapped = new ArrayList<String>();
	ArrayList<String> MNs_mapped = new ArrayList<String>();

	public void map(HashMap<String, String> hmTransfac, ArrayList<String> pnames, ArrayList<String> uniProtIDs, ArrayList<ArrayList<String>> domainholder, ArrayList<String> species, ArrayList<ArrayList<String>> pwms, ArrayList<String> MNs) {
		
		int count = 0;
		int success = 0;
		int empty = 0;
		
		for (int i = 0 ; i < domainholder.size(); i++) {
			
			ArrayList<String> tmp = domainholder.get(i);
			
			//System.out.print(tmp);
			//System.out.print(pnames.get(i)+";");
						
			if (tmp.size() > 1) {
				
				String max = "";
				String hit = "";
				
				for (int j = 0 ; j < tmp.size(); j++) {
					
					String entry = tmp.get(j);
					
					if (entry.isEmpty()) {
						
						empty++;
						
					} else {
						
						if (hmTransfac.containsKey(entry)) {
							
							hit = hmTransfac.get(entry);
							if (hit.length() > max.length()) {
								max = hit;
							}
						}
					}
				}
				
				if (max != "") {
					success++;
					
					Pnames_mapped.add(pnames.get(i));
					UniProtIDs_mapped.add(uniProtIDs.get(i));
					Domainholder_mapped.add(tmp);
					Species_mapped.add(species.get(i));
					Pwms_mapped.add(pwms.get(i));
					Transfac_mapped.add(max);
					//System.out.println(max);
					
					//MN
					MNs_mapped.add(MNs.get(i));
					
				} else {
					
					count++;
					System.out.println("No mapping found for "+tmp);
					
				}

				
				
			} else {
				
				String entry = tmp.get(0); 

				if (entry.isEmpty()) {
					
					empty++;
					
				} else {
					
					if (hmTransfac.containsKey(entry)) {
						
						success++;
						
						Pnames_mapped.add(pnames.get(i));
						UniProtIDs_mapped.add(uniProtIDs.get(i));
						Domainholder_mapped.add(domainholder.get(i));
						Species_mapped.add(species.get(i));
						Pwms_mapped.add(pwms.get(i));
						Transfac_mapped.add(hmTransfac.get(entry));
						
						// MN
						MNs_mapped.add(MNs.get(i));
						
					
					} else {
						
					count++;
					System.out.println("No mapping found for: "+entry);
							
					}
						
				}
					
			}
				
		}
		System.out.println(count+" lookups returned null.");
		System.out.println(empty+" without domain-annotation.");
		System.out.println(success+" entries successfully mapped.");
	}
		


	
	public ArrayList<String> getPnames_mapped() {
		return Pnames_mapped;
	}

	public ArrayList<String> getUniProtIDs_mapped() {
		return UniProtIDs_mapped;
	}

	public ArrayList<ArrayList<String>> getDomainholder_mapped() {
		return Domainholder_mapped;
	}

	public ArrayList<String> getSpecies_mapped() {
		return Species_mapped;
	}

	public ArrayList<ArrayList<String>> getPwms_mapped() {
		return Pwms_mapped;
	}

	public ArrayList<String> getTransfac_mapped() {
		return Transfac_mapped;
	}
	
	// MN
	
	public ArrayList<String> getMNs_mapped() {
		return MNs_mapped;
	}
	
}
