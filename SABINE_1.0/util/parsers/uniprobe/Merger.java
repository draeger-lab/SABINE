package util.parsers.uniprobe;
/*
 * ===============================================
 * (C) Florian Topf, University of Tuebingen, 2010
 * ===============================================
 */

import java.util.ArrayList;


public class Merger {
	
	ArrayList<String> Pnames_merged = new ArrayList<String>();
	ArrayList<String> UniProtIDs_merged = new ArrayList<String>();
	ArrayList<ArrayList<String>> Domainholder_merged = new ArrayList<ArrayList<String>>();
	ArrayList<String> Species_merged = new ArrayList<String>();
	ArrayList<ArrayList<String>> Pwms_merged = new ArrayList<ArrayList<String>>();
	ArrayList<String> Transfac_merged = new ArrayList<String>();
	ArrayList<String> Sequences_merged = new ArrayList<String>();
	ArrayList<ArrayList<String>> DNAbinding_merged = new ArrayList<ArrayList<String>>();
	ArrayList<String> MNs_merged = new ArrayList<String>();
	
	
	public void merge(ArrayList<String> pnames, ArrayList<String> uniProtIDs, ArrayList<ArrayList<String>> domainholder, ArrayList<String> species, ArrayList<ArrayList<String>> pwms, ArrayList<String> transfac, ArrayList<String> rawUniprotIDs, ArrayList<String> rawSequences, ArrayList<ArrayList<String>> rawDNAbinding, ArrayList<String> MNs) {
		
		mapstat(uniProtIDs, rawUniprotIDs);
		
		int count = 0;
		int err = 0;
		
		for (int i = 0; i < rawUniprotIDs.size(); i++) {
			
			if (uniProtIDs.contains(rawUniprotIDs.get(i))) {
				
				for (int j = 0; j < uniProtIDs.size(); j++) {
					
					if (rawUniprotIDs.get(i).equalsIgnoreCase(uniProtIDs.get(j))) {
						
						
						Pnames_merged.add(pnames.get(j));
						UniProtIDs_merged.add(uniProtIDs.get(j));
						Domainholder_merged.add(domainholder.get(j));
						Species_merged.add(species.get(j));
						Pwms_merged.add(pwms.get(j));
						Transfac_merged.add(transfac.get(j));
						Sequences_merged.add(rawSequences.get(i));
						DNAbinding_merged.add(rawDNAbinding.get(i));
						
						// MN
						MNs_merged.add(MNs.get(j));
						
						count++;
						
						break;
						
					}
					
				}
				
			} else {
				err++;
				System.out.println("Sequence/DNA-Binding, but no further information: "+uniProtIDs.get(i));
			}
			
		}
		System.out.println(err+" entries without further information.");
		
		System.out.println(count+" entries merged.");
		
	}
	
	public void mapstat(ArrayList<String> uniProtIDs, ArrayList<String> rawUniprotIDs) {
		
		int count = 0;
		
		for (int i = 0; i < uniProtIDs.size(); i++) {
			
			if (rawUniprotIDs.contains(uniProtIDs.get(i)) != true) {
				
				count++;
				System.out.println("No Sequence/DNA-Binding for "+uniProtIDs.get(i) +" found.");
				
			}
		}
		System.out.println(count+" without corresponding Sequence/DNA-Binding.");
	}
	
	
	public ArrayList<String> getPnames_merged() {
		return Pnames_merged;
	}
	public ArrayList<String> getUniProtIDs_merged() {
		return UniProtIDs_merged;
	}
	public ArrayList<ArrayList<String>> getDomainholder_merged() {
		return Domainholder_merged;
	}
	public ArrayList<String> getSpecies_merged() {
		return Species_merged;
	}
	public ArrayList<ArrayList<String>> getPwms_merged() {
		return Pwms_merged;
	}
	public ArrayList<String> getTransfac_merged() {
		return Transfac_merged;
	}
	public ArrayList<String> getSequences_merged() {
		return Sequences_merged;
	}
	public ArrayList<ArrayList<String>> getDNAbinding_merged() {
		return DNAbinding_merged;
	}
	
	// MN
	
	public ArrayList<String> getMNs_merged() {
		return MNs_merged;
	}
	

}
