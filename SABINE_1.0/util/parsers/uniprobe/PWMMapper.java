package util.parsers.uniprobe;
/*
 * ===============================================
 * (C) Florian Topf, University of Tuebingen, 2010
 * ===============================================
 */

import java.util.ArrayList;


public class PWMMapper {
	
	ArrayList<String> Pnames_mapped = new ArrayList<String>();
	ArrayList<String> UniProtIDs_mapped = new ArrayList<String>();
	ArrayList<ArrayList<String>> Domainholder_mapped = new ArrayList<ArrayList<String>>();
	ArrayList<String> Species_mapped = new ArrayList<String>();
	ArrayList<ArrayList<String>> Pwms_mapped = new ArrayList<ArrayList<String>>();
	
	//MN
	ArrayList<String> MNs_mapped = new ArrayList<String>();
	
	
	public void map(ArrayList<String> Pnames, ArrayList<String> UniProtIDs, ArrayList<ArrayList<String>> Domainholder, ArrayList<String> Species, ArrayList<String> PWMnames, ArrayList<ArrayList<String>> Pwms, ArrayList<String> MNs) {
		
		mapstat(Pnames, PWMnames);
		
		int count = 0;
		int err = 0;
		
		for (int i = 0; i < PWMnames.size(); i++) {
			
			if (Pnames.contains(PWMnames.get(i))) {
				
				for (int j = 0; j < Pnames.size(); j++) {
					
					if (PWMnames.get(i).equalsIgnoreCase(Pnames.get(j))) {
						
						Pnames_mapped.add(Pnames.get(j));
						UniProtIDs_mapped.add(UniProtIDs.get(j));
						Domainholder_mapped.add(Domainholder.get(j));
						Species_mapped.add(Species.get(j));
						Pwms_mapped.add(Pwms.get(i));
						
						//MN
						MNs_mapped.add(MNs.get(i));

						
						count++;
						
						break;
						
					}
					
				}
				
			} else {
				err++;
				System.out.println("PWM, but no further information: "+PWMnames.get(i));
			}
			
		}
		System.out.println(err+" PWM(s) without further information.");
		System.out.println(count+" entries mapped.");
		
	}
	
	public void mapstat(ArrayList<String> Pnames, ArrayList<String> PWMnames) {
		
		int count = 0;
		
		for (int i = 0; i < Pnames.size(); i++) {
			
			if (PWMnames.contains(Pnames.get(i)) != true) {
				
				count++;
				System.out.println("No PWM for "+Pnames.get(i) +" found.");
				
			}
		}
		System.out.println(count+" without corresponding PWM.");
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
	
	// MN
	public ArrayList<String> getMNs_mapped() {
		return MNs_mapped;
	}


}
