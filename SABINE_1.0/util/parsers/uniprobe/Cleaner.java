package util.parsers.uniprobe;
/*
 * ===============================================
 * (C) Florian Topf, University of Tuebingen, 2010
 * ===============================================
 */
import java.util.ArrayList;

/*
 * Cleans UniProbe Factors from the unsupported 
 */
public class Cleaner {
	
	ArrayList<String> Pnames_dub = new ArrayList<String>();	
	ArrayList<String> UniProtIDs_dub = new ArrayList<String>();
	ArrayList<ArrayList<String>> Domainholder_dub = new ArrayList<ArrayList<String>>();
	ArrayList<String> Species_dub = new ArrayList<String>();
	
	// rebuild without duplicates
	public void removedubs() {
		
		for (int i = 0; i < UniProbeFactorsParser.Pnames.size(); i++) {
			
				if ((Pnames_dub.contains(UniProbeFactorsParser.Pnames.get(i))) && (UniProtIDs_dub.contains(UniProbeFactorsParser.UniProtIDs.get(i))) && (Domainholder_dub.contains(UniProbeFactorsParser.Domainholder.get(i))) && (Species_dub.contains(UniProbeFactorsParser.Species.get(i)))  == true ) {
					
					//System.out.println(UniProbeFactorsParser.Pnames.get(i)+" "+UniProbeFactorsParser.UniProtIDs.get(i)+" "+UniProbeFactorsParser.Domainholder.get(i)+" "+UniProbeFactorsParser.Species.get(i));
					System.out.println("Duplicate found: "+UniProbeFactorsParser.Pnames.get(i));
					
				} else {
										
					Pnames_dub.add(UniProbeFactorsParser.Pnames.get(i));
					UniProtIDs_dub.add(UniProbeFactorsParser.UniProtIDs.get(i));
					Domainholder_dub.add(UniProbeFactorsParser.Domainholder.get(i));
					Species_dub.add(UniProbeFactorsParser.Species.get(i));
					
				}
	    	  	    	      	  
		}

    }

	
	ArrayList<String> Pnames_spec = new ArrayList<String>();	
	ArrayList<String> UniProtIDs_spec = new ArrayList<String>();
	ArrayList<ArrayList<String>> Domainholder_spec = new ArrayList<ArrayList<String>>();
	ArrayList<String> Species_spec = new ArrayList<String>();
	
	// rebuild without unsupported species
	public void removespec() {
		
		int count = 0;
		
		for (int i = 0; i < Species_dub.size(); i++) {
			
			for (int j = 0; j < OrganismListParser.Organism.size(); j++) {
				
				if (Species_dub.get(i).equalsIgnoreCase(OrganismListParser.Organism.get(j))) {
					
					count++;
					
					Pnames_spec.add(Pnames_dub.get(i));
					UniProtIDs_spec.add(UniProtIDs_dub.get(i));
					Domainholder_spec.add(Domainholder_dub.get(i));
					Species_spec.add(Species_dub.get(i));
					
					break;
				}
	    	  
			}
	      
	    }
		System.out.println(Species_dub.size()-count+" unsupported species found.");
	}
	
	
	ArrayList<String> Pnames_cln = new ArrayList<String>();	
	ArrayList<String> UniProtIDs_cln = new ArrayList<String>();
	ArrayList<ArrayList<String>> Domainholder_cln = new ArrayList<ArrayList<String>>();
	ArrayList<String> Species_cln = new ArrayList<String>();
	
	// rebuild without missing UniProtIDs
	public void removeids() {
		
		for (int i = 0; i < UniProtIDs_spec.size(); i++) {
			
			if (UniProtIDs_spec.get(i).isEmpty() != true) {
				
				//System.out.println(UniProtIDs_spec.get(i));
				
				Pnames_cln.add(Pnames_spec.get(i));
				UniProtIDs_cln.add(UniProtIDs_spec.get(i));
				Domainholder_cln.add(Domainholder_spec.get(i));
				Species_cln.add(Species_spec.get(i));

			} else {
				
				System.out.println("Missing UniProtID: "+Pnames_spec.get(i));
				
			}
	      
	    }
		
	}
	
	public void cleanup() {

		removedubs();
		removespec();
		removeids();
		
	}
	
	public ArrayList<String> getPnames() {
		return Pnames_cln;
	}
	
	public ArrayList<String> getUniProtIDs() {
		return UniProtIDs_cln;
	}
	
	public ArrayList<ArrayList<String>> getDomainholder() {
		return Domainholder_cln;
	}
	
	public ArrayList<String> getSpecies() {
		return Species_cln;
	}
	
}
