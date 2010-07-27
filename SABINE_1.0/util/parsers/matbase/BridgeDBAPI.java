package util.parsers.matbase;
/*
 * ===============================================
 * (C) Florian Topf, University of Tuebingen, 2010
 * ===============================================
 */

import java.util.ArrayList;
import java.util.Set;

import org.bridgedb.BridgeDb;
import org.bridgedb.DataSource;
import org.bridgedb.IDMapper;
import org.bridgedb.IDMapperException;
import org.bridgedb.Xref;

/*
 * BridgeDB API
 */

public class BridgeDBAPI {
	
	// connect to mapper and gets destinations
	private Set<Xref> domapping(String connectString, String id, DataSource src, DataSource dest) throws IDMapperException {
		
		// connect to the mapper
		IDMapper mapper = BridgeDb.connect(connectString);
		Xref srcRef = new Xref (id, src);
		
		return mapper.mapID(srcRef, dest);
		
	}
	
	
	// maps entrezgeneID to uniprotID using synergizer web service
	private Set<Xref> synergize(String id, String spec) {
		
		 Set<Xref> destSet = null;
		
		try {
			Class.forName("org.bridgedb.webservice.synergizer.IDMapperSynergizer");
			destSet = domapping("idmapper-synergizer:authority=ncbi&species="+spec, id,
					DataSource.getByFullName("entrezgene"),
					DataSource.getByFullName("uniprot"));
		} catch (Exception e) {
			e.printStackTrace();
		}
		return destSet;
	}
	
	
	// transform result from Set<Xref> to ArrayList
	public ArrayList<String> retrieve(String id, String spec) {
		ArrayList<String> results = new ArrayList<String>();
		
		Set<Xref> destSet = synergize(id, spec);
		
		for (Xref destRef : destSet) {
			results.add(destRef.getId());
		}
		
		return results;
	}
	
}
