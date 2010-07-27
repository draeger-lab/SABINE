package util.parsers.matbase;
/*
 * ===============================================
 * (C) Florian Topf, University of Tuebingen, 2010
 * ===============================================
 */

import java.io.BufferedReader;
import java.io.ByteArrayInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.net.MalformedURLException;
import java.net.URL;
import java.util.ArrayList;

/*
 * simple wget for entrez entries
 */

public class EntrezGeneWget {
	
	// fetch data of url
	private String fetchEntry(String url) {
		String line;
		String output = "";
		BufferedReader r;
		
		try {
			r = new BufferedReader(new InputStreamReader(new URL(url).openStream()));
			while ((line = r.readLine()) != null) {
				output = output.concat(line+"\n");
			}
		} catch (MalformedURLException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}

		return output;
	}

	// retrieve uniprot_ids
	public ArrayList<String> retrieveIds(String id) {
		// url of efetch service as txt
		String url = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=gene&id="+id+"&retmode=text";
		//System.out.println(url);
		String result = fetchEntry(url);
		
		ArrayList<String> destArr = new ArrayList<String>();
		
		String line;
		
		try {
			// convert string to bufferedreader
			InputStream is = new ByteArrayInputStream(result.getBytes("UTF-8"));
			InputStreamReader isr = new InputStreamReader(is);
			BufferedReader br = new BufferedReader(isr);
			
			line = br.readLine();
			
			// skip until section reached
			while ((line = br.readLine()) != null) {
				if (line.contains("Related Sequences")) {
					break;
				}
			}
			
			// parse ids
			while ((line = br.readLine()) != null) {
				if (line.contains("UniProtKB/Swiss-Prot") || (line.contains("UniProtKB/TrEMBL"))) {
					line = br.readLine();
					String[] tmp = line.trim().split("\"");
					
					if (!destArr.contains(tmp[1])) {
						destArr.add(tmp[1]);
					}
				}
			}
		}	 
		catch(IOException e) {
			e.printStackTrace();
		}
		
		return destArr;
	}
}

