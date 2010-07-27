package util.parsers.matbase;
/*
 * ===============================================
 * (C) Florian Topf, University of Tuebingen, 2010
 * ===============================================
 */

import org.apache.commons.httpclient.*;
import org.apache.commons.httpclient.methods.*;

/*
 * UniProt API
 */

public class UniProtAPI {
	
	private static final String base = "http://www.uniprot.org";
	
	
	// submits query to uniprot and retrieves result
	private String run(String tool, NameValuePair[] params) throws Exception {
		
		String result = null;
		
		HttpClient client = new HttpClient();
		
		String location = base + '/' + tool + '/';
		HttpMethod method = new PostMethod(location);
		((PostMethod) method).addParameters(params);
		method.setFollowRedirects(false);
		int status = client.executeMethod(method);
		if (status == HttpStatus.SC_MOVED_TEMPORARILY) {
			location = method.getResponseHeader("Location").getValue();
			method.releaseConnection();
			method = new GetMethod(location);
			status = client.executeMethod(method);
		}
		
		while (true) {
			int wait = 0;
			Header header = method.getResponseHeader("Retry-After");
			
			if (header != null) {
				wait = Integer.valueOf(header.getValue());
			}
			if (wait == 0) {
				break;
			}
			
			//log.info("Waiting (" + wait + ")...");
			Thread.sleep(wait * 1000);
			method.releaseConnection();
			method = new GetMethod(location);
			status = client.executeMethod(method);
		}
		
		// get result as stream
		if (status == HttpStatus.SC_OK) {
			result = (method.getResponseBodyAsString());
		}
		
		method.releaseConnection();
		
		return result;
		
	}
	
	
	// get idmapping
	public String getIDs(String from) {
		String result = null;	
		
		try {
			result = run("mapping", new NameValuePair[] {
					new NameValuePair("from", "P_ENTREZGENEID"),
					new NameValuePair("to", "ACC"),
					new NameValuePair("format", "tab"),
					new NameValuePair("query", from),
			});
		} catch (Exception e) {
			e.printStackTrace();
		}
		
		return result;
	}
	

	// get sequences
	public String getSeqs(String query) {
		String result = null;	
		
		try {
			result = run("batch", new NameValuePair[] {
					new NameValuePair("format", "fasta"),
					new NameValuePair("query", query),
			});
		} catch (Exception e) {
			e.printStackTrace();
		}
		
		return result;
		
	}
	
	
	// get binding
	public String getBinds(String query) {
		String result = null;	
		
		try {
			result = run("batch", new NameValuePair[] {
					new NameValuePair("format", "gff"),
					new NameValuePair("query", query),
			});
		} catch (Exception e) {
			e.printStackTrace();
		}
		
		return result;
		
	}
	
}


