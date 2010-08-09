package util.parsers.matbase;
/*
 * ===============================================
 * (C) Florian Topf, University of Tuebingen, 2010
 * ===============================================
 */

import java.sql.ResultSet;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.HashMap;

/*
 * Converter MySQLResultSets -> Java
 */

public class MySQLResultSetRefiner {
	
	// process rs into arraylist<string>
	public ArrayList<String> processStr(ResultSet rs, String token){
		
		// move cursor before first row
		try {
			rs.beforeFirst();
		} catch (SQLException e1) {
			e1.printStackTrace();
		}
		
		ArrayList<String> sqlArray = new ArrayList<String>();
		
		System.out.println("Processing data with token '" + token + "' ...");
		
		try {
			while (rs.next()) {
				String tmp = (rs.getString(token));
				//System.out.println(tmp);
				sqlArray.add(tmp);
			}
		} catch (SQLException e) {
			e.printStackTrace();
		}
		
		System.out.println("Processed " + sqlArray.size() + " entries.");
		
		return sqlArray;
	}

	
	// process rs into arraylist<integer>
	public ArrayList<Integer> processInt(ResultSet rs, String token){
		
		// move cursor before first row
		try {
			rs.beforeFirst();
		} catch (SQLException e1) {
			e1.printStackTrace();
		}
		
		ArrayList<Integer> sqlArray = new ArrayList<Integer>();
		
		System.out.println("Processing data with token '" + token + "' ...");
		
		try {
			while (rs.next()) {
				int tmp = (rs.getInt(token));
				//System.out.println(tmp);
				sqlArray.add(tmp);
			}
		} catch (SQLException e) {
			e.printStackTrace();
		}
		
		System.out.println("Processed " + sqlArray.size() + " entries.");
		
		return sqlArray;
	}
	
	
	// process rs into hashmap
	public HashMap<Integer, String> processHM(ResultSet rs, String tok1, String tok2) {
		
		// move cursor before first row
		try {
			rs.beforeFirst();
		} catch (SQLException e1) {
			e1.printStackTrace();
		}
		
		HashMap<Integer,String> sqlMap = new HashMap<Integer,String>();
		
		System.out.println("Processing data into HashMap...");
		
		try {
			while (rs.next()) {
				int key = (rs.getInt(tok1));
				String val = (rs.getString(tok2));
				//System.out.println(key + " -> " + val);
				sqlMap.put(key, val);
			}
		} catch (SQLException e) {
			e.printStackTrace();
		}
		
		System.out.println("Processed " + sqlMap.size() + " entries.");
		
		return sqlMap;
	}
	
	
	// process rs of matrix-data (dirty) includes absolute to relative values conversion
	public ArrayList<ArrayList<String>> processMat(ResultSet rs) {
		
		// current matrix
		ArrayList<String> current = new ArrayList<String>();
		
		// bunch of matrices
		ArrayList<ArrayList<String>> bunch = new ArrayList<ArrayList<String>>();
		
		// init
		int count = -1;
		
		String astr = "";
		String cstr = "";
		String gstr = "";
		String tstr = "";
		
		try {
			while (rs.next()) {
				
				int pos = (rs.getInt("position"));
				
				if (pos > count) {
					
					count = pos;
					
					float a = (rs.getFloat("a"));
					float c = (rs.getFloat("c"));
					float g = (rs.getFloat("g"));
					float t = (rs.getFloat("t"));
					
					// absolute to relative values conversion
					float sum = a+c+g+t;
					
					a/=sum;
					c/=sum;
					g/=sum;
					t/=sum;
					
					astr = astr.concat(a + "   ");
					cstr = cstr.concat(c + "   ");
					gstr = gstr.concat(g + "   ");
					tstr = tstr.concat(t + "   ");
				} 
				// more than one matrix
				else {
					
					// save current
					current.add(0, astr.trim());
					current.add(1, cstr.trim());
					current.add(2, gstr.trim());
					current.add(3, tstr.trim());
					
					// temporal
					ArrayList<String> tmp = new ArrayList<String>();
					
					tmp.add(current.get(0));
					tmp.add(current.get(1));
					tmp.add(current.get(2));
					tmp.add(current.get(3));
					
					bunch.add(tmp);
					
					// clear current
					current.clear();
					
					astr = "";
					cstr = "";
					gstr = "";
					tstr = "";
					
					// begin new
					count = pos;
														
					float a = (rs.getFloat("a"));
					float c = (rs.getFloat("c"));
					float g = (rs.getFloat("g"));
					float t = (rs.getFloat("t"));
					
					float sum = a+c+g+t;
					
					a/=sum;
					c/=sum;
					g/=sum;
					t/=sum;
					
					astr = astr.concat(a + "   ");
					cstr = cstr.concat(c + "   ");
					gstr = gstr.concat(g + "   ");
					tstr = tstr.concat(t + "   ");
				}
				
			}

		} catch (SQLException e) {
			e.printStackTrace();
		}
		
		// save last
		current.add(0, astr.trim());
		current.add(1, cstr.trim());
		current.add(2, gstr.trim());
		current.add(3, tstr.trim());
		
		bunch.add(current);
				
		return bunch;
	}

	
	// process rs for domain-data
	public ArrayList<String> processDom(ResultSet rs){
		
		// move cursor before first row
		try {
			rs.beforeFirst();
		} catch (SQLException e1) {
			e1.printStackTrace();
		}
		
		ArrayList<String> sqlArray = new ArrayList<String>();
		
		try {
			while (rs.next()) {
				String tmp = (rs.getString("name"));
				
				// only add new domains
				if (sqlArray.contains(tmp) != true) {
					sqlArray.add(tmp);
				}
			}
		} catch (SQLException e) {
			e.printStackTrace();
		}
		
		return sqlArray;
	}
	
}
