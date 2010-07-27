package util.parsers.matbase;
/*
 * ===============================================
 * (C) Florian Topf, University of Tuebingen, 2010
 * ===============================================
 */

import java.sql.*;

/*
 * MySQL Database API
 */

public class MySQLAPI{
	
	// establish connection
	public Connection connect(String sDbUrl, String sUsr, String sPwd) {
		
		Connection cn = null;
	    		
		try {
			Class.forName("com.mysql.jdbc.Driver");
			System.out.println ("Connecting to " + sDbUrl);
			cn = DriverManager.getConnection(sDbUrl, sUsr, sPwd);
			System.out.println ("Database connection established.");
		
		} catch (Exception e) {
			e.printStackTrace();
			System.exit(1);
		}
	
	return cn;

	}
		
	// handle query
	public ResultSet query(Connection cn, String sQuery) {
		
		Statement  st = null;
	    ResultSet  rs = null;
	    
	    try {
			st = cn.createStatement();
			    
			System.out.println("Retrieving available data...");
		
			rs = st.executeQuery(sQuery);
		
	    } catch (SQLException e) {
			e.printStackTrace();
		}
		
		return rs;
	}
	
}