package util.parsers.uniprobe;
/*
 * ===============================================
 * (C) Florian Topf, University of Tuebingen, 2010
 * ===============================================
 */

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Enumeration;
import java.util.zip.ZipEntry;
import java.util.zip.ZipException;
import java.util.zip.ZipFile;


public class ZipHelper {

	ZipFile zipFile = null;

	// return array with zip-entries
	public ArrayList<ZipEntry> getZipEntries(File in) {
		
		// try open archive
		try {
			zipFile = new ZipFile(in);
		} catch (ZipException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		ArrayList<ZipEntry> entries = new ArrayList<ZipEntry>();
		
		Enumeration enu = zipFile.entries();
				

		while (enu.hasMoreElements()) {
	
			ZipEntry ze = (ZipEntry) enu.nextElement();
			
			// filter directories and readme.txt
			if (ze.isDirectory() == false && ze.getName().contains("readme") == false) {
				
				entries.add(ze);
				
			}
	
		}
		
		return entries;
	}
	
}
