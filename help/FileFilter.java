/*
    SABINE predicts binding specificities of transcription factors.
    Copyright (C) 2009 ZBIT, University of Tübingen, Johannes Eichner

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

package help;

import java.io.File;
import java.io.FilenameFilter;

public class FileFilter {

	// regular expressions to match ncbi files	
	private String format = "none";
	private String directory = "none";
	
	public String getFormat() {
		return format;
	}
	
	public void setFormat(String f) {
		format = f;
	}
	
	public String getDirectory() {
		return directory;
	}
	
	public void setDirectory(String d) {
		directory = d;
	}
	
	private class FormatFileFilter implements FilenameFilter {
		public boolean accept(File dir, String name) {
			return name.matches(format);
		}
	}
	
	// public function which returns lists of names of files
	public String[] listFiles() {
		
		if(format.equals("none")) {
			System.out.println("Please specify file format. (regular expression)");
			return new String[] {"crap"};
		}
		
		if(directory.equals("none")) {
			System.out.println("Please specify path to directory.");
			return new String[] {"crap"};
		}
		
	  	File path = new File(System.getProperty("user.dir") + File.separator + directory);
	  	return path.list(this.new FormatFileFilter());
	}
	
	public String[] listFilesFullPath() {
		
		if(format.equals("none")) {
			System.out.println("Please specify file format. (regular expression)");
			return new String[] {"crap"};
		}
		
		if(directory.equals("none")) {
			System.out.println("Please specify path to directory.");
			return new String[] {"crap"};
		}
		
	  	File path = new File(directory);
	  	return path.list(this.new FormatFileFilter());
	}
	
}

