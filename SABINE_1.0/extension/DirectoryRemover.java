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

package extension;

import help.FileFilter;

import java.io.File;

public class DirectoryRemover {
	
	
	public static boolean removeDirectory(String directory) {
		
		boolean abs_path = false;
		if (directory.startsWith("/")) {
			abs_path = true;
		}
		
		File dir = new File(directory);
		
		if (dir.exists()) {
			
			FileFilter fileFilter = new FileFilter();
			fileFilter.setFormat(".*");
			fileFilter.setDirectory(directory);
			
			String[] files;
			if (abs_path) 
				files = fileFilter.listFilesFullPath();
			else 
				files = fileFilter.listFiles();
			
			
			for (int i=0; i<files.length; i++) {
				
				// directory
				if (new File(directory + files[i]).isDirectory()) {
					
					String sub_directory = directory + files[i] +"/";
					
					fileFilter.setDirectory(sub_directory);
					
					String[] contained_files;
					if (abs_path)
						contained_files = fileFilter.listFilesFullPath();
					else
						contained_files = fileFilter.listFiles();
					
					// empty directory
					if (contained_files.length == 0) {
						
						new File(sub_directory).delete();
						
						if ( new File(sub_directory).exists()) {
							System.out.println("Empty Directory : " + files[i] + " not deleted.");
						}
					}
					
					// full directory
					else {
						removeDirectory(sub_directory);
						new File(sub_directory).delete();
						
						if ( new File(sub_directory).exists()) {
							System.out.println("Full Directory  : " + files[i] + " not deleted.");
						}
					}
				}
				
				// file
				else {
					new File(directory + files[i]).delete();
					if ( new File(directory + files[i]).exists()) {
						System.out.println("File            : " + files[i] + " not deleted.");
					}
					
				}
			}
			
			fileFilter.setDirectory(directory);
			
			if (abs_path)
				return fileFilter.listFilesFullPath().length == 0;
			else
				return fileFilter.listFiles().length == 0;
		}
		else {
			System.out.println("Directory \"" + directory + "\"not found.");
			return false;
		}
	}
	
	public static void main(String[] args) {
		
		boolean succ = removeDirectory("tmp/");
		
		if (succ) {
			System.out.println("Files removed successfully");
		}
		else {
			System.out.println("Files were not removed");
		}
	}

}

