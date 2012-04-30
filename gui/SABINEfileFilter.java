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

package gui;
import java.io.File;
import javax.swing.filechooser.FileFilter;

/**
 * 
 */

/**
 * @author <a href="mailto:johannes.eichner@uni-tuebingen.de">Johannes Eichner</a>
 * 
 */
public class SABINEfileFilter extends FileFilter implements java.io.FileFilter {

	private boolean txt, png;

	public static final short PNG = 0;
	public static final short TXT = 1;
	
	/**
	 * 
	 */
	public SABINEfileFilter(short s) {
		
		if (s == PNG) {
			png = true;
			txt = false;
		}
		else {
			png = false;
			txt = true;
		}
	}
	
	public boolean isImageFilter() {
		return png;
	}
	public boolean isTextFilter() {
		return txt;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see javax.swing.filechooser.FileFilter#accept(java.io.File)
	 */
	@Override
	public boolean accept(File f) {
		return f.isDirectory()
				|| (f.getName().toLowerCase().endsWith(".txt") && txt)
				|| (f.getName().toLowerCase().endsWith(".png") && png);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see javax.swing.filechooser.FileFilter#getDescription()
	 */
	@Override
	public String getDescription() {
		StringBuffer descr = new StringBuffer();
		if (png) {
			descr.append("Image files (*.png)\n");
		}
		if (txt)
			descr.append("Text files (*.txt)\n");
		return descr.toString();
	}

}

