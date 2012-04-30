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
	import java.awt.Component;
	import java.awt.Container;
	import java.awt.GridBagConstraints;
	import java.awt.GridBagLayout;

	/** Stellt Methoden zur Verfuegung, die die Arbeit mit
	  * einem LayoutManager vereinfachen.
	  * 
	  * @author Andreas Draeger
	  */
	public class LayoutHelper {

	  /** Methode, die uns beim Anornden von Elementen im GridBagLayout hilft.
	   * @param cont
	   * @param gbl
	   * @param c
	   * @param x
	   * @param y
	   * @param width
	   * @param height
	   * @param weightx
	   * @param weighty
	   */
	 public static void addComponent(
	         Container cont,
	         GridBagLayout gbl,
	         Component c,
	         int x, int y,
	         int width, int height,
	         double weightx, double weighty )
	 {
	   GridBagConstraints gbc = new GridBagConstraints();
	   gbc.fill               = GridBagConstraints.BOTH;
	   gbc.gridx              = x; 
	   gbc.gridy              = y;
	   gbc.gridwidth          = width; 
	   gbc.gridheight         = height;
	   gbc.weightx            = weightx; 
	   gbc.weighty            = weighty;
	   gbl.setConstraints(c, gbc);
	   cont.add(c);
	 }

}

