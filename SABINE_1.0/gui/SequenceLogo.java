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

import java.awt.Color;
import java.awt.Dimension;
import java.awt.Graphics2D;
import java.awt.Paint;
import java.awt.RenderingHints;
import java.awt.image.BufferedImage;

import javax.swing.Box;
import javax.swing.BoxLayout;
import javax.swing.JComponent;
import javax.swing.JPanel;

import org.biojava.bio.dist.Distribution;
import org.biojava.bio.dist.SimpleDistribution;
import org.biojava.bio.gui.DNAStyle;
import org.biojava.bio.gui.DistributionLogo;
import org.biojava.bio.gui.TextLogoPainter;
import org.biojava.bio.seq.DNATools;
import org.biojava.bio.symbol.Symbol;


/** Utility for producing a sequence logo from a distribution
*/
public class SequenceLogo {

    /** Get an array of these to feed into the drawSequenceLogo function. */
    public static DistributionLogo getLogo( double[] acgt ) {
        return getLogo( makeDistribution( acgt ) );
    }

    public static Distribution makeDistribution( double[] acgt ) {
        Symbol[] syms = { DNATools.a(), DNATools.c(), DNATools.g(), DNATools.t() };
        Distribution dist = new SimpleDistribution( DNATools.getDNA() ); 
        try {
            for (int i = 0; i < 4; ++i) {
                dist.setWeight( syms[i], acgt[i] );
            }
        } catch (Exception biojavaCheckedException) {
            throw new RuntimeException( biojavaCheckedException );
        }
        return dist;
    }

    public static DistributionLogo getLogo( Distribution dist ) {
        DistributionLogo dl = new DistributionLogo();
        try {
            dl.setDistribution( dist );
        } catch (Exception biojavaCheckedException) {
            throw new RuntimeException( biojavaCheckedException );
        }
        return dl;
    }

    private static JComponent makeBlank( Dimension size ) {
        JPanel panel = new JPanel();
        panel.setPreferredSize( size );
        panel.setSize( size );
        panel.setBackground( Color.white );
        return panel;
    }

    /** Draw a sequence logo bitmap image to a file. */
    public static BufferedImage drawSequenceLogo( DistributionLogo[] logos, int scale) throws Exception {
        
    	// set number of blank positions at begin and end of logo
    	int leadingBlanks = 0;
    	int trailingBlanks = 0;
    	
    	int imageWidth = scale * (logos.length + leadingBlanks + trailingBlanks);
        int imageHeight = scale * 2;
        
        BufferedImage image = new BufferedImage( imageWidth, imageHeight, BufferedImage.TYPE_INT_RGB ); 
        Graphics2D gr = (Graphics2D)image.getGraphics(); 
        gr.setPaint( Color.black );
        
        Box hBox = new Box( BoxLayout.X_AXIS );
        Dimension imageSize = new Dimension( imageWidth, imageHeight );
        hBox.setSize( imageSize );
        hBox.setPreferredSize( imageSize );
        
        Dimension logoSize = new Dimension( scale, imageHeight );
        DNAStyle style = new DNAStyle();
        Symbol[] syms = { DNATools.a(), DNATools.c(), DNATools.g(), DNATools.t() };

        RenderingHints hints = new RenderingHints(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
        
        for (Symbol sym : syms) {
            // remove black outline
            Paint paint = style.fillPaint( sym );
            if (Color.yellow.equals( paint )) {
                //paint = Color.yellow.darker();
                paint = new Color( 250, 240, 0 );
            }
            style.setOutlinePaint( sym, paint );
        }
        
        for (int i = 0; i < leadingBlanks; ++i) {
            hBox.add( makeBlank( logoSize ) );
        }
        for (int i = 0; i < logos.length; ++i) {
            logos[i].setStyle( style );
            logos[i].setPreferredSize( logoSize );
            logos[i].setBackground( Color.pink );
            logos[i].setRenderingHints( hints );
            logos[i].setLogoPainter( new TextLogoPainter() );
            hBox.add( logos[i] );
        }
        for (int i = 0; i < trailingBlanks; ++i) {
            hBox.add( makeBlank( logoSize ) );
        }
        gr.setPaint( Color.white );
        gr.fillRect( 0, 0, imageWidth, imageHeight );
        hBox.setDoubleBuffered( false );
        hBox.addNotify();
        hBox.validate();
        hBox.print( gr );

        return image;
    }
    
    public static void main(String[] args) {
    	
    	int logoSize = 6;
    	double A, C, G, T;
    	DistributionLogo[] logoData = new DistributionLogo[logoSize];
    	for (int i=0; i<logoSize; i++) {
    		
    		A = 0; //Math.random();
    		C = 0.05; //Math.random();
    		G = 0.95; //Math.random();
    		T = 0; //Math.random();
    		
    		/*
    		sum = A + C + G + T;
    		
    		A = A / sum;
    		C = C / sum;
    		G = G / sum;
    		T = T / sum;
    		*/
    		logoData[i] = getLogo(new double[] {A, C, G, T});
    	}
    	
    	try {
			drawSequenceLogo(logoData, 100);
		} catch (Exception e) {
			e.printStackTrace();
		}
    }
}


