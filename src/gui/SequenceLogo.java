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
import java.io.File;
import java.util.ArrayList;
import java.util.StringTokenizer;

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

import extension.MatrixFileParser;


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
    
    public static void printLogo(String tf_name, ArrayList<String> tf_matrix) {
    	
    	if (tf_matrix == null) {
    		MatrixFileParser matrix_parser = new MatrixFileParser();
    		matrix_parser.readMatrices("/rahome/eichner/Desktop/predicted_pfms.txt");
    		//matrix_parser.readMatrices("/rahome/eichner/projects/sabine/src/SABINE_1.0/trainingsets_public/trainingset_all_classes.fbps");
    		tf_matrix = matrix_parser.obtainMatrix(tf_name);
    	}
    	
    	DistributionLogo[] logoData = new DistributionLogo[tf_matrix.size()];
    	double A, C, G, T;
    	double sum;
    	StringTokenizer strtok;
    	
    	for (int i=0; i<tf_matrix.size(); i++) {
    		
    		strtok = new StringTokenizer(tf_matrix.get(i));
    		
    		strtok.nextToken();								// Index
    		A = Double.parseDouble(strtok.nextToken());
    		C = Double.parseDouble(strtok.nextToken());
    		G = Double.parseDouble(strtok.nextToken());
    		T = Double.parseDouble(strtok.nextToken());

    		sum = A + C + G + T;
    		
    		A = A / sum;
    		C = C / sum;
    		G = G / sum;
    		T = T / sum;
    		
    		logoData[i] = getLogo(new double[] {A, C, G, T});
    	}

    	try {
    		javax.imageio.ImageIO.write(SequenceLogo.drawSequenceLogo(logoData, 100), "png", new File("/rahome/eichner/Desktop/logos/" + tf_name + "_predicted_logo.png"));
			
			
			
		} catch (Exception e) {
			e.printStackTrace();
		}
    }
    
    public static void main(String[] args) {
    	
    	String[] relevant_tfs =  new String[] {"T00140", "T00949", "T01517", "T01599", "T02381", "T03428", "T03548", "T03552", "T03911", "T04372", "T08313", "T09328", "T09645", "T10291"};
		
    	
    	for (int i=0; i<relevant_tfs.length; i++) {
    		printLogo(relevant_tfs[i], null);
    	}
    	/*
    	printLogo("T00505", null);
    	printLogo("T01006", null);
    	printLogo("T01005", null);
    	printLogo("T01009", null);
    	printLogo("T01783", null);
    	*/
    	

    }
}


