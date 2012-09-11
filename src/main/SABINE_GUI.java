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

package main;

import gui.JHelpBrowser;
import gui.LayoutHelper;
import gui.MessageProcessor;
import gui.ReversedSpinnerListModel;
import gui.SABINEfileFilter;
import gui.SequenceLogo;
import help.RandomStringGenerator;
import help.TimeStampGenerator;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.GridBagLayout;
import java.awt.GridLayout;
import java.awt.Image;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.WindowEvent;
import java.awt.event.WindowListener;
import java.awt.image.BufferedImage;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.OutputStream;
import java.io.PrintStream;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Properties;
import java.util.StringTokenizer;
import java.util.Vector;

import javax.swing.BorderFactory;
import javax.swing.ImageIcon;
import javax.swing.JButton;
import javax.swing.JComboBox;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JList;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JProgressBar;
import javax.swing.JScrollPane;
import javax.swing.JSpinner;
import javax.swing.JTable;
import javax.swing.JTextArea;
import javax.swing.SpinnerModel;
import javax.swing.SpinnerNumberModel;
import javax.swing.ToolTipManager;
import javax.swing.UIManager;
import javax.swing.UIManager.LookAndFeelInfo;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;
import javax.swing.event.DocumentEvent;
import javax.swing.event.DocumentListener;
import javax.swing.event.ListSelectionEvent;
import javax.swing.event.ListSelectionListener;

import org.biojava.bio.gui.DistributionLogo;
import org.biojava.bio.gui.sequence.GUITools;
import org.biojava.bio.seq.ProteinTools;
import org.biojava.bio.symbol.IllegalSymbolException;

import resources.Resource;
import de.zbit.util.progressbar.gui.ProgressBarSwing;
import extension.DirectoryRemover;

public class SABINE_GUI extends JFrame implements ActionListener, ChangeListener,
		DocumentListener, ListSelectionListener, GUIListener, WindowListener {

	/**
	 * Generated serial version id for text storage.
	 */
	private static final long serialVersionUID = -2345389446488533671L;

	// private JTextField name_field;
	private JPanel mainPanel;
	
	private JTextArea sequenceArea;

	private JComboBox superClassSelect;

	private JComboBox organismSelect;

	private JSpinner domLeftSpin;

	private JSpinner domRightSpin;

	private Vector<String> domList;

	private JList domSelected;

	private JSpinner bestMatchThreshold;

	private JSpinner maxNumBestMatches;

	private JSpinner outlierFilterThreshold;

	private JTextArea output;
	
	private JScrollPane outputScroller;
	
	private JProgressBar progressBar;
  
  private JPanel outputPanel;

	private MessageProcessor msg;

	private Thread currRun;

	// Buttons
	private JButton runButton;

	private JButton resetButton;
	
	private JButton demoButton;

	private JButton addButton;

	private JButton deleteButton;
	
	private JButton helpButton;
	
	private JButton enlargeButton;
	
	private DistributionLogo[] logoArray;
	
	private Object[][] best_matches;
	
	private Object[][] pred_matrix;
	
	private ArrayList<String> predicted_pfm;

	private JPanel resultsPanel;

	private JTable bestMatchTable;
	
	private int title_width;
	
	boolean seq_added = false;

	boolean dom_added = false;
	
	private String base_dir;

	final String[] superclasses = new String[] { "Basic Domain", "Zinc Finger",
			"Helix-Turn-Helix", "Beta Scaffold", "Other" };

	private JTable predMatrixTable;

	private JButton matrixSaveButton;

	private JButton bestMatchSaveButton;

	public double high_conf_bmt = 0.0;
	public double medium_conf_bmt = 0.0;
	public double low_conf_bmt = 0.0;
	
	static {
		initLaF(SABINE_Main.appName);
	}
	
  /**
   * 
   * @return
   */
  public static boolean isMacOSX() {
    return (System.getProperty("mrj.version") != null)
        || (System.getProperty("os.name").toLowerCase().indexOf("mac") != -1);
  }
  
  /**
   * Initializes the look and feel.
   */
  public static void initLaF() {
    // 15 s for tooltips to be displayed
    ToolTipManager.sharedInstance().setDismissDelay(15000);
    try { 
      UIManager.setLookAndFeel(new javax.swing.plaf.metal.MetalLookAndFeel());
      String osName = System.getProperty("os.name");
      if (osName.equals("Linux") || osName.equals("FreeBSD")) {
        UIManager.setLookAndFeel("com.sun.java.swing.plaf.gtk.GTKLookAndFeel");
        // UIManager.setLookAndFeel("com.sun.java.swing.plaf.motif.MotifLookAndFeel");
      } else if (isMacOSX()) {
        UIManager.setLookAndFeel("ch.randelshofer.quaqua.QuaquaLookAndFeel");
      } else if (osName.contains("Windows")) {
        UIManager.setLookAndFeel("com.sun.java.swing.plaf.windows.WindowsLookAndFeel");
      } else {
        UIManager.setLookAndFeel(UIManager.getSystemLookAndFeelClassName());
      }
    } catch (Throwable e) {
      try {
        UIManager.setLookAndFeel(UIManager.getSystemLookAndFeelClassName());
      } catch (Throwable e1) {
        // If Nimbus is not available, you can set the GUI to another look
        // and feel.
        // Native look and feel for Windows, MacOS X. GTK look and
        // feel for Linux, FreeBSD
        try {
          for (LookAndFeelInfo info : UIManager.getInstalledLookAndFeels()) {
            if ("Nimbus".equals(info.getName())) {
              UIManager.setLookAndFeel(info.getClassName());
              break;
            }
          }
        } catch (Throwable exc) {
          // ignore.
        }
      }
    }
  }
  
  /**
   * Initializes the look and feel. This method is only useful when the calling
   * class contains the main method of your program and is also derived from this
   * class ({@link GUITools}).
   * 
   * @param title
   *        the title to be displayed in the xDock in case of MacOS. Note that
   *        this title can only be displayed if this method belongs to the class
   *        that has the main method (or is in a super class of it), i.e., in
   *        order to make use of this method in a proper way, you have to extend
   *        this {@link GUITools} and to put the main method into your derived
   *        class. From this main method you then have to call this method.
   *        Otherwise, this title will not have any effect.
   */
  public static void initLaF(String title) {
    if (isMacOSX()) {
      Properties p = System.getProperties();
      // also use -Xdock:name="Some title" -Xdock:icon=path/to/icon on command line
      p.setProperty("apple.awt.graphics.EnableQ2DX", "true");
      p.setProperty("apple.laf.useScreenMenuBar", "true");
      p.setProperty("com.apple.macos.smallTabs", "true");
      p.setProperty("com.apple.macos.useScreenMenuBar", "true");
      if ((title != null) && (title.length() > 0)) {
        p.setProperty("com.apple.mrj.application.apple.menu.about.name", title);
      }
      p.setProperty("com.apple.mrj.application.growbox.intrudes", "false");
      p.setProperty("com.apple.mrj.application.live-resize", "true");
    }
    initLaF();
  }
	

	public SABINE_GUI() {
		super(SABINE_Main.appName);

		/*
		 * adapt JFrame to operating system
		 */

		mainPanel = new JPanel();
		GridBagLayout gbl = new GridBagLayout();
		mainPanel.setLayout(gbl);

		// read species list
		ArrayList<String> species_list = new ArrayList<String>();
		
		try {
			BufferedReader br = new BufferedReader(new InputStreamReader(Resource.class.getResourceAsStream("txt/organism_list.txt")));

			String line;
			
			while ((line = br.readLine()) != null) {
				species_list.add(line);
			}
		}
		catch(IOException ioe) {
			ioe.printStackTrace();
		}
		String[] species = species_list.toArray(new String[] {});

		
		/*
		 *  SABINE Logo
		 */
		
		JPanel titlePanel = new JPanel();
		try {
			
			double scaling_factor = 0.7;
			BufferedImage logoImage = javax.imageio.ImageIO.read(Resource.class.getResource("img/sabine_logo.png"));
			Image titleImage = logoImage.getScaledInstance((int) (logoImage.getWidth()* scaling_factor), (int) (logoImage.getHeight()*scaling_factor), Image.SCALE_SMOOTH);
			JLabel titleLabel = new JLabel(new ImageIcon(titleImage));
			title_width = (int) (logoImage.getWidth()* scaling_factor);
			
			titlePanel.add(titleLabel);
			titlePanel.setBackground(Color.white);
			
			ImageIcon imageIcon = new ImageIcon(Resource.class.getResource("img/sabine_icon.png"));
			this.setIconImage(imageIcon.getImage());
			
		} catch (IOException e) {
			e.printStackTrace();
		}

		
		/*
		 *  Organism and Superclass
		 */
		
		// add combo box for species selection
		organismSelect = new JComboBox(species);
		organismSelect.setSelectedItem("Homo Sapiens");
		organismSelect.setBackground(Color.white);

		// add combo box for superclass selection
		superClassSelect = new JComboBox(superclasses);
		superClassSelect.addItem("Auto-detect");
		superClassSelect.setBackground(Color.white);
		
		JPanel firstPanel = new JPanel();
		firstPanel.setBorder(BorderFactory.createTitledBorder(" Select Organism and Superclass "));
		firstPanel.setLayout(new GridBagLayout());
		
		LayoutHelper.addComponent(firstPanel, (GridBagLayout) firstPanel.getLayout(), new JPanel(), 0, 0, 1, 1, 0, 0);
		LayoutHelper.addComponent(firstPanel, (GridBagLayout) firstPanel.getLayout(), new JLabel("Organism"), 1, 0, 1, 1, 0, 0);
		LayoutHelper.addComponent(firstPanel, (GridBagLayout) firstPanel.getLayout(), new JPanel(), 2, 0, 1, 1, 0, 0);
		LayoutHelper.addComponent(firstPanel, (GridBagLayout) firstPanel.getLayout(), organismSelect, 3, 0, 1, 1, 0, 0);
		LayoutHelper.addComponent(firstPanel, (GridBagLayout) firstPanel.getLayout(), new JPanel(), 4, 0, 1, 1, 0, 0);
		
		LayoutHelper.addComponent(firstPanel, (GridBagLayout) firstPanel.getLayout(), new JPanel(), 0, 1, 2, 2, 0, 0);
		
		LayoutHelper.addComponent(firstPanel, (GridBagLayout) firstPanel.getLayout(), new JLabel("Superclass"), 1, 3, 1, 1, 0, 0);
		LayoutHelper.addComponent(firstPanel, (GridBagLayout) firstPanel.getLayout(), superClassSelect, 3, 3, 1, 1, 0, 0);
		
		//LayoutHelper.addComponent(firstPanel, (GridBagLayout) firstPanel.getLayout(), new JPanel(), 0, 4, 2, 2, 0, 0);
		
		LayoutHelper.addComponent(mainPanel, gbl, firstPanel, 0, 0, 1, 1, 0, 0);
		
		
		/*
		 *  Protein sequence
		 */
		
		// add text panel for protein sequence
		JPanel seqPanel = new JPanel();
		seqPanel.setLayout(new GridBagLayout());
		seqPanel.setBorder(BorderFactory.createTitledBorder(" Protein sequence "));
		
		sequenceArea = new JTextArea(1, 20);
		sequenceArea.getDocument().putProperty("name", "seq");
		sequenceArea.getDocument().addDocumentListener(this);
		
		
		JScrollPane scroller = new JScrollPane(sequenceArea,
											   JScrollPane.VERTICAL_SCROLLBAR_AS_NEEDED,
											   JScrollPane.HORIZONTAL_SCROLLBAR_AS_NEEDED);
		scroller.setBorder(BorderFactory.createLoweredBevelBorder());
		scroller.setPreferredSize(new Dimension(200, 100));
		
		LayoutHelper.addComponent(seqPanel, (GridBagLayout) seqPanel.getLayout(), new JPanel(), 0, 0, 1, 1, 0, 0);
		LayoutHelper.addComponent(seqPanel, (GridBagLayout) seqPanel.getLayout(), scroller, 1, 0, 1, 1, 1, 0);
		LayoutHelper.addComponent(seqPanel, (GridBagLayout) seqPanel.getLayout(), new JPanel(), 2, 0, 1, 1, 0, 0);
		
		LayoutHelper.addComponent(mainPanel, gbl, seqPanel, 0, 1, 1, 1, 0, 0);
		
		
		/*
		 *  DNA-binding domains
		 */
		
		domLeftSpin = new JSpinner(new SpinnerNumberModel(1, 1, 10000, 1));
		domLeftSpin.setName("domLeftSpin");
		domLeftSpin.addChangeListener(this);
		domLeftSpin.setEnabled(false);

		domRightSpin = new JSpinner(new SpinnerNumberModel(1, 1, 10000, 1));
		domRightSpin.setName("domRightSpin");
		domRightSpin.addChangeListener(this);
		domRightSpin.setEnabled(false);

		addButton = new JButton("  Add ");
		addButton.setName("add");
		addButton.addActionListener(this); 
		addButton.setEnabled(false);

		deleteButton = new JButton("Delete");
		deleteButton.setName("delete");
		deleteButton.addActionListener(this);
		deleteButton.setEnabled(false);
		
		domList = new Vector<String>();
		domSelected = new JList(domList);
		domSelected.addListSelectionListener(this);
		JScrollPane listScroller = new JScrollPane(domSelected);
		listScroller.setPreferredSize(new Dimension(200, 70));
		
		JPanel buttonPanel = new JPanel();
		buttonPanel.setLayout(new GridBagLayout());
		LayoutHelper.addComponent(buttonPanel, (GridBagLayout) buttonPanel.getLayout(), addButton, 0, 0, 1, 1, 0, 0);
		LayoutHelper.addComponent(buttonPanel, (GridBagLayout) buttonPanel.getLayout(), new JPanel(), 0, 1, 1, 1, 0, 0.1);
		LayoutHelper.addComponent(buttonPanel, (GridBagLayout) buttonPanel.getLayout(), deleteButton, 0, 2, 1, 1, 0, 0);
		LayoutHelper.addComponent(buttonPanel, (GridBagLayout) buttonPanel.getLayout(), new JPanel(), 0, 3, 1, 1, 0, 0.1);
		
		JPanel domainsPanel = new JPanel();
		domainsPanel.setBorder(BorderFactory.createTitledBorder(" Domains "));
		GridBagLayout domLayout = new GridBagLayout();
		domainsPanel.setLayout(domLayout);
		
		
		LayoutHelper.addComponent(domainsPanel, domLayout, new JLabel("       Start"), 0, 0, 1, 1, 0, 0);
		LayoutHelper.addComponent(domainsPanel, domLayout, new JPanel(), 1, 0, 1, 1, 0, 0);
		LayoutHelper.addComponent(domainsPanel, domLayout, domLeftSpin, 2, 0, 1, 1, 0, 0);
		LayoutHelper.addComponent(domainsPanel, domLayout, new JPanel(), 3, 0, 1, 1, 0, 0);
		LayoutHelper.addComponent(domainsPanel, domLayout, new JLabel("End"), 4, 0, 1, 1, 0, 0);
		LayoutHelper.addComponent(domainsPanel, domLayout, new JPanel(), 5, 0, 1, 1, 0, 0);
		LayoutHelper.addComponent(domainsPanel, domLayout, domRightSpin, 6, 0, 1, 1, 0, 0);
		LayoutHelper.addComponent(domainsPanel, domLayout, new JPanel(), 7, 0, 1, 1, 0, 0);
		
		LayoutHelper.addComponent(domainsPanel, domLayout, new JPanel(), 0, 1, 8, 1, 0, 0);
		LayoutHelper.addComponent(domainsPanel, domLayout, new JLabel("<html>Selected<br>domains</html>"), 0, 2, 1, 1, 0, 0);
		LayoutHelper.addComponent(domainsPanel, domLayout, listScroller, 2, 2, 5, 1, 0, 0);
		LayoutHelper.addComponent(domainsPanel, domLayout, buttonPanel, 8, 0, 1, 3, 0, 0);


		LayoutHelper.addComponent(mainPanel, gbl, domainsPanel, 0, 2, 1, 1, 0, 0);

		
		/*
		 * Parameters
		 */
		
		JPanel parameterPanel = new JPanel();
		parameterPanel.setBorder(BorderFactory.createTitledBorder(" Parameters "));
		parameterPanel.setLayout(new GridBagLayout());
		
		ArrayList<String> bmt_list = new ArrayList<String>();
		bmt_list.add("Dynamic");
		for (int l=95; l>=0; l=l-5) {
			bmt_list.add("0," + l);
		}
		
		bestMatchThreshold = new JSpinner(new ReversedSpinnerListModel(bmt_list.toArray()));
		maxNumBestMatches = new JSpinner(new SpinnerNumberModel(5, 1, 5, 1));
		outlierFilterThreshold = new JSpinner(new SpinnerNumberModel(0.5, 0.1, 1, 0.1));

		// add JSpinner for "Best Match Threshold"
		LayoutHelper.addComponent(parameterPanel, (GridBagLayout) parameterPanel.getLayout(), 
								  new JLabel("Best match threshold"), 0, 0, 1, 1, 0, 0);
		
		LayoutHelper.addComponent(parameterPanel, (GridBagLayout) parameterPanel.getLayout(), new JPanel(), 1, 0, 1, 1, 0, 0);
		LayoutHelper.addComponent(parameterPanel, (GridBagLayout) parameterPanel.getLayout(), bestMatchThreshold, 2, 0, 1, 1, 0, 0);
		LayoutHelper.addComponent(parameterPanel, (GridBagLayout) parameterPanel.getLayout(), new JPanel(), 0, 1, 2, 1, 0, 0);
	

		
		// add JSpinner for "Max. Number of Best Matches"
		LayoutHelper.addComponent(parameterPanel, (GridBagLayout) parameterPanel.getLayout(), 
								  new JLabel("Maximum number of best matches"), 0, 2, 1, 1, 0, 0);
		
		LayoutHelper.addComponent(parameterPanel, (GridBagLayout) parameterPanel.getLayout(), new JPanel(), 1, 2, 1, 1, 0, 0);
		LayoutHelper.addComponent(parameterPanel, (GridBagLayout) parameterPanel.getLayout(), maxNumBestMatches, 2, 2, 1, 1, 0, 0);
		LayoutHelper.addComponent(parameterPanel, (GridBagLayout) parameterPanel.getLayout(), new JPanel(), 0, 3, 2, 1, 0, 0);
		
		// add JSpinner for "Outlier filter threshold"
		LayoutHelper.addComponent(parameterPanel, (GridBagLayout) parameterPanel.getLayout(), 
								  new JLabel("Outlier filter threshold"), 0, 4, 1, 1, 0, 0);
		
		LayoutHelper.addComponent(parameterPanel, (GridBagLayout) parameterPanel.getLayout(), new JPanel(), 1, 4, 1, 1, 0, 0);
		LayoutHelper.addComponent(parameterPanel, (GridBagLayout) parameterPanel.getLayout(), outlierFilterThreshold, 2, 4, 1, 1, 0, 0);
		LayoutHelper.addComponent(parameterPanel, (GridBagLayout) parameterPanel.getLayout(), new JPanel(), 0, 5, 2, 1, 0, 0);
		

		LayoutHelper.addComponent(mainPanel, gbl, parameterPanel, 0, 3, 1, 1, 0, 0);
		
		
		/*
		 *  Run-Button and Stop-Button
		 */

		runButton = new JButton(" Run ");
		runButton.setEnabled(false);
		runButton.setName("run");
		runButton.addActionListener(this);
		
		resetButton = new JButton("Reset");
		resetButton.setName("reset");
		resetButton.addActionListener(this);
		
		demoButton = new JButton("Demo");
		demoButton.setName("demo");
		demoButton.addActionListener(this);
		
		helpButton = new JButton("Help");
		helpButton.setName("help");
		helpButton.addActionListener(this);
		
		JPanel runPanel = new JPanel();
		runPanel.setLayout(new GridBagLayout());
		
		LayoutHelper.addComponent(runPanel, (GridBagLayout) runPanel.getLayout(), new JPanel(), 0, 0, 5, 1, 0, 0);
		LayoutHelper.addComponent(runPanel, (GridBagLayout) runPanel.getLayout(), new JPanel(), 0, 1, 1, 1, 0, 0);
		LayoutHelper.addComponent(runPanel, (GridBagLayout) runPanel.getLayout(), runButton, 1, 1, 1, 1, 0, 0);
		LayoutHelper.addComponent(runPanel, (GridBagLayout) runPanel.getLayout(), new JPanel(), 2, 1, 1, 1, 0, 0);
		LayoutHelper.addComponent(runPanel, (GridBagLayout) runPanel.getLayout(), resetButton, 3, 1, 1, 1, 0, 0);
		LayoutHelper.addComponent(runPanel, (GridBagLayout) runPanel.getLayout(), new JPanel(), 4, 1, 1, 1, 0, 0);
		LayoutHelper.addComponent(runPanel, (GridBagLayout) runPanel.getLayout(), demoButton, 5, 1, 1, 1, 0, 0);
		LayoutHelper.addComponent(runPanel, (GridBagLayout) runPanel.getLayout(), new JPanel(), 6, 1, 1, 1, 0, 0);
		LayoutHelper.addComponent(runPanel, (GridBagLayout) runPanel.getLayout(), helpButton, 7, 1, 1, 1, 0, 0);
		LayoutHelper.addComponent(runPanel, (GridBagLayout) runPanel.getLayout(), new JPanel(), 8, 1, 1, 1, 0, 0);
		LayoutHelper.addComponent(runPanel, (GridBagLayout) runPanel.getLayout(), new JPanel(), 0, 2, 9, 1, 0, 0);
		
		LayoutHelper.addComponent(mainPanel, gbl, runPanel, 0, 4, 1, 1, 0, 0);
		

		/*
		 *  Text area for output
		 */

		output = new JTextArea(25, 75);
		output.getDocument().putProperty("name", "output");
		output.getDocument().addDocumentListener(this);


		outputScroller = new JScrollPane(output,JScrollPane.VERTICAL_SCROLLBAR_AS_NEEDED,
				JScrollPane.HORIZONTAL_SCROLLBAR_AS_NEEDED);

		progressBar = new JProgressBar();
		//progressBar.setMinimumSize(new java.awt.Dimension(100,20));

		outputPanel = new JPanel(new BorderLayout());
		outputPanel.add(progressBar, BorderLayout.NORTH);
		outputPanel.add(outputScroller, BorderLayout.CENTER);
		outputPanel.setVisible(false);


		//LayoutHelper.addComponent(mainPanel, (GridBagLayout) mainPanel.getLayout(), outputScroller, 0, 5, 1, 1, 0, 0);

		msg = new MessageProcessor(new PrintStream(new OutputStream() {
			@Override
			public void write(byte[] b, int off, int len) throws IOException {
				output.append(new String(b, off, len));
			}

			@Override
			public void write(byte[] b) throws IOException {
				this.write(b, 0, b.length);
			}

			@Override
			public void write(int b) throws IOException {
				this.write(new byte[] { (byte) b });
			}

		}), new PrintStream(new OutputStream() {
			@Override
			public void write(byte[] b, int off, int len) throws IOException {
				Font f = output.getFont();
				output.setFont(new Font("Arial", Font.BOLD, 12));
				output.append(new String(b, off, len));
				output.setFont(f);
			}

			@Override
			public void write(byte[] b) throws IOException {
				this.write(b, 0, b.length);
			}

			@Override
			public void write(int b) throws IOException {
				this.write(new byte[] { (byte) b });
			}
		}));
		

		/*
		 *  add JPanel to JFrame
		 */ 
		this.getContentPane().setLayout(new BorderLayout());
		this.getContentPane().add(titlePanel, BorderLayout.NORTH);
		this.getContentPane().add(mainPanel, BorderLayout.CENTER);
		this.getContentPane().add(outputPanel, BorderLayout.SOUTH);
		
		this.setDefaultCloseOperation(EXIT_ON_CLOSE);
		this.setName("mainframe");
		this.addWindowListener(this);
		pack();
		
		setLocationRelativeTo(null);
		setVisible(true);
		setResizable(false);
	}
	
	
	public void writeInputFile(String filename, String organism,
			String superclass, String seq, Vector<String> domains) {

		// obtain superclass

		ArrayList<String> classes = new ArrayList<String>();
		for (int i = 0; i < superclasses.length; i++) {
			classes.add(superclasses[i]);
		}
		String class_id = null;
		if (classes.contains(superclass)) {
			class_id = (classes.indexOf(superclass) + 1) % 5 + ".0.0.0.0.";
		}
		BufferedWriter bw;
		try {
			bw = new BufferedWriter(new FileWriter(new File(filename)));
			bw.write("NA  QueryTF\n");
			bw.write("XX\n");
			bw.write("SP  " + organism + "\n");
			bw.write("XX\n");
			if (class_id != null) {
				bw.write("CL  " + class_id + "\n");
				bw.write("XX\n");
			}
			
			// write sequence
			int SEQLINELENGTH = 60;

			for (int i = 0; i < (seq.length() / SEQLINELENGTH); i++) {
				bw.write("S1  ");
				bw.write(seq.toUpperCase(), i * SEQLINELENGTH, SEQLINELENGTH);
				bw.write("\n");
			}

			if (seq.length() - (seq.length() / SEQLINELENGTH) * SEQLINELENGTH > 0) {
				bw.write("S1  ");
				bw.write(seq.toUpperCase(), (seq.length() / SEQLINELENGTH)
						* SEQLINELENGTH, seq.length()
						- (seq.length() / SEQLINELENGTH) * SEQLINELENGTH);
				bw.write("\n");
			}
			bw.write("XX\n");

			// write domains
			String[] split;
			for (int i = 0; i < domains.size(); i++) {
				split = domains.get(i).split("\\s\\p{Punct}\\s"); // " - "
				bw.write("FT  POS  " + split[0] + "  "
						+ new StringTokenizer(split[1]).nextToken() + "\n");
			}
			bw.write("XX\n");
			bw.write("//\n");
			bw.write("XX\n");

			bw.flush();
			bw.close();

		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	@SuppressWarnings("deprecation")
	public void actionPerformed(ActionEvent e) {

		if (e.getSource() instanceof JButton) {
			JButton b = (JButton) e.getSource();

			if (b.getName().equals("run")) {

				// Show a message on windows
				boolean isWindows = (System.getProperty("os.name").toLowerCase().contains("windows"));
		        if (isWindows) {
		          JOptionPane.showMessageDialog(this, "Sorry, " + SABINE_Main.appName + " is using third party" +
		              " libraries that are not available for windows.\n\nPlease run this application" +
		              " on a UNIX system.", SABINE_Main.appName, JOptionPane.WARNING_MESSAGE);
		          return;
		        }
		        else if (isMacOSX()) {
		          // TODO: Show something?
	            //JOptionPane.showMessageDialog(this, appName + " has been designed for LINUX." +
	              //" You are using a MAC", appName, JOptionPane.INFORMATION_MESSAGE);
		        }
		        
	        
				/*
				 * generate temporary base directory and write input file
				 */
			    
				if (base_dir == null) {
				
					int randnum = RandomStringGenerator.randomNumber(3);
					TimeStampGenerator time_gen = new TimeStampGenerator();
					String time_stamp = time_gen.getTimeStamp().replace(':', '.');
					
					// Hack: fixed base directory
					base_dir = SABINE_Main.createBaseDir();
					//base_dir = "tmp/11.09.2012_10:11_13/";
					File base_dir_path = new File(base_dir);
					try {
					  base_dir_path.mkdirs();
					} catch (Exception ex){}
	
					if (!base_dir_path.exists()) {
						System.err.println("Invalid base directory. Aborting.");
						System.err.println("Base directory: " + base_dir + "\n");
						System.exit(1);
					}
				}
				
				String seq = sequenceArea.getText().toUpperCase().replaceAll(
						"\\n", "").replaceAll("\\s", "");
				
				String input_file = base_dir + "infile.tmp";
				
				writeInputFile(input_file, 
						(String) organismSelect.getSelectedItem(), 
						(String) superClassSelect.getSelectedItem(), 
						seq, domList);

				/*
				 * launch SABINE
				 */

				// read parameters
				double bmt;
				boolean dyn_bmt;
				if (bestMatchThreshold.getValue().equals("Dynamic")) {
					bmt = 0.95;
					dyn_bmt = true;
				}
				else {
					bmt = Double.parseDouble("" + bestMatchThreshold.getValue());
					dyn_bmt = false;
				}
				int mnb = Integer.parseInt(maxNumBestMatches.getValue().toString());
				double oft = Double.parseDouble(outlierFilterThreshold.getValue().toString());

				// create directories for temporary files
				SABINE_Main.createTempDirectories(base_dir);

				// reassign output and error stream
				System.setOut(msg.getOutStream());
				System.setErr(msg.getErrorStream());
			

				// run SABINE on generated input file
				SABINE_Runner runner = new SABINE_Runner(base_dir, bmt, dyn_bmt, mnb, oft, this);
				runner.setProgressBar(new ProgressBarSwing(progressBar));
				currRun = new Thread(runner);
				currRun.start();
				
				// get threshold values for dynamic BMT
				double[] thresholds = FBPPredictor.getThresholdValues();
				high_conf_bmt = thresholds[0];
				medium_conf_bmt = thresholds[1];
				low_conf_bmt = thresholds[2];
				
				
				// disable GUI elements and enable stop button
				sequenceArea.setEnabled(false);
				superClassSelect.setEnabled(false);
				organismSelect.setEnabled(false);
				domLeftSpin.setEnabled(false);
				domRightSpin.setEnabled(false);
				domSelected.setEnabled(false);
				bestMatchThreshold.setEnabled(false);
				maxNumBestMatches.setEnabled(false);
				outlierFilterThreshold.setEnabled(false);
				runButton.setEnabled(false);
				addButton.setEnabled(false);
				deleteButton.setEnabled(false);

				// XXX We're hiding the main panel here.
				mainPanel.setVisible(false);
				outputPanel.setVisible(true);
				output.setText("");
				validate();
				pack();
				
			} else if (b.getName().equals("reset")) {

				if (currRun != null) currRun.stop();
				
				// activate GUI elements
				sequenceArea.setEnabled(true);
				superClassSelect.setEnabled(true);
				organismSelect.setEnabled(true);
				domLeftSpin.setEnabled(true);
				domRightSpin.setEnabled(true);
				domSelected.setEnabled(true);
				bestMatchThreshold.setEnabled(true);
				maxNumBestMatches.setEnabled(true);
				outlierFilterThreshold.setEnabled(true);
				runButton.setEnabled(false);
				addButton.setEnabled(true);

				outputPanel.setVisible(false);
				
				// reset 
				organismSelect.setSelectedItem("Homo Sapiens");
				superClassSelect.setSelectedIndex(0);
				sequenceArea.setText("");
				domLeftSpin.getModel().setValue(1);
				domRightSpin.getModel().setValue(1);
				domList.clear();
				domSelected.setListData(domList);
				bestMatchThreshold.getModel().setValue("Dynamic");
				maxNumBestMatches.getModel().setValue(5);
				outlierFilterThreshold.getModel().setValue(0.5);
				output.setText("");
				validate();
				pack();
			
			} else if (b.getName().equals("demo")) {
				
				organismSelect.setSelectedItem("Xenopus Laevis");
				superClassSelect.setSelectedIndex(2);
				
				String demoSeq = 	"MLLERVRTGTQKSSDMCGYTGSPEIPQCAGCNQHIVDRFILKVLDRHWHSKCLKCNDCQI\n" +
								  	"QLAEKCFSRGDSVYCKDDFFKRFGTKCAACQQGIPPTQVVRRAQEFVYHLHCFACIVCKR\n" +
								  	"QLATGDEFYLMEDSRLVCKADYETAKQREAESTAKRPRTTITAKQLETLKNAYNNSPKPA\n" +
									"RHVREQLSSETGLDMRVVQVWFQNRRAKEKRLKKDAGRQRWGQYFRNMKRSRGNSKSDKD\n" +
									"SIQEEGPDSDAEVSFTDEPSMSEMNHSNGIYNSLNDSSPVLGRQAGSNGPFSLEHGGIPT\n" +
									"QDQYHNLRSNSPYGIPQSPASLQSMPGHQSLLSNLAFPDTGLGIIGQGGQGVAPTMRVIG\n" +
									"VNGPSSDLSTGSSGGYPDFPVSPASWLDEVDHTQF\n";
				
				int demoDomStart = 155;
				int demoDomStop = 211;
				
				domLeftSpin.getModel().setValue(demoDomStart);
				domRightSpin.getModel().setValue(demoDomStop);
				
				sequenceArea.setText(demoSeq);
				sequenceArea.setCaretPosition(0);
				
				demoSeq = demoSeq.replaceAll("\\n", "").replaceAll("\\s", "");
				domList.clear();
				domList.add(demoDomStart + " - " + demoDomStop + " :  " + demoSeq.substring(demoDomStart - 1, demoDomStop));
				domSelected.setListData(domList);
				
				domLeftSpin.getModel().setValue(demoDomStart);
				domRightSpin.getModel().setValue(demoDomStop);
				
				bestMatchThreshold.getModel().setValue("Dynamic");
				maxNumBestMatches.getModel().setValue(5);
				outlierFilterThreshold.getModel().setValue(0.5);
				
				
				// enable Run-Button and Reset-Button
				runButton.setEnabled(true);
				
			} else if (b.getName().equals("help")) {
				
				JHelpBrowser browser = new JHelpBrowser(this,"Help");


				browser.setSize(800, 600);
				browser.setLocationRelativeTo(null);
				browser.validate();
				browser.setVisible(true);
				
			} else if (b.getName().equals("saveBM")) {
				
				JFileChooser fileChooser = new JFileChooser();
				
				SABINEfileFilter txt_filter = new SABINEfileFilter(SABINEfileFilter.TXT);
				fileChooser.addChoosableFileFilter(txt_filter);
				fileChooser.setAcceptAllFileFilterUsed(false);
				
				int rVal = fileChooser.showSaveDialog(this);
				
				if (rVal == JFileChooser.APPROVE_OPTION) {
			    	
			    	File file = fileChooser.getSelectedFile();
			    	SABINEfileFilter filter = (SABINEfileFilter) fileChooser.getFileFilter();
			    	
			    	// save best matches
			    	if (filter.isTextFilter()) {
			    		
			    		if (! file.getName().toLowerCase().endsWith(".txt")) {
			    			file = new File(file.getAbsolutePath() + ".txt");
			    		}
			    		try {
			    			BufferedWriter bw = new BufferedWriter(new FileWriter(file));
			    			bw.write("Factor\tPFM similarity\n");
			    			for (int i=0; i<best_matches.length; i++) {
			    				bw.write(best_matches[i][0] + "\t" + best_matches[i][1] + "\n");
			    			}
			    			bw.flush();
			    			bw.close();	
			    		} catch (IOException e1) {
							e1.printStackTrace();
						}
			    	}
				}
				
			} else if (b.getName().equals("savePFM")) {
				
				JFileChooser fileChooser = new JFileChooser();
				
				SABINEfileFilter txt_filter = new SABINEfileFilter(SABINEfileFilter.TXT);
				fileChooser.addChoosableFileFilter(txt_filter);
				fileChooser.setAcceptAllFileFilterUsed(false);
				
				int rVal = fileChooser.showSaveDialog(this);
				
				if (rVal == JFileChooser.APPROVE_OPTION) {
			    	
			    	File file = fileChooser.getSelectedFile();
			    	SABINEfileFilter filter = (SABINEfileFilter) fileChooser.getFileFilter();
			    	
			    	// save predicted PFM to text file (STAMP format)
			    	if (filter.isTextFilter()) {
			    		
			    		if (! file.getName().toLowerCase().endsWith(".txt")) {
			    			file = new File(file.getAbsolutePath() + ".txt");
			    		}
			    		try {
			    			BufferedWriter bw = new BufferedWriter(new FileWriter(file));
			    			for (int i=0; i<predicted_pfm.size(); i++) {
			    				bw.write(predicted_pfm.get(i) + "\n");
			    			}
			    			bw.flush();
			    			bw.close();	
			    		} catch (IOException e1) {
							e1.printStackTrace();
						}
			    	}
				}
				
			} else if (b.getName().equals("saveLogo")) {
				
				JFileChooser fileChooser = new JFileChooser();
				
				SABINEfileFilter png_filter = new SABINEfileFilter(SABINEfileFilter.PNG);
				
				fileChooser.addChoosableFileFilter(png_filter);
				fileChooser.setAcceptAllFileFilterUsed(false);
				
				int rVal = fileChooser.showSaveDialog(this);
				
			    if (rVal == JFileChooser.APPROVE_OPTION) {
			    	
			    	File file = fileChooser.getSelectedFile();
			    	SABINEfileFilter filter = (SABINEfileFilter) fileChooser.getFileFilter();
			    	
			    	// save sequence logo of predicted PMF as image file 
			    	if (filter.isImageFilter()) {
			    	
			    		if (! file.getName().toLowerCase().endsWith(".png")) {
			    			file = new File(file.getAbsolutePath() + ".png");
			    		}
			    		try {
							javax.imageio.ImageIO.write(SequenceLogo.drawSequenceLogo(logoArray, 100), "png", file);
						} catch (Exception e1) {
							e1.printStackTrace();
						}
			    	}
			    }
 
			} else if (b.getName().equals("add")) {

				int domStart = Integer.parseInt(domLeftSpin.getValue()
						.toString());
				int domStop = Integer.parseInt(domRightSpin.getValue()
						.toString());
				String seq = sequenceArea.getText().replaceAll("\\s", "")
						.replaceAll("\\n", "");

				String curr_dom = domStart + " - " + domStop + " :  "
						+ seq.substring(domStart - 1, domStop);

				if (!domList.contains(curr_dom)) {
					domList.add(curr_dom);
					domSelected.setListData(domList);
					dom_added = true;
					if (seq_added)
						runButton.setEnabled(true);
					validate();
				}
			} else if (b.getName().equals("delete")) {
				int indices[] = domSelected.getSelectedIndices();
				for (int i = indices.length - 1; i > -1; domList
						.remove(indices[i--]))
					;

				if (domList.isEmpty()) {
					dom_added = false;
					runButton.setEnabled(false);
				}
				domSelected.setListData(domList);
				validate();
				
			} else if (b.getName().equals("back")) {
				this.getContentPane().remove(resultsPanel);
				this.getContentPane().add(mainPanel, BorderLayout.CENTER);
				this.getContentPane().add(outputPanel, BorderLayout.SOUTH);
				outputPanel.setVisible(false);
				validate();
				pack();
				
			} else if (b.getName().equals("enlarge")) {
				
				BufferedImage image = null;
				
				try {
					image = SequenceLogo.drawSequenceLogo(logoArray, 100);
				} catch (Exception e1) {
					e1.printStackTrace();
				}
			        
				// create panel for sequence logo 
				
				JLabel logo = new JLabel(new ImageIcon(image));
				JPanel logoHelpPanel = new JPanel(new GridLayout(1,1));
				logoHelpPanel.add(logo);
				logoHelpPanel.setBackground(Color.white);
				JScrollPane logoScroller = new JScrollPane(logoHelpPanel, JScrollPane.VERTICAL_SCROLLBAR_NEVER, 
																          JScrollPane.HORIZONTAL_SCROLLBAR_AS_NEEDED);
				
				JFrame logoFrame = new JFrame();
		        logoFrame.setName("logoframe");

				logoFrame.setTitle("Sequence Logo");
				logoFrame.getContentPane().add(logoScroller);
				logoFrame.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
				logoFrame.addWindowListener(this);
				enlargeButton.setEnabled(false);
				logoFrame.validate();
				logoFrame.pack();
				logoFrame.setSize(Math.min(logoFrame.getWidth(), 800), logoFrame.getHeight()+30);
				logoFrame.setLocationRelativeTo(null);
		        logoFrame.setVisible(true);
			}
		}
	}

	public static void main(String... args) {
		new SABINE_GUI();
	}

	public void stateChanged(ChangeEvent e) {

		SpinnerModel model_left = domLeftSpin.getModel();
		SpinnerModel model_right = domRightSpin.getModel();

		int val_left = Integer.parseInt(model_left.getValue().toString());
		int val_right = Integer.parseInt(model_right.getValue().toString());

		int seq_length = sequenceArea.getText().replaceAll("\\s", "")
				.replaceAll("\\n", "").length();

		if (e.getSource() instanceof JSpinner) {

			JSpinner spinner = (JSpinner) e.getSource();

			if (spinner.getName().equals("domLeftSpin") && val_left > val_right) {
				model_right.setValue(model_left.getValue());
			}
			if (spinner.getName().equals("domRightSpin") && val_left > val_right) {
				model_left.setValue(model_right.getValue());
			}
			if (spinner.getName().equals("domRightSpin") && val_right > seq_length) {
				model_right.setValue(Math.max(seq_length,1));
			}
		}
	}

	
	public void updateDomainList(String seq) {

		StringTokenizer strtok;
		boolean[] invalid_domains = new boolean[domList.size()];
		int start, stop;

		for (int i = 0; i < domList.size(); i++) {

			strtok = new StringTokenizer(domList.get(i));
			start = Integer.parseInt(strtok.nextToken());
			strtok.nextToken();
			stop = Integer.parseInt(strtok.nextToken());

			if (stop > seq.length()) {
				invalid_domains[i] = true;
			} else {
				domList.set(i, start + " - " + stop + " :  "
						+ seq.substring(start - 1, stop));
			}
		}

		for (int i = domList.size() - 1; i >= 0; i--) {
			if (invalid_domains[i])
				domList.remove(i);
		}
		if (domList.isEmpty()) {
			runButton.setEnabled(false);
		}
		
		domSelected.setListData(domList);
		validate();
	}

	public void changedUpdate(DocumentEvent e) {
		System.out.println("Change");
	}

	
	public void insertUpdate(DocumentEvent e) {
		
		if (e.getDocument().getProperty("name").equals("output")) {
			
			output.setCaretPosition(output.getDocument().getLength()); 
			
		} else if (e.getDocument().getProperty("name").equals("seq")) {
			
			String seq = sequenceArea.getText().replaceAll("\\s", "").replaceAll("\\n", "");
		
			try {
				ProteinTools.createProtein(seq);
			} catch (IllegalSymbolException ex) {
				JOptionPane.showMessageDialog(this, "Invalid protein sequence.", ex.getClass().getSimpleName(),
						JOptionPane.ERROR_MESSAGE);
			}
			
			if (! seq.isEmpty()) {
				seq_added = true;
				addButton.setEnabled(true);
				domLeftSpin.setEnabled(true);
				domRightSpin.setEnabled(true);
			
				if (dom_added) {
					runButton.setEnabled(true);
				}
				updateDomainList(seq);
			}
		}
	}

	public void removeUpdate(DocumentEvent e) {
		
		if (e.getDocument().getProperty("name").equals("seq")) {
		
			String seq = sequenceArea.getText().replaceAll("\\s", "").replaceAll("\\n", "");
			
			if (seq.isEmpty()) {
				seq_added = false;
				addButton.setEnabled(false);
				domLeftSpin.setEnabled(false);
				domRightSpin.setEnabled(false);
				runButton.setEnabled(false);
			}
			
			updateDomainList(seq);
	
			// update limit for right JSpinner
			int val_right = Integer.parseInt(domRightSpin.getModel().getValue().toString());
			if (val_right > seq.length()) {
				domRightSpin.getModel().setValue(Math.max(seq.length(), 1));
			}
		}
	}

	public void valueChanged(ListSelectionEvent e) {
		if (domSelected.getSelectedIndex() == -1) {
			deleteButton.setEnabled(false);

		} else {
			deleteButton.setEnabled(true);
		}
	}
	
	public boolean readOutfile(String outfile) {
		
		try {
			BufferedReader br = new BufferedReader(new FileReader(new File(outfile)));
			
			StringTokenizer strtok;
			String line;
			int pos;
			double A, C, G, T;
			ArrayList<DistributionLogo> logoData = new ArrayList<DistributionLogo>();
			predicted_pfm = new ArrayList<String>();
			
			ArrayList<Object[]> bestMatches = new ArrayList<Object[]>();
			ArrayList<Object[]> predMatrix = new ArrayList<Object[]>();
			
			// skip best matches
			while (!(line = br.readLine()).startsWith("BM")); 
			
			while (line.startsWith("BM")) {
				
				if (line.startsWith("BM  none")) {
					return false;
				}
				
				DecimalFormat fmt = new DecimalFormat();
				fmt.setMaximumFractionDigits(4);
				fmt.setMinimumFractionDigits(4);
				
				strtok = new StringTokenizer(line.substring(4));
				bestMatches.add(new Object[]{ strtok.nextToken(), fmt.format(Double.parseDouble(strtok.nextToken()))});

				line = br.readLine();
			}
			line = br.readLine();	// XX
			best_matches = bestMatches.toArray(new Object[][] {});
			
			while (line.startsWith("MA")) {
				
				predicted_pfm.add(line.substring(4));
				strtok = new StringTokenizer(line.substring(4));
				
				pos = Integer.parseInt(strtok.nextToken()); 
				A = Double.parseDouble(strtok.nextToken()) / 100;
				C = Double.parseDouble(strtok.nextToken()) / 100;
				G = Double.parseDouble(strtok.nextToken()) / 100;
				T = Double.parseDouble(strtok.nextToken()) / 100;
				
				predMatrix.add(new Object[] {pos, A, C, G, T});
				logoData.add(SequenceLogo.getLogo(new double[] {A, C, G, T}));

				line = br.readLine();
			}
			pred_matrix = predMatrix.toArray(new Object[][] {});
			
			logoArray = logoData.toArray(new DistributionLogo[] {});

		} catch (IOException e) {
			e.printStackTrace();
		}
		return true;
	}

	public void done() {
		
		sequenceArea.setEnabled(true);
		superClassSelect.setEnabled(true);
		organismSelect.setEnabled(true);
		domLeftSpin.setEnabled(true);
		domRightSpin.setEnabled(true);
		domSelected.setEnabled(true);
		bestMatchThreshold.setEnabled(true);
		maxNumBestMatches.setEnabled(true);
		outlierFilterThreshold.setEnabled(true);
		runButton.setEnabled(true);
		addButton.setEnabled(true);
		mainPanel.setVisible(true);

		validate();
		pack();
		
		boolean pfm_transferred = readOutfile(base_dir + "prediction.out");
		
		if (pfm_transferred) {
		  this.getContentPane().remove(outputPanel);
		  this.getContentPane().remove(mainPanel);
		  showResults();	
		}
	}

	public void showResults() {
		
		/*
		 *  add panel for predicted best matches
		 */
		
		bestMatchTable = new JTable(best_matches, new String[] {"Factor" , "PFM similarity"});
		bestMatchTable.setPreferredScrollableViewportSize(new Dimension(400, 50));
		bestMatchTable.setFillsViewportHeight(true);

		JScrollPane bestMatchScroller = new JScrollPane(bestMatchTable, JScrollPane.VERTICAL_SCROLLBAR_AS_NEEDED,	
																		JScrollPane.HORIZONTAL_SCROLLBAR_AS_NEEDED);
		
		bestMatchSaveButton = new JButton("Save", UIManager.getIcon("FileView.floppyDriveIcon"));
		bestMatchSaveButton.setName("saveBM");
		bestMatchSaveButton.addActionListener(this);
		
		JPanel bestMatchPanel = new JPanel();
		bestMatchPanel.setBorder(BorderFactory.createTitledBorder(" Best Matches "));
		bestMatchPanel.setLayout(new GridBagLayout());
		LayoutHelper.addComponent(bestMatchPanel, (GridBagLayout) bestMatchPanel.getLayout(), bestMatchScroller, 0, 0, 2, 1, 0, 0);
		LayoutHelper.addComponent(bestMatchPanel, (GridBagLayout) bestMatchPanel.getLayout(), new JPanel(), 0, 1, 1, 1, 1, 0);
		LayoutHelper.addComponent(bestMatchPanel, (GridBagLayout) bestMatchPanel.getLayout(), bestMatchSaveButton, 1, 1, 1, 1, 0, 0);

		
		/*
		 *  add panel for predicted matrix
		 */		
		
		predMatrixTable = new JTable(pred_matrix, new String[] {"Pos." , "A", "C", "G", "T"});
		predMatrixTable.setPreferredScrollableViewportSize(new Dimension(400, 145));
		predMatrixTable.setFillsViewportHeight(true);
		JScrollPane predMatrixScroller = new JScrollPane(predMatrixTable, JScrollPane.VERTICAL_SCROLLBAR_AS_NEEDED,	
									      								  JScrollPane.HORIZONTAL_SCROLLBAR_AS_NEEDED);
		
		matrixSaveButton = new JButton("Save", UIManager.getIcon("FileView.floppyDriveIcon"));
		matrixSaveButton.setName("savePFM");
		matrixSaveButton.addActionListener(this);
		
		JPanel predMatrixPanel = new JPanel();
		predMatrixPanel.setBorder(BorderFactory.createTitledBorder(" Transferred PFM "));
		predMatrixPanel.setLayout(new GridBagLayout());
		
		// confidence sign
		double lowest_score = Double.parseDouble(("" + best_matches[1][best_matches[1].length-1]).replace(',', '.'));
		
		String conf_img_path;
		if (lowest_score > high_conf_bmt)
			conf_img_path = "img/high_conf.png";
		else if (lowest_score > medium_conf_bmt) {
			conf_img_path = "img/medium_conf.png";
			System.out.println("Medium confidence");
		}
		else {
			conf_img_path = "img/low_conf.png";
			System.out.println("Low confidence");
		}
		JLabel confLabel = new JLabel();
		try {
			
			double scaling_factor = 0.35;
			BufferedImage confImage = javax.imageio.ImageIO.read(Resource.class.getResource(conf_img_path));
			Image confsignImage = confImage.getScaledInstance((int) (confImage.getWidth()* scaling_factor), (int) (confImage.getHeight()*scaling_factor), Image.SCALE_SMOOTH);
			confLabel = new JLabel(new ImageIcon(confsignImage));
			
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		LayoutHelper.addComponent(predMatrixPanel, (GridBagLayout) predMatrixPanel.getLayout(), predMatrixScroller, 0, 0, 3, 1, 0, 0);
		LayoutHelper.addComponent(predMatrixPanel, (GridBagLayout) predMatrixPanel.getLayout(), confLabel, 0, 1, 1, 1, 0, 0);
		LayoutHelper.addComponent(predMatrixPanel, (GridBagLayout) predMatrixPanel.getLayout(), new JPanel(), 1, 1, 1, 1, 1, 0);
		LayoutHelper.addComponent(predMatrixPanel, (GridBagLayout) predMatrixPanel.getLayout(), matrixSaveButton, 2, 1, 1, 1, 0, 0);
		
		resultsPanel = new JPanel();
		resultsPanel.setLayout(new GridBagLayout());
		

		
		/*
		 *  add panel for sequence logo
		 */		
		
		JPanel logoPanel = new JPanel();
		BufferedImage image = null;
		
		try {
			image = SequenceLogo.drawSequenceLogo(logoArray, 35);
		} catch (Exception e1) {
			e1.printStackTrace();
		}
	        
		// create panel for sequence logo 
		logoPanel.setBorder(BorderFactory.createTitledBorder(" Sequence Logo "));
        logoPanel.setLayout(new GridBagLayout());
        logoPanel.setMinimumSize(new Dimension(title_width, 200));

        JLabel logo = new JLabel(new ImageIcon(image));
        JPanel logoHelpPanel = new JPanel(new GridLayout(1,1));
        logoHelpPanel.add(logo);
        logoHelpPanel.setBackground(Color.white);
        
        JScrollPane logoScrollPane = new JScrollPane(logoHelpPanel, 
        											 JScrollPane.VERTICAL_SCROLLBAR_AS_NEEDED, 
        											 JScrollPane.HORIZONTAL_SCROLLBAR_AS_NEEDED);
        
        // increase height of JScrollPane, if scrollbar is needed
        if (image.getWidth() > title_width)
        	logoScrollPane.setPreferredSize(new Dimension(title_width, image.getHeight() + 45));
        else {
        	logoScrollPane.setPreferredSize(new Dimension(title_width, image.getHeight() + 25));
        }
        logoScrollPane.setBorder(BorderFactory.createLoweredBevelBorder());
		
        
        enlargeButton = new JButton("Enlarge");
        enlargeButton.setName("enlarge");
        enlargeButton.addActionListener(this);
        
        JButton logoSaveButton = new JButton("Save", UIManager.getIcon("FileView.floppyDriveIcon"));
        logoSaveButton.setName("saveLogo");
        logoSaveButton.addActionListener(this);
        
        LayoutHelper.addComponent(logoPanel, (GridBagLayout) logoPanel.getLayout(), logoScrollPane, 0, 0, 4, 1, 0, 0);
        LayoutHelper.addComponent(logoPanel, (GridBagLayout) logoPanel.getLayout(), new JPanel(), 0, 1, 1, 1, 1, 0);
        LayoutHelper.addComponent(logoPanel, (GridBagLayout) logoPanel.getLayout(), enlargeButton, 1, 1, 1, 1, 0, 0);        
        LayoutHelper.addComponent(logoPanel, (GridBagLayout) logoPanel.getLayout(), new JPanel(), 2, 1, 1, 1, 0, 0);
        LayoutHelper.addComponent(logoPanel, (GridBagLayout) logoPanel.getLayout(), logoSaveButton, 3, 1, 1, 1, 0, 0);        	

		
        JButton backButton = new JButton("Back");
        backButton.setName("back");
        backButton.addActionListener(this);
        
        JPanel backPanel = new JPanel();
		backPanel.add(backButton);

		
		LayoutHelper.addComponent(resultsPanel, (GridBagLayout) resultsPanel.getLayout(), bestMatchPanel, 0, 0, 1, 1, 0, 0);
		LayoutHelper.addComponent(resultsPanel, (GridBagLayout) resultsPanel.getLayout(), predMatrixPanel, 0, 1, 1, 1, 0, 0);
		LayoutHelper.addComponent(resultsPanel, (GridBagLayout) resultsPanel.getLayout(), logoPanel, 0, 2, 1, 1, 0, 0);
		LayoutHelper.addComponent(resultsPanel, (GridBagLayout) resultsPanel.getLayout(), backPanel, 0, 3, 1, 1, 0, 0);
		
		this.getContentPane().add(resultsPanel, BorderLayout.CENTER);
		validate();
		pack();
	}
	
	/*
	 * (non-Javadoc)
	 * @see java.lang.Object#finalize()
	 */
	@Override
	protected void finalize() throws Throwable {
		super.finalize();
		DirectoryRemover.removeDirectory(base_dir);
	}
	
	public void windowClosed(WindowEvent e) {
		if (e.getSource() instanceof JFrame) {
			JFrame frame = (JFrame) e.getSource();
			
			if (frame.getName().equals("logoframe")) {
				enlargeButton.setEnabled(true);
			}
		}
	}
	
	public void windowClosing(WindowEvent e) {
	
		if (e.getSource() instanceof JFrame) {
			JFrame frame = (JFrame) e.getSource();
			
			if (frame.getName().equals("mainframe") && base_dir != null) {
				
				// CHANGED
				//DirectoryRemover.removeDirectory("tmp/");
			}
		}
	}
	public void windowActivated(WindowEvent e) {}	
	public void windowDeactivated(WindowEvent e) {}
	public void windowDeiconified(WindowEvent e) {}
	public void windowIconified(WindowEvent e) {}
	public void windowOpened(WindowEvent e) {}
}

