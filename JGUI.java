
import extension.DirectoryRemover;
import gui.JHelpBrowser;
import gui.LayoutHelper;
import gui.MessageProcessor;
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
import java.io.OutputStream;
import java.io.PrintStream;
import java.text.DecimalFormat;
import java.util.ArrayList;
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
import javax.swing.JScrollPane;
import javax.swing.JSpinner;
import javax.swing.JTable;
import javax.swing.JTextArea;
import javax.swing.SpinnerModel;
import javax.swing.SpinnerNumberModel;
import javax.swing.UIManager;
import javax.swing.UnsupportedLookAndFeelException;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;
import javax.swing.event.DocumentEvent;
import javax.swing.event.DocumentListener;
import javax.swing.event.ListSelectionEvent;
import javax.swing.event.ListSelectionListener;

import org.biojava.bio.gui.DistributionLogo;
import org.biojava.bio.seq.ProteinTools;
import org.biojava.bio.symbol.IllegalSymbolException;
import org.biojava.bio.symbol.SymbolList;

import resources.Resource;

public class JGUI extends JFrame implements ActionListener, ChangeListener,
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
	
	private PrintStream stdout;
	
	boolean seq_added = false;

	boolean dom_added = false;
	
	private String base_dir;

	final String[] superclasses = new String[] { "Basic Domain", "Zinc Finger",
			"Helix-Turn-Helix", "Beta Scaffold", "Other" };

	final String[] species = new String[] { "APLYSIA CALIFORNICA",
			"LOLIGO OPALESCENS", "ARTEMIA FRANCISCANA", "METAPENAEUS ENSIS",
			"PORCELLIO SCABER", "AEDES AEGYPTI", "CHIRONOMUS TENTANS",
			"DROSOPHILA ANANASSAE", "DROSOPHILA MELANOGASTER",
			"DROSOPHILA PSEUDOOBSCURA", "DROSOPHILA FUNEBRIS",
			"DROSOPHILA VIRILIS", "DROSOPHILA", "LUCILIA CUPRINA",
			"SARCOPHAGA PEREGRINA", "BOMBYX MORI", "MANDUCA SEXTA",
			"GALLERIA MELLONELLA", "HELIOTHIS VIRESCENS", "JUNONIA COENIA",
			"CHORISTONEURA FUMIFERANA", "TRIBOLIUM CASTANEUM",
			"ACHETA DOMESTICUS", "SCHISTOCERCA AMERICANA",
			"SCHISTOCERCA GREGARIA", "ONCOPELTUS FASCIATUS",
			"THERMOBIA DOMESTICA", "CUPIENNIUS SALEI", "HELOBDELLA ROBUSTA",
			"HELOBDELLA TRISERIALIS", "HIRUDO MEDICINALIS", "NEREIS VIRENS",
			"LINEUS SANGUINEUS", "LINGULA ANATINA", "CNEMIDOPHORUS UNIPARENS",
			"POEPHILA GUTTATA TAENIOPYGIA GUTTATA", "SERINUS CANARIA",
			"STURNUS VULGARIS", "COTURNIX COTURNIX",
			"COTURNIX COTURNIX JAPONICA", "GALLUS GALLUS",
			"MELEAGRIS GALLOPAVO", "TRACHEMYS SCRIPTA", "CALLITHRIX JACCHUS",
			"SAGUINUS OEDIPUS", "SAIMIRI SCIUREUS",
			"SAIMIRI BOLIVIENSIS BOLIVIENSIS", "AOTUS NANCYMAAE",
			"ATELES BELZEBUTH CHAMEK", "HOMO SAPIENS", "PAN TROGLODYTES",
			"HYLOBATES LAR", "CERCOPITHECUS AETHIOPS", "MACACA FASCICULARIS",
			"PAPIO HAMADRYAS", "EULEMUR FULVUS COLLARIS", "CRICETULUS GRISEUS",
			"CRICETULUS LONGICAUDATUS", "CRICETULUS SP.",
			"MESOCRICETUS AURATUS", "MUS DOMESTICUS TORINO STRAIN",
			"MUS MUSCULUS", "MUS SPEC.", "RATTUS NORVEGICUS", "RATTUS RATTUS",
			"MARMOTA MONAX", "SPERMOPHILUS BEECHEYI",
			"TAMIAS ASIATICUS SIBIRICUS", "CAVIA PORCELLUS",
			"ORYCTOLAGUS CUNICULUS", "TUPAIA GLIS BELANGERI", "BOS TAURUS",
			"BOVINE PAPILLOMA VIRUS TYPE 1", "OVIS ARIES", "SUS SCROFA",
			"CROCUTA CROCUTA", "FELIS SILVESTRIS", "FELIS SILVESTRIS CATUS",
			"CANIS FAMILIARIS CANIS CANIS", "EQUUS SPEC.", "EQUUS CABALLUS",
			"SMINTHOPSIS MACROURA", "RANA CATESBEIANA", "RANA RUGOSA",
			"XENOPUS BOREALIS", "XENOPUS LAEVIS", "XENOPUS",
			"XENOPUS TROPICALIS", "AMBYSTOMA MEXICANUM",
			"NOTOPHTHALMUS VIRIDESCENS", "PAGRUS MAJOR CHRYSOPHRYS MAJOR",
			"SPARUS AURATA", "MICROPOGONIAS UNDULATUS", "OREOCHROMIS AUREUS",
			"OREOCHROMIS NILOTICUS TILAPIA NILOTICA", "FUGU RUBRIPES",
			"FUNDULUS HETEROCLITUS", "ORYZIAS LATIPES",
			"HIPPOGLOSSUS HIPPOGLOSSUS", "PARALICHTHYS OLIVACEUS",
			"ONCORHYNCHUS MYKISS", "SALMO SALAR", "CARASSIUS AURATUS",
			"CYPRINUS CARPIO", "BRACHYDANIO RERIO", "ICTALURUS PUNCTATUS",
			"ANGUILLA JAPONICA", "HETERODONTUS FRANCISCI",
			"PETROMYZON MARINUS", "BRANCHIOSTOMA FLORIDAE",
			"CIONA INTESTINALIS", "PHALLUSIA MAMMILATA", "HALOCYNTHIA RORETZI",
			"ASTERIAS VULGARIS", "ECHINOIDEA", "HELIOCIDARIS ERYTHROGRAMMA",
			"PARACENTROTUS LIVIDUS", "STRONGYLOCENTROTUS PURPURATUS",
			"LYTECHINUS VARIEGATUS", "PTYCHODERA FLAVA",
			"CAENORHABDITIS BRIGGSAE", "CAENORHABDITIS ELEGANS",
			"CAENORHABDITIS VULGARENSIS", "ONCHOCERCA VOLVULUS",
			"PRISTIONCHUS PACIFICUS", "SCHISTOSOMA MANSONI",
			"CHLOROHYDRA VIRIDISSIMA", "HYDRA LITTORALIS",
			"HYDRA MAGNIPAPILLATA", "HYDRA VULGARIS", "ELEUTHERIA DICHOTOMA",
			"PODOCORYNE CARNEA", "CHRYSAORA QUINQUECIRRHA",
			"EPHYDATIA FLUVIATILIS", "GEODIA CYDONIUM", "FUSARIUM SOLANI",
			"MAGNAPORTHE GRISEA", "NEUROSPORA CRASSA", "ASPERGILLUS FLAVUS",
			"ASPERGILLUS NIGER", "ASPERGILLUS ORYZAE", "PENICILLIUM URTICAE",
			"PENICILLIUM CHRYSOGENUM", "EMERICELLA NIDULANS",
			"COCHLIOBOLUS CARBONUM", "KLUYVEROMYCES LACTIS",
			"SACCHAROMYCES CEREVISIAE", "CANDIDA ALBICANS YEAST",
			"CANDIDA GLABRATA", "SCHIZOSACCHAROMYCES POMBE", "USTILAGO MAYDIS",
			"ANTIRRHINUM MAJUS", "CRATEROSTIGMA PLANTAGINEUM",
			"PERILLA FRUTESCENS", "SYRINGA VULGARIS", "CAPSICUM ANNUUM",
			"LYCOPERSICON ESCULENTUM", "LYCOPERSICON PERUVIANUM",
			"SOLANUM CHACOENSE", "SOLANUM DULCAMARA", "SOLANUM TUBEROSUM",
			"NICOTIANA SP.", "NICOTIANA SYLVESTRIS", "NICOTIANA TABACUM",
			"PETUNIA HYBRIDA", "PETUNIA INTEGRIFOLIA", "IPOMOEA BATATAS",
			"CATHARANTHUS ROSEUS", "PETROSELINUM CRISPUM P. HORTENSE",
			"PIMPINELLA BRACHYCARPA", "DAUCUS CAROTA", "PANAX GINSENG",
			"FLAVERIA BIDENTIS", "FLAVERIA TRINERVIA CLUSTERD YELLOWTOPS",
			"HELIANTHUS ANNUUS", "GERBERA HYBRIDA", "BRASSICA NAPUS",
			"BRASSICA OLERACEA VAR. BOTRYTIS",
			"BRASSICA OLERACEA VAR. ITALICA", "BRASSICA OLERACEA",
			"ARABIDOPSIS THALIANA", "JONOPSIDIUM ACAULE", "RAPHANUS SATIVUS",
			"SINAPIS ALBA", "GOSSYPIUM HIRSUTUM", "FRAGARIA ANANASSA",
			"MALUS DOMESTICA", "HUMULUS LUPULUS", "HEVEA BRASILIENSIS",
			"JATROPHA CURCAS", "POPULUS TRICHOCARPA", "BETULA PENDULA",
			"CORYLUS AVELLANA", "GLYCINE MAX", "PHASEOLUS VULGARIS",
			"VIGNA RADIATA PHASEOLUS AUREUS", "CICER ARIETINUM",
			"LOTUS JAPONICUS", "MEDICAGO SATIVA", "MEDICAGO TRUNCATULA",
			"PISUM SATIVUM", "CUCUMIS SATIVUS", "CUCURBITA MAXIMA",
			"EUCALYPTUS GLOBULUS SUBSP. BICOSTATA",
			"EUCALYPTUS GLOBULUS SUBSP. GLOBULUS", "EUCALYPTUS GRANDIS",
			"EUCALYPTUS GUNNII", "DIANTHUS CARYOPHYLLUS", "SILENE LATIFOLIA",
			"RUMEX ACETOSA", "VITIS VINIFERA", "RANUNCULUS BULBOSUS",
			"RANUNCULUS FICARIA", "CALTHA PALUSTRIS", "DELPHINIUM AJACIS",
			"DICENTRA EXIMIA", "PAPAVER NUDICAULE", "PACHYSANDRA TERMINALIS",
			"ARANDA DEBORAH", "PHALAENOPSIS SP.",
			"DENDROBIUM GREX MADAME THONG-IN", "ONCIDIUM CV. GOWER RAMSEY",
			"ASPARAGUS OFFICINALIS GARDEN ASPARAGUS", "HYACINTHUS ORIENTALIS",
			"LOLIUM PERENNE", "LOLIUM TEMULENTUM", "FESTUCA ARUNDINACEA",
			"HORDEUM VULGARE", "HORDEUM VULGARE L.", "SECALE CEREALE",
			"TRITICUM AESTIVUM", "AVENA FATUA", "ORYZA SATIVA",
			"COIX LACRYMA-JOBI", "SORGHUM BICOLOR", "ZEA MAYS",
			"PENNISETUM GLAUCUM", "ERAGROSTIS TEF", "LILIUM LONGIFLORUM",
			"PEPEROMIA HIRTA", "PIPER MAGNIFICUM", "ASARUM EUROPAEUM",
			"LIRIODENDRON TULIPIFERA", "MICHELIA FIGO", "PICEA ABIES",
			"PINUS RADIATA", "PINUS RESINOSA", "PINUS TAEDA",
			"CRYPTOMERIA JAPONICA", "CERATOPTERIS PTERIDOIDES",
			"CERATOPTERIS RICHARDII", "ADIANTUM RADDIANUM",
			"OPHIOGLOSSUM PEDUNCULOSUM", "LYCOPODIUM ANNOTINUM",
			"LYCOPODIUM SPEC.", "MARCHANTIA POLYMORPHA",
			"PHYSCOMITRELLA PATENS", "CHLAMYDOMONAS REINHARDTII",
			"ACANTHAMOEBA CASTELLANII", "DICTYOSTELIUM DISCOIDEUM",
			"ESCHERICHIA COLI", "EPSTEIN-BARR VIRUS", "HUMAN HERPESVIRUS 8",
			"HERPES SIMPLEX VIRUS TYPE 1", "ADENOVIRUS", "ADENOVIRUS TYPE 2",
			"HUMAN PAPILLOMA VIRUS TYPE 16", "SIMIAN VIRUS 40",
			"VACCINIA VIRUS", "AVIAN MUSCULOAPONEUROTIC FIBROSARCOMA VIRUS",
			"AVIAN MYELOCYTOMATOSIS VIRUS CMII",
			"AVIAN MYELOCYTOMATOSIS VIRUS", "AVIAN RETROVIRUS",
			"AVIAN SARCOMA VIRUS 31", "AVIAN SARCOMA VIRUS 17",
			"AVIAN MYELOBLASTOSIS VIRUS", "AVIAN MYELOCYTOMATOSIS VIRUS MC29",
			"AVIAN ERYTHROBLASTOSIS VIRUS", "FBJ MURINE OSTEOSARCOMA VIRUS",
			"FBR MURINE OSTEOSARCOMA VIRUS", "FELINE LEUKEMIA PROVIRUS FTT",
			"RETICULOENDOTHELIOSIS VIRUS STRAIN T",
			"HUMAN IMMUNODEFICIENCY VIRUS TYPE 1", "HUMAN HEPATITIS B VIRUS",
			"HEPATITIS C VIRUS GENOTYPE 1A ISOLATE H", "MINUTE VIRUS OF MICE" };

	private JTable predMatrixTable;

	private JButton matrixSaveButton;

	private JButton bestMatchSaveButton;



	static {
		try {
			UIManager.setLookAndFeel(UIManager.getSystemLookAndFeelClassName());
		} catch (ClassNotFoundException exc) {
			JOptionPane.showMessageDialog(null, "<html>" + exc.getMessage()
					+ "</html>", exc.getClass().getName(),
					JOptionPane.WARNING_MESSAGE);
			exc.printStackTrace();
		} catch (InstantiationException exc) {
			JOptionPane.showMessageDialog(null, "<html>" + exc.getMessage()
					+ "</html>", exc.getClass().getName(),
					JOptionPane.WARNING_MESSAGE);
			exc.printStackTrace();
		} catch (IllegalAccessException exc) {
			JOptionPane.showMessageDialog(null, "<html>" + exc.getMessage()
					+ "</html>", exc.getClass().getName(),
					JOptionPane.WARNING_MESSAGE);
			exc.printStackTrace();
		} catch (UnsupportedLookAndFeelException exc) {
			JOptionPane.showMessageDialog(null, "<html>" + exc.getMessage()
					+ "</html>", exc.getClass().getName(),
					JOptionPane.WARNING_MESSAGE);
			exc.printStackTrace();
		}

	}

	public JGUI() {
		super("SABINE 1.0");

		/*
		 * adapt JFrame to operating system
		 */

		mainPanel = new JPanel();
		GridBagLayout gbl = new GridBagLayout();
		mainPanel.setLayout(gbl);

		// rewrite species list
		String[] species_new = new String[species.length];
		StringTokenizer strtok;
		boolean first = true;

		for (int i = 0; i < species.length; i++) {

			species_new[i] = "";
			first = true;
			strtok = new StringTokenizer(species[i]);
			String partial_name;

			while (strtok.hasMoreTokens()) {
				partial_name = strtok.nextToken();
				partial_name = partial_name.substring(0, 1).toUpperCase()
						+ partial_name.substring(1).toLowerCase();

				if (first) {
					species_new[i] += partial_name;
					first = false;
				} else {
					species_new[i] += " " + partial_name;
				}
			}
		}
		
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
		organismSelect = new JComboBox(species_new);
		organismSelect.setSelectedItem("Homo Sapiens");
		organismSelect.setBackground(Color.white);

		// add combo box for superclass selection
		superClassSelect = new JComboBox(superclasses);
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

		bestMatchThreshold = new JSpinner(new SpinnerNumberModel(0.95, 0, 1.0, 0.05));
		maxNumBestMatches = new JSpinner(new SpinnerNumberModel(5, 1, 10, 1));
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
		
		output = new JTextArea(4, 20);
		output.getDocument().putProperty("name", "output");
		output.getDocument().addDocumentListener(this);
		
		
		outputScroller = new JScrollPane(output,JScrollPane.VERTICAL_SCROLLBAR_AS_NEEDED,
				 								JScrollPane.HORIZONTAL_SCROLLBAR_AS_NEEDED);
		outputScroller.setVisible(false);
		
		LayoutHelper.addComponent(mainPanel, (GridBagLayout) mainPanel.getLayout(), outputScroller, 0, 5, 1, 1, 0, 0);
		
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
		this.setDefaultCloseOperation(EXIT_ON_CLOSE);
		this.setName("mainframe");
		this.addWindowListener(this);
		pack();
		setLocationRelativeTo(null);
		setVisible(true);
	}
	
	
	public void writeInputFile(String filename, String organism,
			String superclass, String seq, Vector<String> domains) {

		// obtain superclass

		ArrayList<String> classes = new ArrayList<String>();
		for (int i = 0; i < superclasses.length; i++) {
			classes.add(superclasses[i]);
		}
		String class_id = (classes.indexOf(superclass) + 1) % 5 + ".0.0.0.0.";

		BufferedWriter bw;
		try {
			bw = new BufferedWriter(new FileWriter(new File(filename)));

			bw.write("NA  QueryTF\n");
			bw.write("XX\n");
			bw.write("SP  " + organism + "\n");
			bw.write("XX\n");
			bw.write("CL  " + class_id + "\n");
			bw.write("XX\n");

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

				/*
				 * generate temporary base directory and write input file
				 */
				
				// CHANGED
				base_dir = "tmp/my_basedir/";
				/*
				if (base_dir == null) {
				
					int randnum = RandomStringGenerator.randomNumber(3);
					TimeStampGenerator time_gen = new TimeStampGenerator();
					String time_stamp = time_gen.getTimeStamp();
	
					base_dir = "tmp/" + time_stamp + "_" + randnum + "/";
					File base_dir_path = new File(base_dir);
					
					
	
					if (!base_dir_path.exists() && !base_dir_path.mkdir()) {
						System.out.println("\nInvalid base directory. Aborting.");
						System.out.println("Base directory: " + base_dir + "\n");
						System.exit(0);
					}
				}
				*/
				String seq = sequenceArea.getText().toUpperCase().replaceAll(
						"\\n", "").replaceAll("\\s", "");
				
				String input_file = base_dir + "infile.tmp";
				
				writeInputFile(input_file, (String) organismSelect
						.getSelectedItem(), (String) superClassSelect
						.getSelectedItem(), seq, domList);

				/*
				 * launch SABINE
				 */

				// read parameters
				double bmt = Double.parseDouble(bestMatchThreshold.getValue()
						.toString());
				int mnb = Integer.parseInt(maxNumBestMatches.getValue()
						.toString());
				double oft = Double.parseDouble(outlierFilterThreshold
						.getValue().toString());

				// create directories for temporary files
				SABINE_Caller dir_creator = new SABINE_Caller();
				dir_creator.createTempDirectories(base_dir);

				// reassign output and error stream
				stdout = System.out;
				System.setOut(msg.getOutStream());
				System.setErr(msg.getErrorStream());
			

				// run SABINE on generated input file
				currRun = new Thread(new SABINE_Runner(base_dir, bmt, mnb, oft, this));
				currRun.start();

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

				outputScroller.setVisible(true);
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

				outputScroller.setVisible(false);
				
				// reset 
				organismSelect.setSelectedItem("Homo Sapiens");
				superClassSelect.setSelectedIndex(0);
				sequenceArea.setText("");
				domLeftSpin.getModel().setValue(1);
				domRightSpin.getModel().setValue(1);
				domList.clear();
				domSelected.setListData(domList);
				bestMatchThreshold.getModel().setValue(0.95);
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
				
				bestMatchThreshold.getModel().setValue(0.95);
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
				LayoutHelper.addComponent(mainPanel, (GridBagLayout) mainPanel.getLayout(), outputScroller, 0, 5, 1, 1, 0, 0);
				outputScroller.setVisible(false);
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
		new JGUI();
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
				
				SymbolList aa = ProteinTools.createProtein(seq);
			} catch (IllegalSymbolException ex) {
				JOptionPane.showMessageDialog(this, "Invalid protein sequence.", ex.getClass().getName().substring(ex.getClass().getName().lastIndexOf(".") + 1),
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

		validate();
		pack();
		
		boolean pfm_transferred = readOutfile(base_dir + "prediction.out");
		
		if (pfm_transferred) {
			mainPanel.remove(outputScroller);
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
		
		LayoutHelper.addComponent(predMatrixPanel, (GridBagLayout) predMatrixPanel.getLayout(), predMatrixScroller, 0, 0, 2, 1, 0, 0);
		LayoutHelper.addComponent(predMatrixPanel, (GridBagLayout) predMatrixPanel.getLayout(), new JPanel(), 0, 1, 1, 1, 1, 0);
		LayoutHelper.addComponent(predMatrixPanel, (GridBagLayout) predMatrixPanel.getLayout(), matrixSaveButton, 1, 1, 1, 1, 0, 0);
		
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

        JLabel logo = new JLabel(new ImageIcon(image));
        JPanel logoHelpPanel = new JPanel(new GridLayout(1,1));
        logoHelpPanel.add(logo);
        logoHelpPanel.setBackground(Color.white);
        
        JScrollPane logoScrollPane = new JScrollPane(logoHelpPanel, 
        											 JScrollPane.VERTICAL_SCROLLBAR_AS_NEEDED, 
        											 JScrollPane.HORIZONTAL_SCROLLBAR_AS_NEEDED);
        
        // increase height of JScrollPane, if scrollbar is needed
        if (image.getWidth() > title_width)
        	logoScrollPane.setPreferredSize(new Dimension(title_width, image.getHeight() + 25));
        else {
        	logoScrollPane.setPreferredSize(new Dimension(title_width, image.getHeight() + 5));
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
