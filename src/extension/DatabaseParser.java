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

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Calendar;
import java.util.PriorityQueue;
import java.util.StringTokenizer;

public class DatabaseParser {

  ArrayList<String> gene_IDs = new ArrayList<String>();						// Ensembl gene IDs
  ArrayList<String> prot_IDs = new ArrayList<String>();						// Ensembl protein IDs
  ArrayList<String[]> uniprot_IDs = new ArrayList<String[]>();				// UniProt IDs
  ArrayList<String> transfac_IDs = new ArrayList<String>();					// TRANSFAC accession numbers
  ArrayList<String> sequences = new ArrayList<String>();						// sequences from Ensembl or UniProt

  ArrayList<ArrayList<String>> domains = new ArrayList<ArrayList<String>>();  // domain assignments from DBD or Interpro
  ArrayList<ArrayList<String>> hmm_IDs = new ArrayList<ArrayList<String>>();	// HMM-IDs from PFAM or Superfamily
  ArrayList<String> superclasses = new ArrayList<String>();					// superclasses from TRANSFAC

  String species_ID = "";
  String species_name = "";


  public void transferDBDFields(DBDParser p) {
    gene_IDs = p.gene_IDs;
    prot_IDs = p.prot_IDs;
    sequences = p.sequences;
    domains = p.domains;
    hmm_IDs = p.hmm_IDs;
  }


  public void transferInterproFields(InterproParser p) {

    for (int i=0; i<p.uniprot_IDs.size(); i++) {
      uniprot_IDs.add(new String[] {p.uniprot_IDs.get(i)});
    }
    sequences = p.sequences;
    domains = p.domains;
    hmm_IDs = p.hmm_IDs;
  }

  public void transferTransfacFields(TransfacParser p) {

    for (int i=0; i<p.crossrefs.size(); i++) {
      uniprot_IDs.add(new String[] {p.crossrefs.get(i)});
    }
    sequences = p.sequences1;
    domains = p.domains;
    hmm_IDs = p.hmm_IDs;
    transfac_IDs = p.tf_names;
  }


  /*
   *  given a list of TRANSFAC superclasses this method returns the consensus superclass
   */

  public static String getConsensusClass(ArrayList<String> classes) {

    // generate list of unique superclasses
    ArrayList<String> unique_classes = new ArrayList<String>();

    for (int i=0; i<classes.size(); i++) {
      if (! unique_classes.contains(classes.get(i))) {
        unique_classes.add(classes.get(i));
      }
    }

    // sort superclasses
    PriorityQueue<String> sorted_classes = new PriorityQueue<String>();
    sorted_classes.addAll(classes);
    classes.clear();
    while(! sorted_classes.isEmpty()) {
      classes.add(sorted_classes.poll());
    }

    // count occurences of the classes and compute consensus
    int[] class_counter = new int[unique_classes.size()];
    int max_val, max_pos;
    max_val = max_pos = 0;
    for (int i=0; i<unique_classes.size(); i++) {
      class_counter[i] = (classes.lastIndexOf(unique_classes.get(i)) - classes.indexOf(unique_classes.get(i))) + 1;
      if (class_counter[i] >= max_val) {
        max_val = class_counter[i];
        max_pos = i;
      }
    }
    return unique_classes.get(max_pos);
  }



  public void getSpeciesID(String infile) {

    if (species_name.length() == 0) {
      System.out.println("Fatal Error. Unable to map species ID to name. Global variable \"species_name\" was not initialized. Aborting.");
      System.exit(0);
    }

    String line, curr_ID;
    StringTokenizer strtok;

    try {
      BufferedReader br = new BufferedReader(new FileReader(new File(infile)));

      while ((line = br.readLine()) != null) {

        strtok = new StringTokenizer(line, "\t");
        curr_ID = strtok.nextToken().trim();

        if (strtok.nextToken().trim().equals(species_name.toUpperCase())) {
          species_ID = curr_ID;
          return;
        }
      }
    }
    catch(IOException ioe) {
      System.out.println(ioe.getMessage());
      System.out.println("IOException occurred while mapping species name to ID.");
    }

    System.out.println("Unable to parse species name \"" + species_ID + "\". Aborting.");
    System.exit(0);
  }



  public void parseTFClassMapping(String infile) {

    ArrayList<String> parsed_HmmIDs = new ArrayList<String>();
    ArrayList<String> parsed_Superclasses = new ArrayList<String>();

    StringTokenizer strtok;
    String line;

    /*
     *  parse TF class mapping from input file
     */

    try {
      BufferedReader br = new BufferedReader(new FileReader(new File(infile)));

      while ((line = br.readLine()) != null) {

        strtok = new StringTokenizer(line, "\t");

        strtok.nextToken(); 									// family name
        strtok.nextToken(); 									// family ID
        parsed_HmmIDs.add(strtok.nextToken().trim());			// HMM ID
        parsed_Superclasses.add(strtok.nextToken().trim());	    // TRANSFAC superclass
      }
      br.close();
    }

    catch(IOException ioe) {
      System.out.println(ioe.getMessage());
      System.out.println("IOException occurred while parsing TF class mapping.");
    }

    /*
     *  map HMM-IDs to TRANSFAC superclasses
     */

    int curr_idx;
    ArrayList<String> curr_superclasses;
    String consensus_class;

    for (int i=0; i<hmm_IDs.size(); i++) {

      // map HMM-IDs to superclasses
      curr_superclasses = new ArrayList<String>();
      for (int j=0; j<hmm_IDs.get(i).size(); j++) {

        curr_idx = parsed_HmmIDs.indexOf(hmm_IDs.get(i).get(j));
        if (curr_idx != -1) {
          curr_superclasses.add(parsed_Superclasses.get(curr_idx));
        }
      }

      // compute consensus superclass
      if (curr_superclasses.isEmpty()) {
        consensus_class = "NA";
      }
      else {
        consensus_class = getConsensusClass(curr_superclasses);
      }
      superclasses.add(consensus_class);

    }
  }


  public void writeInputFile(String outfile, int[] relevant_factors) {

    String curr_seq;
    int SEQLINELENGTH = 60;
    int i;

    try {
      BufferedWriter bw = new BufferedWriter(new FileWriter(new File(outfile)));

      for (int x=0; x<relevant_factors.length; x++) {

        i = relevant_factors[x];

        /*
         *    write name and species
         */

        if (! prot_IDs.isEmpty()) {
          bw.write("NA  " + prot_IDs.get(i) + "\tgene:" + gene_IDs.get(i) + "\n");
        }
        else if (! transfac_IDs.isEmpty()) {
          bw.write("NA  " + transfac_IDs.get(i) + "\n");
        }
        else {
          bw.write("NA  " + uniprot_IDs.get(i)[0] + "\n");
        }
        bw.write("XX\n" +
            "SP  " + species_name + "\n" +
            "XX\n");

        /*
         *    write crossref. to UniProt
         */

        if (!prot_IDs.isEmpty() || !transfac_IDs.isEmpty()) {
          bw.write("RF  " + uniprot_IDs.get(i)[0]);
          if (uniprot_IDs.get(i).length > 1) {
            for (int j=1; i<uniprot_IDs.get(i).length; j++) {
              bw.write(", " + uniprot_IDs.get(i)[j]);
            }
          }
          bw.write("\n" +
              "XX\n");
        }

        /*
         *    write TRANSFAC superclass
         */

        bw.write("CL  " + superclasses.get(i) + ".1.1.1.1.\n" +
            "XX\n");

        /*
         *    write sequence
         */

        curr_seq = sequences.get(i);

        for(int j=0; j<(curr_seq.length()/SEQLINELENGTH); j++) {

          bw.write("S1  ");
          bw.write(curr_seq.toUpperCase(), j*SEQLINELENGTH, SEQLINELENGTH);
          bw.write("\n");

        }

        if((curr_seq.length()-((curr_seq.length()/SEQLINELENGTH)*SEQLINELENGTH)) > 0) {

          bw.write("S1  ");
          bw.write(curr_seq.toUpperCase(), (curr_seq.length()/SEQLINELENGTH)*SEQLINELENGTH, curr_seq.length()-((curr_seq.length()/SEQLINELENGTH)*SEQLINELENGTH));
          bw.write("\n");

        }

        bw.write("XX\n");

        /*
         *    write domains
         */

        for (int j=0; j<domains.get(i).size(); j++) {
          bw.write("FT  " + hmm_IDs.get(i).get(j) + "\t" + domains.get(i).get(j) + "\n");
        }

        bw.write("XX\n" +
            "//\n" +
            "XX\n");
      }

      bw.flush();
      bw.close();
    }

    catch(IOException ioe) {
      System.out.println(ioe.getMessage());
      System.out.println("IOException occurred while writing input file for SABINE.");
    }
  }

  public void filterDomains() {

    boolean[] empty_entrys = new boolean[domains.size()];

    for (int i=0; i<domains.size(); i++) {

      if (domains.get(i).isEmpty()) {
        empty_entrys[i] = true;
      }
    }

    for (int i=domains.size()-1; i>=0; i--) {

      if (empty_entrys[i]) {

        if (! prot_IDs.isEmpty()) {
          gene_IDs.remove(i);
          prot_IDs.remove(i);
        }
        else if (! transfac_IDs.isEmpty()) {
          transfac_IDs.remove(i);
          uniprot_IDs.remove(i);
        }
        else {
          uniprot_IDs.remove(i);
        }
        sequences.remove(i);
        domains.remove(i);
        hmm_IDs.remove(i);
      }
    }
  }



  public void filterSuperclasses() {

    if (superclasses.size() == 0) {
      System.out.println("Fatal Error. Unable to filter entrys. Global variable \"superclasses\" was not initialized. Aborting.");
      System.exit(0);
    }

    boolean[] NA_entrys = new boolean[superclasses.size()];

    for (int i=0; i<superclasses.size(); i++) {

      if (superclasses.get(i).equals("NA")) {
        NA_entrys[i] = true;
      }
    }

    for (int i=superclasses.size()-1; i>=0; i--) {

      if (NA_entrys[i] == true) {

        if (! prot_IDs.isEmpty()) {
          gene_IDs.remove(i);
          prot_IDs.remove(i);
        }
        if (! transfac_IDs.isEmpty()) {
          transfac_IDs.remove(i);
        }
        sequences.remove(i);
        domains.remove(i);
        hmm_IDs.remove(i);
        superclasses.remove(i);
        uniprot_IDs.remove(i);
      }
    }
  }

  public void parseGeneIDMapping(String infile) {

    if (gene_IDs.size() == 0) {
      System.out.println("Fatal Error. Unable to map gene IDs. Global variable \"gene_IDs\" was not initialized. Aborting.");
      System.exit(0);
    }

    String line, ensembl_ID, uniprot_ID;
    StringTokenizer strtok;
    int line_counter = 0;

    try {
      BufferedReader br = new BufferedReader(new FileReader(new File(infile)));


      while ((line = br.readLine()) != null) {

        strtok = new StringTokenizer(line);

        ensembl_ID = strtok.nextToken().trim();        	// parse Ensembl gene-ID

        // assert correct Ensembl gene-ID
        if (! ensembl_ID.equals(gene_IDs.get(line_counter++))) {
          System.out.println("Parse Error. Ensembl IDs are inconsistent. Aborting");
          System.out.println("Line " + line_counter + ":\n" +
              "  Ensembl-ID in text file:\t" + ensembl_ID + "\n" +
              "  Ensembl-ID in memory:\t" + gene_IDs.get(line_counter-1));
          System.exit(0);
        }

        // check if Ensembl identifiers could be successfully mapped to UniProt
        if (strtok.hasMoreTokens()) {
          strtok.nextToken().trim();						// skip UniProt name(s)
          uniprot_ID = strtok.nextToken().trim();			// parse UniProt ID(s)

          uniprot_IDs.add(uniprot_ID.split(","));         // add UniProt ID(s)
        }
        else {
          uniprot_IDs.add(new String[] {"NA"});
        }
      }
      br.close();

    }
    catch(IOException ioe) {
      System.out.println(ioe.getMessage());
      System.out.println("IOException occurred while mapping Ensembl gene IDs to UniProt IDs.");
    }
  }


  public void writeToFile(String protID_outfile,
    String seq_outfile,
    String domain_outfile,
    String hmmID_outfile,
    String geneID_outfile) {

    try {
      BufferedWriter bw_prot = new BufferedWriter(new FileWriter(new File(protID_outfile)));
      BufferedWriter bw_seq = new BufferedWriter(new FileWriter(new File(seq_outfile)));
      BufferedWriter bw_dom = new BufferedWriter(new FileWriter(new File(domain_outfile)));
      BufferedWriter bw_hmm = new BufferedWriter(new FileWriter(new File(hmmID_outfile)));

      BufferedWriter bw_gene = null;
      if (geneID_outfile != null) {
        bw_gene = new BufferedWriter(new FileWriter(new File(geneID_outfile)));
      }


      /*
       *  write domains and HMM-IDs
       */

      for (int i=0; i<domains.size(); i++) {
        for (int j=0; j<domains.get(i).size(); j++) {

          bw_dom.write(domains.get(i).get(j) + "\n");
          bw_hmm.write(hmm_IDs.get(i).get(j) + "\n");
        }
        bw_dom.write("\n");
        bw_hmm.write("\n");
      }

      bw_dom.flush();
      bw_dom.close();
      bw_hmm.flush();
      bw_hmm.close();

      /*
       *  write gene IDs, protein IDs and sequences
       */

      for (int i=0; i<hmm_IDs.size(); i++) {

        if (! prot_IDs.isEmpty()) {
          bw_prot.write(prot_IDs.get(i) + "\n");
        }
        else {
          bw_prot.write(uniprot_IDs.get(i)[0] + "\n");
        }
        bw_seq.write(sequences.get(i) + "\n");
        if ((bw_gene != null) && !gene_IDs.isEmpty()) {
          bw_gene.write(gene_IDs.get(i) + "\n");
        }
        if ((bw_gene != null) && !transfac_IDs.isEmpty()) {
          bw_gene.write(transfac_IDs.get(i) + "\n");
        }
      }
      bw_prot.flush();
      bw_prot.close();
      bw_seq.flush();
      bw_seq.close();

      if (bw_gene != null) {
        bw_gene.flush();
        bw_gene.close();
      }
    }

    catch(IOException ioe) {
      System.out.println(ioe.getMessage());
      System.out.println("IOException occurred while writing variables to text files.");
    }
  }

  public void readFromFile(String protID_infile,
    String seq_infile,
    String domain_infile,
    String hmmID_infile,
    String geneID_infile) {

    // clear variables if necessary
    if (! hmm_IDs.isEmpty()) {
      domains.clear(); hmm_IDs.clear(); gene_IDs.clear(); transfac_IDs.clear();
      prot_IDs.clear(); sequences.clear(); uniprot_IDs.clear();
    }

    String line, line_hmm, line_dom;

    /*
     *  read domains and HMM IDs
     */

    ArrayList<String> curr_domains = new ArrayList<String>();
    ArrayList<String> curr_hmm_IDs = new ArrayList<String>();

    try {
      BufferedReader br_dom = new BufferedReader(new FileReader(new File(domain_infile)));
      BufferedReader br_hmm = new BufferedReader(new FileReader(new File(hmmID_infile)));
      BufferedReader br_prot = new BufferedReader(new FileReader(new File(protID_infile)));
      BufferedReader br_seq = new BufferedReader(new FileReader(new File(seq_infile)));

      BufferedReader br_gene = null;
      if (geneID_infile != null) {
        br_gene = new BufferedReader(new FileReader(new File(geneID_infile)));
      }

      /*
       *  read domains and HMM-IDs
       */

      while (((line_dom = br_dom.readLine()) != null) && ((line_hmm = br_hmm.readLine()) != null)) {

        if (line_dom.length() == 0) {

          domains.add(curr_domains);
          hmm_IDs.add(curr_hmm_IDs);
          curr_domains = new ArrayList<String>();
          curr_hmm_IDs = new ArrayList<String>();
        }
        else {
          curr_domains.add(line_dom);
          curr_hmm_IDs.add(line_hmm);
        }
      }
      br_dom.close();
      br_hmm.close();


      /*
       *  read gene IDs, protein IDs and sequences
       */


      if ((br_gene != null) && !geneID_infile.contains("transfac")) {
        while ((line = br_prot.readLine()) != null) {
          prot_IDs.add(line);
        }
        while ((line = br_gene.readLine()) != null) {
          gene_IDs.add(line);
        }
      }
      else if ((br_gene != null) && geneID_infile.contains("transfac")) {
        while ((line = br_prot.readLine()) != null) {
          uniprot_IDs.add(new String[] {line});
        }
        while ((line = br_gene.readLine()) != null) {
          transfac_IDs.add(line);
        }
      }
      else {
        while ((line = br_prot.readLine()) != null) {
          uniprot_IDs.add(new String[] {line});
        }
      }
      while ((line = br_seq.readLine()) != null) {
        sequences.add(line);
      }

      br_prot.close();
      br_seq.close();
      if (br_gene != null) {
        br_gene.close();
      }
    }

    catch(IOException ioe) {
      System.out.println(ioe.getMessage());
      System.out.println("IOException occurred while reading variables from text files.");
    }
  }


  public void filterTransfacFlatfile(String infile, String outfile, boolean use_transfac_AC) {

    TransfacParser tf_parser = new TransfacParser();
    tf_parser.parseFactors(infile);

    /*
     *  filter parsed factors according to species, ref. to UniProt and superclass
     */

    int num_factors = tf_parser.tf_names.size();
    boolean[] relevant_factors = new boolean[num_factors];

    ArrayList<String> all_uniprot_IDs = new ArrayList<String>();
    ArrayList<Integer> index_map = new ArrayList<Integer>();

    if (! use_transfac_AC) {

      for (int i=0; i<uniprot_IDs.size(); i++) {
        for (int j=0; j<uniprot_IDs.get(i).length; j++) {
          all_uniprot_IDs.add(uniprot_IDs.get(i)[j]);
          index_map.add(i);
        }
      }

      for (int i=0; i<num_factors; i++) {

        if (species_name.toUpperCase().equals(tf_parser.species.get(i).toUpperCase()) &&
            all_uniprot_IDs.contains(tf_parser.crossrefs.get(i)) &&
            ! tf_parser.classes.get(i).equals("NA"))  {

          relevant_factors[i] = true;
        }
      }
    }
    else {

      for (int i=0; i<num_factors; i++) {

        if (species_name.toUpperCase().equals(tf_parser.species.get(i).toUpperCase()) &&
            transfac_IDs.contains(tf_parser.tf_names.get(i)) &&
            ! tf_parser.classes.get(i).equals("NA"))  {

          relevant_factors[i] = true;
        }
      }
    }


    for (int i=num_factors-1; i>=0; i--) {

      if (! relevant_factors[i]) {

        tf_parser.tf_names.remove(i);
        tf_parser.species.remove(i);
        tf_parser.crossrefs.remove(i);
        tf_parser.classes.remove(i);
        tf_parser.sequences1.remove(i);
        tf_parser.sequences2.remove(i);
        tf_parser.domains.remove(i);
        tf_parser.pfm_names.remove(i);
        tf_parser.pfms.remove(i);
      }
    }

    /*
     *  filter multiple entries
     */

    relevant_factors = new boolean[tf_parser.tf_names.size()];
    ArrayList<String> unique_IDs = new ArrayList<String>();

    if (use_transfac_AC) {
      for (int i=0; i<tf_parser.tf_names.size(); i++) {

        if (! unique_IDs.contains(tf_parser.tf_names.get(i))) {
          relevant_factors[i] = true;
          unique_IDs.add(tf_parser.tf_names.get(i));
        }
      }
    }
    else if (! prot_IDs.isEmpty() && ! gene_IDs.isEmpty()) {
      ArrayList<String> unique_ensemble_IDs = new ArrayList<String>();
      int id_idx;
      String curr_prot_ID;

      for (int i=0; i<tf_parser.crossrefs.size(); i++) {

        if (! unique_IDs.contains(tf_parser.crossrefs.get(i))) {
          id_idx = all_uniprot_IDs.indexOf(tf_parser.crossrefs.get(i));
          curr_prot_ID = prot_IDs.get(index_map.get(id_idx).intValue());

          if (! unique_ensemble_IDs.contains(curr_prot_ID)) {
            relevant_factors[i] = true;
            unique_IDs.add(tf_parser.tf_names.get(i));
            unique_ensemble_IDs.add(curr_prot_ID);
          }
        }
      }
    }
    else {

      for (int i=0; i<tf_parser.crossrefs.size(); i++) {

        if (! unique_IDs.contains(tf_parser.crossrefs.get(i))) {
          relevant_factors[i] = true;
          unique_IDs.add(tf_parser.crossrefs.get(i));
        }
      }
    }

    for (int i=relevant_factors.length-1; i>=0; i--) {

      if (! relevant_factors[i]) {

        tf_parser.tf_names.remove(i);
        tf_parser.species.remove(i);
        tf_parser.crossrefs.remove(i);
        tf_parser.classes.remove(i);
        tf_parser.sequences1.remove(i);
        tf_parser.sequences2.remove(i);
        tf_parser.domains.remove(i);
        tf_parser.pfm_names.remove(i);
        tf_parser.pfms.remove(i);
      }
    }

    /*
     *  filter outlier PFMs and merge them using STAMP
     */

    PFMFormatConverter pfm_converter = new PFMFormatConverter();

    ArrayList<ArrayList<String>> stamp_pfms = new ArrayList<ArrayList<String>>();
    ArrayList<String> merged_pfm = new ArrayList<String>();
    ArrayList<String[]> merged_transfac_pfm;

    for (int i=0; i<tf_parser.tf_names.size(); i++) {

      // more than one PFM ?
      if (tf_parser.pfms.get(i).size() > 1) {

        stamp_pfms = pfm_converter.convertAllTransfacToSTAMP(tf_parser.pfms.get(i));
        merged_pfm = pfm_converter.mergePFMs(stamp_pfms);

        merged_transfac_pfm = new ArrayList<String[]>();
        merged_transfac_pfm.add(pfm_converter.convertSTAMPToTransfac(merged_pfm));
        tf_parser.pfms.set(i, merged_transfac_pfm);
      }
    }

    // adjust matrix IDs
    tf_parser.pfm_names = pfm_converter.mergePFMnames(tf_parser.pfm_names);

    /*
     *  write filtered TRANSFAC flatfile
     */

    System.out.println("  Matrices found for " + tf_parser.tf_names.size() +
      " / " + uniprot_IDs.size() + " proteins.");

    tf_parser.writeLabelFile(outfile);
  }


  public void extractTestSet(String infile, String outfile, boolean use_transfac_AC) {

    ArrayList<String> all_uniprot_IDs = new ArrayList<String>();
    ArrayList<Integer> index_map = new ArrayList<Integer>();

    if (! use_transfac_AC) {
      for (int i=0; i<uniprot_IDs.size(); i++) {
        for (int j=0; j<uniprot_IDs.get(i).length; j++) {
          all_uniprot_IDs.add(uniprot_IDs.get(i)[j]);
          index_map.add(i);
        }
      }
    }

    try {
      BufferedReader br = new BufferedReader(new FileReader(new File(infile)));

      String line, curr_ID;
      int curr_idx;
      ArrayList<Integer> testset_idx = new ArrayList<Integer>();

      while ((line = br.readLine()) != null) {

        if (! use_transfac_AC) {

          if (line.startsWith("RF")) {
            curr_ID = line.substring(4).trim();
            curr_idx = index_map.get(all_uniprot_IDs.indexOf(curr_ID));
            testset_idx.add(curr_idx);

          }
        }
        else {

          if (line.startsWith("NA")) {
            curr_ID = line.substring(4).trim();
            curr_idx = transfac_IDs.indexOf(curr_ID);
            testset_idx.add(curr_idx);
          }
        }
      }

      /*
       *  write test set to output file
       */

      int[] print_perm = new int[testset_idx.size()];
      for (int i=0; i<print_perm.length; i++) {
        print_perm[i] = testset_idx.get(i);
      }

      writeInputFile(outfile, print_perm);
    }
    catch(IOException ioe) {
      System.out.println(ioe.getMessage());
      System.out.println("IOException occurred while extracting test set.");
    }
  }


  public void runDatabaseParser(String organism, String database, String base_dir) {

    String curr_dir = System.getProperty("user.dir") + "/";


    /*
     *  set organism
     */

    if      (organism.equals("human")) {
      species_name = "Homo sapiens";
    } else if (organism.equals("mouse")) {
      species_name = "Mus musculus";
    } else if (organism.equals("rat")) {
      species_name = "Rattus norvegicus";
    } else {
      usage();
    }


    /*
     *  set base directory
     */

    if (base_dir == null) {

      DecimalFormat fmt = new DecimalFormat();
      fmt.setMaximumIntegerDigits(2);
      fmt.setMinimumIntegerDigits(2);

      // get current time and date
      Calendar cal = Calendar.getInstance ();
      String curr_date = (fmt.format(cal.get(Calendar.DAY_OF_MONTH)) + "." +
          fmt.format((cal.get(Calendar.MONTH) + 1)) + "." +
          cal.get(Calendar.YEAR));
      String curr_time = (fmt.format(cal.get(Calendar.HOUR_OF_DAY)) + ":" +
          fmt.format(cal.get(Calendar.MINUTE)));

      if (database.equals("uniprot")) {
        base_dir = curr_dir + "results/uniprot/";
      } else if (database.equals("ensembl")) {
        base_dir = curr_dir + "results/ensembl/";
      } else if (database.equals("transfac")) {
        base_dir = curr_dir + "results/transfac/";
      }

      base_dir += organism + "/" + curr_date + "_" + curr_time + "/";

      if (! new File(base_dir).mkdir()) {
        System.out.println("\nInvalid base directory. Aborting.");
        System.out.println("Base directory: " + base_dir + "\n");
        System.exit(0);
      }
    }

    String tmp_dir;
    if (! base_dir.endsWith("/")) {
      base_dir += "/";
    }

    tmp_dir = base_dir + "tmp/";

    if (! new File(tmp_dir).exists() && ! new File(tmp_dir).mkdir()) {
      System.out.println("\nInvalid base directory. Aborting.");
      System.out.println("Base directory: " + base_dir + "\n");
      System.exit(0);
    }


    // set paths
    String map_dir = curr_dir + "mappings/";
    String ensembl_dir = curr_dir + "data/ensembl/";
    String uniprot_dir = curr_dir + "data/uniprot/";
    String interpro_dir = curr_dir + "data/interpro/";
    String transfac_dir = curr_dir + "data/transfac/";
    tmp_dir += "/";

    String seq_filename, dom_filename;
    String geneID_filename = null;
    String protID_filename = null;

    // get species ID
    getSpeciesID(map_dir + "SpeciesMapping.txt");

    boolean parse_transfac = false;
    String database_seq = "UniProt  (www.uniprot.org)";
    String database_dom = "InterPro (www.ebi.ac.uk/interpro/)";

    if (database.equals("ensembl")) {
      database_seq = "Ensembl (www.ensembl.org)";
      database_dom = "DBD     (www.transcriptionfactor.org)";
    }
    if (database.equals("transfac")) {
      database_seq = "TRANSFAC (www.gene-regulation.com)";
      database_dom = null;
      parse_transfac = true;
    }


    System.out.println();
    System.out.println("  ---------------------------------------------------------------------");
    System.out.println("  Database Parser for SABINE - StandAlone BINding specificity Estimator");
    System.out.println("  ---------------------------------------------------------------------");
    System.out.println("\n");
    System.out.println("  Organism:  " + species_name);
    System.out.println("  Databases: " + database_seq);
    if (database_dom != null) {
      System.out.println("             " + database_dom);
    }
    System.out.println();
    System.out.println("  Base directory: " + base_dir);
    System.out.println("\n");


    if (database.equals("uniprot")) {

      UniprotParser uniprot_parser = new UniprotParser();
      uniprot_parser.species_name = species_name;


      seq_filename = uniprot_dir + "uniprot-keyword-805-AND-keyword-DNA-binding-238.fasta";
      System.out.println("  Parsing UniProt sequences.\n    (file: " + seq_filename + ")\n");
      uniprot_parser.parseSequences(seq_filename);
      System.out.println("  " + uniprot_parser.sequences.size() + " sequences parsed.\n");

      InterproParser interpro_parser = new InterproParser(uniprot_parser.uniprot_IDs,
        uniprot_parser.sequences);

      dom_filename = interpro_dir + "match_complete.xml";
      System.out.println("  Parsing InterPro domain assignments.\n    (file: " + dom_filename + ")\n");
      interpro_parser.parseDomains(dom_filename,
        tmp_dir + "parsed_domains.txt");

      // set file name for gene-IDs
      geneID_filename = null;
      protID_filename = tmp_dir + "uniprot_IDs.tmp";

      // transfer fields to DatabaseParser object
      transferInterproFields(interpro_parser);
    }

    else if (database.equals("ensembl")) {


      DBDParser dbd_parser =  new DBDParser();

      // set species name and ID
      dbd_parser.species_name = species_name;
      dbd_parser.species_ID = species_ID;

      // parse datasets and save variables to disk
      seq_filename = ensembl_dir + dbd_parser.species_ID + ".fa";
      dom_filename = ensembl_dir + dbd_parser.species_ID + ".tf.ass";
      System.out.println("  Parsing Ensembl sequences.\n    (file: " + seq_filename + ")\n");
      dbd_parser.parseProtIDs(dom_filename);
      dbd_parser.parseSequences(seq_filename);
      System.out.println("  " + dbd_parser.sequences.size() + " sequences parsed.\n");


      System.out.println("  Parsing DBD domain assignments.\n    (file: " + dom_filename + ")\n");
      dbd_parser.parseDomains(dom_filename);

      // set file name for gene-IDs
      geneID_filename = tmp_dir + "ensembl_gene_IDs.tmp";
      protID_filename = tmp_dir + "ensembl_prot_IDs.tmp";

      // transfer fields to DatabaseParser object
      transferDBDFields(dbd_parser);
    }

    else if (database.equals("transfac")) {

      TransfacParser tf_parser = new TransfacParser();

      //seq_filename = "/data/dat0/transfac/2008.3/dat/factor.dat";
      seq_filename = "data/transfac/factor.dat";
      System.out.println("  Parsing TRANSFAC sequences and domains.\n    (file: " + seq_filename + ")\n");
      tf_parser.parseSequencesAndDomains(seq_filename, species_name);
      int[] dom_counter = tf_parser.countDomains();
      System.out.println("  " + tf_parser.sequences1.size() + " sequences and " +
          dom_counter[0] + " PFAM, " +
          dom_counter[1] + " PROSITE, " +
          dom_counter[2] + " SMART domains parsed.\n");


      System.out.println("  Parsing TRANSFAC domain assignments.\n" +
          "    (files: " + interpro_dir + "interpro.xml\n" +
          "            " + map_dir + "TFclassMapping.txt)\n");
      tf_parser.filterDomains(interpro_dir + "interpro.xml",
        map_dir + "TFclassMapping.txt");

      // set file name for gene-IDs
      geneID_filename = tmp_dir + "transfac_IDs.tmp";
      protID_filename = tmp_dir + "uniprot_IDs.tmp";

      // transfer fields to DatabaseParser object
      transferTransfacFields(tf_parser);
    }

    int old_num_proteins = domains.size();
    filterDomains();
    System.out.println("  Domain assignment found for " + domains.size() + " / "
        + old_num_proteins + " proteins.\n");

    writeToFile(protID_filename, tmp_dir + "seqs.tmp", tmp_dir + "domains.tmp",
      tmp_dir + "hmm_IDs.tmp", geneID_filename);

    readFromFile(protID_filename, tmp_dir + "seqs.tmp", tmp_dir + "domains.tmp",
      tmp_dir + "hmm_IDs.tmp", geneID_filename);

    if (!gene_IDs.isEmpty()) {
      System.out.println("  Mapping Ensembl gene identifiers to UniProt identifiers.\n" +
          "    (file: " + map_dir + species_ID + "_IDmapping.txt)\n");
      parseGeneIDMapping(map_dir + species_ID + "_IDmapping.txt");
      int succ_parsed_IDs = 0;
      for (int i=0; i<uniprot_IDs.size(); i++) {
        if (!uniprot_IDs.get(i)[0].equals("NA")) {
          succ_parsed_IDs++;
        }
      }
      System.out.println("  UniProt identifier found for " + succ_parsed_IDs + " / " +gene_IDs.size() + " proteins.\n");
    }

    System.out.println("  Mapping PFAM and SUPERFAMILY identifiers to TRANSFAC superclasses.\n" +
        "    (file: " + map_dir + "TFclassMapping.txt)\n");
    parseTFClassMapping(map_dir + "TFclassMapping.txt");

    int num_proteins =domains.size();
    filterSuperclasses();
    System.out.println("  TRANSFAC superclass found for " +domains.size() + " / "
        + num_proteins + " proteins.\n");

    // obtain set of labeled factors from TRANSFAC
    String matrix_file = transfac_dir + "transpall_interleaved_classes.out";
    String label_file = base_dir + species_ID + "_testset.labels";
    System.out.println("  Parsing matrices from TRANSFAC and JASPAR.\n    (file: " + matrix_file + ")\n");
    filterTransfacFlatfile(matrix_file, label_file, parse_transfac);
    System.out.println("\n  Writing label file which contains the obtained matrices.\n"
        + "    (file: " + label_file + ")\n");


    // write input file for Standalone predictor
    String outfile_test_set = base_dir +species_ID + "_testset.input";
    System.out.println("  Writing input file for SABINE which contains the evaluation set.\n"
        + "    (file: " + outfile_test_set + ")\n");
    extractTestSet(base_dir + species_ID + "_testset.labels",
      outfile_test_set, parse_transfac);


    // write input file for Standalone predictor
    int[] all_factors = new int[sequences.size()];
    for (int i=0; i<all_factors.length; i++) {
      all_factors[i] = i;
    }

    String outfile_all_factors = base_dir + species_ID + "_all_factors.input";
    System.out.println("  Writing input file for SABINE which contains all factors.\n" +
        "    (file: " + outfile_all_factors + ")\n");

    writeInputFile(outfile_all_factors, all_factors);
  }


  public static void main(String[] args) {

    DatabaseParser parser = new DatabaseParser();

    // default options
    String organism = "human";
    String database = "uniprot";
    String base_dir = null;

    /*
     * parse options
     */

    if ((args.length == 1) && (args[0].equals("-help") || args[0].equals("--help"))) {
      parser.usage();
    }

    for(int i=0; i<(args.length-1); i+=2) {

      if(args[i].equals("-o")) { organism			   = args[i+1]; 	continue; }
      if(args[i].equals("-d")) { database 	       = args[i+1]; 	continue; }
      if(args[i].equals("-b")) { base_dir	   		   = args[i+1]; 	continue; }

      if( !args[i].equals("-o") && !args[i].equals("-d") && !args[i].equals("-b") ) {
        System.out.println("\n  Invalid argument: " + args[i]);
        parser.usage();
      }
    }

    /*
     *  check input
     */

    if (! (database.equals("uniprot") ||
        database.equals("ensembl") ||
        database.equals("transfac"))) {

      System.out.println("\n  Invalid argument. Choose between \"uniprot\", \"ensembl\" and \"transfac\".");
      parser.usage();
    }

    if (! (organism.equals("human") ||
        organism.equals("mouse") ||
        organism.equals("rat"))) {

      System.out.println("\n  Invalid argument. Choose between \"human\", \"mouse\" and \"rat\".");
      parser.usage();
    }

    /*
     *  call main function
     */

    parser.runDatabaseParser(organism, database, base_dir);
  }


  private void usage() {

    System.out.println();
    System.out.println("  ---------------------------------------------------------------------");
    System.out.println("  Database Parser for SABINE - StandAlone BINding specificity Estimator");
    System.out.println("  ---------------------------------------------------------------------");
    System.out.println("\n");
    System.out.println("  Usage   : dbparser [OPTIONS]\n");
    System.out.println("  OPTIONS : -o <organism>        (e.g. human, mouse, rat)           default = human");
    System.out.println("            -d <database>        (e.g. uniprot, ensembl, transfac)  default = uniprot");
    System.out.println("            -b <base_directory>  (directory to save the results) \n\n");

    System.exit(0);

  }

}



