
  -------------------------------------------------
  SABINE - StandAlone BINding specificity Estimator
  -------------------------------------------------
  (version 1.2, Copyright (C) 2012 Johannes Eichner)


  Contents:
 _________________________
  - License
  - Installation
  - Manual
  - Format specification
  - Website and questions
____________________________________________________________________________________________________________________________  

  -------
  License
  -------

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 3 of the License, or (at
  your option) any later version.

  This program is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
  General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, see <http://www.gnu.org/licenses/>.

____________________________________________________________________________________________________________________________

  ------------
  Installation
  ------------

  To extract the gzipped tar archive of SABINE obtained from our downloads section use the command:

    tar -xzf sabine.tar.gz

  For convenience, a shell script was implemented to simplify the installation of SABINE. 
  Start the installation script with the following command:

    install.sh

  The script will install SABINE and all required third-party software packages and libraries on your system.
  To test if the program is working properly, you can validate your installation using the command:

    sabine.sh -check-install

  
  Note: Please ensure that the tcsh shell is installed on your system, as it is required by the tool PSIPRED, 
        which is employed by SABINE to predict secondary structures. 
  ____________________________________________________________________________________________________________________________


  ------
  Manual
  ------

  DESCRIPTION:  The tool predicts the binding specificity of a transcription factor (TF) in terms of a position frequency 
		matrix (PFM), given the amino acid sequence, DNA-binding domains, superclass and species of the TF of 
		interest.  
		The given factor (query factor) is compared to factors of the same superclass (training factors) 
		for which experimentally determined PFMs are available. Based on diverse alignment- and kernel-based 
		similarity measures comparing the DNA-binding domains of the query factor and the training factors, the 
		PFM-similarity of each training factor to the query factor is estimated using a support-vector-machine-based 
		regression model. 
		The tool reports a set of best matches, i.e. training factors for which a PFM-similarity greater than a 
		predefined similarity threshold (see OPTIONS: -s) was predicted. The PFMs of those best matches are in turn 
		filtered to remove outliers (see OPTIONS: -o) and then merged using STAMP to generate the predicted PFM. 
		The maximal number of PFMs that shall be merged can be specified by the user (see OPTIONS: -m). 

  INPUT:  	SABINE input file 
  			contains the query factors in the SABINE input file format (see SABINE Format Specification)
  		
  OUTPUT:	SABINE output file
  	  		contains the predicted best matches and PFMs
			the output filename can be specified by the user (see OPTIONS: -f)
  
  USAGE:	sabine.sh <input_filename> [OPTIONS]

  OPTIONS : 	-s <similarity_threshold> (min. FBP-similarity of best matches)      default = 0.95
            	-m <max_num_best_matches> (max. number of best matches)              default = 5
            	-o <outlier_filter_param> (max. deviation of a single best match)    default = 0.5
				-b <base_dir>             (directory that contains temporary files)  
            	-f <output_filename>      (file to save the results)                 default = <input_filename>.out
            	-v <verbose_mode>         (write status to standard output)          default = y (yes)


  A full documentation including a tutorial is available at the supplementary website of SABINE:
  http://www.cogsys.cs.uni-tuebingen.de/software/SABINE/ 

____________________________________________________________________________________________________________________________

  --------------------
  Format specification
  --------------------


  SABINE input file:
  __________________________________________________________________________________________
  NA  Identifier
  XX
  SP  Organism
  XX
  RF  Reference to UniProt (optional)
  XX
  CL  Classification (class acc. no. or decimal classification no. as in TRANSFAC)
  XX
  S1  Amino acid sequence
  XX
  S2  Alternative amino acid sequence (optional)
  XX
  FT  DNA-binding domain (domain ID   start position   end position)
  XX
  //
  XX


  SABINE output file:
  __________________________________________________________________________________________
  NA  Identifier
  XX
  BM  Best match (transcription factor ID   PFM similarity score) 
  XX
  MA  A   C   G   T   rows:	positions within the aligned sequences
  MA		  first column: position index
  MA               columns: relative frequencies of A, C, G, T residues
  MA           last column: consensus sequence in IUPAC code
  XX
  //
  XX

  ____________________________________________________________________________________________________________________________

  ---------------------
  Website and questions
  ---------------------

  To obtain more detailed information about SABINE, see the website: http://www.cogsys.cs.uni-tuebingen.de/software/SABINE/ 
  If you have any further questions, please contact me by e-mail:    johannes.eichner@uni-tuebingen.de

