# SABINE
<img align="right" src="src/resources/img/sabine_logo.png"/>

**StandAlone BINding specificity Estimator**

[![License (GPL version 3)](https://img.shields.io/badge/license-GPLv3.0-blue.svg?style=plastic)](http://opensource.org/licenses/GPL-3.0)
[![Latest version](https://img.shields.io/badge/Latest_version-1.2-brightgreen.svg?style=plastic)](https://github.com/draeger-lab/SABINE/releases/)
[![DOI](http://img.shields.io/badge/DOI-10.1371%20%2F%20journal.pone.0082238-blue.svg?style=plastic)](http://dx.doi.org/10.1371%2Fjournal.pone.0082238)

*Authors:* [Johannes Eichner](https://github.com/jeichner), [Adrian Schröder](http://www.cogsys.cs.uni-tuebingen.de/mitarb/schroeder/), [Andreas Dräger](https://github.com/draeger/), Jonas Eichner, [André Hennig](https://github.com/AndreHennig), Florian Topf, Dierk Wanke, [Klaus Harter](http://www.zmbp.uni-tuebingen.de/plant-physiol/research-groups/harter.html), [Andreas Zell](https://github.com/ZellTuebingen)

____________________________________________________________________________________________________________________________  

### Short description
SABINE is a tool to predict the binding specificity of a transcription factor (TF), given its amino acid sequence, species, structural superclass and DNA-binding domains. For convenience, the superclass and DNA-binding domains of a given TF can be predicted based on sequence homology with TFs in the training set of SABINE. Alternatively, the tool [TFpredict](https://github.com/draeger-lab/TFpredict), which predicts all structural characteristics of TFs required by SABINE, can be employed in an additional preprocessing step. SABINE compares a given factor to a predefined set of TFs of the same superclass for which experimentally confirmed position frequency matrices (PFM) are available. Based on various features capturing evolutionary, structural and physicochemical similarity of the DNA-binding domains, the PFM similarity is predicted utilizing support vector regression. The TFs with highest PFM similarity to the factor of interest is reported, and their PFMs are merged using STAMP to generate the predicted consensus PFM.

## Publications

Article citations are **critical** for us to be able to continue support for SABINE.  If you use SABINE and you publish papers about work that uses SABINE, we ask that you **please cite the SABINE papers**.

<dl>
  <dt>Original method paper:</dt>
  <dd>Adrian Schröder, Johannes Eichner, Jochen Supper, Jonas Eichner, Dierk Wanke, Carsten Henneges, and Andreas Zell. <a href="http://dx.doi.org/10.1371%2Fjournal.pone.0013876">Predicting DNA-Binding Specificities of Eukaryotic Transcription Factors</a>. PLoS ONE , 5(11):e13876, November 2010. [ <a href="http://dx.doi.org/10.1371%2Fjournal.pone.0013876">DOI</a> | <a href="http://journals.plos.org/plosone/article/file?id=10.1371/journal.pone.0013876&type=printable">PDF</a> ]</dd>
  <dt>Research Article:</dt>
  <dd>Johannes Eichner, Florian Topf, Andreas Dräger, Clemens Wrzodek, Dierk Wanke, and Andreas Zell. <a href="http://dx.doi.org/10.1371%2Fjournal.pone.0082238">TFpredict and SABINE: Sequence-Based Prediction of Structural and Functional Characteristics of Transcription Factors</a>. PLoS ONE, 8(12):e82238, December 2013.
  [ <a href="http://dx.doi.org/10.1371/journal.pone.0082238">DOI</a> | <a href="http://www.plosone.org/article/fetchObject.action?uri=info%3Adoi%2F10.1371%2Fjournal.pone.0082238&representation=PDF">PDF</a> ]  
  </dd>
</dl>

## Table of Contents:
  - [Installation](#installation)
  - [Manual](#manual)
  - [Format specification](#format-specification)
  - [License](#license)
  - [Integrated data](#integrated-data)
  - [Contact](#contact)
____________________________________________________________________________________________________________________________  

## Installation

To extract the gzipped tar archive of SABINE obtained from our downloads section use the command:

    tar -xzf sabine.tar.gz

For convenience, a shell script was implemented to simplify the installation of SABINE. Start the installation script with the following command:

    install.sh

The script will install SABINE and all required third-party software packages and libraries on your system. To test if the program is working properly, you can validate your installation using the command:

    sabine.sh -check-install

  
**Note:** Please ensure that the tcsh shell is installed on your system, as it is required by the tool `PSIPRED`, which is employed by SABINE to predict secondary structures. 

## Manual

**DESCRIPTION:**
The tool predicts the binding specificity of a transcription factor (TF) in terms of a position frequency matrix (PFM), given the amino acid sequence, DNA-binding domains, superclass and species of the TF of interest.
              
The given factor (query factor) is compared to factors of the same superclass (training factors) for which experimentally determined PFMs are available. Based on diverse alignment- and kernel-based similarity measures comparing the DNA-binding domains of the query factor and the training factors, the PFM-similarity of each training factor to the query factor is estimated using a support-vector-machine-based regression model. 
              
The tool reports a set of best matches, i.e., training factors for which a PFM-similarity greater than a predefined similarity threshold (see OPTIONS: `-s`) was predicted. The PFMs of those best matches are in turn filtered to remove outliers (see OPTIONS: `-o`) and then merged using `STAMP` to generate the predicted PFM. The maximal number of PFMs that shall be merged can be specified by the user (see OPTIONS: `-m`). 

**INPUT:**    SABINE input file contains the query factors in the SABINE input file format (see SABINE Format Specification)

**OUTPUT:**   SABINE output file
* contains the predicted best matches and PFMs
* the output filename can be specified by the user (see **OPTIONS:** `-f`)

**USAGE:**

    sabine.sh <input_filename> [OPTIONS]

**OPTIONS:**

| Option                     | Meaning                                  | Default Value                    |
|----------------------------|------------------------------------------|----------------------------------|
|`-s <similarity_threshold>` | min. FBP-similarity of best matches      | `default = 0.95`                 |
|`-m <max_num_best_matches>` | max. number of best matches              | `default = 5`                    |
|`-o <outlier_filter_param>` | max. deviation of a single best match    | `default = 0.5`                  |
| `-b <base_dir>`            | directory that contains temporary files  |                                  |
|`-f <output_filename>`      | file to save the results                 | `default = <input_filename>.out` |
|`-v <verbose_mode>`         | write status to standard output          | `default = y (yes)`              |

A more detailed description of the progam including usage examples can be found below in the [Users' Guide](#users-guide)

## Format specification

**SABINE input file:**
```
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
```

**SABINE output file:**
```
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
```

## License
<img align="right" src="https://www.gnu.org/graphics/gplv3-127x51.png"/>

This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program; if not, see <http://www.gnu.org/licenses/>.

### Included third-party software

* The Local Alignment Kernel by Hiroto Saigo et al. was integrated in SABINE as a sequence feature to measure the similarity of the DNA-binding domains. Website: http://sunflower.kuicr.kyoto-u.ac.jp/~hiroto/project/homology.html
* We included the Mismatch Kernel by Christina Leslie et al. as an additional feature which incorporates the domain similarity of transcription factors. Website: http://cbio.mskcc.org/leslielab/software/string-kernels
* To capture the structural similarity of TF binding domains we applied the PSIPRED secondary structure prediction method by David Jones *et al.* Website: http://bioinf.cs.ucl.ac.uk/psipred/
* We integrated the Motif statistic software suite (MoSta) by Utz Pape *et al.* to measure the similarity of PFMs. Website: http://mosta.molgen.mpg.de/
* In order to generate familial binding profiles by merging PFMs, we applied the tool STAMP by Shaun Mahony *et al.* Website: http://www.benoslab.pitt.edu/stamp/
* For training and evaluation of the support vector regression models we used the libSVM implementation provided by Chih-Chung Chang *et al.* Website: http://www.csie.ntu.edu.tw/~cjlin/libsvm/
* We computed pairwise sequence alignment scores and generated sequence logos using the open source framework BioJava by Richard Holland *et al.* Website: http://biojava.org/wiki/Main_Page

### Integrated data

The basis of our supervised machine learning based approach to predicting DNA-binding specificities of TFs was the generation of a non-redundant training data set. We restricted the sources to databases providing experimentally validated DNA-binding specificity information in terms of DNA-binding sites, consensus sequences or PFMs. Besides intergrating large databases spanning the whole eucaryotic kingdom, we extracted data from diverse smaller databases whose content is specific to particular organisms. An overview of our data sources can be found in the table below.

| Database  | URL                                                              |
|-----------|------------------------------------------------------------------|
| TRANSFAC  | http://www.biobase-international.com/pages/index.php?id=transfac |
| JASPAR    | http://jaspar.cgb.ki.se/                                         |
| YEASTRACT | http://www.yeastract.com/                                        |
| SCPD      | http://rulai.cshl.edu/SCPD/                                      |
| AGRIS     | http://arabidopsis.med.ohio-state.edu/                           |
| FLYREG    | http://www.flyreg.org/                                           |

## Users' Guide

### Introduction
The key features of gene-regulatory networks are the interconnections between specific transcription factors (TF) and cis-regulatory elements of the DNA. These connections function as an interface between signaling pathways and the regulation of gene expression. After decades of intensive research effort only a small fraction of these connections is known. A major and labour-intensive part of this effort is the characterization of the DNA-binding specificity of TFs whose DNA-binding domain enables the specific recognition of short DNA motives in the promoter region of their proximal target genes. 

To further increase our knowledge about the specific interactions between TFs and cis-regulatory elements, we designed an algorithm which allows for predicting the binding specificity of transcription factors with high accuracy. Apparently, the protein sequence and structure of the DNA-binding domain determines its function, which is the molecular recognition and binding to a defined set of DNA motives. Thus we approach to estimate the binding specificity of TFs based on structural, physicochemical and properties of their DNA-binding domains. Employing support vector regression we estimate the similarity of the binding specificities of two TFs based on diverse features incorporating domain sequence similarity, secondary structure and phylogenetic distance. This approach provides a quantitative measure for the functional similarity of two factors and enables the transfer of DNA-binding specificity data with low error.

### Example screenshots

Will be added soon.

## Contact

If you have any further questions, please contact [Johannes Eichner](https://github.com/jeichner).
