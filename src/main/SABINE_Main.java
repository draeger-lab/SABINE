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

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.Calendar;
import java.util.StringTokenizer;

import org.apache.commons.lang3.time.StopWatch;

import extension.PredictionEvaluator;
import model.ModelGenerator;


public class SABINE_Main {


  double best_match_threshold = 0.95;

  int max_number_of_best_matches = 5;

  double outlier_filter_threshold = 0.5;

  boolean silent = false;

  boolean stopTime = false;
  StopWatch stopwatch = new StopWatch();

  boolean dynamic_threshold = true;

  public final static String version = "1.2";
  public final static String appName = "SABINE " + version;


  public static String createBaseDir() {

    String curr_dir = System.getProperty("user.dir") + "/";

    DecimalFormat fmt = new DecimalFormat();
    fmt.setMaximumIntegerDigits(2);
    fmt.setMinimumIntegerDigits(2);

    // get current time and date
    Calendar cal = Calendar.getInstance ();
    String curr_date = (fmt.format(cal.get(Calendar.DAY_OF_MONTH)) + "." +
        fmt.format((cal.get(Calendar.MONTH) + 1)) + "." +
        cal.get(Calendar.YEAR));
    String curr_time = (fmt.format(cal.get(Calendar.HOUR_OF_DAY)) + ":" +
        fmt.format(cal.get(Calendar.MINUTE)) + "_" +
        fmt.format(cal.get(Calendar.MILLISECOND)));

    String base_dir = curr_dir + "tmp/" + curr_date + "_" + curr_time + "/";

    if (!(new File(base_dir).mkdirs())) {
      System.out.println("\nInvalid base directory. Aborting.");
      System.out.println("Base directory: " + base_dir + "\n");
      System.exit(0);
    }
    return(base_dir);
  }

  public static void createTempDirectories(String base_dir) {

    String relpairs_dir, allpairs_dir, libsvm_dir, mmkernel_dir, matlign_dir, mosta_dir, psipred_dir, stamp_dir;

    allpairs_dir = base_dir + "allpairs/";
    relpairs_dir = base_dir + "relevantpairs/";
    libsvm_dir = base_dir + "libsvmfiles/";
    mmkernel_dir = base_dir + "mismatchkernel/";
    matlign_dir = base_dir + "matlign/";
    mosta_dir = base_dir + "mosta/";
    stamp_dir = base_dir + "stamp/";
    psipred_dir = base_dir + "psipred/";


    if ((! new File(allpairs_dir).exists() && ! new File(allpairs_dir).mkdirs()) |
        (! new File(relpairs_dir).exists() && ! new File(relpairs_dir).mkdirs()) |
        (! new File(libsvm_dir).exists() && ! new File(libsvm_dir).mkdirs()) |
        (! new File(mmkernel_dir).exists() && ! new File(mmkernel_dir).mkdirs()) |
        (! new File(matlign_dir).exists() && ! new File(matlign_dir).mkdirs()) |
        (! new File(mosta_dir).exists() && ! new File(mosta_dir).mkdirs()) |
        (! new File(stamp_dir).exists() && ! new File(stamp_dir).mkdirs()) |
        (! new File(psipred_dir).exists() && ! new File(psipred_dir).mkdirs())) {

      System.out.println("\nInvalid base directory. Aborting.");
      System.out.println("Base directory: " + base_dir + "\n");
      System.exit(0);
    }
  }

  /*
   *  print copyright message
   */

  public static void printCopyright() {

    System.out.println("\n-------------------------------------------------");
    System.out.println("SABINE - StandAlone BINding specificity Estimator");
    System.out.println("-------------------------------------------------");
    System.out.println("(version " + version + ")\n");
    System.out.println("Copyright (C) 2009-2013 Center for Bioinformatics T\u00fcbingen (ZBIT),");
    System.out.println("University of T\u00fcbingen, Johannes Eichner.\n");
    System.out.println("This program comes with ABSOLUTELY NO WARRANTY.");
    System.out.println("This is free software, and you are welcome to redistribute it under certain conditions.");
    System.out.print("For details see: ");
    System.out.println("http://www.gnu.org/licenses/gpl-3.0.html\n");
    System.out.println("Third-party software used by this program:");
    System.out.println("  LIBSVM. Copyright (C) 2000-2009 Chih-Chung Chang and Chih-Jen Lin. All rights reserved.");
    System.out.println("  BioJava. Copyright (C) 2008 Richard C. G. Holland. All rights reserved.");
    System.out.println("  Mismatch Kernel. Copyright (C) 2002 Christina Leslie. All rights reserved.");
    System.out.println("  Local Alignment Kernel. Copyright (C) 2005 Jean-Philippe Vert and Hiroto Saigo. All rights reserved.");
    System.out.println("  PSIPRED. Copyright (C) 2000 David T. Jones. All rights reserved.");
    System.out.println("  MoSta. Copyright (C) 2007 Utz J. Pape. All rights reserved.");
    System.out.println("  STAMP. Copyright (C) 2007 Shaun Mahony and Takis Benos. All rights reserved.");
    System.out.println();
  }


  public void launch_SABINE(String infile, String outfile, String verbose_option, String base_dir, String train_dir, String model_dir) {

    /*
     *  set base directory
     */

    if (base_dir == null) {
      base_dir = createBaseDir();
    }
    if (! base_dir.endsWith("/")) {
      base_dir += "/";
    }

    if (base_dir.startsWith("~")) {
      base_dir = System.getProperty("user.home") + base_dir.substring(1);
    }

    createTempDirectories(base_dir);

    try{

      BufferedReader br = new BufferedReader(new FileReader(new File(infile)));
      BufferedWriter bw = new BufferedWriter(new FileWriter(new File(outfile)));

      BufferedReader br_tmp;
      BufferedWriter bw_tmp;

      String line;
      String curr_name;
      StringTokenizer strtok;

      while ((line = br.readLine()) != null) {


        // write current input file
        bw_tmp = new BufferedWriter(new FileWriter(new File(base_dir + "infile.tmp")));

        if (!line.startsWith("NA")) {
          System.out.println("Parse Error. Line is expected to start with \"NA\".");
          System.exit(0);
        }

        strtok = new StringTokenizer(line);
        strtok.nextToken();					// NA
        curr_name = strtok.nextToken();		// name of current TF



        while ((line != null) && !line.startsWith("//")) {    // write current input file and go to next TF

          bw_tmp.write(line + "\n");
          line = br.readLine();
        }
        br.readLine(); 						// XX

        bw_tmp.flush();
        bw_tmp.close();

        // run SABINE on current input file
        String[] arguments = {base_dir + "infile.tmp",
          "-s", "" + best_match_threshold,
          "-m", "" + max_number_of_best_matches,
          "-o", "" + outlier_filter_threshold,
          "-v", verbose_option,
          "-b", base_dir,
          "-t", train_dir,
          "-c", model_dir,
          "-d", "" + dynamic_threshold};

        FBPPredictor.main(arguments);

        // write output file
        bw.write("NA  " + curr_name + "\nXX\n");
        br_tmp = new BufferedReader(new FileReader(new File(base_dir + "prediction.out")));

        while ((line = br_tmp.readLine()) != null) {
          bw.write(line + "\n");
        }
        br_tmp.close();
        bw.write("//\nXX\n");

      }
      br.close();
      bw.flush();
      bw.close();

    }
    catch(IOException ioe) {
      System.out.println(ioe.getMessage());
      System.out.println("IOException occurred while calling SABINE.");
    }
  }


  /*
   *  computes the prediction rate depending on the best match score threshold
   */

  public static void main(String[] args) {

    SABINE_Main caller = new SABINE_Main();

    // run Galaxy Mode
    if ((args.length > 1) && args[0].equals("--galaxy")) {
      String[] newArgs = new String[args.length-1];
      for (int i=1; i<args.length; i++) {
        newArgs[i-1] = args[i];
      }
      SABINE_Galaxy.main(newArgs);

      // run Model Training Mode
    } else if ((args.length >= 1) && args[0].equals("--train")) {
      String[] newArgs = new String[args.length-1];
      for (int i=1; i<args.length; i++) {
        newArgs[i-1] = args[i];
      }
      ModelGenerator.main(newArgs);

      // run Evaluation Mode
    } else if ((args.length >= 1) && args[0].equals("--eval")) {
      String[] newArgs = new String[args.length-1];
      for (int i=1; i<args.length; i++) {
        newArgs[i-1] = args[i];
      }
      PredictionEvaluator.main(newArgs);

      // run Installation Validation Mode
    } else if ((args.length == 1) && (args[0].equals("-check-install") || args[0].equals("--check-install"))) {
      SABINE_Main.printCopyright();
      SABINE_Validator validator = new SABINE_Validator();
      validator.verifyInstallation();

      // run Stand-Alone Mode
    } else {

      SABINE_Main.printCopyright();

      if ((args.length == 0) || ((args.length == 1) && (args[0].equals("-help") || args[0].equals("--help")))) {
        usage();
      }
      if ((args.length == 1) && (args[0].equals("-gui") || args[0].equals("--gui"))) {
        new SABINE_GUI();
        return;
      }

      String infile = args[0];
      String outfile = infile + ".out";
      String verbose_option = "y";
      String stopwatch_option = "n";
      String base_dir = null;
      String train_dir = FBPPredictor.public_trainingset;
      String model_dir = FBPPredictor.defaultModelDir;

      if (infile.endsWith(".input")) {
        outfile = infile.substring(0,infile.length()-5) + "output";
      }


      for(int i=1; i<(args.length-1); i+=2) {

        if(args[i].equals("-s")) { caller.best_match_threshold	  		= Double.parseDouble (args[i+1]);
        caller.dynamic_threshold 			= false;							continue;}
        if(args[i].equals("-m")) { caller.max_number_of_best_matches 	= Integer.parseInt   (args[i+1]); 	continue; }
        if(args[i].equals("-o")) { caller.outlier_filter_threshold		= Double.parseDouble (args[i+1]); 	continue; }
        if(args[i].equals("-f")) { outfile	   		   					= args[i+1]; 						continue; }
        if(args[i].equals("-b")) { base_dir		 						= args[i+1]; 						continue; }
        if(args[i].equals("-t")) { train_dir		 					= args[i+1]; 						continue; }
        if(args[i].equals("-c")) { model_dir							= args[i+1]; 						continue; }
        if(args[i].equals("-v")) { verbose_option	   					= args[i+1]; 						continue; }
        if(args[i].equals("-w")) { stopwatch_option						= args[i+1]; 						continue; }


        if( !args[i].equals("-s") && !args[i].equals("-m") && !args[i].equals("-t") && !args[i].equals("-b") &&
            !args[i].equals("-o") && !args[i].equals("-f") && !args[i].equals("-c") && !args[i].equals("-v") && !args[i].equals("-w")) {

          System.out.println("\n  Invalid argument: " + args[i]);
          usage();
        }
      }

      if (verbose_option.equals("n") || verbose_option.equals("no") || verbose_option.equals("h")) {
        caller.silent = true;
      }
      if (stopwatch_option.equals("y") || stopwatch_option.equals("yes")) {
        caller.stopTime = true;
        caller.stopwatch.start();
      }

      if (! train_dir.endsWith("/")) {
        train_dir += "/";
      }

      /*
       *  call main function
       */

      caller.launch_SABINE(infile, outfile, verbose_option, base_dir, train_dir, model_dir);

      String curr_dir = System.getProperty("user.dir") + "/";

      if (! outfile.startsWith("/")) {
        outfile = curr_dir + outfile;
      }

      if (! caller.silent) {
        System.out.println("\nOutput file: " + outfile + "\n");
      }

      if (caller.stopTime) {
        System.out.println("Time elapsed: " + caller.stopwatch.toString());
      }
    }
  }


  public static void usage() {

    System.out.println("  Usage   : sabine <input_filename> [OPTIONS]\n");
    System.out.println("  OPTIONS : -s <similarity_threshold> (min. FBP-similarity of best matches)     default = 0.95");
    System.out.println("            -m <max_num_best_matches> (max. number of best matches)             default = 5");
    System.out.println("            -o <outlier_filter_param> (max. deviation of a single best match)   default = 0.5");
    System.out.println("            -b <base_dir>             (directory that contains temporary files)");
    //System.out.println("            -t <training_set_dir>  	(directory that contains training sets)        default = data/trainingsets");
    //System.out.println("            -c <model_dir>  	        (directory that contains custom-built models)  default = data/models");
    System.out.println("            -f <output_filename>      (file to save the results)                default = <input_filename>.out");
    System.out.println("            -v <verbose_mode>         (write status to standard output)         default = y (yes)\n\n");
    System.exit(0);

  }
}

