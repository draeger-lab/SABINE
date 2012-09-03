package cv;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Arrays;

import org.apache.commons.cli.BasicParser;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.OptionGroup;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;

public class CVToolRunner {

	private static String old_logfile = "";

	/**
	 * 
	 * @param args
	 */
	@SuppressWarnings("static-access")
	public static void main(String[] args) {
		CVToolRunner cvtool = new CVToolRunner();
		cvtool.runDirectoryMode(args);
	}

	@SuppressWarnings("static-access")
	private void runDirectoryMode(String args[]) {
		CommandLine lvCmd = null;
		HelpFormatter lvFormater = new HelpFormatter();
		CommandLineParser lvParser = new BasicParser();

		Options options = new Options();

		OptionGroup optionGroup = new OptionGroup();
		optionGroup.setRequired(true);
		optionGroup.addOption(OptionBuilder.withLongOpt("directory")
				.isRequired(true).withDescription(
						"directory with LIBSVM compatible files").hasArg(true)
				.create("d"));

		optionGroup.addOption(OptionBuilder.withLongOpt("file")
				.isRequired(true).withDescription("LIBSVM compatible file")
				.hasArg(true).create("f"));

		Option optHelp = new Option("h", "help", false, "shows this help");

		Option optLogCstart = OptionBuilder.withLongOpt("log2c_start")
				.withDescription("log2C parameter grid start point")
				.isRequired(false).hasArg(true).create("cs");

		Option optLogCend = OptionBuilder.withLongOpt("log2c_end")
				.withDescription("log2C parameter grid end point").isRequired(
						false).hasArg(true).create("ce");

		Option optLogCincr = OptionBuilder.withLongOpt("log2c_incr")
				.withDescription("log2C parameter grid point increment")
				.isRequired(false).hasArg(true).create("ci");

		Option optLogEstart = OptionBuilder.withLongOpt("log2e_start")
				.withDescription("log2e parameter grid start point")
				.isRequired(false).hasArg(true).create("es");

		Option optLogEend = OptionBuilder.withLongOpt("log2e_end")
				.withDescription("log2e parameter grid end point").isRequired(
						false).hasArg(true).create("ee");

		Option optLogEincr = OptionBuilder.withLongOpt("log2e_incr")
				.withDescription("log2e parameter grid point increment")
				.isRequired(false).hasArg(true).create("ei");

		Option optLogGstart = OptionBuilder.withLongOpt("log2g_start")
				.withDescription("log2g parameter grid start point")
				.isRequired(false).hasArg(true).create("gs");

		Option optLogGend = OptionBuilder.withLongOpt("log2g_end")
				.withDescription("log2g parameter grid end point").isRequired(
						false).hasArg(true).create("ge");

		Option optLogGincr = OptionBuilder.withLongOpt("log2g_incr")
				.withDescription("log2g parameter grid point increment")
				.isRequired(false).hasArg(true).create("gi");

		Option optRuns = OptionBuilder.withLongOpt("runs").withDescription(
				"number of multiruns for each cross-validation").isRequired(
				false).hasArg(true).create("r");

		Option optFolds = OptionBuilder.withLongOpt("kfolds").withDescription(
				"number of folds (k) in k-fold cross-validation").isRequired(
				false).hasArg(true).create("k");

		Option optSVMCache = OptionBuilder.withLongOpt("cache")
				.withDescription("cache in MByte for LIBSVM").isRequired(false)
				.hasArg(true).create("m");

		Option optRegressionThreshold = OptionBuilder.withLongOpt("rthreshold")
				.withDescription("maximum number of labels for classification")
				.isRequired(false).hasArg(true).create("th");

		Option optGridInfo = OptionBuilder.withLongOpt("grid").withDescription(
				"displays infos of grid sizes").isRequired(false).create("g");

		options.addOption(optHelp);
		options.addOptionGroup(optionGroup);
		options.addOption(optLogCend);
		options.addOption(optLogCincr);
		options.addOption(optLogCstart);
		options.addOption(optLogEend);
		options.addOption(optLogEincr);
		options.addOption(optLogEstart);
		options.addOption(optLogGend);
		options.addOption(optLogGincr);
		options.addOption(optLogGstart);
		options.addOption(optRuns);
		options.addOption(optFolds);
		options.addOption(optSVMCache);
		options.addOption(optGridInfo);
		options.addOption(optRegressionThreshold);

		/**
		 * parse command line
		 */
		try {
			lvCmd = lvParser.parse(options, args);

			if (lvCmd.hasOption('h')) {
				lvFormater.printHelp("", options);
				return;
			}

		} catch (ParseException pvException) {
			lvFormater.printHelp("CVTool", options);
			System.out.println("Parse error: " + pvException.getMessage());
			return;
		}

		/**
		 * If required, override parameters
		 */
		String dir = null;
		// -f,--file <arg> LIBSVM compatible file
		String file = null;
		dir = lvCmd.getOptionValue("d");
		file = lvCmd.getOptionValue("f");
		// -ce,--log2c_end <arg> log2C parameter grid end point
		if (lvCmd.hasOption("ce")) {
			CVGlobalSettings.log2c_end = new Integer(lvCmd.getOptionValue("ce"));
		}
		// -ci,--log2c_incr <arg> log2C parameter grid point increment
		if (lvCmd.hasOption("ci")) {
			CVGlobalSettings.log2c_incr = new Integer(lvCmd
					.getOptionValue("ci"));
			if (CVGlobalSettings.log2c_incr == 0){
				System.err.println("Increment of log2C grid is zero!");
				System.exit(1);
			}
		}
		// -cs,--log2c_start <arg> log2C parameter grid start point
		if (lvCmd.hasOption("cs")) {
			CVGlobalSettings.log2c_start = new Integer(lvCmd
					.getOptionValue("cs"));
		}
		if (lvCmd.hasOption("ee")) {
			// -ee,--log2e_end <arg> log2e parameter grid end point
			CVGlobalSettings.log2e_end = new Integer(lvCmd.getOptionValue("ee"));
		}
		// -ei,--log2e_incr <arg> log2e parameter grid point increment
		if (lvCmd.hasOption("ei")) {
			CVGlobalSettings.log2e_incr = new Integer(lvCmd
					.getOptionValue("ei"));
			if (CVGlobalSettings.log2e_incr == 0){
				System.err.println("Increment of log2e grid is zero!");
				System.exit(1);
			}
		}
		// -es,--log2e_start <arg> log2e parameter grid start point
		if (lvCmd.hasOption("es")) {
			CVGlobalSettings.log2e_start = new Integer(lvCmd
					.getOptionValue("es"));
		}
		// -ge,--log2g_end <arg> log2g parameter grid end point
		if (lvCmd.hasOption("ge")) {
			CVGlobalSettings.log2g_end = new Integer(lvCmd.getOptionValue("ge"));
		}
		// -gi,--log2g_incr <arg> log2g parameter grid point increment
		if (lvCmd.hasOption("gi")) {
			CVGlobalSettings.log2g_incr = new Integer(lvCmd
					.getOptionValue("gi"));
			if (CVGlobalSettings.log2g_incr == 0){
				System.err.println("Increment of log2g grid is zero!");
				System.exit(1);
			}

		}
		// -gs,--log2g_start <arg> log2g parameter grid start point
		if (lvCmd.hasOption("gs")) {
			CVGlobalSettings.log2g_start = new Integer(lvCmd
					.getOptionValue("gs"));
		}
		// -k,--kfolds <arg> number of folds (k) in k-fold cross-validation
		if (lvCmd.hasOption("k")) {
			CVGlobalSettings.folds = new Integer(lvCmd.getOptionValue("k"));
		}
		// -m,--cache <arg> cache in MByte for LIBSVM
		if (lvCmd.hasOption("m")) {
			CVGlobalSettings.svm_cache = new Integer(lvCmd.getOptionValue("m"));
		}
		// -r,--runs <arg> number of multiruns for each cross-validation
		if (lvCmd.hasOption("r")) {
			CVGlobalSettings.runs = new Integer(lvCmd.getOptionValue("r"));
		}
		// classification threshold
		if (lvCmd.hasOption("th")) {
			CVGlobalSettings.regression_threshold = new Integer(lvCmd
					.getOptionValue("th"));
		}
		if (lvCmd.hasOption("g")) {
			GridSizeTester.main(args);
		}

		if (file != null) {
			String str_file = null;
			try {
				str_file = (new File(file)).getCanonicalPath();
			} catch (IOException e) {
				e.printStackTrace();
			}
			String str_filename = (new File(file)).getName();
			String logfile = str_filename + ".dat";
			FileWriter fw = null;
			try {
				fw = new FileWriter(logfile);
			} catch (IOException e1) {
				e1.printStackTrace();
			}
			if (str_filename.endsWith(".matrix")) {
				runCV(fw, "", str_file, 0);
			}
			if (str_filename.endsWith(".att")) {
				runCV(fw, "", str_file, 1);
			}
			try {
				fw.close();
			} catch (IOException e) {
				e.printStackTrace();
			}
		}

		if (dir != null) {
			DirectoryManager dm = new DirectoryManager();
			String[] str_files = dm.getDirectoryContent(new File(dir));
			if (str_files == null) {
				System.err.println(new File(dir).getPath());
				System.err.println("No files found, check directory!");
				System.exit(1);
			}

			Arrays.sort(str_files);
			FileWriter fw = null;

			for (int i = 0; i < str_files.length; i++) {

				/**
				 * use online defined LIBSVM compatible files
				 */
				if (str_files[i].endsWith(".matrix")|| str_files[i].endsWith(".dat")) {
					/**
					 * determine *.dat file name
					 */
					String logfile = null;
					try {
						logfile = str_files[i].substring(0, str_files[i].indexOf(".lp"))+ ".dat";
					} catch (Exception e) {
						continue;
					}

					if (!old_logfile.equals(logfile) && logfile != null) {
						old_logfile = logfile;
						try {
							if (fw != null) {
								fw.close();
							}
							fw = new FileWriter(logfile);
						} catch (IOException e1) {
							e1.printStackTrace();
						}
					}
					if(fw != null){
					if (str_files[i].endsWith(".matrix")) {
						runCV(fw, dir, str_files[i], 0);
					}

					if (str_files[i].endsWith(".att")) {
						runCV(fw, dir, str_files[i], 1);
					}}
					System.out.println("");
				}
			}
			try {
				if (fw != null) {
					fw.close();
				}
			} catch (IOException e) {
				e.printStackTrace();
				System.exit(1);
			}
		}
	}

	private static void runCV(FileWriter fw, String dir, String lsvmfile,
			int vmode) {
		/**
		 * determine problem
		 */
		boolean svreg = false;
		boolean isempty = false;
		File f = new File(dir + "//" + lsvmfile);
		try {
			isempty = CVHelper.isEmptyFile(f);
			if (!isempty) {
				svreg = CVHelper.isARegressionProblem(f);
			}
		} catch (IOException e) {
			e.printStackTrace();
		}

		if (!isempty) {
			if (!svreg) {
				// numerical kernel or precomputed kernel?
				if (Integer.valueOf(vmode) == 1) {
					SVClassificationNumerical.performClassification(lsvmfile,
							dir, fw);
				} else {
					SVClassification.performClassification(lsvmfile, dir, fw);
				}
			} else {
				if (Integer.valueOf(vmode) == 1) {
					SVEpsilonRegressionNumerical.performEpsilonRegression(
							lsvmfile, dir, fw);
				} else {
					SVEpsilonRegression.performEpsilonRegression(lsvmfile, dir,
							fw);
				}
			}
		}
	}
}
