package proBoundTools;

import java.io.File;
import java.net.URISyntaxException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.security.CodeSource;
import java.util.ArrayList;

import org.apache.commons.cli.*;

public class App {
	
	public static void main(String[] args) {
		
		//Determines if internal commands are available 
		boolean iv=false;


		//Gets the schema file.
		String generalSchemaFile = null;
		CodeSource codeSource = App.class.getProtectionDomain().getCodeSource();
		try {
			File jarFile = new File(codeSource.getLocation().toURI().getPath());
			String jarDir = jarFile.getParentFile().getPath();
			generalSchemaFile = jarDir+"/../../config/schema.general.json";
		} catch (URISyntaxException e1) {
			e1.printStackTrace();
			System.exit(1);
		}

		if(!Files.exists(Paths.get(generalSchemaFile))) {
			System.err.println("The general schema file '"+generalSchemaFile+"' cannot be found: ");
			System.exit(1);
		}
		
		///////////////////////////
		// Parses the arguments. //
		///////////////////////////
		CommandLineParser parser = new DefaultParser();
		Options options          = new Options();
		Option commandOption     = new Option("c", "command", true, "Command to be executed (required).");
		Option helpOption        = new Option( "help",    "prints help" );
		Option verboseOption     = new Option( "verbose", "be extra verbose" );
		options.addOption(commandOption);
		options.addOption(helpOption);
		options.addOption(verboseOption);
		HelpFormatter formatter  = new HelpFormatter();
		CommandLine cmd;
		try {
			cmd = parser.parse(options, args);
		} catch (ParseException e) {
			System.err.println(e.getMessage());
			System.out.println("");
			printUsage(new Toolbox(generalSchemaFile, false), formatter, options, iv);
			System.exit(1);
			return;
		}
		//Retrieves the arguments
		boolean verbose      = cmd.hasOption("verbose");
		boolean help         = cmd.hasOption("help");
		String commandString = cmd.getOptionValue("command");
		
		//Parses the help
		if(help) {
			printUsage(new Toolbox(generalSchemaFile, false), formatter, options, iv);
			System.exit(0);
		}
		
		//Parses the command string
		ArrayList<String> functionNames        = new ArrayList<String>();
		ArrayList<ArrayList<String>> arguments = new ArrayList<ArrayList<String>>();
		boolean commandOK                      = commandString==null ? false : parseCommandString(commandString, functionNames, arguments);
		if((!commandOK) || functionNames.size()==0) {
			System.err.println("ERROR: No valid command specified (using -c).");
			System.out.println("");
			printUsage(new Toolbox(generalSchemaFile, false), formatter, options, iv);
			System.exit(1);
		}
			
		if(verbose) {
			System.out.println("Starting ProBoundTools.");
			System.out.println("=======================");
			System.out.println("command = '"+commandString+"'");
		}
		
		////////////////////////////
		// Executes API commands. //
		////////////////////////////
		Toolbox t = new Toolbox(generalSchemaFile, verbose);
		for(int iF=0; iF<functionNames.size(); iF++) {
			if(verbose) {
				System.out.println("");
				//Formating the current function calls.
				System.out.print(functionNames.get(iF) + "(");
				for(int iA=0; iA<arguments.get(iF).size(); iA++) {
					System.out.print(arguments.get(iF).get(iA));
					if(iA < arguments.get(iF).size()-1)
						System.out.print(", ");
				}
				System.out.println(")");
				System.out.println("---------------------------------------------");
			}
			
			
			//Calls functions.
			String func           = functionNames.get(iF);
			ArrayList<String> arg = arguments.get(iF);
			
			if(     func.equals("inputTableTSV"))                t.inputTableTSV(arg);           
			else if(func.equals("inputFastq"))                   t.inputFastq(arg);              
			else if(func.equals("inputTXT"))                     t.inputTXT(arg);                
			else if(func.equals("loadScoringModel"))             t.loadScoringModel(arg);        
			else if(func.equals("loadFitLine"))                  t.loadFitLine(arg);            
			else if(func.equals("loadMotifCentralModel"))        t.loadMotifCentralModel(arg);
			else if(func.equals("writeModel"))                   t.writeModel(arg);             
			else if(func.equals("setMismatchGauge"))             t.setMismatchGauge(arg);        
			else if(func.equals("buildConsensusModel"))          t.buildConsensusModel(arg);    
			else if(func.equals("addNScoring"))                  t.addNScoring(arg);
			else if(func.equals("setAlphabet"))                  t.setAlphabet(arg);
			else if(func.equals("guessAlphabet"))                t.guessAlphabet(arg);
			else if(func.equals("removeBindingMode"))            t.removeBindingMode(arg);       
			else if(func.equals("selectBindingMode"))            t.selectBindingMode(arg);       
			else if(func.equals("removeInteraction"))            t.removeInteraction(arg);       
			else if(func.equals("selectAndCleanBindingMode"))    t.selectAndCleanBindingMode(arg);
			else if(func.equals("computeMonoInformation"))       t.computeMonoInformation(arg); 
			else if(func.equals("affinitySum"))                  t.affinitySum(arg);            
			else if(func.equals("bindingModeScores"))            t.bindingModeScores(arg);      
			else if(func.equals("kMerCount"))                    t.kMerCount(arg);               
			else if(func.equals("kMerAverage"))                  t.kMerAverage(arg);             
			else if(func.equals("predictCountTable"))            t.predictCountTable(arg);       
			else if(func.equals("positionBaseCount"))            t.positionBaseCount(arg);       
			else {
				System.out.println("Invalid function: "+func+"\n");
				printUsage(t, formatter, options, iv);
				System.exit(1);
			}
		}
		
		/////////////////////
		// END THE PROGRAM //
		/////////////////////
		t.close();
		
		if(verbose)
			System.out.println("\n=> END OF PROGRAM. <=");
		System.exit(0);
	}


	public static void error(String s) {
		System.err.println("ERROR: "+s);
		System.exit(1);
	}

	public static void warning(String s) {
		System.err.println("WARNING: "+s);
	}

	public static void printUsage(Toolbox t, HelpFormatter formatter, Options options, boolean iv) {
		
		System.out.println("OVERVIEW:");
		System.out.println("=========");
		formatter.printHelp("utility-name", options);
		System.out.println("");

		
		System.out.println("COMMAND STRING:");
		System.out.println("===============");
		System.out.println("The command string consists of a sequence of function calls separated by '.':");
		System.out.println("function1(var1,...,varN).function2(var1,...,varN)...");
		System.out.println("Each function calls as a numer of aguments that can be enclosed by ' or \".");
		System.out.println("");
		System.out.println("Functions for loading sequence data:");
		System.out.println("------------------------------------");
		t.inputTableTSV(   (ArrayList<String>)null);
		t.inputFastq(      (ArrayList<String>)null);
		t.inputTXT(        (ArrayList<String>)null);
		System.out.println("");
		System.out.println("Functions for loading binding models:");
		System.out.println("-------------------------------------");
		t.loadScoringModel(          (ArrayList<String>)null);
		t.loadFitLine(               (ArrayList<String>)null);
		t.loadMotifCentralModel(     (ArrayList<String>)null);
		
		System.out.println("");
		System.out.println("Functions for manipulating binding models:");
		System.out.println("------------------------------------------");
		t.writeModel(                    (ArrayList<String>)null);
		t.setMismatchGauge(              (ArrayList<String>)null);
		t.buildConsensusModel(           (ArrayList<String>)null);
		t.addNScoring(                   (ArrayList<String>)null);
		t.setAlphabet(                   (ArrayList<String>)null);
		t.guessAlphabet(                 (ArrayList<String>)null);
		t.removeBindingMode(             (ArrayList<String>)null);
		t.selectBindingMode(             (ArrayList<String>)null);
		t.removeInteraction(             (ArrayList<String>)null);
		t.selectAndCleanBindingMode(     (ArrayList<String>)null);
		t.computeMonoInformation(        (ArrayList<String>)null);
		System.out.println("");
		System.out.println("Functions for scoring sequences using binding modes:");
		System.out.println("----------------------------------------------------");
		t.affinitySum(              (ArrayList<String>)null);
		t.bindingModeScores(        (ArrayList<String>)null);
		t.predictCountTable(        (ArrayList<String>)null); 
		
		System.out.println("");
		System.out.println("Functions for charecterizing the input sequences/tables:");
		System.out.println("--------------------------------------------------------");
		t.kMerCount(                (ArrayList<String>)null);
		t.kMerAverage(              (ArrayList<String>)null);
		t.positionBaseCount(        (ArrayList<String>)null);
		
		System.out.println("");
		System.out.println("EXAMPLES:");
		System.out.println("=========");
		System.out.println("Command for loading the scoring model with ID=1000 from motifcentral.org and scoring sequences:");
		System.out.println("loadScoringModelDB(1000).inputTXT(sequenceFile.txt).bindingModeScores(outFile.tsv)");
		System.out.println("");
		System.out.println("Command for extracting the first binding mode from a ProBound fit and saving it as a scoring model:");
		System.out.println("loadFitLine(fit.json).setMismatchGauge().selectAndCleanBindingMode(1).writeModel(scoringModel.json)");
		System.out.println("");
	}
	
	//Parses the command string into a list of function calls and arguments.
	public static boolean parseCommandString(String commandString, ArrayList<String> functionNames, ArrayList<ArrayList<String>> arguments) {

		////////////////////////////////
		// Parses the command string. //
		////////////////////////////////
		// Parsing states:
		// 0: reads function name
		// 1: reads arguments.
		// 2: State after ")", waiting for ".".
		int parseState      = 0; 
		String currentF     = "";
		ArrayList<String> a = null;
		String currentArg   = "";
		boolean inQuote     = false;
		Character quoteSign = null;
		for(int i=0; i<commandString.length(); i++) {
			Character c = commandString.charAt(i);
			
			//Checks of we we currently are in the beginning or end of a "" or ''
			if(c.equals('"') || c.equals('\'')) {
				if(inQuote == false) {
					inQuote = true;
					quoteSign = c;
					continue;
				} else {
					if(quoteSign == c) {
						inQuote = false;
						continue;
					}
				}
			}
			
			//Removes blank spaces unless inQuote=true
			if(!inQuote && c.equals(' '))
				continue;
			
			//Transitions between parsing states, saves strings.
			if(parseState==0) {
				if(!inQuote && c.equals('(')) {
					functionNames.add(currentF);
					parseState = 1;
					a = new ArrayList<String>();
					currentArg = "";
					continue;
				} else {
					currentF += c;
					continue;
				}
			} else if(parseState==1) {
				if(!inQuote && c.equals(',')) {
					a.add(currentArg);
					currentArg = "";
					continue;
				} else if(!inQuote && c.equals(')')) {
					if(!currentArg.equals("")) {
						a.add(currentArg);
					}
					arguments.add(a);
					parseState = 2;
					continue;
				} else {
					currentArg += c;
					continue;
				}

			} else if(parseState==2) {
				if(!inQuote && c.equals('.')) { 
					parseState = 0;
					currentF   = "";
					continue;
				} else {
					System.err.println("ERROR: Did expected '.' on position "+i+" of argument. Found '"+c+"'.");
				}
			} 
		}
		//Checking so the parser did not end in the argument-reading state. 
		if(inQuote) {
			System.err.println("ERROR: When parsing command string, a substring enclosed by the character "+ quoteSign +" did not end.");
			System.out.println("");
			return false;
		}
		if(parseState == 1) {
			System.err.println("ERROR: The last function in the command ("+currentF+") did not complete.");
			System.out.println("");
			return false;
		}
		return true;
	}

}