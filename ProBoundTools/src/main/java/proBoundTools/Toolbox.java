package proBoundTools;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.net.MalformedURLException;
import java.net.URL;
import java.net.URLConnection;
import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.ResultSet;
import java.sql.Statement;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Random;
import java.util.stream.Collectors;
import java.util.zip.GZIPInputStream;

import javax.management.RuntimeErrorException;

import org.json.*;

import Jama.EigenvalueDecomposition;
import Jama.Matrix;
import base.Array;
import modelOptimization.CombinedLikelihood;
import modelOptimization.LikelihoodOptimizer;
import proBoundTools.JSONModel;
import proBoundTools.Misc;
import sequenceTools.*;
import modelComponents.*;

public class Toolbox {

	//Mode parameters
	JSONObject oModel;
	ArrayList<Integer> nFrames;
	ArrayList<Integer> iBindingMode;
	int nModes;
	ArrayList<Integer> size, dInt;
	ArrayList<SlidingWindow> sw;
	int nPredictors;
	boolean mismatchGauge;

	//Information about the sequence input stream.
	private String sequenceFile;         // Sequence input file
	private String inputFormat;          // Indicates if the input is "table" or "fasta"
	private BufferedReader br;           // Input file stream.
		
	//Current entry in the sequence file
	private int stringColumn;           // Index of string column.
	public String sequence;             // The currently loaded
	private int firstNeededCharacter;   // Index of the last processed base in the string.
								        //This should be increased by the scoring function as the 
	                                    // some bases are not needed. It should be decreased as the next
	                                    // substring is read.
	double[] tableData;                 // Current table data.
	boolean verbose;
	
	
	//ProBound objects
	MultiRoundData dataTable;
	EnrichmentModel enrichmentModel;
	ArrayList<BindingMode> bindingModes;
	ArrayList<BindingModeInteraction> interactions;
	ArrayList<ModelComponent> componentList;
	CountTable countTable;
	
	
	//List of k-mers
	int k;
	ArrayList<String> kMers, kMersRC;           //List of k-mers and reverse-complement in the k-mer
	Map<String,Integer> kMerIndex, kMerIndexRC; //Lookup to map k-mer to index.
	int nKmers;
	
	//Alphabet
	LongSequence.SequenceClass sc;
	public String letterComplement, sortedLetters;
	public String letterOrder;
	int nMono, nDi;
	
	String generalSchemaFile;

	// Constructor
	//////////////
	Toolbox(String generalSchemaFileIn, boolean verboseIn) {
		verbose           = verboseIn;
		mismatchGauge     = false;
		generalSchemaFile = generalSchemaFileIn;
	}
	

	public static void error(String s) {
		System.err.println("ERROR: "+s);
		System.exit(1);
	}

	public static void warning(String s) {
		System.err.println("WARNING: "+s);
	}
	
	private void setupAlphabet(String letterComplementIn, String letterOrderIn) {
		letterOrder      = letterOrderIn;
		letterComplement = letterComplementIn;
		sc               = new LongSequence.SequenceClass(letterComplement);
		nMono            = letterOrder.length(); 
		nDi              = nMono*nMono;
		kMers            = null;
		kMersRC          = null;
		
		//Checks if the alphabet is "ACGT" or "ACGU":
        char tempArray[] = letterOrder.toCharArray(); 
        Arrays.sort(tempArray);
        sortedLetters    = new String(tempArray);
		
	}
	
	//Creates lists of k-mers and lookup dictionaries for indices.
	private void createKMerList(int kIn) {
		k = kIn;
		//Creates list of k-mers
		kMers       = new ArrayList<String>();
		kMersRC     = new ArrayList<String>();
		kMers.add("");
		kMersRC.add("");
		char[] nucl   = letterOrder.toCharArray();
		for(int iAdd=0; iAdd<k; iAdd++) {
			ArrayList<String> tempList   = new ArrayList<String>();
			ArrayList<String> tempListRC = new ArrayList<String>();
			for(int iKmer = 0; iKmer<kMers.size(); iKmer++) {
				String ts   = kMers.get(iKmer);
				String tsRC = kMersRC.get(iKmer);
				for(int iNucl=0; iNucl<nucl.length; iNucl++) {
					tempList.add(  ts  +nucl[iNucl]);
					tempListRC.add(nucl[nMono-1-iNucl]+tsRC);
				}
			}
			kMers    = tempList;
			kMersRC  = tempListRC;
		}
		
		//Creates k-mer lookup.
		kMerIndex   = new HashMap<String,Integer>();
		kMerIndexRC = new HashMap<String,Integer>();
		for(int iK=0; iK<kMers.size(); iK++) {
			kMerIndex.put(  kMers.get(  iK), iK);
			kMerIndexRC.put(kMersRC.get(iK), iK);
		}
		nKmers = kMers.size();
				
	}
	
	//Functions for reading from the sequence file
	///////////////////////////////////////////////

	private void openSequenceFile(String seqF) {
		

		if(letterComplement==null)
			setupAlphabet("C-G,A-T", "ACGT");
		
		sequenceFile         = seqF;
		firstNeededCharacter = 0;
		sequence             = "";
		
		String[] d = sequenceFile.split("\\.");
		if(d[d.length-1].equals("gz")) {
			if(verbose)
				System.out.println("> Opening gzip file "+seqF);
			try {
				br          = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(sequenceFile))));
			} catch(Exception e) {
			    System.err.println("Error: Could not open gz input sequence file "+sequenceFile);
			    e.printStackTrace();
			}
			
		} else {
			if(verbose)
				System.out.println("> Opening text file "+seqF);
			try{
				br          = new BufferedReader(new InputStreamReader(new FileInputStream(sequenceFile)));
			} catch(Exception e) {
			    System.err.println("Error: Could not open input sequence file "+sequenceFile);
			    e.printStackTrace();
			}
		}
		
	}
	
	/**
	 * Sets the input to be a fastq file
	 * @param {string} faFile Path to input fasta file.
	 */
	public void inputFasta(String faFile) {
		
		inputFormat  = "fasta";
		stringColumn = 0;
		tableData    = null;
		this.openSequenceFile(faFile);

	}

	/**
	 * Sets the input to be a specific fastq file 
	 * @param {string} fqFile Path to input fastq file.
	 */
	public void inputFastq(ArrayList<String> args) {
		
		if(args==null) {
			System.out.println("inputFastq(fqFile)                - Sets the input data stream to be the fastq file 'fqFile' containing sequences of uniform length.");
			return;
		}
		
		int nArgs = args.size(); 
		if(nArgs == 1) {
			String fqFile = args.get(0);
			inputFastq(fqFile);
		} else {
			error("The function inputFastq(fqFile) must have exactly one arguments.");
		}
	}
	
	public void inputFastq(String fqFile) {
		
		inputFormat  = "fastq";
		stringColumn = 0;
		tableData    = null;
		this.openSequenceFile(fqFile);
		
	}
	
	/**
	 * Sets the input to be a specific TSV file 
	 * @param {string} txtFile Path to input fastq file.
	 */
	public void inputTXT(ArrayList<String> args) {
		if(args==null) {
			System.out.println("inputTXT(txtFile)                 - Sets the input data stream to be the text file 'txtFile' containing sequences of uniform length.");
			return;
		}
		int nArgs = args.size();

		if(nArgs == 1) {
			String txtFile = args.get(0);
			inputTXT(txtFile);
		} else {
			error("The function inputTXT(txtFile) must have exactly one arguments.");
		}
		
	}
	
	public void inputTXT(String txtFile) {

		inputFormat  = "txt";
		stringColumn = 0;
		tableData    = null;
		this.openSequenceFile(txtFile);
		
	}
	
	/**
	 * Sets the input to be a specific TSV file 
	 * @param {string} tsvFile Path to input TSV file.
	 * @param {string} iString index of string column (first=0).
	 */
	
	public void inputTableTSV(ArrayList<String> args) {
		
		if(args==null) {
			System.out.println("inputTableTSV(tsvFile, iString=0) - Sets the input data stream to be the TSV table file 'tsvFile' with sequences (of uniform length) in column 'iString'.");
			return;
		}
		
		int nArgs = args.size();
		if(nArgs == 1 || nArgs==2) {
			String tsvFile = args.get(0);
			int iString    = 0;
			if(nArgs == 2)
				iString    = Integer.parseInt(args.get(1));
			inputTableTSV(tsvFile, iString);
		} else {
			error("The function inputTableTSV(tsvFile, iString=0) must have one or two arguments.");
		}
	}
	
	public void inputTableTSV(String tsvFile) {
		inputTableTSV(tsvFile, 0);
	}
	
	public void inputTableTSV(String tsvFile, int iString) {

		inputFormat  = "table";
		stringColumn = iString;
		tableData    = null;
		this.openSequenceFile(tsvFile);
	}
	
	
    // Code for reading strings
    ///////////////////////////
	/**
	 * Gets the next substring of the current sequence and returns true. Returns false if no more sequence can be added.
	 */
	public boolean nextSubstring() {
		boolean out;
		if(inputFormat.equals("fastq"))
			out = nextSubstring_fastq();

		else if(inputFormat.equals("fasta"))
			out = nextSubstring_fasta();
		
		else if(inputFormat.equals("fasta"))
			out = nextSubstring_txt();
		
		else if(inputFormat.equals("txt"))
			out = nextSubstring_txt();
			
		else if(inputFormat.equals("table"))
			out = nextSubstring_table();
		
		else 
		    throw new RuntimeErrorException(null, "Error: Invalid input format.");
		
		return out;
	}
	
	public boolean nextSubstring_fasta() {
		String line = null;
		try {
			//Finds the next substring in the current fasta element
			while(true) {
				line = br.readLine();
				if(line == null)
					break;
				
				if(line.length() == 0)
					continue;
				char firstChar = line.charAt(0);
				if( firstChar == '>' ) {
					//A new FASTA element was reached before finding a sequence.
					line = null;
					break;
				}
					
				else
					break;
			}
		}  catch(Exception e) {
		    System.err.println("Error: Could not read from the file "+sequenceFile);
		    e.printStackTrace();
		}
		
		//Trims of the already-used characters.
		sequence             = sequence.substring(firstNeededCharacter);
		firstNeededCharacter = 0;
		
		if(line==null) {
			return false;
			
		} else {
			//Adds the new characters.
			sequence             = sequence + line;
			return true;
		}
	}

	public boolean nextSubstring_fastq() {
		String line = null;
		try {
			//Reads the next sequence in the fasta file
			while(true) {
				line = br.readLine();
				
				if(line == null)
					//End of file
					break;
				
				if(line.length() == 0)
					//Empty line, continues
					continue;
				char firstChar = line.charAt(0);

				if(firstChar == '+' ) {
					//Sequence has already been read. 
					line = null;
					break;
					
				} else if ( firstChar == '@' || firstChar == '#' ) {
					//Ignore lines starting with @ and #.
					continue;
					
				} else
					break;
			}
		}  catch(Exception e) {
		    System.err.println("Error: Could not read from the file "+sequenceFile);
		    e.printStackTrace();
		}
		
		//Trims of the already-used characters.
		sequence             = sequence.substring(firstNeededCharacter);
		firstNeededCharacter = 0;


		if(line==null) {
			if(verbose)
				System.out.println("End of entry!");
			return false;
			
		} else {
			sequence         = sequence + line;
			return true;
		}
	}

	public boolean nextSubstring_txt() {
		//In a TXT file, each line contains an independent sequence, so there is no next substring 
		return false;
	}
	
	public boolean nextSubstring_table() {
		//In a table, each line contains an independent sequence, so there is no next substring 
		return false;
	}
	
	/**
	 * Starts reading the next sequence and loads the table data.
	 */
	public boolean nextSequence() {

		boolean out;
		if(inputFormat.equals("fasta"))
			 out = nextSequence_fasta();
		else if(inputFormat.equals("fastq"))
			out =  nextSequence_fastq();
		else if(inputFormat.equals("table"))
			out =  nextSequence_table();
		else if(inputFormat.equals("txt"))
			out =  nextSequence_txt();
		else 
		    throw new RuntimeErrorException(null, "Error: Invalid input format.");

		if(sortedLetters == "ACGT") 
			sequence = sequence.replace("U", "T");
		else if (sortedLetters == "ACGU")
			sequence = sequence.replace("T", "U");
		
		return out;
	}
	
	private boolean nextSequence_fastq() {
		sequence             = "";
		firstNeededCharacter = 0;
		return nextSubstring_fastq();
	}
	
	private boolean nextSequence_fasta() {
		sequence             = "";
		firstNeededCharacter = 0;
		return nextSubstring_fastq();
	}
	
	private boolean nextSequence_txt() {
		
		String line = null;
		try {
			line = br.readLine();
		}  catch(Exception e) {
		    System.err.println("Error: Could not read from the file "+sequenceFile);
		    e.printStackTrace();
		}
		
		if(line == null)
			return false;
		
		//Saves the sequence
		sequence             = line;
		firstNeededCharacter = 0;
		
		return true;
	}
	
	private boolean nextSequence_table() {
		
		String line = null;
		try {
			line = br.readLine();
		}  catch(Exception e) {
		    System.err.println("Error: Could not read from the file "+sequenceFile);
		    e.printStackTrace();
		}
		
		if(line == null)
			return false;
		
		//Makes sure that the columns have the correct number of columns.
		String[] d = line.split("\t");
		if(tableData == null) {
			if(d.length < 2) 
			    throw new RuntimeErrorException(null, "Error: A table must contain at least two columns.");
			tableData = new double[d.length-1];
		} else { 
			if(d.length-1 != tableData.length)
				throw new RuntimeErrorException(null, "Error: All rows must contain "+tableData.length+" columns.");
		}
		
		//Saves the sequence
		sequence             = d[stringColumn];
		firstNeededCharacter = 0;
		
		//Copies the data.
		int iWrite = 0;
		for(int iCol=0; iCol<d.length; iCol++) {
			if(iCol == stringColumn)
				continue;
			tableData[iWrite] = Double.parseDouble(d[iCol]);
			iWrite++;
		}
		
		return true;
	}

	/**
	 * Closes the sequence input file
	 */
	public void close() {
		
		if(br != null) {
			try {
				br.close();
			} catch(Exception e) {
			    System.err.println("Error: Could not close the file "+sequenceFile);
			    e.printStackTrace();
			}
		}
		
		
		inputFormat          = null;
		sequenceFile         = null;
		stringColumn         = 0;
		tableData            = null;
		firstNeededCharacter = 0;
		sequence             = null;
		
	}
	
	private void restartFile() {
		try {
			br.close();
			this.openSequenceFile(sequenceFile);
		} catch(Exception e) {
		    System.err.println("Error: Could not close and reopen the file "+sequenceFile);
		    e.printStackTrace();
		}
	}
	
	// Code for loading a scoring model JSON file
	/////////////////////////////////////////////
	/**
	 * Loads a scoring model JSON file..
	 * @param {string} Path to JSON file.
	 */
	public void loadScoringModel(ArrayList<String> args) {
		if(args==null) {
			System.out.println("loadScoringModel(modelJSON)       - Loads a scoring model from a JSON file.");
			return;
		}
		int nArgs = args.size();

		if(nArgs==1) {
			String filePath = args.get(0);
			loadScoringModel(filePath, generalSchemaFile);
		} else {
			error("The function loadScoringModel(modelJSON) must have exactly one arguments.");
		}
		
	}
	
	public void loadScoringModel(String jsonPath, String generalSchemaFile) {
		
		
		if(verbose)
			System.out.println("> Loading "+jsonPath);
		//Only using the last model in the file.
		
		oModel = JSONModel.loadJSONFile(jsonPath);
		
		JSONModel.validateSchemaFile_O(generalSchemaFile, oModel);
		//Adds model fitting constraints if they do not already exist
		JSONModel.addEmptyModelFittingConstraints(generalSchemaFile, oModel);
		
		loadOModel();
	}
	
	private void loadOModel() {
		
		//Sets the alphabet.
		if(oModel.has("modelSettings")) {
			JSONObject oSett = oModel.getJSONObject("modelSettings");
			letterComplement = oSett.has("letterComplement") ? oSett.getString("letterComplement") : "C-G,A-T";
			letterOrder      = oSett.has("letterOrder")      ? oSett.getString("letterOrder")      : CombinedLikelihood.alphabetDefToOrderedLetters(letterComplement);
			
			//Saves the alphabet
			if(!oSett.has("letterComplement"))
				oSett.put("letterComplement", letterComplement);
			if(!oSett.has("letterOrder"))
				oSett.put("letterOrder", letterOrder);
		} else {
			letterComplement = "C-G,A-T";
			letterOrder      = "ACGT";
		}
		setupAlphabet(letterComplement, letterOrder);
		
		JSONArray aBindingModeSett  = oModel.getJSONObject("modelSettings").getJSONArray("bindingModes");
		JSONArray aBindingModeCoef  = oModel.getJSONObject("coefficients").getJSONArray("bindingModes");
		nModes                      = aBindingModeSett.length();
		size                        = new ArrayList<Integer>();
		dInt                        = new ArrayList<Integer>();
		mismatchGauge               = false;
		//Keeping track of the number of frames.
		int iCurrent=0;
		int currentNFrames=0;
		iBindingMode                = new ArrayList<Integer>();
		iBindingMode.add(0);
		nFrames                     = new ArrayList<Integer>();
		
		//Sliding window
		sw                          = new ArrayList<SlidingWindow>();
		
		//Loops over binding modes and fills in information
		for(int iMode=0; iMode<nModes; iMode++) {
			JSONObject oBindingModeSett = aBindingModeSett.getJSONObject(iMode);
			JSONObject oBindingModeCoef = aBindingModeCoef.getJSONObject(iMode);
			int currentSize             = oBindingModeSett.getInt("size");
			int currentDInt             = Math.min(oBindingModeSett.getInt("dinucleotideDistance"), currentSize-1);
			size.add(currentSize);
			dInt.add(currentDInt);
			
			if(currentSize>0) {
				//currentNFrames               = 2*(L-currentSize+1);
				ArrayList<ArrayList<double[]>> aal = new ArrayList<ArrayList<double[]>>();
				aal.add(new ArrayList<double[]>());
				
				JSONArray aMonoBetas         = oBindingModeCoef.getJSONArray("mononucleotide");
				int nBetas                   = nMono*currentSize;
				double[] tempBetas           = new double[nBetas];
				for(int iBeta=0; iBeta<nBetas; iBeta++)
					tempBetas[iBeta]         = aMonoBetas.getDouble(iBeta);
				aal.get(0).add(tempBetas);

				if( currentDInt>0 ) {
					aal.add(new ArrayList<double[]>());
					for(int id=0; id< currentDInt; id++) {
						JSONArray aDiBetas   = oBindingModeCoef.getJSONArray("dinucleotide").getJSONArray(id);
						nBetas               = nDi*(currentSize-id-1);
						tempBetas            = new double[nBetas];
						for(int iBeta=0; iBeta<nBetas; iBeta++)
							tempBetas[iBeta] = aDiBetas.getDouble(iBeta); 
						aal.get(1).add(tempBetas);
					}
				}
				
				sw.add(new SlidingWindow("", "", sc, aal, letterOrder));
				
			} else {
				currentNFrames = 1;
				sw.add(null);
			}
			nFrames.add(currentNFrames);
			iCurrent += currentNFrames; 
			iBindingMode.add(iCurrent);
				
		}
		
		//Total number of binding frames across all modes:
		nPredictors = iBindingMode.get(iBindingMode.size()-1);
	}
	
	// Code for a fit line models
	/////////////////////////////
	/**
	 * Loads a model from a JSON file (one fit per line). If multiple models are contained in the file, only the last is retained.
	 * @param {string} Path to JSON file.
	 */
	public void loadFitLine(ArrayList<String> args) {
		
		if(args==null) {
			System.out.println("loadFitLine(fitJSON,iLine=-1)     - Loads the model on line 'iLine' from a ProBound output file 'fitJSON'. By default the line with the smallest -log(likelihood) is retrived.");
			return;
		}
		int nArgs = args.size();

		if(nArgs==1 || nArgs==2) {
			String filePath = args.get(0);
			int iLine       = -1;
			if(nArgs==2)
				iLine       = Integer.parseInt(args.get(1));
			loadFitLine(filePath, generalSchemaFile, iLine);
		} else {
			error("The function loadFitLine(filePathLine, iLine=-1) must at least one arguments.");
		}
	}
	
	public void loadFitLine(String jsonPath, String generalSchemaFile, int iLine) {
		if(verbose)
			System.out.println("> Loading "+jsonPath);
		//Only using the last model in the file.
		
		oModel = JSONModel.loadJSONLine(jsonPath, iLine);

		JSONModel.validateSchemaFile_O(generalSchemaFile, oModel);
		//Adds model fitting constraints if they do not already exist
		JSONModel.addEmptyModelFittingConstraints(generalSchemaFile, oModel);

		loadOModel();
		
		
	}


	
    public void writeModel(ArrayList<String> args) {
    	if(args==null) {
    		System.out.println("writeModel(outFile)               - Writes the current model to a JSON file 'outFile'.");
    		return;
    	}
    	int nArgs = args.size();
		if(nArgs == 1) {
			String outFile = args.get(0);
			writeModel(outFile);
		} else {
			error("The function writeModel(outFile) must have exactly one arguments.");
		}
    }
    
    public void writeModel(String outFile) {
		LikelihoodOptimizer.writeCompactJSONModel(oModel, outFile, false);
    }
    
    //Functions for fixing the gauge
    ////////////////////////////////
	/**
	 * Sets the mismatch gauge and absorbs the affinity of the highest-scoring sequence in the activities. 
	 */
    public void setMismatchGauge(ArrayList<String> args) {
    	if(args==null) {
    		System.out.println("setMismatchGauge()                - Imposes the mismatch gauge on the binding modes, meaning the top sequence has score zero.");
    		return;
    	}
    	
    	int nArgs = args.size();
		if(nArgs == 0) {
			setMismatchGauge();
		} else {
			error("The function setMismatchGauge() must have exactly zero arguments.");
		}
    }
    
    public void setMismatchGauge() {
    	
    	//1.  Fixes the scoring matrices.
    	for(int iBM=0; iBM<nModes; iBM++) {
    		//1.0 Read all coefficients.
    		//Reads coefficients for the current binding mode
    		JSONObject oBMSetting       = oModel.getJSONObject("modelSettings").getJSONArray("bindingModes").getJSONObject(iBM);
    		JSONObject oBMCoeff         = oModel.getJSONObject("coefficients").getJSONArray("bindingModes").getJSONObject(iBM);
    		int size                    = oBMSetting.getInt("size");
    		
    		if(size>0) {
        		//Mononucleotide
        		double[] monoBetas          = ModelComponent.readFromJSON_d(oBMCoeff.getJSONArray("mononucleotide"));
        		//Dinucleotide
        		ArrayList<double[]> diBetas;
        		if(oBMCoeff.has("dinucleotide"))
        			diBetas                 = ModelComponent.readFromJSON_Ad(oBMCoeff.getJSONArray("dinucleotide"));
        		else
        			diBetas                 = new ArrayList<double[]>();
        		
        		//Info about current binding mode.
            	int nNuc   = monoBetas.length / size;
            	int nDi    = nNuc * nNuc;

            	
            	//1.1  Identify the highest-affinity sequence for all binding modes.
        		int[] topSeq = null;
        		if(dInt.get(iBM) > 1)
        			topSeq = sampleTopSequence(size, monoBetas, diBetas);
        		else
        			topSeq = sampleTopSequence(size, monoBetas, diBetas);
        		
            	//1.2 Compute the energy of the highest affinity sequence.
        		double topScore = scoreSeq(topSeq, monoBetas, diBetas, size);
        		//System.out.println("topScore = "+topScore);
        		
            	//1.3 Compute the energy of all single-base mismatches relative the highest-affinity sequence. Save as new PWM.
        		int[] tempSeq              = new int[size];
        		double[] monoBetaMismatch = new double[monoBetas.length];
        		for(int i=0; i<size; i++)
        			tempSeq[i] = topSeq[i];
        		for(int x=0; x<size; x++) {
        			for(int n=0; n<nNuc; n++) {
        				int old = tempSeq[x];
        				tempSeq[x] = n;
        				monoBetaMismatch[x*nNuc+n] = scoreSeq(tempSeq, monoBetas, diBetas, size) - topScore;
        				tempSeq[x] = old;
        			}
        		}
        		oBMCoeff.put("mononucleotide", ModelComponent.JSONArrayConvert_d(monoBetaMismatch));
        				
            	//1.4 Compute the energy of all double-base mismatches relative the highest-affinity sequence. Save as new interactions.
        		ArrayList<double[]> diBetaMismatch = new ArrayList<double[]>();
        		for(int d=1; d<=diBetas.size(); d++) {
        			double[] newDiBeta = new double[diBetas.get(d-1).length];
        			for(int x=0; x<size-d; x++) {
        				for(int n1=0; n1<nNuc; n1++) {
        					for(int n2=0; n2<nNuc; n2++) {
        						//Updates the sequence
        						int old1     = tempSeq[x];
        						int old2     = tempSeq[x+d];
        						tempSeq[x]   = n1;
        						tempSeq[x+d] = n2;
        						
        						//Computes the dinucleotide interaction.
        						newDiBeta[x*nDi+nNuc*n1+n2] = 
        								scoreSeq(tempSeq, monoBetas, diBetas, size) 
        								- (topScore + monoBetaMismatch[x*nNuc+n1] + monoBetaMismatch[(x+d)*nNuc+n2]);
        						
        						//Reverts to o the old sequence
        						tempSeq[x]   = old1;
        						tempSeq[x+d] = old2;
        					}
        				}
        			}
        			diBetaMismatch.add(newDiBeta);
        		}
        		oBMCoeff.put("dinucleotide",   ModelComponent.JSONArrayConvert_Ad(diBetaMismatch));
        		
        		
        		//1.5 Shifts the position bias.  
        		double[] meanPB = null;
        		if(oBMSetting.getBoolean("positionBias")) {
        			if(!oBMCoeff.has("positionBias") || oBMCoeff.getJSONArray("positionBias").length()==0) 
        				throw new java.lang.RuntimeException("ERROR: positionBias=true but no coefficients were specified.");

        			JSONArray aPB = oBMCoeff.getJSONArray("positionBias");
        			meanPB        = new double[aPB.length()];
        			for(int iExp=0; iExp<aPB.length(); iExp++) {
            			//1.1 Computes the log(mean(exp(positionbias)) for each experiment
            			double expSum = 0;
            			int nVal      = 0;
    					for(int s=0; s<2; s++) {
    						for(int x=0; x<aPB.getJSONArray(iExp).getJSONArray(s).length(); x++) {
        						expSum += Math.exp(aPB.getJSONArray(iExp).getJSONArray(s).getDouble(x));
    							nVal   += 1;
        					}
    					}
    					
    					//Saves the log-mean position bias.
    					meanPB[iExp] = Math.log(expSum/nVal);
        				
        				
            			//1.2 Shifts the relevant activities, removes the position bias.
    					for(int s=0; s<2; s++)
    						for(int x=0; x<aPB.getJSONArray(iExp).getJSONArray(s).length(); x++)
    							aPB.getJSONArray(iExp).getJSONArray(s).put(x, aPB.getJSONArray(iExp).getJSONArray(s).getDouble(x) - meanPB[iExp]);
    					
        			}
        		} 
        		
        		//1.6 Shifts the highest-scoring value to the activities.
        		//1.6.1 Reads activities, if any exits, otherwise adds a single zero.
        		JSONArray act = null;
        		if(oBMCoeff.has("activity")) {
        			act = oBMCoeff.getJSONArray("activity");
        		} else {
        			act = new JSONArray();
        			JSONArray temp = new JSONArray();
        			act.put(temp);
        			temp.put(0);
        		}
        		if(meanPB!=null && meanPB.length>0 && meanPB.length != act.length())
        			throw new RuntimeErrorException(null, "ERROR: The position bias differs across experiments the number of experiments. However, the length of the position bias and the activity vectors do not match.");

        		
        		//1.6.2 Shifts the activities
        		for(int iExp=0; iExp<act.length(); iExp++) {
        			JSONArray temp = act.getJSONArray(iExp);

        			double activityShift = topScore;
    				if(meanPB!=null)
    					activityShift += meanPB[iExp%(meanPB.length)];

        			for(int iRound=0; iRound<temp.length(); iRound++) {
        				temp.put(iRound, temp.getDouble(iRound) + activityShift);
        			}
        		}
        		//1.6.3 Sets the activity.
        		oBMCoeff.put("activity", act);
        		
        		//1.7. Updates the activities of the interactions.
        		JSONArray aInt         = oModel.getJSONObject("modelSettings").getJSONArray("bindingModeInteractions");
        		for(int iInt=0; iInt<aInt.length(); iInt++) {
        			JSONObject oInt    = aInt.getJSONObject(iInt);
        			
        			//Reads activities.
        			JSONArray intAct   = null;
            		if(oBMCoeff.has("activity")) {
            			intAct         = oBMCoeff.getJSONArray("activity");
            		} else {
            			intAct         = new JSONArray();
            			JSONArray temp = new JSONArray();
            			intAct.put(temp);
            			temp.put(0);
            		}
            		
            		//Shifts the activities.
            		for(int iI=0; iI<2; iI++) {
            			if(iBM == oInt.getJSONArray("bindingModes").getInt(iI)) {
                    		for(int iExp=0; iExp<intAct.length(); iExp++) {
                    			JSONArray temp = intAct.getJSONArray(iExp);
                    			
                				double activityShift = topScore;
                				if(meanPB!=null)
                					activityShift += meanPB[iExp%(meanPB.length)];

                    			for(int iRound=0; iRound<temp.length(); iRound++) 
                    				temp.put(iRound, temp.getDouble(iRound)+activityShift);
                    		}
            			}
            		}
        			
            		//Saves the activities.
            		oInt.put("activity", intAct);
        		}
    		}
    		
    		//1.8. Shift the maximum position-bias score into the activities.
    	}
    	//3. Fix the interactions
		//TODO: Shift the interactions into the activities.
		//3.1 For each interaction, shift the top value to the interaction-activity. 
		//3.2 Subtract the top value from all interactions.
    	if(oModel.has("modelSettings") &&
    			oModel.getJSONObject("modelSettings").has("bindingModeInteractions") &&
    			oModel.getJSONObject("modelSettings").getJSONArray("bindingModeInteractions").length()>0)
    				System.err.println("WARNING: The current model has binding mode interactions, but fixing of interactions is not implemented for setMismatchGauge().");
    	
    	loadOModel();
    	mismatchGauge = true;
    	
    }
    
    public void removeBindingMode(ArrayList<String> args) {

    	if(args==null) {
    		System.out.println("removeBindingMode(iBM)            - Removes the binding mode with index 'iBM' from the current model.");
    		return;
    	}
    	
    	int nArgs = args.size();
		if(nArgs == 1) {
			int iBM = Integer.parseInt(args.get(0));
			removeBindingMode(iBM);
		} else {
			error("The function removeBindingMode(iBM) must have exactly one arguments.");
		}
	
    }
    
    
    public void removeBindingMode(int iBM) {

    	CombinedLikelihood.removeBindingMode_JSON(oModel, iBM);
    	loadOModel();
    	
    }
    
    public void selectBindingMode(ArrayList<String> args) {
    	
    	if(args==null) {
    		System.out.println("selectBindingMode(iBM)            - Removes all binding modes except that with index 'iBM' from the current model.");
    		return;
    	}

    	int nArgs = args.size();
		if(nArgs == 1) {
			int iBM = Integer.parseInt(args.get(0));
			selectBindingMode(iBM);
		} else {
			error("The function selectBindingMode(iBM) must have exactly one arguments.");
		}

		
		
    }
    
    public void selectBindingMode(int iBM) {

    	int nModes = oModel.getJSONObject("modelSettings").getJSONArray("bindingModes").length();
    	for(int iBMRemove=nModes-1; iBMRemove>-1; iBMRemove--) {
    		if(iBMRemove == iBM)
    			continue;
   			CombinedLikelihood.removeBindingMode_JSON(oModel, iBMRemove);
        		
    	}
    	loadOModel();
    	
    }

    
    public void selectAndCleanBindingMode(ArrayList<String> args) {
    	
    	if(args==null) {
    			System.out.println("selectAndCleanBindingMode(iMode)  - Selects a binding mode and removes all other binding modes, interactions, enrichment models etc.");
    			return;
    	}
    	
    	int nArgs = args.size();
		if(nArgs == 1) {
			int iMode = Integer.parseInt(args.get(0));
			selectAndCleanBindingMode(iMode);
		} else {
			error("The function selectAndCleanBindingMode(iMode) must have exactly one arguments.");
		}
    }
    
    public void selectAndCleanBindingMode(int iMode) {

    	JSONObject oSett = oModel.getJSONObject("modelSettings");
		JSONObject oCoef = oModel.getJSONObject("coefficients");

		//1: Selects the identified binding mode
    	selectBindingMode(iMode);

    	//2: Removes binding mode interactions.
    	if(oSett.has("bindingModeInteractions"))
    		oSett.remove("bindingModeInteractions");
    	if(oCoef.has("bindingModeInteractions"))
    		oCoef.remove("bindingModeInteractions");

    	//3: Sets the activity to zero.
    	JSONObject oBMCoef = oCoef.getJSONArray("bindingModes").getJSONObject(0);
    	JSONObject oBMSett = oSett.getJSONArray("bindingModes").getJSONObject(0);
    	oBMCoef.put("activity", ModelComponent.JSONArrayConvert_dd(ModelComponent.zero_dd(1, 1)));

    	//4: Removes position bias:
    	oBMSett.put("positionBias", false );
    	oBMCoef.put("positionBias", new JSONArray());

    	//5: Clean up model and loads it:
    	/////////////////////////////////
    	//5.1. Cleans up the scoring model
    	cleanUpScoringModel();
    	//5.2. Rebuilds the sliding windows etc.
    	loadOModel();
    	
    }

    public void removeInteraction(ArrayList<String> args) {
    	
    	if(args==null) {
    		System.out.println("removeInteraction(iInt)           - Removes the binding mode interaction with index 'iInt'.");
    		return;
    	}
    	
    	int nArgs = args.size();
		if(nArgs == 1) {
			int iInt = Integer.parseInt(args.get(0));
			removeInteraction(iInt);
		} else {
			error("The function removeInteraction(iInt) must have exactly one arguments.");
		}
		
    }
    
    public void removeInteraction(int iInt) {
    	
    	CombinedLikelihood.removeBindingModeInteraction_JSON(oModel, iInt);
    	loadOModel();

    	
    }
    
	/**
	 * Computes the highest-affinity sequence of a binding mode using sampling.
	 * @param {int} Index iBM of binding mode in model.
	 */
    private int[] sampleTopSequence(int size, double[] monoBetas, ArrayList<double[]> diBetas) {

    	//Returns an empty vector if the size is zero.
    	if(size==0)
    		return new int[0];
    	
    	//0. Defines meta parameters
    	int nRepeats = 40; //Number of monte-carlo runs.

    	//1.  Re-scales scoring matrix
    	//1.1 Compute an upper bound of the maximum range of values
    	int nNucl   = monoBetas.length / size;
    	int nDi     = nNucl * nNucl;
    	int dInt    = diBetas.size();
    	double r    = 0;
    	for(int x1=0; x1<size; x1++) {
    		r      += valueRange(monoBetas, nNucl*x1, nNucl*(x1+1));
    		int x2H = Math.min(size, x1+dInt+1);
    		for(int x2=x1+1; x2<x2H; x2++)
    			r  += valueRange(diBetas.get(x2-x1-1), nDi*x1,nDi*(x1+1));
    	}
    		
    	//1.2 Re-scales the matrices to the range of values is less than or equal to one.
    	double[] scaledMono          = ModelComponent.clone_d(monoBetas);
    	ArrayList<double[]> scaledDi = ModelComponent.clone_Ad(diBetas);
    	for(int i=0; i<scaledMono.length; i++) scaledMono[i] /= r;
    	for(int i1=0; i1<scaledDi.size(); i1++)
    		for(int i2=0; i2<scaledDi.get(i1).length; i2++)
    			scaledDi.get(i1)[i2] /= r;
    	
    	//2.  Repeatedly performs Monte-Carlo simulations 
    	int[] maxSeq        = null;
    	double maxScore     = 0;
    	for(int iRep=0; iRep<nRepeats; iRep++) {
    		int[] newSeq    = singleMonteCarloRun(size, scaledMono, scaledDi);
    		double newScore = scoreSeq(newSeq, scaledMono, scaledDi, size);
    		if(maxSeq==null || newScore>maxScore) {
    			maxSeq      = newSeq;
    			maxScore    = newScore;
    		}
    	}
    	
    	return maxSeq;
    }
    
    //Computes the difference between the maximum and minimum value in a vector. 
    private double valueRange(double[] v, int iMin, int iMax) {
    	double min = v[iMin], max=v[iMin];
    	for(int i=iMin+1; i<iMax; i++) {
    		if(v[i]>max)
    			max=v[i];
    		if(v[i]<min)
    			min=v[i];
    	}
    	return max-min;
    }
    
    private int[] singleMonteCarloRun(int size, double[] monoBetas, ArrayList<double[]> diBetas) {
    	
    	//Meta parameters.
    	int nEFolds       = 20;
    	int editsPerEFold = 10;
    	
    	//Properties of the scoring model. 
    	int dInt   = diBetas.size();
    	int nNuc   = monoBetas.length / size;
    	int nDi    = nNuc * nNuc;
    	
    	//Generates random sequence.
    	int[] seq  = new int[size];
		Random generator	= new Random();
    	for(int x=0; x<seq.length; x++)
    		seq[x] = generator.nextInt(nNuc);
    	
    	//Generates a list of edits.
    	ArrayList<int[]> edits = new ArrayList<int[]>();
    	for(int x=0; x<size; x++) {
    		for(int n=0; n<nNuc; n++) {
    			int[] e = new int[2];
    			e[0]    = x;
    			e[1]    = n;
    			edits.add(e);
    		}
    	}
    	

    	for(int iFold=0; iFold<nEFolds; iFold++) {
    		double foldFactor = Math.exp(iFold);
    		for(int iEditRound=0; iEditRound<editsPerEFold; iEditRound++) {
    			//Permute the edits.
    			Collections.shuffle(edits);
    			
    			for(int iEdit=0; iEdit<edits.size(); iEdit++) {
    				int[] e = edits.get(iEdit);
    				int eX  = e[0];
    				int eN  = e[1];
    				
    				//Only tests edits where the edit changes the sequence.
    				if(seq[eX] == eN)
    					continue;
    				
    				//Computes the change in energy.
    				double deltaE = 0;
    				int xMin = Math.max(0,    eX-dInt);
    				int xMax = Math.min(size, eX+dInt+1);
    				for(int x=xMin; x<xMax; x++) {
    					if(x<eX) {
    						double[] db = diBetas.get(eX-x-1);
    						deltaE     += db[ nDi*x  + nNuc*seq[x]  + eN      ] 
    									- db[ nDi*x  + nNuc*seq[x]  + seq[eX] ];
    					} else if(x>eX) {
    						double[] db = diBetas.get(x-eX-1);
    						deltaE     += db[ nDi*eX + nNuc*eN      + seq[x] ] 
    									- db[ nDi*eX + nNuc*seq[eX] + seq[x] ];
    					} else {
    						deltaE += monoBetas[eX*nNuc + eN     ] 
    								- monoBetas[eX*nNuc + seq[eX]];
    					}
    				}
    				
    				//Saves the edit.
    				if(deltaE>0 || generator.nextDouble() < Math.exp(deltaE * foldFactor)) {
    					//System.out.println("deltaE = "+deltaE);
    					seq[eX] = eN;
    				} 
    			}
    		}
    	}
    	
    	return seq;
    	
    }
    
	/**
	 * Scores a sequence
	 * @param {int} Index iBM of binding mode in model.
	 */
    private double scoreSeq(int[] seq, double[] monoBetas, ArrayList<double[]> diBetas, int size) {
    	
    	double out = 0;
    	int dInt   = diBetas.size();
    	int nNuc   = monoBetas.length / size;
    	int nDi    = nNuc * nNuc;
    	
    	for(int x1=0; x1<size; x1++) {
    		out    += monoBetas[nNuc*x1 + seq[x1]];
    		int x2H = Math.min(size, x1+dInt+1);
    		for(int x2=x1+1; x2<x2H; x2++)
    			out += diBetas.get(x2-x1-1)[nDi*x1+nNuc*seq[x1]+seq[x2]];
    		
    	}
    	return out;
    	
    }
    
    
	/**
	 * Builds a consensus model by integrating all experiments/rounds.	
	 */
    public void buildConsensusModel(ArrayList<String> args) {
    	if(args == null) {
    		System.out.println("buildConsensusModel()             - Processess a ProBound-fit to remove experiment-specific parameters such as experiment/round-specific activities, position bias, flank length and the enrichment model.");
    		return;
    	}
    	int nArgs = args.size();

		if(nArgs == 0) {
			buildConsensusModel();
		} else {
			error("The function buildConsensusModel() must have exactly zero arguments.");
		}
    }

	public void buildConsensusModel() {
	
    	JSONObject oCoef = oModel.getJSONObject("coefficients");
    	JSONObject oSett = oModel.getJSONObject("modelSettings");

    	setMismatchGauge();
    	
    	//0. Reads all binding-mode activities and binding-mode interaction activities into nested structures.
    	int nExp                              = 0;
    	ArrayList<Integer> nActivity          = null;             // Stores the number of activities in each experiment.
    	ArrayList<ArrayList<double[]>> act    = null; // For each bm and exp, stores the activities. 
    	ArrayList<ArrayList<double[]>> intAct = null; // For each int and exp, stores the activities.
    	
    	//0.1 Saves the activities for the binding modes.
    	for(int iBM=0; iBM<nModes; iBM++) {
    		JSONObject oBMCoeff = oCoef.getJSONArray("bindingModes").getJSONObject(iBM);
    		
    		ArrayList<double[]> newAct = null;
    		if(oBMCoeff.has("activity") || oBMCoeff.getJSONArray("activity").length()==0) {
    			newAct    = ModelComponent.readFromJSON_Ad(oBMCoeff.getJSONArray("activity"));
    		} else {
    			newAct        = new ArrayList<double[]>();
    			newAct.add(new double[1]);
    		}

			if(act==null) {
    			//Adds the number of activities in each experiment for the first binding mode...
    			nExp                      = newAct.size();
    			act                       = new ArrayList<ArrayList<double[]>>();
				nActivity                 = new ArrayList<Integer>();
				for(int iExp=0; iExp<nExp; iExp++)
					nActivity.add(newAct.get(iExp).length);

			} else {
				//.. or checks so that the correct number of activities are given for the second.
				if(nExp != newAct.size())
					throw new RuntimeErrorException(null, "Error: All binding modes must have the same number of experiments, but this is not the case.");
				for(int iExp=0; iExp<nExp; iExp++)
					if( nActivity.get(iExp) != newAct.get(iExp).length )
						throw new RuntimeErrorException(null, "Error: All binding modes must have the same number of activities, but this is not the case.");
			}
			act.add(newAct);
			
    	}
    	

    	//0.2 Saves the interaction activities  
    	if(oCoef.has("bindingModeInteractions")) {
    		JSONArray aIntCoeff   = oCoef.getJSONArray("bindingModeInteractions");
    		JSONArray aIntSetting = oSett.getJSONArray("bindingModeInteractions");

        	int nInt = aIntCoeff.length();
        	for(int iInt=0; iInt<nInt; iInt++) {
        		JSONObject oIntCoeff          = aIntCoeff.getJSONObject(iInt);
    			JSONObject oIntSetting        = aIntSetting.getJSONObject(iInt);
    			ArrayList<double[]> newIntAct = null;
    			if(oIntCoeff.has("activity") || oIntCoeff.getJSONArray("activity").length()==0) {
    				newIntAct = ModelComponent.readFromJSON_Ad(oIntCoeff.getJSONArray("activity"));
    			} else {
    				newIntAct = new ArrayList<double[]>();
    				newIntAct.add(new double[1]);
    			}
    			
    			if(nExp != newIntAct.size())
    				throw new RuntimeErrorException(null, "Error: All binding-mode interactions must have the same number of experiments, but this is not the case.");
    			
				for(int iExp=0; iExp<nExp; iExp++)
					if( nActivity.get(iExp) != newIntAct.get(iExp).length )
						throw new RuntimeErrorException(null, "Error: All binding-mode interactions must have the same number of activities, but this is not the case.");

        		if(intAct==null) {
        			intAct = new ArrayList<ArrayList<double[]>>();
        		}
        		
        		int[] newBMs = new int[2];
        		newBMs[0] = oIntSetting.getJSONArray("bindingModes").getInt(0);
        		
        		intAct.add(newIntAct);
        	}
    	} 
    	
    	
    	//1. Remove the position bias (the position bias should be exp-mean-zeroed after setMismatchGauge(), so it is simply removed here).
    	for(int iBM=0; iBM<nModes; iBM++) {
    		JSONObject oPBCoef   = oCoef.getJSONArray("bindingModes").getJSONObject(iBM);
    		JSONObject oPBSett   = oSett.getJSONArray("bindingModes").getJSONObject(iBM);
    		if(oPBCoef.has("positionBias")) {
    			oPBCoef.remove("positionBias");
    			oPBSett.put("positionBias", false);
    		}
    	}
    	
    	//2. Collapse positionMatrix into spacingVector.
    	if(intAct!=null) {
    		System.err.println("WARNING: Gauge-fixing and averaging of binding-mode interations have note been implemented. Removes binding mode interactions.");
    		if(oCoef.has("bindingModeInteractions")) {
    			oCoef.remove("bindingModeInteractions");
    			oCoef.put("bindingModeInteractions", new JSONArray());
    		}
    		if(oSett.has("bindingModeInteractions")) {
    			oSett.remove("bindingModeInteractions");
    			oSett.put("bindingModeInteractions", new JSONArray());
    		}
    		intAct=null;
    		/*
        	for(int iInt=0; iInt<intAct.size(); iInt++) {

        		JSONObject oIntCoef = oCoef.getJSONArray("bindingModeInteractions").getJSONObject(iInt);
        		JSONObject oIntSett = oSett.getJSONArray("bindingModeInteractions").getJSONObject(iInt);
        		
            	//2.0 The gauge fixing should rescale the interaction alphas so that the max imum value is one.
    			
        		//2.1 If experiment as spacing vector, use default functions to convert it into a position matrix.
        		//Implement these functions in BindingModeInteractions class. Use to read seeding for spacing.
        		if(oIntCoef.has("spacingVector")) {
        			JSONArray aSV = oIntCoef.getJSONArray("spacingVector");
        			// Implement
            	}

        		ArrayList<ArrayList<double[][]>> averagedInteractions = null;
        		if(oIntSett.has("positionMatrix")) {
            		//2.2 If there is a position matrix, compute the aligned, averaged interaction across experiment 
        			
        			//Implement:
        			// 1. Convert interaction betas into nested ArrayList/double[]. Exponentiate to compute alphas.
        			// 2. Create a 'target' packing matrix that can fit all interaction components.
        			// 3. Make the target packing 'homogenous' using:
        			//      public static void makeInteractionPackingHomogenous(ArrayList<ArrayList<int[][]>> interactionPacking)
        			// 4. Create packing matrices for all sub-matrices.
        			// 5. Loop through the matrices and fill in the correct index from in the 'target' packing matrix.
        			// 6. Use 'removePackingGaps' to remove packing gaps in the packing.
        			// 7. Use packGradient_A(interactions, packing) to sum over Exp[interaction] and
        			//    use packGradient_A(ones, packing) to compute the sumer of summed DOFs. Use this to compute the Exp-Mean
        			// 		- Consider renaming packGradient_{O,A} and packParameters_{O,A} to addUsingPacking_{O,A} and setusingPacking_{O,A}
        			// 8. Unpack to get homogenous averaged interaction MATRIX.

        			//Potentially useful functions:
        			//	ArrayList<ArrayList<ArrayList<int[][]>>> interactionRange = range_AAAdd(interactionBetas, iCurr+1);
        			//	Functions for manipulating positionMatrix object:
        			//	maskInteractionMatrix(ArrayList<ArrayList<double[][]>> interactionAlphas)
        			//	maskInteractionPackingMatrix(ArrayList<ArrayList<int[][]>> interactionPacking)
        			//	makeInteractionPackingHomogenous(ArrayList<ArrayList<int[][]>> interactionPacking)
        			//	packModel_component(JSONObject packing, int iFirst)

        			oIntCoef.remove("positionMatrix");
        			oIntSett.put("positionBias", false);
        		}
        		
        		//2.3 Convert the homogenous positionMatrix into a 'spacing vector'.
        		// - This should be implemented in the class 'bindingModeInteraction'.

        	}*/
    	}
    	
    	//3. For each experiment and round, shift the binding-mode activities so that the largest is zero.
    	for(int iExp=0; iExp<nActivity.size(); iExp++) {
    		for(int iRound=0; iRound<nActivity.get(iExp); iRound++) {
    			//3.1 Finds the maximum activity.
    			double maxAct = act.get(0).get(iExp)[iRound];
    			for(int iBM=1; iBM<act.size(); iBM++) 
    				maxAct = Math.max(maxAct, act.get(iBM).get(iExp)[iRound]);

    			//Subtracts the maximum activity from all binding mode and binding-mode interaction activities:
    			for(int iBM=0; iBM<act.size(); iBM++) 
    				act.get(iBM).get(iExp)[iRound] -= maxAct;
    			if(intAct!=null)
    				for(int iInt=0; iInt<intAct.size(); iInt++)
    					intAct.get(iInt).get(iExp)[iRound] -= maxAct;
    			
    		}
    	}
    			
    	//4. For each binding mode, compute the mean exp-activity across the experiments/round, ignoring round zero.
    	for(int iBM=0; iBM<act.size(); iBM++) {
    		
    		
    		//4.1 Sums over activities, ignoring round zero.
    		int nVal      = 0;
    		double expSum = 0;
        	for(int iExp=0; iExp<nActivity.size(); iExp++) {

        		
        		for(int iRound=0; iRound<nActivity.get(iExp); iRound++) {
        			//Ignores the first round of each SELEX experiment
        			if(		iRound==0 &&
        					oSett.has("enrichmentModel") && 
        					oSett.getJSONArray("enrichmentModel").getJSONObject(iExp).getString("modelType").equals("SELEX")) {
        				
        				if(nActivity.get(iExp)==1)
            				throw new RuntimeErrorException(null, "Error: All SELEX experiments must have at least two binding-mode activities recorded, but this is not the casea.");
        				else
        					continue;
        				
        			}
        			
        			expSum += Math.exp(act.get(iBM).get(iExp)[iRound]);
        			nVal   += 1;
        			
        		}
        	}
        	
        	//4.2 Sets a new value in the binding model.
        	JSONArray newAct = new JSONArray();
        	JSONArray temp = new JSONArray();
        	newAct.put(temp);
        	temp.put(Math.log(expSum/nVal));
        	oCoef.getJSONArray("bindingModes").getJSONObject(iBM).put("activity", newAct);
        	
    	}
    	
    	
    	//5. Take the mean exp-interaction-activity of experiment/rounds, ignoring round zero.
/*    	if(intAct!=null) {
    		for(int iInt=0; iInt<intAct.size(); iInt++) {
        		//5.1 Sums over activities, ignoring round zero.
        		int nVal      = 0;
        		double expSum = 0;

        		for(int iExp=0; iExp<nActivity.size(); iExp++) {

        			for(int iRound=0; iRound<nActivity.get(iExp); iRound++) {
        				//Ignores the first round of each SELEX experiment
        				if(		iRound==0 &&
        						oSett.has("enrichmentModel") && 
        						oSett.getJSONArray("enrichmentModel").getJSONObject(iExp).getString("modelType").equals("SELEX")) {
        					continue;

        				}

        				expSum += Math.exp(intAct.get(iInt).get(iExp)[iRound]);
        				nVal   += 1;
        			}
        		}
        		
            	//5.2 Sets a new value in the binding model.
            	JSONArray newAct = new JSONArray();
            	JSONArray temp = new JSONArray();
            	newAct.put(temp);
            	temp.put(Math.log(expSum/nVal));
            	oCoef.getJSONArray("bindingModeInteractions").getJSONObject(iInt).put("activity", newAct);
        	}
    	}*/

    	//6) Transforms k=0,1 binding modes into a scoring matrix with beta=0.
    	//6.1) Determines the size of the scoring matrix: finds the smallest matrix with k>1
    	int newSize = -1, newFlankLength = 0;
    	for(int iBM=0; iBM<nModes; iBM++) {
    		JSONObject oBMSett  = oSett.getJSONArray("bindingModes").getJSONObject(iBM);
    		int kCurr = oBMSett.getInt("size");
    		if(kCurr>1) {
				if(kCurr<newSize || newSize==-1) {
					newSize        = kCurr;
					newFlankLength = oBMSett.getInt("flankLength");
				} else if(kCurr==newSize) {
					newFlankLength = Math.max(newFlankLength, oBMSett.getInt("flankLength"));
				}
    		}
    	}
    	
    	//6.2) Transforms the binding modes (only possible if enrichment models and count tables are defined).
    	if(newSize!=-1 && oSett.has("countTable") && oSett.has("enrichmentModel")) {
    		JSONArray aTableSett = oSett.getJSONArray("countTable");
    		JSONArray aEnrSett   = oSett.getJSONArray("enrichmentModel");

        	for(int iBM=0; iBM<nModes; iBM++) {
        		
        		JSONObject oBMSett   = oSett.getJSONArray("bindingModes").getJSONObject(iBM);
        		JSONObject oBMCoef   = oCoef.getJSONArray("bindingModes").getJSONObject(iBM);
        		
            	int oldSize          = oBMSett.getInt("size");
            	int oldFlankLength   = oBMSett.getInt("flankLength");

        		double betaShift     = 0;
        		
        		// 6.2.1 Skips binding modes larger than 1pb
        		if(oldSize>1)
        			continue;

        		// 6.2.2 Computes the number of windows for the current binding modes.
        		//Creates a list of the number of binding windows.
        		ArrayList<Double> newFrames_A = new ArrayList<Double>();
        		ArrayList<Double> oldFrames_A = new ArrayList<Double>();

        		//Computes the average variableRegionLength in the experiments where this binding modes is present. 
        		for(int iTable=0; iTable<aEnrSett.length(); iTable++) {
        			JSONObject oTableSett     = aTableSett.getJSONObject(iTable);
        			JSONObject oEnrSett       = aEnrSett.getJSONObject(  iTable);

        			//Checks if the current binding mode is included in the current experiment
        			boolean bmIncluded = false;
        			JSONArray inclBMs  = oEnrSett.getJSONArray("bindingModes");
        			if(inclBMs.length()==1&&inclBMs.getInt(0)==-1)
        				bmIncluded = true;
        			else 
        				for(int iIncl=0; iIncl<inclBMs.length(); iIncl++)
        					if(inclBMs.getInt(iIncl)==iBM)
        						bmIncluded = true;
        			if(!bmIncluded) {
        				continue;
        			}

        			//Saves the length of the variable region.
        			int varLen = oTableSett.getInt("variableRegionLength");
        			oldFrames_A.add((double) (oldSize == 0 ? 1 : ( 2 * oldFlankLength + varLen - oldSize + 1 ) * (oBMSett.getBoolean("singleStrand") ? 1 : 2)) );
        			newFrames_A.add((double) (                   ( 2 * newFlankLength + varLen - newSize + 1 ) *                                           2 ) );
        		}
        		
        		//Skips the current binding mode if there is no experiment for which the number of binding windows is defined.
    			if(newFrames_A.size()==0) {
            		continue;
    			}
    			
    			//Computes the average number of frames across all experiments for which it is defined.
        		double oldNFrames    = 0;
    			for(double ofi : oldFrames_A)
    				oldNFrames += ofi;
    			oldNFrames /= oldFrames_A.size();

    			double newNFrames    = 0; 
    			for(double nfi : newFrames_A)
    				newNFrames += nfi;
    			newNFrames /= newFrames_A.size();

        		//6.2.3 Corrects for the change in the number of binding frames.
    			
        		betaShift     += Math.log(oldNFrames/newNFrames);
        		
        		//6.2.4 If size=-1, computes the expected value of the bias term.
        		if(oldSize==1) {
        			betaShift += Math.log(Array.mean(ModelComponent.exp_d(ModelComponent.readFromJSON_d(oBMCoef.getJSONArray("mononucleotide")))));
        		}
        		
        		
        		
        		//6.2.5 Updates the binding mode.
        		oBMSett.put("size",           newSize);
        		oBMSett.put("flankLength",    newFlankLength);
        		oBMSett.put("singleStrand",   false);
        		oBMCoef.put("mononucleotide", ModelComponent.JSONArrayConvert_d(new double[nMono*newSize]));
        		
            	
        		oBMCoef.getJSONArray("activity").getJSONArray(0).put(0, 
        				oBMCoef.getJSONArray("activity").getJSONArray(0).getDouble(0) + betaShift);
            	

        	}    		
    	}

    	//7. Cleans up the scoring model
    	cleanUpScoringModel();

		
		//8. Rebuilds the sliding windows etc.
		loadOModel();
		
    }    
	
	// Cleans up the scoring model JSON Object (removes elements specific to fits.)
	public void cleanUpScoringModel() {
		
    	JSONObject oCoef = oModel.getJSONObject("coefficients");
    	JSONObject oSett = oModel.getJSONObject("modelSettings");

		// 1) Puts the flank length to zero. 
    	for(int iBM=0; iBM<nModes; iBM++) 
    		oModel.getJSONObject("modelSettings").getJSONArray("bindingModes").getJSONObject(iBM).put("flankLength", 0);
    	
    	// 2) Removes the count table and enrichment model.
    	if(oCoef.has("enrichmentModel"))
    		oCoef.remove("enrichmentModel");
    	if(oCoef.has("countTable"))
    		oCoef.remove("countTable");

    	if(oSett.has("enrichmentModel"))
    		oSett.remove("enrichmentModel");
    	if(oSett.has("countTable"))
    		oSett.remove("countTable");

    	JSONObject oConstr = oModel.getJSONObject("modelFittingConstraints");
    	if(oConstr.has("enrichmentModel"))
    		oConstr.remove("enrichmentModel");
    	if(oConstr.has("countTable"))
    		oConstr.remove("countTable");
    	
    	// 3) Variation descriptions 
    	for(int iBM=0; iBM<nModes; iBM++) {
    		JSONObject oBMSett  = oModel.getJSONObject("modelSettings").getJSONArray("bindingModes").getJSONObject(iBM);
    		String[] removeKeys = {"variationsOptimized", "variationDescription", "variationName", "freezingLevel", "includeComponent"};
    		for(String key: removeKeys) 
    			if(oBMSett.has(key))
    				oBMSett.remove(key);
    	}
    	
    	// 4) Removes optimizer settings
		if(oModel.has("optimizerSetting"))
			oModel.remove("optimizerSetting");
		
		// 5) Removes model-fitting constraints:
		if(oModel.has("modelFittingConstraints"))
			oModel.remove("modelFittingConstraints");
		
		// 6) Removes metadata:
		if(oModel.has("metadata"))
			oModel.remove("metadata");

	}
    
	// Counts kMers and writes the counts to a TSV file
	///////////////////////////////////////////////////
	//If table numbers are non-zero present, they are assumed to be counts of the sequence.
	//If no table numbers are present, we the each window has weight 1
	public void kMerCount(ArrayList<String> args) {
		if(args==null) {
			System.out.println("kMerCount(outTSV, k=8)                                       - Tabulates the occurrences of k-mer with length 'k' in the input sequence stream and stores the result in 'outTSV'");
			return;
		}
		int nArgs = args.size();
		//DoStuff

		if(nArgs>0 && nArgs<=3) {
			String outTSV = args.get(0);
			int k = 8;
			boolean forwardOnly = false;
			if(nArgs>= 2)
				k = Integer.parseInt(args.get(1));
			if(nArgs >= 3)
				forwardOnly = Boolean.parseBoolean(args.get(2));
			
			kMerCount(outTSV, k, forwardOnly);
		} else {
			error("The function kMerCount(outTSV, k=8) must have exactly one arguments.");
		}
	}
	
	
	public void kMerCount(String outTSV, int kIn, boolean forwardOnly) {
		
		if(letterComplement==null)
			setupAlphabet("C-G,A-T", "ACGT");
		
		if(kIn!=k)
			createKMerList(kIn);
		
		//Loops over sequences
		double[][] kMerCounts = null; //This is allocated once we know what columns are present. 
		long nKmers       = kMers.size();
		while(nextSequence()) {
			do {
				if(kMerCounts==null) {
					if(tableData != null)
						kMerCounts = new double[tableData.length][kMers.size()];
					else
						kMerCounts = new double[1][kMers.size()];
				}
				
				if(forwardOnly) {
					
					for(int x=0; x<sequence.length()-k+1; x++) {
						String subSeq        = sequence.substring(x,x+k);
						if( kMerIndex.containsKey(subSeq) ) {
							if(tableData==null)
								//Simply sums the number of windows if weights were given in the table
								kMerCounts[0][kMerIndex.get(  subSeq)]++;
							else
								//Otherwise adds the value in the table 
								for(int iCol=0; iCol<tableData.length; iCol++) 
									kMerCounts[iCol][kMerIndex.get(  subSeq)] += tableData[iCol];
						}
						firstNeededCharacter = sequence.length()-k+1;
			    	}
					
				} else {
					
					for(int x=0; x<sequence.length()-k+1; x++) {
						String subSeq        = sequence.substring(x,x+k);
						if( kMerIndex.containsKey(subSeq) ) {
							if(tableData==null) {
								//Simply sums the number of windows if weights were given in the table
								kMerCounts[0][kMerIndex.get(  subSeq)]++;
								kMerCounts[0][kMerIndexRC.get(subSeq)]++;
							} else {
								//Otherwise adds the value in the table 
								for(int iCol=0; iCol<tableData.length; iCol++) {
									kMerCounts[iCol][kMerIndex.get(  subSeq)] += tableData[iCol];
									kMerCounts[iCol][kMerIndexRC.get(subSeq)] += tableData[iCol];
								}
							}
						}
						firstNeededCharacter = sequence.length()-k+1;
			    	}
					
				}
				
				
			} while(nextSubstring());
		}

		//Writes output file.
		if(verbose)
			System.out.println("> Writes k-mer counts "+outTSV);
		try {
			BufferedWriter writer = new BufferedWriter(new FileWriter(outTSV, false));
			
			for(int iKmer=0; iKmer<nKmers; iKmer++) {
				writer.write(kMers.get(iKmer));
				for(int iCol=0; iCol<kMerCounts.length; iCol++)
					writer.write("\t"+kMerCounts[iCol][iKmer]);
				writer.write("\n");

			}
			writer.close();
		} catch(Exception e) {
		    System.err.println("Error: Cannot write to output file.");
		    e.printStackTrace();
		} 
		
		restartFile();
	}
	
	
	
	public void createProBoundObjects(ArrayList<String> modifications) {
		
		// 1) Reads the sequences and creates a ArrayList<LongSeequence> and a zero
		if(verbose)
			System.out.println("> Reading sequences.");

		if(sequenceFile == null)
			throw new RuntimeErrorException(null, "Error: No sequence file has been specified.");
		
		//Creates a count table object using the input stream.
		ArrayList<int[]> counts        = new ArrayList<int[]>();
		ArrayList<LongSequence> probes = new ArrayList<LongSequence>();
		Integer sequenceLength         = null;
		while(nextSequence()) {
			//Reads the full sequence.
			String fullSeq = "";
			do {
				fullSeq += sequence;
				firstNeededCharacter = sequence.length();
			} while(nextSubstring());

			if(sequenceLength==null) 
				sequenceLength = fullSeq.length();
			else if(sequenceLength!=fullSeq.length())
				throw new RuntimeErrorException(null, "Error: All sequences must have the same length. The first sequence has L="
						+sequenceLength+" but the current sequence ("+fullSeq+") has length "+fullSeq.length()+".");

			
			probes.add(sc.build(fullSeq));
			counts.add(new int[1]);
		}
		dataTable              = new MultiRoundData("", "", counts, probes, new int[1], sc);
		if(verbose)
			System.out.println(">> Number of sequences: "+probes.size());
		
		// 2) Create a config/model
		// 2.1) Clone settings
		if(verbose)
			System.out.println("> Checking the model.");

		JSONObject configClone = ModelComponent.clone_JSON_O(oModel);
    	JSONObject oCoef       = configClone.getJSONObject("coefficients");
    	JSONObject oSett       = configClone.getJSONObject("modelSettings");
    	
    	//Gets the alphabet
		letterComplement      = oSett.has("letterComplement") ? oSett.getString("letterComplement")       : "C-G,A-T";
		letterOrder           = oSett.has("letterOrder")      ? oSett.getString("letterOrder") : CombinedLikelihood.alphabetDefToOrderedLetters(letterComplement);
		setupAlphabet(letterComplement, letterOrder);
		
		
		// 2.2) Checks the model
		// 2.2.1) Checks the binding modes.
		JSONArray aBM = oCoef.getJSONArray("bindingModes");
		for(int iBM=0; iBM<aBM.length(); iBM++) {
			JSONObject oBM = aBM.getJSONObject(iBM);
			if(!oBM.has("activity")) 
				throw new RuntimeErrorException(null, "Error: Binding mode "+iBM+" does not have activities defined. Run buildConsensusModel().");
			else if(oBM.getJSONArray("activity").length()!=1)
				throw new RuntimeErrorException(null, "Error: Binding mode "+iBM+" has activities specified for multiple experiments. Run buildConsensusModel().");
			else if(oBM.getJSONArray("activity").getJSONArray(0).length() != 1) 
				throw new RuntimeErrorException(null, "Error: Binding mode "+iBM+" has multiple activities for the first experiments. Run buildConsensusModel().");
			if(oBM.has("positionBias") && oBM.getJSONArray("positionBias").length()>0)
				throw new RuntimeErrorException(null, "Error: Binding mode "+iBM+" has position bias. Run buildConsensusModel().");
		}
		// 2.2.2) Checks the binding mode interactions
		JSONArray aInt = oCoef.getJSONArray("bindingModeInteractions");
		for(int iInt=0; iInt<aInt.length(); iInt++) {
			JSONObject oInt = aInt.getJSONObject(iInt);
			if(!oInt.has("activity")) 
				throw new RuntimeErrorException(null, "Error: Binding mode interaction "+iInt+" does not have activities defined. Run buildConsensusModel().");
			else if(oInt.getJSONArray("activity").length()!=1)
				throw new RuntimeErrorException(null, "Error: Binding mode interaction "+iInt+" has activities specified for multiple experiments. Run buildConsensusModel().");
			else if(oInt.getJSONArray("activity").getJSONArray(0).length() != 1) 
				throw new RuntimeErrorException(null, "Error: Binding mode interaction "+iInt+" has multiple activities for the first experiments. Run buildConsensusModel().");
			if(oInt.has("positionMatrix"))
				throw new RuntimeErrorException(null, "Error: Binding mode interaction "+iInt+" has a positionMatrix. Run buildConsensusModel().");
		}
		// 2.2.3) Checks so the enrichment model and count tables has been removed:
    	if(oCoef.has("enrichmentModel") || oSett.has("enrichmentModel"))
			throw new RuntimeErrorException(null, "Error: The model has coefficients/settings for the enrichment model. Run buildConsensusModel().");
    	if(oCoef.has("countTable")      || oSett.has("countTable"))
			throw new RuntimeErrorException(null, "Error: The model has coefficients/settings for the count table. Run buildConsensusModel().");

    	
		// 2.3) Create a bare-bones enrichment model and count table
    	// 2.3.1) CountTable
    	
    	// 2.3.1.1) Coefficients
    	JSONObject oCount = new JSONObject();
    	JSONArray aCount  = new JSONArray();
    	oCount.put("h",                    ModelComponent.JSONArrayConvert_d(new double[1]) );
    	aCount.put(oCount);
    	oCoef.put("countTable",            aCount);
    	
    	// 2.3.1.2) Settings
    	oCount            = new JSONObject();
    	oCount.put("countTableFile",       "N/A");
    	oCount.put("inputFileType",        "N/A");
    	oCount.put("rightFlank",           "");
    	oCount.put("leftFlank",            "");
    	oCount.put("nColumns",             1);
    	oCount.put("variableRegionLength", probes.get(0).getLength());
    	aCount  = new JSONArray();
    	aCount.put(oCount);
    	oSett.put("countTable", aCount);
    	
    	// 2.3.1) Enrichment model
    	
    	// 2.3.2.1) Coefficients
    	JSONObject oEnr = new JSONObject();
    	JSONArray aEnr  = new JSONArray();
    	oEnr.put("rho",   ModelComponent.JSONArrayConvert_i(new int[1]));
    	oEnr.getJSONArray("rho").put(0, 1);
    	oEnr.put("gamma", ModelComponent.JSONArrayConvert_i(new int[1]));
    	aEnr.put(oEnr);
    	oCoef.put("enrichmentModel", aEnr);
		
    	// 2.3.2.2) Settings
    	oEnr            = new JSONObject();
    	int[] tempNeg   = new int[1];
    	tempNeg[0]      = -1;
    	oEnr.put("modelType",               "RhoGamma");
    	JSONArray aMod = new JSONArray();
    	if(modifications!=null)
    		for(String mod : modifications)
    			aMod.put(mod);
    	oEnr.put("modifications",           aMod);
    	oEnr.put("bindingModes",            ModelComponent.JSONArrayConvert_i(tempNeg));
    	oEnr.put("bindingModeInteractions", ModelComponent.JSONArrayConvert_i(tempNeg));
    	aEnr  = new JSONArray();
    	aEnr.put(oEnr);
    	oSett.put("enrichmentModel", aEnr);

		// 3) Create ProBound objects
		if(verbose)
			System.out.println("> Builds model-component objects.");
    	componentList                                  = new ArrayList<ModelComponent>();
		// 3.1) CounTable
		countTable                                     = new CountTable(configClone, 0, sc, letterComplement, letterOrder, false);
		countTable.setNThreads(1);           //Sets the number of threads
		countTable.loadDataTable(dataTable); //Loads the data
		componentList.add(countTable);
		
		// 3.2) Add binding modes.
		bindingModes                                   = new ArrayList<BindingMode>();
		int nBindingModes                              = configClone.getJSONObject("modelSettings").getJSONArray("bindingModes").length();
		for(int iBM=0; iBM<nBindingModes; iBM++) {
			BindingMode newBM = new BindingMode(configClone, iBM, sc, letterComplement, letterOrder);
			newBM.addCountTable(countTable);
			newBM.setComponentInclusion(true);
			bindingModes.add(   newBM);
			componentList.add(  newBM);
		}
		
		// 3.3) Build binding model interactions.
		interactions = new ArrayList<BindingModeInteraction>();
		int nInteractions                              = configClone.getJSONObject("modelSettings").getJSONArray("bindingModeInteractions").length();
		for(int iInt=0; iInt<nInteractions; iInt++) {
			BindingModeInteraction newInt = new BindingModeInteraction(configClone, iInt, bindingModes);
			newInt.addCountTable(countTable);
			interactions.add(    newInt);
			componentList.add(   newInt);
		}
		
    	// 3.4). Build the enrichment models.
		enrichmentModel                  = EnrichmentModel.buildEnrichmentModel(configClone, 0, countTable, bindingModes, interactions);
		enrichmentModel.includeComponent = true;
		countTable.setEnrichmentModel(enrichmentModel);
		componentList.add(            enrichmentModel);

		// 4) Sets the coefficients in the model
		if(verbose)
			System.out.println("> Sets the coeffficients.");
		ModelComponent.readFromJSON(configClone, "coefficients", componentList);
		
	}
	
	public void affinitySum(ArrayList<String> args) {
		if(args==null) {
			System.out.println("affinitySum(outFile, modifications='')                                - Computes the (activity-weighted) affinity sums given a stream of sequences");
			return;
		}
		int nArgs = args.size();

		if(nArgs==1 || nArgs==2) {
			String outFile = args.get(0);
			ArrayList<String> modifications = new ArrayList<String>();
			if(nArgs == 2 && args.get(1).length()>0) 
				for(String s : args.get(1).split(","))
					modifications.add(s);
			affinitySum(outFile, modifications);
		} else {
			error("The function affinitySum(outFile, modifications='') must have exactly one argument.");
		}
		
	}
	
	public void affinitySum(String outFile, ArrayList<String> modifications) {
		
		//TODO: Currently the non-specific binding is added once per probe. 
		

		// 1) Creates proBound Objects.
		createProBoundObjects(modifications);
		
		// 2) Loop over sequences in MultiRoundData object, score and write.
		if(verbose)
			System.out.println("> Scores the sequences..");
		countTable.writeAlphaTable(outFile);

	}

	public void bindingModeScores(ArrayList<String> args) {
		if(args==null) {
			System.out.println("bindingModeScores(outFile, sum/max/mean/profile, modifications)       - Scores each sequence in the input stream using all loaded binding modes. The output file 'outFile' contains, for each probe, etiher the affinity-sum (Defualt), the maximum score, the mean score, or the score at each offset and strand.");
			return;
		}
		int nArgs = args.size();
		//DoStuff
		

		if(1 <= nArgs && nArgs <=3) {
			String outFile                  = args.get(0);
			String outFormat                = nArgs>=2 ? args.get(1) : "sum";
			ArrayList<String> modifications = new ArrayList<String>();
			if(nArgs == 3 && args.get(2).length()>0) 
				for(String s : args.get(2).split(","))
					modifications.add(s);
			bindingModeScores(outFile, outFormat, modifications);
		} else {
			error("The function bindingModeScores(outFile, sum/max/mean/profile, modifications) must have one, two, or three arguments.");
		}
	}
	
	
	public void bindingModeScores(String outFile, String outFormat, ArrayList<String> modifications) {

		//1.0 Creates proBound Objects.
		createProBoundObjects(modifications);

		//Compute binding modes scores
		countTable.writeBindingModeAlphas(outFile, outFormat);
	}
	
	public void kMerAverage(ArrayList<String> args) {
		if(args==null) {
			System.out.println("kMerAverage(outTable, k=8, forwardOnly=False, pseudocount=0) - Given an input stream that is a table with one sequence column and additional numeric columns, the mean numeric value is computed for each k-mer.");
			return;
		}
		int nArgs = args.size();

		if(nArgs>=1 && nArgs<=4) {
			String outTable = args.get(0);
			int k               = 8;
			boolean forwardOnly = false;
			int pseudocount     = 0;
			if(nArgs >= 2)
				k           = Integer.parseInt(    args.get(1));
			if(nArgs >= 3)
				forwardOnly = Boolean.parseBoolean(args.get(2));
			if(nArgs >= 4)
				pseudocount = Integer.parseInt(    args.get(3));
			kMerAverage(outTable, k, forwardOnly, pseudocount);
		} else {
			error("The function kMerAverage(outTable, k=8, forwardOnly=False, pseudocount=0) must have between one and four arguments.");
		}
	}
	
	
	public void kMerAverage(String outTable, int kIn, boolean forwardOnly, int pseudocount) {
		
		if(kIn!=k)
			createKMerList(kIn);
		
		//int iCurrent                = 0;
		//int currentNFrames;
		iBindingMode                = new ArrayList<Integer>();
		iBindingMode.add(0);
		nFrames                     = new ArrayList<Integer>();
		int L                       = 0;
		int nCols                   = 0;
		int nKMerFrames             = 0;
		double[][] binnedSignal     = null;
		double[] kMerCounts         = null;
		int[] iKmer                 = null;
		int iSeq                    = 0;

		while(nextSequence()) {

			if(tableData==null || tableData.length==0)
				throw new RuntimeErrorException(null, "Error: Cannot run kMerPositionBias since there is no columnd data.");
				
			//Reads all substrings to build a full sequence 
			while(nextSubstring());
			iSeq++;

			if(iSeq==1) {

				// Sets up the variables at the first sequence when its length is known. 
				L     = sequence.length();
				nCols = tableData.length;

				//Counts the k-mers. 
				if(forwardOnly)
					nKMerFrames      = (L-k+1);
				else
					nKMerFrames      = 2*(L-k+1);

				if(verbose)
					System.out.println("Number of k-mers: "+nKmers);

				binnedSignal     = new double[nKmers][nCols];
				kMerCounts       = new double[nKmers];
				iKmer            = new int[nKMerFrames];

			} else {
				if(L != sequence.length())
					throw new RuntimeErrorException(null, "Error: All sequences in the table must have the same length.");
			}


			//Records the indices of the k-mers that are present
			for(int x=0; x<L-k+1; x++) {
				iKmer[      x]       = kMerIndex.get(  sequence.substring(x,x+k));
				if(!forwardOnly) {
					iKmer[x+L-k+1]   = kMerIndexRC.get(sequence.substring(x,x+k));
				}
			}

			//Computes the k-mer frequency and k-mer-binned signal.
			for(int iK=0; iK<nKMerFrames; iK++) {
				for(int iCol=0; iCol<nCols; iCol++) 
					binnedSignal[iKmer[iK]][iCol] += tableData[iCol];
				kMerCounts[ iKmer[iK]] += 1;
			}
		}
		
		//Computes the column average across k-mers
		double[] colAverage = new double[nCols];
		double nKTot = 0;
		for(int iK=0; iK<nKmers; iK++) { 
			nKTot += kMerCounts[iK];
			for(int iCol=0; iCol<nCols; iCol++)
				colAverage[iCol] += binnedSignal[iK][iCol];
		}
		for(int iCol=0; iCol<nCols; iCol++)
			colAverage[iCol] /= nKTot;
		
		
		//Writes the predicted data
		if(verbose)
			System.out.println("> Writes k-mer table averages to "+outTable);
		try {
			BufferedWriter writer = new BufferedWriter(new FileWriter(outTable, false));
			
			for(int iK=0; iK<nKmers; iK++) {
				int iData = 0;
				for(int iCol=0; iCol<nCols+1; iCol++) {
					if(iCol==stringColumn) {
						writer.write(kMers.get(iK));
					} else {
						writer.write(""+((binnedSignal[iK][iData] + pseudocount * colAverage[iData]) / (kMerCounts[iK]+pseudocount)));
						iData++;
					}
					
					if(iCol==nCols)
						writer.write("\n");
					else
						writer.write("\t");
				}
			}
			
			writer.close();
			
		} catch(Exception e) {
		    System.err.println("Error: Cannot write to output file.");
		    e.printStackTrace();
		} 
		
		restartFile();
	}
	

	
	
	public void kMerMedian(String outTable, int kIn) {
		
		if(letterComplement==null || letterOrder ==null) {
			letterOrder = "ACGT";
			letterComplement = "A-T,C-G";
			setupAlphabet(letterComplement, letterOrder);
		}
		
		if(kIn!=k)
			createKMerList(kIn);
		
		iBindingMode                = new ArrayList<Integer>();
		iBindingMode.add(0);
		nFrames                     = new ArrayList<Integer>();
		int L                       = 0;
		int nCols                   = 0;
		ArrayList<ArrayList<ArrayList<Double>>> bins = new ArrayList<ArrayList<ArrayList<Double>>>();
		int iSeq                    = 0;

		while(nextSequence()) {

			if(tableData==null || tableData.length==0)
				throw new RuntimeErrorException(null, "Error: Cannot run kMerPositionBias since there is no columnd data.");
				
			//Reads all substrings to build a full sequence 
			while(nextSubstring());
			iSeq++;

			if(iSeq==1) {

				// Sets up the variables at the first sequence when its length is known. 
				L     = sequence.length();
				nCols = tableData.length;

				//Counts the k-mers. 
				if(verbose)
					System.out.println("Number of k-mers: "+nKmers);
				
				//Allocates lists for binning signal.
				for(int i = 0; i<nKmers; i++) {
					ArrayList<ArrayList<Double>> temp = new ArrayList<ArrayList<Double>>();
					bins.add(temp);
					for(int iCol=0; iCol<nCols; iCol++)
						temp.add(new ArrayList<Double>());
				}
						
			} else {
				if(L != sequence.length())
					throw new RuntimeErrorException(null, "Error: All sequences in the table must have the same length.");
			}

			//Identifies the k-mer indices and records the column value for all relevant values. 
			for(int x=0; x<L-k+1; x++) {
				int iKmerF = kMerIndex.get(  sequence.substring(x,x+k));
				int iKMerR = kMerIndexRC.get(sequence.substring(x,x+k));
				for(int iCol=0; iCol<nCols; iCol++) {
					bins.get(iKmerF).get(iCol).add(tableData[iCol]);
					bins.get(iKMerR).get(iCol).add(tableData[iCol]);
				}
			}
		}
		
		//Writes the predicted data
		if(verbose)
			System.out.println("> Writes k-mer table averages to "+outTable);
		try {
			BufferedWriter writer = new BufferedWriter(new FileWriter(outTable, false));
			
			for(int iK=0; iK<nKmers; iK++) {
				for(int iCol=0; iCol<nCols+1; iCol++) {
					if(iCol==stringColumn)
						writer.write(kMers.get(iK));
					else {
						if(bins.get(iK).get(iCol).size()>0)
							writer.write(""+Misc.median(bins.get(iK).get(iCol)));
						else 
							writer.write("NA");
					}
					if(iCol==nCols)
						writer.write("\n");
					else
						writer.write("\t");
				}
			}
			
			writer.close();
			
		} catch(Exception e) {
		    System.err.println("Error: Cannot write to output file.");
		    e.printStackTrace();
		} 
		
		restartFile();
	}
	
	public void predictCountTable(ArrayList<String> args) {
		if(args==null) {
			System.out.println("predictCountTable(predictedTable, predictCountTable=null, iTable=0)   - Predicts the count table with index 'iTable' using a full, loaded ProBound model. A different count table than that specified by the ProBound model can be predicted using 'predictCountTable'.");
			return;
		}
		int nArgs = args.size();

		if(nArgs==1 || nArgs==2 || nArgs==3) {
			String predictedTable  = args.get(0);
			String inputCountTable = null;
			int iTable             = 0;
			if(nArgs>=2) {
				inputCountTable    = args.get(1);
				if(inputCountTable.toLowerCase().equals("null"))
					inputCountTable=null;
			}
			if(nArgs==3)
				iTable             = Integer.parseInt(args.get(2));

			predictCountTable(predictedTable, inputCountTable, iTable);
		} else {
			error("The function predictCountTable(predictedTable, predictCountTable=null, iTable=0) must have exactly one or two arguments.");
		}
	}
	
	
	public void predictCountTable(String predictedTable, String tableFile, int iTable) {
		String proBoundDir = System.getenv("PROBOUND_DIR");
		String generalSchemaFile  = proBoundDir + "/config/" + "schema.general.json";
		//String configSchemaFile   = proBoundDir + "/config/" + "schema.config.json";
		
        if(verbose) {
        	System.out.println("predictedTable = "+predictedTable);
        	System.out.println("tableFile      = "+tableFile);
        	System.out.println("iTable         = "+iTable);
        	
        }
        if(tableFile!=null)
        	oModel.getJSONObject("modelSettings").getJSONArray("countTable").getJSONObject(iTable).put("countTableFile", tableFile);
        
        JSONModel.validateSchemaFile_O(generalSchemaFile, oModel);
        
        //Creates the likelihood model.
        CombinedLikelihood l = new CombinedLikelihood(oModel);
        l.setStateJSON(oModel);
        
        if(verbose) {
        	System.out.println("modelSettings: "+oModel.getJSONObject("modelSettings").toString());
        	
        	System.out.println("coefficients:");
        	JSONModel.printJSONState(l.getJSONState(), "coefficients");
        }

        //Includes all components.
        for(int iComp=0; iComp<l.componentList.size(); iComp++)
        	l.componentList.get(iComp).setComponentInclusion(true);

        l.tableModels.get(iTable).writePCTable(predictedTable);
        
	}

	// Counts the nucleotides across positions.
	///////////////////////////////////////////////////
	//If table numbers are non-zero present, they are assumed to be counts of the sequence.
	//If no table numbers are present, we the each window has weight 1
	public void positionBaseCount(ArrayList<String> args) {
		if(args==null) {
			System.out.println("positionBaseCount(countFile)                                 - Counts the letters at each position in the input data stream (weighted by column values, if available).");
			return;
		}
		int nArgs = args.size();
		

		if(nArgs==1 ) {
			String countFile  = args.get(0);
			positionBaseCount(countFile);
		} else {
			error("The function positionBaseCount(countFile) must have exactly one arguments.");
		}
	}
	
	
	public void positionBaseCount(String outTSV) {
		
		if(letterComplement==null)
			setupAlphabet("C-G,A-T", "ACGT");
		
		//Loops over sequences
		double[][][] positionBaseTensor      = null; //This is allocated once we know what columns are present.
		int nBases                           = letterOrder.length();
		int L                                = sequence.length();
		HashMap<Character,Integer> alphabet  = sc.getALPHABET_LOOKUP();
		HashMap<Integer,Character> alphabetR = sc.getALPHABET_LOOKUP_REVERSE();

		while(nextSequence()) {
			do {
				if(positionBaseTensor==null) {
					L = sequence.length();
					if(tableData != null)
						positionBaseTensor = new double[tableData.length][L][nBases];
					else
						positionBaseTensor = new double[               1][L][nBases];
				}
				
				for(int x=0; x<L; x++) {
//					System.out.println(alphabet);
					int iBase = alphabet.get(sequence.charAt(x));
					
					if(tableData==null)
						//Simply sums the number of windows if weights were given in the table
						positionBaseTensor[       0][x][iBase]++;
					else
						//Otherwise adds the value in the table 
						for(int iCol=0; iCol<tableData.length; iCol++) 
							positionBaseTensor[iCol][x][iBase] += tableData[iCol];
					
					firstNeededCharacter = sequence.length();
		    	}
				
			} while(nextSubstring());
		}

		
		//Writes output file.
		if(verbose)
			System.out.println("> Writes k-mer counts "+outTSV);
		try {
			BufferedWriter writer = new BufferedWriter(new FileWriter(outTSV, false));
			int nCol              = positionBaseTensor.length;
			for(int x=0; x<L; x++) {
				//Writes first columns
				for(int n=0; n<nBases; n++) {
					if(n>0)
						writer.write(',');
					writer.write(alphabetR.get(n)+":"+x);
				}
				
				for(int iCol=0; iCol<nCol; iCol++) {
					writer.write('\t');
					for(int iBase=0; iBase<nBases; iBase++) {
						if(iBase>0)
							writer.write(',');
						writer.write(""+positionBaseTensor[iCol][x][iBase]);
					}
			
				}
				writer.write('\n');
			}
			writer.close();
		} catch(Exception e) {
		    System.err.println("Error: Cannot write to output file.");
		    e.printStackTrace();
		} 
		
		restartFile();
	}
	

	public void guessAlphabet(ArrayList<String> args) {
		if(args==null) {
			System.out.println("guessAlphabet()                   - Guesses the alphabet.");
			return;
		}
		int nArgs = args.size();

		if(nArgs==0) {
			guessAlphabet();
		} else {
			error("The function guessAlphabet() must have exactly zero arguments.");
		}
	
		
	}
	
	
	public void guessAlphabet() {
		
		HashSet<Character> chars = new HashSet<Character>();
		HashSet<Character> acgt  = new HashSet<Character>();
		acgt.add('A');
		acgt.add('C');
		acgt.add('G');
		acgt.add('T');
		while(nextSequence()) {
			do {
				int L = sequence.length();
				for(int x=0; x<L; x++) {
				chars.add(sequence.charAt(x));
				}
				
			} while(nextSubstring());
		}
		if(chars.size()==4 && chars.contains('A') && chars.contains('C') && chars.contains('G') && chars.contains('T')) {
			setupAlphabet("A-A,C-C,G-G,T-T","ACGT");
		} else {
			String lo = "";
			String lc = "";
			for(Character c : chars) {
				lo += c;
				if(lc.length()>0)
					lc += ",";
				lc += c + "-" + c;
			}
			setupAlphabet(lc, lo);

		}
	}
	


	
	public void computeMonoInformation(ArrayList<String> args) {
		if(args==null) {
			System.out.println("computeMonoInformation()          - Computes the information content of each binding mode using the mononucleotide coefficients.");
			return;
		}
		
		int nArgs = args.size();

		if(nArgs == 0) {
			computeMonoInformation();
		} else {
			error("The function computeMonoInformation() must have exactly one arguments.");
		}
		
	}
	
	public void computeMonoInformation() {
		//Imposes mismatch gauge
		if(!mismatchGauge)
			setMismatchGauge();
		
        JSONArray bmSett = oModel.getJSONObject("modelSettings").getJSONArray("bindingModes");
        JSONArray bmCoef = oModel.getJSONObject("coefficients" ).getJSONArray("bindingModes");
        double[] out     = new double[bmSett.length()];
        for(int iBM=0; iBM<bmSett.length(); iBM++) {
        	double[] monoBetas = ModelComponent.readFromJSON_d(bmCoef.getJSONObject(iBM).getJSONArray("mononucleotide"));
        	int size           = bmSett.getJSONObject(iBM).getInt("size");
        	out[iBM]           = BindingMode.computeMonoInformation(monoBetas, size);
        }
        System.out.println(Misc.formatVector_d(out));
	}
	

	public void addNScoring(ArrayList<String> args) {
		if(args==null) {
			System.out.println("addNScoring()                     - Given a binding model and an alphabet that doesn't contain 'N', adds 'N' to denote wild cards and extends the binding model so 'N' corresponds to the average affinity.");
			return;
		}
		int nArgs = args.size();
		
		if(nArgs==0 ) {
			addNScoring();
		} else {
			error("The function addNScoring() must have exactly zero arguments.");
		}
	
		
	}
	//Gets n-scoring.
	public void addNScoring() {
		//Only add N-scoring if it doesn't allready exist in alphabet
		if(letterOrder.indexOf("N")==-1) {
			//Saves the old alphabet.
			String oldAlphabet = letterOrder;
			char[] oldChar     = oldAlphabet.toCharArray();
			int oldNChar       = oldAlphabet.length();
			int oldNDiChar     = oldNChar * oldNChar;

			//1. Update oModel.modelSettings.alphabet, saves parameters for the new.
			setupAlphabet(letterComplement+",N-N", letterOrder+"N");
			String newAlphabet = letterOrder;
			oModel.getJSONObject("modelSettings").put("letterComplement", letterComplement);
			oModel.getJSONObject("modelSettings").put("letterOrder",      letterOrder);
			int newNChar       = newAlphabet.length();
			int newNDiChar     = newNChar*newNChar;
			HashMap<Integer,Integer> oldToNew = new HashMap<Integer,Integer>();
			for(int i=0; i<oldChar.length; i++) 
				oldToNew.put(i, newAlphabet.indexOf(oldChar[i]));
			int indexNNew      = letterOrder.indexOf("N");

			//1.  Loops over binding modes
			for(int iBM=0; iBM<nModes; iBM++) {
				//1.0 Read all coefficients.
				//Reads coefficients for the current binding mode
				JSONObject oBMSetting           = oModel.getJSONObject("modelSettings").getJSONArray("bindingModes").getJSONObject(iBM);
				JSONObject oBMCoeff             = oModel.getJSONObject("coefficients").getJSONArray("bindingModes").getJSONObject(iBM);
				int size                        = oBMSetting.getInt("size");

				if(size>0) {
					//Mononucleotide
					double[] monoBetas          = ModelComponent.readFromJSON_d(oBMCoeff.getJSONArray("mononucleotide"));
					//Dinucleotide
					ArrayList<double[]> diBetas;
					if(oBMCoeff.has("dinucleotide"))
						diBetas                 = ModelComponent.readFromJSON_Ad(oBMCoeff.getJSONArray("dinucleotide"));
					else
						diBetas                 = new ArrayList<double[]>();

						//1.1  Identify the highest-affinity sequence for all binding modes.
						int[] topSeq                = null;
						if(dInt.get(iBM) > 1)
							topSeq                  = sampleTopSequence(size, monoBetas, diBetas);
						else
							topSeq                  = sampleTopSequence(size, monoBetas, diBetas);

						//1.2 Compute the energy of the highest affinity sequence.
						double topScore = scoreSeq(topSeq, monoBetas, diBetas, size);
						int[] tempSeq               = new int[size];
						for(int i=0; i<size; i++)
							tempSeq[i]              = topSeq[i];

						//1.3 Compute the energy of all single-base mismatches relative the highest-affinity sequence. Save as new PWM.
						//1.3.1 Computes the impact of consensus mismatches
						double[] monoBetaMismatch   = new double[size * newNChar];
						for(int x=0; x<size; x++) {
							double expMismatchSum   = 0;
							for(int iNucInOld=0; iNucInOld<oldNChar; iNucInOld++) {
								int old                          = tempSeq[x];
								tempSeq[x]                       = iNucInOld;
								double deltaMismatch             = scoreSeq(tempSeq, monoBetas, diBetas, size) - topScore;
								expMismatchSum                  += Math.exp(deltaMismatch);
								int iInNewMonoBeta               = x*newNChar + oldToNew.get(iNucInOld);
								monoBetaMismatch[iInNewMonoBeta] = deltaMismatch;
								tempSeq[x]                       = old;
							}
							monoBetaMismatch[x*newNChar + indexNNew] = Math.log(expMismatchSum/oldNChar);
						}
						oBMCoeff.put("mononucleotide", ModelComponent.JSONArrayConvert_d(monoBetaMismatch));

						//1.3.2 Shifts mononucleotide matrix so that it gives the correct top score.
						JSONArray aMono = oBMCoeff.getJSONArray("mononucleotide");
						for(int iMono = 0; iMono<aMono.length(); iMono++)
							aMono.put(iMono, aMono.getDouble(iMono) + topScore/size);

						//1.4 Compute the energy of all double-base mismatches relative the highest-affinity sequence. Save as new interactions.
						ArrayList<double[]> diBetaMismatch         = new ArrayList<double[]>();
						for(int d=1; d<=diBetas.size(); d++) {
							double[] newDiBeta                     = new double[(size-d) * newNDiChar];
							double expBetaSum                      = 0;
							double[] expBetaSum_P1                 = new double[oldNChar];
							double[] expBetaSum_P2                 = new double[oldNChar];

							for(int x=0; x<size-d; x++) {
								//Resets the running sums.
								expBetaSum = 0;
								for(int i=0; i<oldNChar; i++) {
									expBetaSum_P1[i] = 0;
									expBetaSum_P2[i] = 0;
								}
								for(int iNuc1InOld=0; iNuc1InOld<oldNChar; iNuc1InOld++) {
									for(int iNuc2InOld=0; iNuc2InOld<oldNChar; iNuc2InOld++) {
										//Updates the sequence
										int old1                   = tempSeq[x];
										int old2                   = tempSeq[x+d];
										tempSeq[x]                 = iNuc1InOld;
										tempSeq[x+d]               = iNuc2InOld;

										//Computes the dinucleotide interaction.
										int iInNewDiBeta           = x*newNDiChar   + oldToNew.get(iNuc1InOld)*newNChar + oldToNew.get(iNuc2InOld);
										int i1InNewMonoBeta        = x*newNChar     + oldToNew.get(iNuc1InOld);
										int i2InNewMonoBeta        = (x+d)*newNChar + oldToNew.get(iNuc2InOld);

										//Saves the mismatch beta
										double deltaMismatch       = scoreSeq(tempSeq, monoBetas, diBetas, size) - topScore;
										newDiBeta[iInNewDiBeta]    = deltaMismatch 
												- monoBetaMismatch[i1InNewMonoBeta] 
														- monoBetaMismatch[i2InNewMonoBeta];

										double expDeltaMismatch    = Math.exp(deltaMismatch);
										expBetaSum                += expDeltaMismatch;
										expBetaSum_P1[iNuc1InOld] += expDeltaMismatch;
										expBetaSum_P2[iNuc2InOld] += expDeltaMismatch;

										//Reverts to o the old sequence
										tempSeq[x]                 = old1;
										tempSeq[x+d]               = old2;

									}
								}


								//Computes the value of beta_N = Log[Mean[Exp[beta]]
								int N1InNewMonoBeta = x*newNChar     + indexNNew;
								int N2InNewMonoBeta = (x+d)*newNChar + indexNNew;

								//beta[N,N]
								newDiBeta[    x*newNDiChar + indexNNew           * newNChar + indexNNew          ] = Math.log(expBetaSum           / oldNDiChar)
										- monoBetaMismatch[N1InNewMonoBeta] 
												- monoBetaMismatch[N2InNewMonoBeta];
								//beta[i,N]
								for(int iNuc1InOld=0; iNuc1InOld<oldNChar; iNuc1InOld++) {
									int i1InNewMonoBeta        = x*newNChar     + oldToNew.get(iNuc1InOld);
									newDiBeta[x*newNDiChar + oldToNew.get(iNuc1InOld) * newNChar + indexNNew          ] = Math.log(expBetaSum_P1[iNuc1InOld] / oldNChar  )
											- monoBetaMismatch[i1InNewMonoBeta] 
													- monoBetaMismatch[N2InNewMonoBeta];
								}
								//beat[N,i]
								for(int iNuc2InOld=0; iNuc2InOld<oldNChar; iNuc2InOld++) {
									int i2InNewMonoBeta        = (x+d)*newNChar + oldToNew.get(iNuc2InOld);
									newDiBeta[x*newNDiChar + indexNNew           * newNChar + oldToNew.get(iNuc2InOld)] = Math.log(expBetaSum_P2[iNuc2InOld] / oldNChar  )
											- monoBetaMismatch[N1InNewMonoBeta] 
													- monoBetaMismatch[i2InNewMonoBeta];
								}

							}
							diBetaMismatch.add(newDiBeta);
						}

						oBMCoeff.put("dinucleotide",   ModelComponent.JSONArrayConvert_Ad(diBetaMismatch));

				}
			}


			




			/////////////////////////////

			//2. Updates the binding modes.
			/*		JSONObject oCoef = oModel.getJSONObject("coefficients");
		JSONArray aBM    = oCoef.getJSONArray("bindingModes");
		for(int iBM=0; iBM<aBM.length(); iBM++) {
			JSONObject oBM          = aBM.getJSONObject(iBM);

			//Saves the mononucleotide betas
			JSONArray aOldMonoBetas = oBM.getJSONArray("mononucleotide");
			int size                = aOldMonoBetas.length()/oldNChar;

			JSONArray aNewMonoBetas = ModelComponent.JSONArrayConvert_d(new double[size*newNChar]);
			for(int x=0; x<size; x++) {

				//Copies old betas, sums exp[beta]
				double expBetaSum = 0;
				for(int iOld=0; iOld<oldNChar; iOld++) {
					double currentBeta  = aOldMonoBetas.getDouble(x*oldNChar + iOld);
					expBetaSum         += Math.exp(currentBeta);
					aNewMonoBetas.put(x*newNChar + oldToNew.get(iOld), currentBeta);
				}


				aNewMonoBetas.put(x*newNChar + indexNNew, Math.log(expBetaSum / oldNChar));
			}
			oBM.put("mononucleotide", aNewMonoBetas);

			//Saves the dinucleotide betas
			JSONArray aaOldDiBetas = oBM.getJSONArray("dinucleotide");
			JSONArray aaNewDiBetas = new JSONArray();
			for(int iD=0; iD<aaOldDiBetas.length(); iD++) {
				JSONArray aNewDiBetas   = ModelComponent.JSONArrayConvert_d(new double[(size-iD-1)*newNDiChar]);
				aaNewDiBetas.put(aNewDiBetas);
				for(int x=0; x<size-iD-1; x++) {
					double expBetaSum = 0;
					double[] expBetaSum_P1 = new double[oldNChar];
					double[] expBetaSum_P2 = new double[oldNChar];

					for(int iOld1=0; iOld1<oldNChar; iOld1++) {
						for(int iOld2=0; iOld2<oldNChar; iOld2++) {
							double currentBeta    = aaOldDiBetas.getJSONArray(iD).getDouble(x*oldNDiChar + iOld1*oldNChar + iOld2);
							double expCurrentBeta = Math.exp(currentBeta);
							expBetaSum           += expCurrentBeta;
							expBetaSum_P1[iOld1] += expCurrentBeta;
							expBetaSum_P2[iOld2] += expCurrentBeta;
							aNewDiBetas.put(x*newNDiChar + oldToNew.get(iOld1)*newNChar + oldToNew.get(iOld2), currentBeta);
						}
					}
					//Computes the value of beta_N = Log[Mean[Exp[beta]]
					//beta[N,N]
					aNewDiBetas.put(x*newNDiChar     + indexNNew           * newNChar + indexNNew,           Math.log(expBetaSum           / oldNDiChar));
					//beta[i,N]
					for(int iOld1=0; iOld1<oldNChar; iOld1++)
						aNewDiBetas.put(x*newNDiChar + oldToNew.get(iOld1) * newNChar + indexNNew,           Math.log(expBetaSum_P1[iOld1] / oldNChar));
					//beat[N,i]
					for(int iOld2=0; iOld2<oldNChar; iOld2++)
						aNewDiBetas.put(x*newNDiChar + indexNNew           * newNChar + oldToNew.get(iOld2), Math.log(expBetaSum_P2[iOld2] / oldNChar));
				}
			}
			oBM.put("dinucleotide", aaNewDiBetas);
		}*/

			//3. Loads the model
			loadOModel();
		}
	}
	
	public void setAlphabet(ArrayList<String> args) {
		if(args==null) {
			System.out.println("setAlphabet(letterOrder,compl)    - Sets the alphabet.");
			return;
		}
		int nArgs = args.size();

		if(nArgs==1) {
			String letterOrder = args.get(0);
			String complement = "";
			for(int i=0; i<letterOrder.length(); i++) {
				if(i>0)
					complement = complement+",";
				complement = complement + letterOrder.charAt(i) + "-" + letterOrder.charAt(i);
			}
			setAlphabet(complement, letterOrder);
		} else if(nArgs==2) {
			String letterOrder = args.get(0);
			String complement  = args.get(1);
			setAlphabet(complement, letterOrder);
		} else {
			error("The function setAlphabet(letterOrder,compl) must have exactly two arguments.");
		}
	
		
	}
	//Gets n-scoring.
	public void setAlphabet(String complement, String letterOrder) {
		setupAlphabet(complement, letterOrder);
	}
	
	//Loads model from the MotifCentral database
	public void loadMotifCentralModel(ArrayList<String> args) {
		if(args==null) {
			System.out.println("loadMotifCentralModel(fit_id)     - Loads a model with from motifcentral.org");
			return;
		}
		int nArgs = args.size();

		if(nArgs==1) {
			int fit_id = Integer.parseInt(args.get(0));
			loadMotifCentralModel(fit_id);
		} else {
			error("The function motifCentralModel(fit_id) must have exactly one argument.");
		}
		
	}
	
	public void loadMotifCentralModel(int fitID) {
		
		String sURL ="https://prod-gateway.motifcentral.org/cellx/api/web/utility/fit/"+String.format("%d", fitID);
		System.out.println(sURL);

	    // Connect to the URL using java's native library
		String scoringJSON = null;
	    URL url            = null;
		try {
			url = new URL(sURL);
		} catch (MalformedURLException e) {
			e.printStackTrace();
			System.exit(1);
		}
		
	    URLConnection request;
		try {
			request     = url.openConnection();
			request.connect();
			scoringJSON = (new BufferedReader(new InputStreamReader((InputStream) request.getContent()))).lines().collect(Collectors.joining("\n"));

		} catch (IOException e) {
			e.printStackTrace();
			throw new java.lang.RuntimeException("ERROR: Error downloading model from motifcentral.org for rfit_id="+fitID);
		}
		
		//Checks so scoringJSON is not empty
    	if(scoringJSON==null || scoringJSON=="")
    		throw new java.lang.RuntimeException("ERROR: No scoring model exists for fit_id="+fitID);

    	//Converts the scoring_json string to a model object
		oModel = new JSONObject(scoringJSON); 
    	
		//Validates.
		JSONModel.validateSchemaFile_O(generalSchemaFile, oModel);
		
		//Adds model fitting constraints if they do not already exist
		JSONModel.addEmptyModelFittingConstraints(generalSchemaFile, oModel);

		//Loads the model.
		loadOModel();
		
	}
}
