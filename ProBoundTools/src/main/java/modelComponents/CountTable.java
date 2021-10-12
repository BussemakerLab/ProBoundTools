package modelComponents;

import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.concurrent.*;

//import org.apache.commons.math3.random.MersenneTwister;
import org.json.*;

import base.Array;
import base.MersenneTwisterFast;
import modelComponents.MultiRoundData;
import proBoundTools.Misc;
import sequenceTools.*;

public class CountTable extends ModelComponent  {
	
	//GOAL: Holds the data, computes the likelihood and its gradient.  
	
	//Variables defining dataset.
	private LongSequence.SequenceClass sc;
	public String letterComplement;
	public String letterOrder;
	int nMono, nDi;
	public MultiRoundData fullTable;
	ArrayList<LongSequence> longProbes;
	public int l, nColumns;
	public int nReads;
	String countTableFile, inputFileType, rightFlank, leftFlank;
	
	MersenneTwisterFast randGenerator;
	//Variables defining the batch table.
	public ArrayList<Integer> batchProbeIndices;
	public ArrayList<int[]> batchProbeCounts;
	public int[] batchCountPerRound;
	public boolean batchFixedLibrarySize;
	ArrayList<int[]> readIndices = null;

	//Round that should be included in the likelihood.
	public int[] modeledColumns;
	
	//Threading information
	private int nThreads;
	private int[][] threadRange;
	private ExecutorService pool;

	//Variables defining the enrichment model. 
	public double[] eta, h;

	//Variables for computing the function and its gradients.
	public double weight;
	double lambda;
	public double functionValue; //Function value
	public JSONObject gradient;  //Combined gradient
	public double functionValueSquared; //Sum over the squared value
	public boolean computeVariance; 
	
	private ArrayList<String> trIn, trOut;
	
	//Object defining the enrichment between the columns
	EnrichmentModel enr;
	
	
	///////////////////////
	// Setting up object //
	///////////////////////
	
	
	//Default constructor (Constructor that doesn't read the data).
	public CountTable(JSONObject config, int iExpIn, LongSequence.SequenceClass scIn, String letterComplementIn, String letterOrderIn, boolean loadData) {
		super("countTable");

		iComp          = iExpIn;
		componentName  = "Count table "+iComp;
		if(verbose)
			System.out.println(">> Creating "+componentName+".");

		sc               = scIn;
		letterComplement = letterComplementIn;
		letterOrder      = letterOrderIn;
		nMono            = letterOrder.length();
		nDi              = nMono * nMono;
		
		maxFreezeLevel = 0;
		readFromJSON_settings(   config);
		readFromJSON_constraints(config);

		seed_component(config);

		//Sets up threading.
		if(config.has("optimizerSetting") && config.getJSONObject("optimizerSetting").has("nThreads"))
			setNThreads(config.getJSONObject("optimizerSetting").getInt("nThreads"));
		else
			setNThreads(1);

		randGenerator = new MersenneTwisterFast();
		batchFixedLibrarySize = false;
		setFreezingLevel(0);
				
		computeVariance     = false;
		
		if(loadData) {
			loadDataTable(Misc.readDataFile(config, sc, iComp));
		}
	
	}
		
	public void loadDataTable(MultiRoundData dataIn) {
		
		//Reads data, builds a naive batch.
		fullTable          = dataIn;
		if(fullTable.longProbes==null)
			fullTable.encodeLongToLongSequence(sc);

		//Transliterates characters
		if(trIn.size()>0) 
			MultiRoundData.transliterate(fullTable, sc, trIn, trOut);

		longProbes         = fullTable.longProbes;
		nReads             = 0;
		for(int r: modeledColumns)
			nReads += (int) fullTable.countPerRound[r];
		
		nextBatchAll();
		
	}
	
	//Sets up threading.
	public void setNThreads(int nThreads) {
		this.nThreads = nThreads;
		pool          = Executors.newFixedThreadPool(nThreads);
		
	}
	
	
	public void setEnrichmentModel(EnrichmentModel enrIn) {
		enr           = enrIn;
	}
	
	//////////////////////////////////////////////
	// Function for implementing ModelComponent //
	//////////////////////////////////////////////
		
	public void allocateParameters() {
		
		h     = new double[nColumns];
		eta   = exp_d(h);
		updateAlphas();
		return;
		
	}
	
	@Override
	public void setComponentFiting(boolean fit) {
		
		if(fit) {
			if(includeComponent)
				fitComponent = true;
			else
				fitComponent = false;
		} else {
			fitComponent = false;
		}
	}
	
	public void updateAlphas() {
		
		if(eta==null)
			eta =new double[nColumns];
		else
			for(int r=0; r<nColumns; r++)
				eta[r] = 0;
		
		for(int r: modeledColumns) 
			eta[r] = Math.exp(h[r]);
			
		return;
	}

	@Override
	public void seed_component(JSONObject config) {

		String coefficientKey = "modelSeeding";
		h=null;
		if(config.has(coefficientKey)) {
			JSONObject oSeed = config.getJSONObject(coefficientKey); 
			if( oSeed.has(componentKey) ) {
				JSONObject oEnr = oSeed.getJSONArray(componentKey).getJSONObject(iComp);
				if(oEnr.has("h"))
					h     = readFromJSON_d(oEnr.getJSONArray("h"));
			}
		}
		
		if(h==null)
			h = new double[nColumns];

		updateAlphas();


	}
	
	public double[] computeExpectedPartitionFunction() {
		
		if(!includeComponent)
			return null;
		
		double[] alphaSum = new double[nColumns];

		//Sums over binding modes.
		for(int iBM=0; iBM<enr.bindingModes.size(); iBM++) {
			BindingMode oBM = enr.bindingModes.get(iBM);
			
			//Only sums over included binding modes
			if(!oBM.includeComponent)
				continue;
			
			//Expectation value of the scoring matrix.
			double alphaSeq = oBM.computeExpectedAlpha();
			
			//Computes the sum of the position bias/interactions.
			double alphaPBSum = (oBM.usePositionBias && oBM.k>0) ? tr_Ad(oBM.positionBiasAlphas.get(iComp)) : oBM.maxFrames.get(iComp);
			
			//Adds contribution to the alpha sum
			for(int iCol=0; iCol<nColumns; iCol++)
				alphaSum[iCol] += alphaSeq * alphaPBSum * oBM.activityAlphas.get(iComp)[iCol] * enr.concentration;
		}
		
		//Sums over interactions
		for(int iInt=0; iInt<enr.interactions.size(); iInt++) {
			BindingModeInteraction oInt = enr.interactions.get(iInt);
			
			//Only sums over included interactions
			if(!oInt.includeComponent)
				continue;
			
			//Expectation value of the scoring matrix.
			double alphaSeq = oInt.b0.computeExpectedAlpha() * oInt.b0.computeExpectedAlpha();

			//Computes the sum of the position bias/interactions.
			double alphaIntSum = tr_AAdd(oInt.interactionAlphas.get(iComp));

			//Adds contribution to the alpha sum
			for(int iCol=0; iCol<nColumns; iCol++)
				alphaSum[iCol] += alphaSeq * alphaIntSum * oInt.activityAlphas.get(iComp)[iCol] * enr.concentration;


		}
		
		return alphaSum;
	}
	
	@Override
	int packModel_component(JSONObject packing, int iFirst) {
		
		String coefficientKey = "packing";
		addEmptyJSON_component_O(packing, coefficientKey, componentKey, iComp);

		int iCurr = iFirst;
		
		if(fitComponent) {

			int[] hRange = ModelComponent.constant_i(h, -1);
			for(int i=0; i<modeledColumns.length; i++) 
				hRange[modeledColumns[i]] = iCurr+i+1;
			iCurr += modeledColumns.length;

			saveToJSON_h_i(packing, coefficientKey, hRange);
		}
		return iCurr;
	}

	@Override
	public void addZeroJSON_component(JSONObject in, String coefficientKey) {
		
		addEmptyJSON_component_O(in, coefficientKey, componentKey, iComp);
		saveToJSON_h_d(in,     coefficientKey, zero_d(h));
		
		return;
	}

	@Override
	void saveToJSON_settings(JSONObject out) {

		String coefficientKey = "modelSettings";
		addEmptyJSON_component_O(out,       coefficientKey, componentKey, iComp);
		JSONObject oCT   =                  out.getJSONObject(coefficientKey).getJSONArray(componentKey).getJSONObject(iComp);
		
		oCT.put("countTableFile",       countTableFile);
		oCT.put("inputFileType",        inputFileType);
		oCT.put("nColumns",             nColumns);
		oCT.put("variableRegionLength", l);
		oCT.put("rightFlank",           rightFlank);
		oCT.put("leftFlank",            leftFlank);
		oCT.put("modeledColumns",       ModelComponent.JSONArrayConvert_i(modeledColumns));
		
		JSONObject oTR = new JSONObject();
		oCT.put("transliterate", oTR);
		oTR.put("in",            trIn);
		oTR.put("out",           trOut);

	}
	
	@Override
	void readFromJSON_settings(JSONObject in) {

		String coefficientKey = "modelSettings";
		JSONObject oCT        = in.getJSONObject(coefficientKey).getJSONArray(componentKey).getJSONObject(iComp);
		countTableFile        = oCT.getString("countTableFile");
		inputFileType         = oCT.getString("inputFileType");
		nColumns              = oCT.getInt("nColumns");
		l                     = oCT.getInt("variableRegionLength");
		rightFlank            = oCT.getString("rightFlank");
		leftFlank             = oCT.getString("leftFlank");
		
		//By default, all columns are modeled.
		modeledColumns = new int[nColumns];
		for(int c=0; c<nColumns; c++)
			modeledColumns[c] = c;
		//Reads the modeled columns from the file if they exist. 
		if(oCT.has("modeledColumns")) {
			JSONArray mc = oCT.getJSONArray("modeledColumns");
			if( !(mc.length()==1&&mc.getInt(0)==-1) )
				modeledColumns = ModelComponent.readFromJSON_i(mc);
		}
		
		trIn  = new ArrayList<String>();
		trOut = new ArrayList<String>();
		if(oCT.has("transliterate")) {
			JSONObject oTR = oCT.getJSONObject("transliterate");
			if(oTR.has("in") && oTR.has("out")) {
				//Gets the in and out JSON arrays
				JSONArray oIn  = (oTR.get("in")  instanceof JSONArray) ? oTR.getJSONArray("in")  : (new JSONArray()).put(oTR.getString("in"));
				JSONArray oOut = (oTR.get("out") instanceof JSONArray) ? oTR.getJSONArray("out") : (new JSONArray()).put(oTR.getString("out"));

				//Checks so 'in' and 'out' have the same number of elements.
				if(oIn.length() != oOut.length())
					throw new IllegalArgumentException("For modelSettings.countTable.transliterate, the length of 'in' does not match the length of 'out.");
				//Loops over elements
				for(int iStr=0; iStr<oIn.length(); iStr++) {
					//Gets the strings
					String sIn  = oIn.getString(iStr);
					String sOut = oOut.getString(iStr);
					//Checks so the in and out strings have the same length 
					if(sIn.length() != sOut.length())
						throw new IllegalArgumentException("For modelSettings.countTable.transliterate, the length of 'in' string "+sIn+" does not match the length of 'out string "+sOut+".");
					//Saves the strings.
					trIn.add(sIn);
					trOut.add(sOut);
				}
			}
		}
	}
	
	@Override
	void saveToJSON_constraints(JSONObject out) {
		
		String coefficientKey = "modelFittingConstraints";
		
		//Creates coefficients if required. 
		if(!out.has(coefficientKey))
			out.append(coefficientKey, new JSONObject());
		JSONObject oCoeff = out.getJSONObject(coefficientKey);
		
		//Creates binding modes if required. 
		if(!oCoeff.has(componentKey))
			oCoeff.append(componentKey, new JSONObject());

	}

	@Override
	void readFromJSON_constraints(JSONObject in) {
		
	}
	
	@Override
	void saveToJSON_parameters(JSONObject in, String coefficientKey) {
		
		addEmptyJSON_component_O(in, coefficientKey, componentKey, iComp);
		saveToJSON_h_d(in,     coefficientKey, h);

	}

	@Override
	void readFromJSON_parameters(JSONObject in, String coefficientKey) {
		JSONObject oEnr = in.getJSONObject(coefficientKey).getJSONArray(componentKey).getJSONObject(iComp);
		h               = readFromJSON_d(oEnr.getJSONArray("h"));
		updateAlphas();
		
	}
	
	//Methods for saving double objects
	void saveToJSON_h_d(JSONObject out, String coefficientKey, double[] h) {
		
		JSONObject oBM = out.getJSONObject(coefficientKey).getJSONArray(componentKey).getJSONObject(iComp);
		oBM.put("h", JSONArrayConvert_d(h));

		return;
	}
	

	//Methods for saving long objects
	void saveToJSON_h_i(JSONObject out, String coefficientKey, int[] h) {
		
		JSONObject oBM = out.getJSONObject(coefficientKey).getJSONArray(componentKey).getJSONObject(iComp);
		oBM.put("h", JSONArrayConvert_i(h));

		return;
	}
	
	////////////////////////////////////////////////////////
	// Function for computing function value and gradient //
	////////////////////////////////////////////////////////
	
	//Creates a new batch table given a list of read indices for each column.
	public void nextBatch(ArrayList<int[]> randomReads, double expectedReadCount) {
		
		//FIRST TIME ONLY: Creates a list of reads: readIndices[iColumn][iRead] = <index of probe>
		if(readIndices==null) {
			readIndices = new ArrayList<int[]>();
			for(int iCol=0; iCol<nColumns; iCol++) {
				int[] newIndices = new int[fullTable.countPerRound[iCol]];
				readIndices.add(newIndices);
				int iRead = 0;
				for(int iProbe=0; iProbe<longProbes.size(); iProbe++) {
					int readCount = fullTable.countTable.get(iProbe)[iCol];
					if(readCount>0) {
						for(int a=0; a<readCount; a++) {
							newIndices[iRead] = iProbe;
							iRead++;
						}
					}
				}
			}
		}
		
		//Variables for the batch table.
		batchProbeIndices  = new ArrayList<Integer>();
		batchProbeCounts   = new ArrayList<int[]>();
		batchCountPerRound = new int[nColumns];
		
		//Map from index of probe in the full table to index of the probe in the batch table.
		HashMap<Integer,Integer> probeToBatchRow  = new HashMap<Integer,Integer>();

		for(int iCol : modeledColumns) {
			//Adds the randomly selected reads to the batch count table.
			
			for(int ir: randomReads.get(iCol)) {
				int iProbe = readIndices.get(iCol)[ir]; //Index of the probe corresponding to the current read.
				
				//If the probe corresponding to the current read does not exist in the batch table, create a new entry.  
				if(!probeToBatchRow.containsKey(iProbe)) {
					batchProbeIndices.add(iProbe);
					batchProbeCounts.add(new int[nColumns]);
					probeToBatchRow.put(iProbe, batchProbeCounts.size()-1);
				}
				
				//Increments the probe count
				batchProbeCounts.get(probeToBatchRow.get(iProbe))[iCol] += 1;
				batchCountPerRound[iCol] ++;
			}
		}
		
		weight        = 1. / expectedReadCount;
		threadSchedule(nThreads);	
	}
	
	//Creates a new batch given a list of read indices (not grouped by column)
	public void nextBatch(int[] randomUnsortedReads, double expectedReadCount) {
		
		int[] cpr = fullTable.countPerRound;
		
		//A temporary list storing the read indices in each column.  
		ArrayList<ArrayList<Integer>> tempColReads = new ArrayList<ArrayList<Integer>>();
		for(int iCol=0; iCol<nColumns; iCol++)
			tempColReads.add(new ArrayList<Integer>());
		
		//Step 2: Sorts the random reads by columns
		for(int iRead: randomUnsortedReads) {
			int iInCol = iRead;
			//for(int iCol=0; iCol<nColumns; iCol++) {
			for(int iCol : modeledColumns) {
				if(iInCol < cpr[iCol]) {
					tempColReads.get(iCol).add(iInCol);
					break;
				} else {
					iInCol -= cpr[iCol];
				}
			}
		}
		
		//Step 3: Converts the array list to a int[];
		//Creates a random list of the read indices for each column  
		ArrayList<int[]> sortedReads = new ArrayList<int[]>();
		for(int iCol=0; iCol<nColumns; iCol++)  {
			ArrayList<Integer> aI = tempColReads.get(iCol);
			int[] tableCol = new int[aI.size()];
			for(int iI=0; iI<aI.size(); iI++)
				tableCol[iI] = aI.get(iI);
			sortedReads.add(tableCol);
		}

		nextBatch(sortedReads, expectedReadCount);
	}
	
	//Creates a new batch using all reads.
	public void nextBatchAll() {
		//Uses all reads in the table
		batchProbeCounts   = fullTable.countTable;
		batchCountPerRound = fullTable.countPerRound;
		batchProbeIndices  = new ArrayList<Integer>();
		for(int i=0; i<batchProbeCounts.size(); i++)
			batchProbeIndices.add(i);
		
		//Schedules the threads.
		threadSchedule(nThreads);
		weight             = 1.0 / nReads;
	}
	
	//Creates a new batch given the total number of desired reads.
	public void nextBatch(int nDataIn) {
		
		if(fullTable.countTable!=null) {

			//Counts the number of reads
			nReads             = 0;
			for(int r: modeledColumns)
				nReads += (int) fullTable.countPerRound[r];
			
			//Determines the whole dataset should be used.
			int nData;
			if(nDataIn>=nReads) //Makes sure that the batch is not lager than the number of probe sequences
				nData = -1;
			else if(nDataIn==0) //Makes sure we use at least one probe.
				nData = 1;
			else
				nData = nDataIn;
			
			if(nData == -1) {
				nextBatchAll();
				return;

			} else {
				
				//int nReadsTot = (int) Array.sum(fullTable.countPerRound);
				if(batchFixedLibrarySize) {
					//CASE 1: Number of reads per column proportional to sequencing depth 
					//////////////////////////////////////////////////////////////////////
					//Counts the total number of reads
					int[] tempCountPerRound  = new int[nColumns];
					//Creates a random list of the read indices for each column  
					ArrayList<int[]> randomReads = new ArrayList<int[]>();
					int nReadsCum = 0; //Cumulative number of reads in the columns
					for(int iCol=0; iCol<nColumns; iCol++) {
						//Only pull reads from the modeled columns.
						if(Array.memberQ(modeledColumns, iCol)) {

							//Step 1: Determines how many reads we should have in each column. 
							int lastNReadsCum        = nReadsCum;
							nReadsCum               += fullTable.countPerRound[iCol];
							
							tempCountPerRound[iCol]  =  ((int) Math.floor(1.0*nData*nReadsCum    /nReads)) 
										              - ((int) Math.floor(1.0*nData*lastNReadsCum/nReads));
							//Step 2: Creates a list of read indices.
							randomReads.add(Misc.randomSample(fullTable.countPerRound[iCol], tempCountPerRound[iCol], randGenerator));
						} else {
							//Add no reads if the columns isn't modeled.
							randomReads.add(new int[0]);
						}
					}
					
					//Creates count table.
					nextBatch(randomReads, nData);
					return;

				} else {
					//CASE 2: Reads drawn randomly regardless of round 
					//////////////////////////////////////////////////
					nextBatch(Misc.randomSample(nReads, nData, randGenerator), nData);
					return;
				}
			}
			
			
		} else {
			nReads        = 0;
			weight        = 0;
			
			//Schedules the threads.
			threadSchedule(nThreads);	
			return;
		}

		
	}
	
	public void writeBatchTable(String outTable) {
		
		//Loops over probes.
		PrintStream originalStream = System.out;
		PrintStream fileStream     = null;
		try {
			fileStream = new PrintStream(new FileOutputStream(outTable, false));
			System.setOut(fileStream);
		} catch (FileNotFoundException e) {
			System.out.println("Cannot create batch table file: "+outTable);
			e.printStackTrace();
			System.exit(1);
		}
		
		for(int iProbe=0; iProbe<batchProbeCounts.size(); iProbe++)
			System.out.println(longProbes.get(batchProbeIndices.get(iProbe)).toString()+"\t"+Misc.formatVector_i(batchProbeCounts.get(iProbe), "\t", "", ""));
			
		//Returns to original stream.
		System.setOut(originalStream);
		fileStream.close();
				
	}

	// Sets up threading
	////////////////////
	public void threadPoolShutdown() {
		pool.shutdown();
		while (!pool.isShutdown()) {}
		return;
	}
	
	private void threadSchedule(int nThreads) {
		
		//Threading setup for cycling the dats
		threadRange		                = new int[nThreads][2];
		int nProbes                     = batchProbeCounts.size();
		
		for(int iThread=0; iThread<nThreads; iThread++ ) {
			threadRange[iThread][0]     = ((int) Math.floor((1.0*(iThread  )/nThreads)*nProbes )); 
			threadRange[iThread][1]     = ((int) Math.floor((1.0*(iThread+1)/nThreads)*nProbes )); 
			
			// Checks if the thread range contains zero reads: threadRange[iThread][0]=threadRange[iThread][1]
			// could imply that the the thread contains zero probes or that it contains all probes. 
			if( nProbes<nThreads && threadRange[iThread][0]==threadRange[iThread][1] ) {
				threadRange[iThread][0] = -1;
				threadRange[iThread][1] = -1;
			}
		}
	}		


	private void computePRI(double[] pRI, double[] kappaRI, double[] deltaKappaRI) {
		
		//Calculates Kappa
		if(enr.cumulativeEnrichment){
			double runningKappaProduct = 1.0;
			for(int r=0; r<nColumns; r++){
				runningKappaProduct *= deltaKappaRI[r];
				kappaRI[r] = runningKappaProduct;
			}
		} else {
			for(int r=0; r<nColumns; r++){
				kappaRI[r] = deltaKappaRI[r];
			}
		}

		//Updates pRI using only the columns specified in modeledColumns.	
		double zProbe = 0.0;
		for(int r=0; r<nColumns; r++)
			pRI[r] = 0;
		for(int r: modeledColumns) {
			pRI[r] = eta[r] * kappaRI[r];
			zProbe += pRI[r];
		}
		
		for(int r: modeledColumns) pRI[r] /= zProbe;

	}
	
	public void writePCTable(String predictedTable) /*throws Exception*/ {
		
		int nRounds = fullTable.countPerRound.length;

		double[] alphaSeq             = new double[enr.nModes];
		double[] alphaInt             = new double[enr.nInteractions];
		double[] alphaRI              = new double[nRounds];
		ArrayList<ArrayList<ArrayList<Double>>> longAlphaList = new ArrayList<ArrayList<ArrayList<Double>>>(); 
		for(int bm=0; bm<enr.nModes; bm++)
			longAlphaList.add(null);

		double[] deltaKappaRI         = new double[nRounds];
		double[] kappaRI              = new double[nRounds];
		double[] pRI                  = new double[nRounds];
		
		ArrayList<SlidingWindow> sw  = new ArrayList<SlidingWindow>();
		for(int bm=0; bm<enr.nModes; bm++)
			sw.add(enr.bindingModes.get(bm).getSlidingWindow(iComp, enr.modifications));

		//Loops over probes.
		PrintStream original	= System.out;
		try {
			System.setOut(new PrintStream(new FileOutputStream(predictedTable, false)));
		} catch (FileNotFoundException e) {
			System.out.println("Cannot create predicted-table file at this location: "+predictedTable);
			e.printStackTrace();
			System.exit(1);
		}

		int nDataPoints		= longProbes.size();
		
		//Loops over all probes (ignoring possible batching and multithreading.)
		for(int iProbe=0; iProbe<nDataPoints; iProbe++) {
			//Updates alpha values.
			enr.computeAlphas(sw, alphaSeq, alphaInt, alphaRI, longAlphaList, longProbes.get(iProbe));
			//Calculates deltaKappa:
			enr.updateDeltaKappaRI(deltaKappaRI, alphaRI, nRounds);
			//Computes pRI.
			computePRI(pRI, kappaRI, deltaKappaRI);
			//Counts the reads in the modeled columns.
			long nRow = 0;
			for(int r:modeledColumns)
				nRow += fullTable.countTable.get(iProbe)[r];
			
			//Prefroms reverse transliteration to get the original sequence.
			String probeSeq = longProbes.get(iProbe).toString();
			for(int iTr=0; iTr<trIn.size(); iTr++) 
				probeSeq = probeSeq.replaceAll(trOut.get(iTr), trIn.get(iTr));
			
			//Multiplies pRI by the number of reads.
			System.out.println(probeSeq + "\t" + Misc.formatVectorE_d(Array.scalarMultiply(pRI, nRow),"\t","","", 10));
		} 
		
		//Returns to original stream.
		System.setOut(original);
		
		return;
	}
	
	
	public void writeAlphaTable(String outAlphaFile) /*throws Exception*/ {
		
		int nRounds                   = fullTable.countPerRound.length;
		double[] alphaSeq             = new double[enr.nModes];
		double[] alphaInt             = new double[enr.nInteractions];
		double[] alphaRI              = new double[nRounds];
		
		ArrayList<ArrayList<ArrayList<Double>>> longAlphaList
		                              = new ArrayList<ArrayList<ArrayList<Double>>>();
		
		for(int bm=0; bm<enr.nModes; bm++)
			longAlphaList.add(null);

		ArrayList<SlidingWindow> sw   = new ArrayList<SlidingWindow>();
		for(int bm=0; bm<enr.nModes; bm++)
			sw.add(enr.bindingModes.get(bm).getSlidingWindow(iComp, enr.modifications));

		//Loops over probes.
		int nSequences		          = longProbes.size();
		PrintStream original          = System.out;
		
		try {
			System.setOut(new PrintStream(new FileOutputStream(outAlphaFile, false)));
		} catch (FileNotFoundException e) {
			System.out.println("Cannot create output table at this location: "+outAlphaFile);
			e.printStackTrace();
			System.exit(1);
		}
		
		//Loops over sequences, computes affinity sums (alphas), and writes to file
		for(int iProbe=0; iProbe<nSequences; iProbe++) {
			LongSequence currentProbeSeq = longProbes.get(iProbe);
			enr.computeAlphas(sw, alphaSeq, alphaInt, alphaRI, longAlphaList, currentProbeSeq);

			//Writes the table
			System.out.println(longProbes.get(iProbe).toString() + "\t" + Misc.formatVector_d(alphaRI, "\t", "", "", 8));

		}
		
		//Returns to original stream.
		System.setOut(original);
		
		return;
	}
	
	public void writeBindingModeAlphas(String outAlphaFile, String outFormat) {
		
		//Creates a list of relevant sliding windows and binding modes
		ArrayList<BindingMode> bmList     = new ArrayList<BindingMode>();
		ArrayList<SlidingWindow> swList   = new ArrayList<SlidingWindow>();
		for(BindingMode oBM: enr.bindingModes) {
			if(oBM.includeComponent && oBM.k>0) {
				bmList.add(oBM);
				swList.add(oBM.getSlidingWindow(iComp, enr.modifications));
			}
		}
		
		//Loops over probes.
		int nSequences		          = longProbes.size();
		PrintStream original          = System.out;
		
		//Opens output file
		if(!outAlphaFile.equals("-")) {
			try {
				System.setOut(new PrintStream(new FileOutputStream(outAlphaFile, false)));
			} catch (FileNotFoundException e) {
				System.out.println("Cannot create output table at this location: "+outAlphaFile);
				e.printStackTrace();
				System.exit(1);
			}
		}
		
		//Loops over sequences, computes affinity sums (alphas), and writes to file
		ArrayList<ArrayList<Double>> longAlphaList;
		for(int iProbe=0; iProbe<nSequences; iProbe++) {
			
			LongSequence currentProbeSeq = longProbes.get(iProbe);
			
			System.out.print(currentProbeSeq.toString());
			
			for(int iBM=0; iBM<bmList.size(); iBM++) {
				
				//Scores the sequence
				BindingMode oBM   = bmList.get(iBM);
				SlidingWindow oSW = swList.get(iBM);
				boolean ss        = oBM.singleStrand;
				
				if(!oBM.swIncludeDi) {
					longAlphaList = oSW.slidePN(currentProbeSeq,                       1 );
				} else {
					longAlphaList = oSW.slidePN(currentProbeSeq, Math.min(oBM.dInt+1, 2) );
				}
				String colString = "";
				int nDigits      = 5;
				if(outFormat.equals("max")) {
					double maxValue  = Array.max(longAlphaList.get(0));
					if(!ss)
						maxValue = Math.max(maxValue,  Array.max(longAlphaList.get(1)));
					colString        = String.format("%."+nDigits+"e", maxValue);

				} else if(outFormat.equals("sum")) {
					double sumValue  = Array.sum(longAlphaList.get(0));
					if(!ss)
							sumValue+= Array.sum(longAlphaList.get(1));
					colString        = String.format("%."+nDigits+"e", sumValue);
					
				} else if(outFormat.equals("mean")) {
					int n = longAlphaList.get(0).size();
					double meanValue =                 (Array.sum(longAlphaList.get(0))/n);
					if(!ss) 
						meanValue = 0.5 * (meanValue + (Array.sum(longAlphaList.get(1))/n));
					colString        = String.format("%."+nDigits+"e", meanValue);

				} else if(outFormat.equals("profile")) {
					Collections.reverse(longAlphaList.get(1));
					colString        = "" +   Misc.formatVectorE_d(                      longAlphaList.get(0),         ",", "", "", 5);
					if(!ss)
						colString   += "\t" + Misc.formatVectorE_d(                      longAlphaList.get(1),         ",", "", "", 5);
					else
						colString   += "\t" + Misc.formatVectorE_d(ModelComponent.zero_d(longAlphaList.get(1).size()), ",", "", "", 5);
					
					
					
				} else {
					System.err.println("ERROR: Invalid output format: '"+outFormat+"'.");
					System.exit(1);
				}
				System.out.print("\t"+colString);
			}
			System.out.println("");
		}
		
		//Returns to original stream.
		System.setOut(original);
		
		return;
		
	}
	
	//Uses the count table to compute a weight:
	// (maxStrandMean / minStrandMean) * maxValue
	public HashMap<Integer,Double> getEmpericalBindingModeActivities() {
		
		HashMap<Integer,Double> out = new HashMap<Integer,Double>();
				
		//Loops over binding modes (that are relevant for this count table)
		ArrayList<ArrayList<Double>> longAlphaList;

		//BindingMode oBM = enr.bindingModes.get(iBMInt);
		for(BindingMode oBM : enr.bindingModes) {

			//Checks so the binding mode is included and has size>0, skips if not.
			if(!oBM.includeComponent || oBM.k<=1)
				continue;

			//Gets the sliding window.
			int tempFL        = oBM.flankLength;
			oBM.flankLength   = 0;
			SlidingWindow oSW = oBM.getSlidingWindow(iComp, enr.modifications);
			oBM.flankLength   = tempFL;

			//Loops over probes
			int countSum = 0;
			double forwardSum=0, reverseSum=0, maxVal=0, activitySum = 0;
			
			//TODO: TEMP, REMOVE:
			double[] vSumFwrd = new double[nColumns], vSumRvrs = new double[nColumns];

			
			for(int iProbe=0; iProbe<longProbes.size(); iProbe++) {
				LongSequence currentProbeSeq = longProbes.get(iProbe);
				
				//Scores using sliding window (that doesn't go into the flanks).
				if(!oBM.swIncludeDi) {
					longAlphaList = oSW.slidePN(currentProbeSeq,                       1 );
				} else {
					longAlphaList = oSW.slidePN(currentProbeSeq, Math.min(oBM.dInt+1, 2) );
				}

				//Computes the max value and the strand sums
				double probeMeanForward = Array.mean(longAlphaList.get(0));
				double probeMeanReverse = Array.mean(longAlphaList.get(1));
				double weight           = 0;
				boolean bmSelected      = false;
				for(int iCol : modeledColumns) {
					//Gets the counts, computes a count-based probe weight 
					int cnt         = fullTable.countTable.get(iProbe)[iCol];
					weight         += cnt;
					
					//Counts the reads and sums the count-weighted activities, ignoring the first 
					//round if cumulativeEnrichment=true:  
					if( !(enr.cumulativeEnrichment&&iCol==0) ) {
						activitySum += cnt * oBM.activityAlphas.get(iComp)[iCol];
						countSum    += cnt;
						bmSelected   = true;
					}
					
					
				}
				
				//For computing the round-wise ratio
				for(int iCol=0; iCol<nColumns; iCol++) {
					vSumFwrd[iCol] += probeMeanForward * fullTable.countTable.get(iProbe)[iCol];
					vSumRvrs[iCol] += probeMeanReverse * fullTable.countTable.get(iProbe)[iCol];
				}
				
				//Computes the appropriate round-and-activity weighted sums
				forwardSum += weight * probeMeanForward;
				reverseSum += weight * probeMeanReverse;
				
				//Computes maximum value (only keeping probes with non-zero counts in a modeled column)
				if(bmSelected) 
						maxVal = Math.max(maxVal, Math.max(Array.max(longAlphaList.get(0)), Array.max(longAlphaList.get(1))));
				
			}
			
			//Ratios for verbose output


			
			//Computes the strand ratio using trend across rounds
			double e1=0, eRatio=0, eRound=0, eRatioRound=0, eRound2=0;
			for(int iCol : modeledColumns) {
				if(fullTable.countPerRound[iCol]>0) {
					double logRatio = Math.log(vSumRvrs[iCol]/vSumFwrd[iCol]);
					
					e1          += fullTable.countPerRound[iCol];
					eRatio      += fullTable.countPerRound[iCol]*logRatio;
					eRound      += fullTable.countPerRound[iCol]*iCol;
					eRatioRound += fullTable.countPerRound[iCol]*iCol*logRatio;
					eRound2     += fullTable.countPerRound[iCol]*iCol*iCol;
				}
			}
			double covRatioRound = (eRatioRound/e1) - (eRatio/e1)*(eRound/e1);
			double varRound      = (eRound2/e1)     - (eRound/e1)*(eRound/e1);
			double slope         = covRatioRound / varRound;
			double modelRatio    = Math.exp(-Math.abs(slope)*(Array.max(modeledColumns)-Array.min(modeledColumns)));

			//Computes the ratio strand ratio by taking read-weighted strand average
			double frwdRevRatio = Math.min(reverseSum/forwardSum, forwardSum/reverseSum);

			//Identifies the most extreme ratio
			double strandRatio = Math.min(modelRatio, frwdRevRatio);
			
			//Deals with None:
			if(Double.isNaN(strandRatio))
				strandRatio = 0;
			
			//Computes the mean position bias (ignoring positions overlapping flanks..
			double meanPB = 1;
			if(oBM.positionBiasAlphas!=null && iComp<oBM.positionBiasAlphas.size()) {
				int fl            = oBM.flankLength;
				double[] pbFwrd   = oBM.positionBiasAlphas.get(iComp).get(0);
				double[] bpRvrs   = oBM.positionBiasAlphas.get(iComp).get(1);
				double pbFwrdMean = Array.mean(Array.copyOfRange(pbFwrd, fl, pbFwrd.length-fl));
				double pbRvrsMean = Array.mean(Array.copyOfRange(bpRvrs, fl, pbFwrd.length-fl));
				meanPB = (pbFwrdMean + pbRvrsMean) / 2;
			}
	
			out.put(oBM.iComp, maxVal * meanPB * (activitySum / countSum) * strandRatio);

			
			if(verbose) {

				double[] vRatio = new double[nColumns];
				for(int iCol=0; iCol<nColumns; iCol++) {
					if(fullTable.countPerRound[iCol]>0) {
						vRatio[iCol] = vSumRvrs[iCol]/vSumFwrd[iCol];
					}
				}
				
				
				//System.out.println("iBM = "+oBM.iComp+", ratio = "+tempRatio);
				System.out.println("iBM = "+oBM.iComp+", ratio1 = "+Misc.formatVector_d(vRatio)+" => "+modelRatio+", ratio2 = "+frwdRevRatio+" => min="+strandRatio);
				System.out.println("         max="+maxVal+", PB="+(meanPB) + ", act="+(activitySum / countSum)+", ratio="+strandRatio);
			}
		}
				
		return out;
	}
	
	static public void printJSONObjectCoefficients(JSONObject model, String coefficientKey, int iExp) {
		
		JSONObject oTable = model.getJSONObject(coefficientKey).getJSONArray("countTable").getJSONObject(iExp);
		System.out.println("h:                 "+Misc.formatVector_d(ModelComponent.readFromJSON_d(oTable.getJSONArray("h"))));
		
	}
}