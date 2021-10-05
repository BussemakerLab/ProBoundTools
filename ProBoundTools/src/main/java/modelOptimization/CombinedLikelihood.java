package modelOptimization;

import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.PrintStream;
import java.sql.Timestamp;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.Set;

import org.json.JSONArray;
import org.json.JSONObject;

import base.MersenneTwisterFast;
import sequenceTools.*;
import modelComponents.MultiRoundData;
import proBoundTools.*;
import modelComponents.*;


public class CombinedLikelihood extends LossFunction{
	
		
	// - Our generalized optimizer class takes a lossFunction as an argument.
	
	// GOAL: This class holds
	// 1) Holds all components: binding modes, binding mode interactions, and experiments.
	// 2) Holds a list of all model components in a format that can be converted between the expanded and the compactified format.
	// 3) 
	
	JSONObject config;
	public ArrayList<BindingMode> bindingModes;             // Holds all binding modes.
	public ArrayList<BindingModeInteraction> interactions;  // Holds all binding mode interactions.
	public ArrayList<EnrichmentModel> enrichmentModels;     // Models of the individual experiments.
	public ArrayList<CountTable> tableModels;               // Models of the individual experiments.
	public ArrayList<ModelComponent> componentList;         // List of all components in "ComponentList" format.
	public ArrayList<ArrayList<ModelComponent>> fittingOrder; // Order in which the binding modes will be added and fitted.
	public JSONObject packing;                              // 
	private boolean valueUpdated=false;
	LongSequence.SequenceClass sc;                   // Object for encoding sequences.
	
	ArrayList<MultiRoundData> datasets = new ArrayList<MultiRoundData>();
	boolean experimentSpecificActivity;
	boolean addBindingModesSequentially;
	
	MersenneTwisterFast randGenerator;

	public double fitStartTime, logLikelihood, logLikelihoodPerRead;
	public long fitSteps, functionCalls, gradientCalls;
	public double dataLoops;
	public boolean verbose;
	public double lambdaL2            = 0;
	public double pseudocount         = 0;
	public double expBound            = 40;
	private double regularizationTerm = 0; //Value of the regularization term.
	public boolean fixedLibrarySize  = false;
	
	public String letterComplement;
	public String letterOrder;
	
	//Variables for saving the trajectory
	boolean printTrajectory = false;
	String trajectoryFile = null;
	

	
	public CombinedLikelihood(JSONObject configIn) {
		

		verbose = configIn.has("optimizerSetting") 
				&& configIn.getJSONObject("optimizerSetting").has("output")
				&& configIn.getJSONObject("optimizerSetting").getJSONObject("output").has("verbose") ?
						configIn.getJSONObject("optimizerSetting").getJSONObject("output").getBoolean("verbose") :
						false;
						
		if(verbose)
			System.out.println(">> Creating CombinedLikelihood object.");
		config           = configIn;
		componentList    = new ArrayList<ModelComponent>();
		
		//Reads the alphabet
		JSONObject oSett = config.getJSONObject("modelSettings");
		letterComplement = oSett.has("letterComplement") ? oSett.getString("letterComplement") : "C-G,A-T";
		letterOrder      = oSett.has("letterOrder")      ? oSett.getString("letterOrder")      : alphabetDefToOrderedLetters(letterComplement);
		if(verbose) {
			System.out.println("Alphabet");
			System.out.println("========");
			System.out.println("Letter Complement: "+letterComplement);
			System.out.println("Letter Order:      "+letterOrder);
			System.out.println("");
		}
		sc = new LongSequence.SequenceClass(letterComplement);
		
		if(configIn.has("optimizerSetting")) {
			JSONObject oOpt   = configIn.getJSONObject("optimizerSetting");
			lambdaL2          = oOpt.has("lambdaL2")         ? oOpt.getDouble("lambdaL2")          : 0;
			pseudocount       = oOpt.has("pseudocount")      ? oOpt.getDouble("pseudocount")       : 0;
			expBound          = oOpt.has("expBound")         ? oOpt.getDouble("expBound")          : 40;
			fixedLibrarySize  = oOpt.has("fixedLibrarySize") ? oOpt.getBoolean("fixedLibrarySize") : false;
			if(verbose) {
				System.out.println("Optimizer settings:");
				System.out.println("===================");
				System.out.println("lambdaL2         = "+lambdaL2);
				System.out.println("pseudocount      = "+pseudocount);
				System.out.println("expBound         = "+expBound);
				System.out.println("fixedLibrarySize = "+fixedLibrarySize);
				System.out.println("");
			}
		}

		randGenerator = new MersenneTwisterFast();
		
		//1. Build the count table models.
		tableModels = new ArrayList<CountTable>();
		int nExperiments = config.getJSONObject("modelSettings").getJSONArray("enrichmentModel").length();
		for(int iExp=0; iExp<nExperiments; iExp++) {
			
			tableModels.add(new CountTable(config, iExp, sc, letterComplement, letterOrder, true));
			tableModels.get(iExp).includeComponent = true;
			tableModels.get(iExp).batchFixedLibrarySize = fixedLibrarySize;
		}
		
		//2. Build binding modes.
		bindingModes     = new ArrayList<BindingMode>();
		int nBindingModes = config.getJSONObject("modelSettings").getJSONArray("bindingModes").length();
		for(int iBM=0; iBM<nBindingModes; iBM++) {
			bindingModes.add(new BindingMode(config, iBM, sc, letterComplement, letterOrder));
			for(int iExp=0; iExp<nExperiments; iExp++) 
				bindingModes.get(iBM).addCountTable(tableModels.get(iExp));
		}
		
		//3. Build binding model interactions.
		interactions     = new ArrayList<BindingModeInteraction>();
		int nInteractions = config.getJSONObject("modelSettings").getJSONArray("bindingModeInteractions").length();
		for(int iInt=0; iInt<nInteractions; iInt++) {
			interactions.add(new BindingModeInteraction(config, iInt, bindingModes));
			//Links binding mode interaction to the experiment
			for(int iExp=0; iExp<nExperiments; iExp++) 
				interactions.get(iInt).addCountTable(tableModels.get(iExp));
		}

		//4. Build the enrichment models.
		enrichmentModels = new ArrayList<EnrichmentModel>();
		for(int iExp=0; iExp<nExperiments; iExp++) {
			enrichmentModels.add(EnrichmentModel.buildEnrichmentModel(config, iExp, tableModels.get(iExp), bindingModes, interactions));
			
			enrichmentModels.get(iExp).includeComponent = true;
			tableModels.get(iExp).setEnrichmentModel(enrichmentModels.get(iExp));
		}
		//5. Compiles list of fitting components.
		componentList    = new ArrayList<ModelComponent>();
		for(int iBM=0;  iBM<bindingModes.size();      iBM++ ) componentList.add(bindingModes.get(iBM));
		for(int iInt=0; iInt<interactions.size();     iInt++) componentList.add(interactions.get(iInt));
		for(int iEnr=0; iEnr<enrichmentModels.size(); iEnr++) componentList.add(enrichmentModels.get(iEnr));
		for(int iTab=0; iTab<tableModels.size(); iTab++)      componentList.add(tableModels.get(iTab));
		
		//6. Seeds model.
		for(int iComp=0; iComp<componentList.size(); iComp++)
			componentList.get(iComp).seed_component(config);
		
		//7. Load parameters
		addBindingModesSequentially = config.getJSONObject("modelFittingConstraints").getBoolean("addBindingModesSequentially");

		//8. Pack models, saves parameters.
		updatePacking();
		parameters   = new double[nParameters];
		lastReported = null;
		ModelComponent.packParameters_O(getJSONState().getJSONObject("coefficients"), packing.getJSONObject("packing"), parameters);
		
		//9. Determines fitting order of the binding modes and the interactions
		if(verbose)
			System.out.println(">> Determining fitting order.");
		fittingOrder = new ArrayList<ArrayList<ModelComponent>>();
		//9.1 Adds binding modes and, as the interactions become possible, binding mode interactions.
		ArrayList<Boolean> bmAdded  = new ArrayList<Boolean>(Collections.nCopies(bindingModes.size(), false));
		ArrayList<Boolean> intAdded = new ArrayList<Boolean>(Collections.nCopies(interactions.size(), false));
		ArrayList<ModelComponent> temp;
		for(int iBM=0; iBM<bindingModes.size(); iBM++) {
			//Adding the new interaction.
			temp = new ArrayList<ModelComponent>();
			temp.add(bindingModes.get(iBM));
			fittingOrder.add(temp);
			bmAdded.set(iBM,true);
			//Checking if any new interaction should be unlocked
			for(int iInt=0; iInt<interactions.size(); iInt++) 
				if(!intAdded.get(iInt) && bmAdded.get(interactions.get(iInt).i0) && bmAdded.get(interactions.get(iInt).i1)) {
					temp = new ArrayList<ModelComponent>();
					temp.add(interactions.get(iInt));
					fittingOrder.add(temp);
					intAdded.set(iInt,true);
				}
		}
		//9.2 Adds the enrichment models (only gives non-trivial variations if bindingSaturation=false and trySaturation=true)
		ArrayList<ModelComponent> expKinTempList = new ArrayList<ModelComponent>();//All ExponentialKineticsModel are fitted at the same time.
		for(EnrichmentModel em: enrichmentModels) {
			if( em instanceof SELEXModel) {
				SELEXModel sm = (SELEXModel) em;
				if(sm.bindingSaturation == false&&sm.trySaturation     == true) {
					temp = new ArrayList<ModelComponent>();
					temp.add(em);
					fittingOrder.add(temp);
				}
			} else if(em instanceof ExponentialKineticsModel) {
				if(em.freezingLevel>0)
					expKinTempList.add(em);
			}
		}
		if(expKinTempList.size()>0)
			fittingOrder.add(expKinTempList);
		
		//10.1 Writes what components are included in what experiments:
		if(verbose) {
			System.out.println("");
			System.out.println(" Summary of experiments ");
			System.out.println(" ====================== ");
			for(int iExp=0; iExp<tableModels.size(); iExp++) {
				System.out.println("Experiment "+iExp+":");
				System.out.println("-------------");
				System.out.println("Count table:      "+tableModels.get(iExp).componentName);
				System.out.println("Enrichment model: "+enrichmentModels.get(iExp).componentName);
				System.out.println("   Concentration: "+enrichmentModels.get(iExp).concentration);
				
				System.out.println("Binding modes: ");
				if(enrichmentModels.get(iExp).bindingModes.size()>0)
					for(int iBM=0; iBM<enrichmentModels.get(iExp).bindingModes.size(); iBM++)
						System.out.println("                  "+enrichmentModels.get(iExp).bindingModes.get(iBM).componentName);
				else
				System.out.println("                  NONE");
				System.out.println("Binding mode interactions: ");
				if(enrichmentModels.get(iExp).interactions.size()>0)
					for(int iInt=0; iInt<enrichmentModels.get(iExp).interactions.size(); iInt++)
						System.out.println("                  "+enrichmentModels.get(iExp).interactions.get(iInt).componentName);
				else
					System.out.println("                  NONE");
				System.out.println("");
			}
		}

		resetMetadataVariables();
		
	}
		
	
	
	public void updatePacking() {
		
		packing              = ModelComponent.packModel(componentList);
		Integer lastIndex    = ModelComponent.maxIndex_O(packing.getJSONObject("packing"));
		nParameters          = lastIndex==null ? 0 : lastIndex.intValue() + 1;
		parameters           = new double[nParameters];
		gradient             = new double[nParameters];
		valueUpdated         = false;
		
		ModelComponent.packParameters_O(
				getJSONState().getJSONObject("coefficients"), 
				packing.getJSONObject("packing"), 
				parameters);
	}
	
	//Creates a new JSON object containing the current state.
	public JSONObject getJSONState() {
		JSONObject out   = ModelComponent.saveToJSON(componentList, "coefficients"); 
		JSONObject oSett = out.getJSONObject("modelSettings"); 
		oSett.put("letterComplement", letterComplement);
		oSett.put("letterOrder",      letterOrder);

		
		JSONObject oOpt = new JSONObject();
		out.put("optimizerSetting", oOpt);
		oOpt.put("lambdaL2",         lambdaL2);
		oOpt.put("pseudocount",      pseudocount);
		oOpt.put("expBound",         expBound);
		oOpt.put("fixedLibrarySize", fixedLibrarySize);
		
		return out;
	}

	//Returns the current state, including only added components
	public void setStateJSON(JSONObject in) {

		valueUpdated    = false;
		
		ModelComponent.readFromJSON(in, "coefficients", componentList);

		return;
	}
	
	public ArrayList<ModelComponent> getNextFitComponent() {
		if(fittingOrder==null)
			return null;
		ArrayList<ModelComponent> out = fittingOrder.get(0);
		if(fittingOrder.size() > 0)
			fittingOrder.remove(0);
		else
			fittingOrder = null;
		return out;
		
	}
	
	// Functions for implementing Loss Function
	///////////////////////////////////////////
	
	// Function for setting parameters
	@Override
	public void lossFunction_setParameters(double[] in) {

		//Saves the raw parameters
		parameters = Array.clone(in);

		//Builds a JSON object with updated parameters.
		JSONObject unpackedParameters = getJSONState();
		ModelComponent.unpackParameters_O(parameters, packing.getJSONObject("packing"), unpackedParameters.getJSONObject("coefficients"));
		
		//Updates the parameters in the model.
		ModelComponent.readFromJSON(unpackedParameters, "coefficients", componentList);
		valueUpdated    = false;
	}
	
	// Function for updating function value
	@Override
	public void lossFunction_updateValue() {
		
		if(valueUpdated)
			return;
		value         = 0;
		logLikelihood = 0;
		
		if(computeVariance)
			valueVariance = 0;

		for(int iTable=0; iTable<tableModels.size(); iTable++) {
			CountTable table         = tableModels.get(iTable);
			if(table.includeComponent) {
				try {
					table.updateValue();
				} catch (Exception e) {
					e.printStackTrace();
				}
				double meanFunctionValue = table.functionValue * table.weight;
				value           += meanFunctionValue;
				logLikelihood   += table.functionValue;
				
				if(computeVariance) {
					//TODO: Q: Is the function value and gradient normalized by the number of reads??
					double meanFunctionValueSquared = table.functionValueSquared * table.weight;
					valueVariance += meanFunctionValueSquared - meanFunctionValue * meanFunctionValue;
				}
			}
		}
		
		//Adds regularization (to value, not to logLikelihoodPerRead and to logLikelihood).
		//TODO: Check what regularization is best.
		JSONObject oReg = getJSONState();

		double L2Sum = ModelComponent.trJSONObject_O(ModelComponent.multiplyJSONObject_O(oReg.getJSONObject("coefficients"), oReg.getJSONObject("coefficients")));
		double regularizationTerm = 0;
		
		//L2 regularization
		regularizationTerm       += lambdaL2*L2Sum;
		
		//Dirichlet 
		for(BindingMode bm: bindingModes)
			if(bm.includeComponent)
				regularizationTerm += bm.getDirichletValue(pseudocount);
		
		//Exponential bounding term: Exp[value-expBound] + Exp[-value-expBound]
		regularizationTerm       += Array.sum(Array.exp(Array.add(     ModelComponent.constant_d(parameters, -expBound), parameters)))
						       	  + Array.sum(Array.exp(Array.subtract(ModelComponent.constant_d(parameters, -expBound), parameters)));
		
		value += regularizationTerm;
		
		logLikelihoodPerRead = value;
		
		functionCalls += 1;
		dataLoops     += 1.0 * this.lossFunction_getBatchSize() / this.lossFunction_getDataSize();
		valueUpdated   = true; 
		

		
		
		
	}

	
	// Function that is called every time the position is updated.
	@Override
	public void lossFunction_reportPosition() {
		

		if(printTrajectory){
			boolean append = true;
			PrintStream original	= System.out;
			
			//Changes output stream.
			if(trajectoryFile != null){
				try {
					System.setOut(new PrintStream(new FileOutputStream(trajectoryFile, append)));
					
				} catch (FileNotFoundException e) {
					System.out.println("Cannot create trajectory file at this "
							+ "location: "+trajectoryFile);
					e.printStackTrace();
				}
			
			}

			System.out.println(Misc.formatVector_d(parameters,",","",""));
			
			//Returns to original stream.
			System.setOut(original);

		}
		 
		lastReported  = Array.clone(parameters);
		fitSteps     += 1;

	}
	
	
	@Override
	
	public void lossFunction_nextBatch(int nDataIn) {
		//TODO: What to do when not all columns are used??
		
		//Counts the total number of reads
		int[] readsInTable = new int[tableModels.size()];
		int nReadsTot      = 0;
		int nTables        = tableModels.size();
		
		for(int iTable=0; iTable<nTables; iTable++) {
			CountTable c = tableModels.get(iTable);
			for(int r : c.modeledColumns) {
				int nReads            = c.fullTable.countPerRound[r];
				readsInTable[iTable] += nReads;
				nReadsTot            += nReads;
			}
		}
		
		//Checks if all reads should be used (nData=-1)
		int nData = nDataIn==-1 || nDataIn>= nReadsTot ? -1 : nDataIn;
		
		if(nData==-1) {
			//Uses all data if nData=-1;
			for(CountTable tm : tableModels)
				tm.nextBatchAll();
		} else {
			
			if(fixedLibrarySize) {

				//CASE 1: Uses an equal fraction of reads in each data set.
				///////////////////////////////////////////////////////////
				//Creates a new batch by requesting a specific number of reads from each table:
				for(CountTable tm : tableModels) {
					tm.batchFixedLibrarySize = true;
					tm.nextBatch((int) (1.0 * nData * (tm.nReads / nReadsTot)));
				}

			} else {
				//CASE 2: Draws reads randomly across all tables.
				//////////////////////////////////////////////////
				//STEP 1: Creates a list of random reads (with index 0...(nReadsTot-1))
				int[] randomReads = Misc.randomSample(nReadsTot, nData, randGenerator);
				
				//STEP 2: Sorts the reads by the tables.
				//Allocates variables
				ArrayList<ArrayList<Integer>> randomReadsInTables = new ArrayList<ArrayList<Integer>>();
				for(int iTable=0; iTable<nTables; iTable++) 
					randomReadsInTables.add(new ArrayList<Integer>());
				//Loops over reads and sorts them into the correct table.
				for(int iR : randomReads) {
					int iRead = iR;
					for(int iTable=0; iTable<nTables; iTable++) {
						if(iRead < readsInTable[iTable]) {
							randomReadsInTables.get(iTable).add(iRead);
							break;
						} else {
							iRead -= readsInTable[iTable];
						}
					}
				}

				//STEP 3: Creates new batches for the count tables with the specified reads
				for(int iTable=0; iTable<nTables; iTable++) {
					double expectedTableCount = 1.0 * nData * readsInTable[iTable] / nReadsTot;
					ArrayList<Integer> aI = randomReadsInTables.get(iTable);
					int[] tableReads = new int[aI.size()];
					for(int iI=0; iI<aI.size(); iI++)
						tableReads[iI] = aI.get(iI);
					tableModels.get(iTable).nextBatch(tableReads, expectedTableCount );
				}
			}
		}
		
		valueUpdated    = false;
				
	}


	//Functions for writing trajectory file
	void startNewTrajectoryFile(String path, JSONObject packing) {

		//Saves parameters.
		printTrajectory = true;
		trajectoryFile  = path;
		
		//Writes header 
		PrintStream original	= System.out;
		boolean append=false;
		if(trajectoryFile != null){
			try {
				System.setOut(new PrintStream(new FileOutputStream(trajectoryFile, append)));
				System.out.println(packing.toString());
			} catch (FileNotFoundException e) {
				System.out.println("Cannot create trajectory file at this "
						+ "location: "+trajectoryFile);
				e.printStackTrace();
			}
		}
		
		//Returns to original stream.
		System.setOut(original);

	}
	
	void stopTrajectoryOutput() {
		
		printTrajectory = false;
		trajectoryFile  = null;
		
	}

	public void resetMetadataVariables() {
		fitStartTime = System.currentTimeMillis();
		fitSteps = 0;
		functionCalls = 0;
		gradientCalls = 0;

	}
	
	
	public JSONObject getMetadataJSON(Integer index, String name) {
		
    	JSONObject metadata = new JSONObject();
    	if(index != null)
    		metadata.put("index", index);
    	metadata.put("timeStamp", (new Timestamp(System.currentTimeMillis())).toString());
    	metadata.put("fitter", "proBound");
    	if(name != null)
    		metadata.put("fitName", name);
    	metadata.put("fitTime", (System.currentTimeMillis()-fitStartTime)/1000);
    	metadata.put("fitSteps", fitSteps);
    	metadata.put("functionCalls", functionCalls);
    	metadata.put("gradientCalls", gradientCalls);
    	metadata.put("dataLoops",     dataLoops);
    	metadata.put("logLikelihood", logLikelihood);
    	metadata.put("logLikelihoodPerRead", logLikelihoodPerRead);
    	metadata.put("regularization", regularizationTerm);
        
        return metadata;
	}

	@Override
	public long lossFunction_getDataSize() {
		long nReads = 0; 
		for(CountTable m : tableModels) 
			nReads += m.nReads;
			
		return nReads;
	}

	public static String alphabetDefToOrderedLetters(String alphabetDef) {
		String[] tempLetters = alphabetDef.replace("-", ",").split(",");
		Set<String> letterSet = new HashSet<String>();
		for(String li: tempLetters)letterSet.add(li);
		String[] tempUnion = letterSet.toArray(new String[letterSet.size()]);
		Arrays.sort(tempUnion);
		return Misc.joinStrings(tempUnion);
	}
	
	
	
    public static void removeBindingMode_JSON(JSONObject oModel, int iBM) {

    	//1. Removes binding mode.
    	String[] mainKeys = {"metadata", "modelError", "modelSeeding", "scoringModel", "modelSettings", "coefficients", "modelFittingConstraints" };
    	for(String k : mainKeys)
        	if( oModel.has(k) && oModel.getJSONObject(k).has("bindingModes")) {
        		JSONArray aBM = oModel.getJSONObject(k).getJSONArray("bindingModes");
        		if( aBM.length()>iBM ) 
        			aBM.remove(iBM);
        		else
        			System.err.println("WARNING: Binding mode "+iBM+" does not exist in "+k+" and cannot be removed.");
        	}
        		
    	//2. Updates enrichment model to remove the binding mode.
    	if( oModel.getJSONObject("modelSettings").has("enrichmentModel") ) {
    		JSONArray aEnr = oModel.getJSONObject("modelSettings").getJSONArray("enrichmentModel");
    		for(int iEnr=0; iEnr<aEnr.length(); iEnr++) {
    			JSONObject oEnr = aEnr.getJSONObject(iEnr);
    			if(oEnr.has("bindingModes")) {
    				JSONArray aBM   = oEnr.getJSONArray("bindingModes");
    				if(aBM.length()>0 && aBM.getInt(0)!=-1){
    					for(int iiBM=aBM.length()-1; iiBM>=0; iiBM--) {
    						int iBMtemp = aBM.getInt(iiBM);
    						if(iBMtemp>iBM)
    							aBM.put(iiBM, iBMtemp-1);
    						else if(iBMtemp==iBM)
    							aBM.remove(iiBM);
    					}
    				}
    			}
    		}
    	}
    		
    	//3. Removes all interactions that involve the removed binding mode, update the index of binding mode after...
    	Set<Integer> iIntRemove = new HashSet<Integer>();
    	if(oModel.getJSONObject("modelSettings").has("bindingModeInteractions")) {
    		JSONArray aInt                = oModel.getJSONObject("modelSettings").getJSONArray("bindingModeInteractions");

        	for(int iInt=0; iInt<aInt.length(); iInt++) {
        		JSONObject oInt           = aInt.getJSONObject(iInt);
        		JSONArray aBM             = oInt.getJSONArray("bindingModes");
        		for(int i=0; i<aBM.length(); i++) {
        			int iBMtemp = aBM.getInt(i);
        			if(iBMtemp > iBM)
        				aBM.put(i, iBMtemp-1);
        			else if(iBMtemp==iBM) 
        				iIntRemove.add(iInt);
        		}
        	}
    	}
    	//Sorts the binding modes that are removed
		Integer[] sortedIIntRemove = iIntRemove.toArray(new Integer[iIntRemove.size()]);
		Arrays.sort(sortedIIntRemove);
		//Removes the interactions
		for(int iInt : sortedIIntRemove) 
			removeBindingModeInteraction_JSON(oModel, iInt);
		
    }
    
    
    static public void removeBindingModeInteraction_JSON(JSONObject oModel, int iInt) {

    	//1. Removes binding mode interaction.
    	String[] mainKeys = {"metadata", "modelError", "modelSeeding", "scoringModel", "modelSettings", "coefficients", "modelFittingConstraints" };
    	for(String k : mainKeys)
        	if( oModel.has(k) && oModel.getJSONObject(k).has("bindingModesInteractions")) {
        		JSONArray aInt = oModel.getJSONObject(k).getJSONArray("bindingModesInteractions");
        		if( aInt.length()>iInt ) 
        			aInt.remove(iInt);
        		else
        			System.err.println("WARNING: Binding mode interaction "+iInt+" does not exist in "+k+" and cannot be removed.");
        	}
        		
    	//2. Updates enrichment model to remove the binding mode interaction.
    	if( oModel.getJSONObject("modelSettings").has("enrichmentModel") ) {
    		JSONArray aEnr = oModel.getJSONObject("modelSettings").getJSONArray("enrichmentModel");
    		for(int iEnr=0; iEnr<aEnr.length(); iEnr++) {
    			JSONObject oEnr = aEnr.getJSONObject(iEnr);
    			if(oEnr.has("bindingModesInteractions")) {
    				JSONArray aInt   = oEnr.getJSONArray("bindingModesInteractions");
    				if(aInt.length()>0 && aInt.getInt(0)!=-1){
    					for(int iiInt=aInt.length()-1; iiInt>=0; iiInt--) {
    						int iBMtemp = aInt.getInt(iiInt);
    						if(iBMtemp>iInt)
    							aInt.put(iiInt, iBMtemp-1);
    						else if(iBMtemp==iInt)
    							aInt.remove(iiInt);
    					}
    				}
    			}
    		}
    	}
    }


	@Override
	public int lossFunction_getBatchSize() {
		int N = 0;
		for(CountTable c : tableModels)
			for(int r : c.modeledColumns)
				N += c.batchCountPerRound[r];
		return N;

	}
    		
		
}

