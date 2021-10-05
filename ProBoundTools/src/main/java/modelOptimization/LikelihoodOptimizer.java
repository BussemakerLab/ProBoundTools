package modelOptimization;


import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.PrintStream;

import org.json.*;



public class LikelihoodOptimizer {
	
	//GOAL: Abstract class containing code:
	// - Running optimization.
    // - Writing trajectory etc.
	// - Knows about what binding modes are updated / optimized.
	// - 
	//
	
	
	// QUESTION: How do we coordinate the fitting strategy?
	// APPROACH:
	//  1) Each model component (binding mode, binding mode interaction, or rho, gamma or h) has:
	//    - boolean included;        //Indicates if this component has been included yet.
	//    - boolean optimized;       //Indicates if the component has been fully optimized yet.
	//    - boolean fit;             //Indicates if the binding mode should be fit.
	// 2) Given this, we use getJSONState() and setStateJSON() to record and update all INCLUDED components.
	// 3) LikelihoodOptimizer creates an ordered indicating when the binding modes should be included.
	//    - Binding modes are ordered in order of ID
	//    - Once the relevant interactions have been included, the interactions are included.
	//    - Potentially we could add one experiment at a time.
	// 4) Each binding component can somehow create a list of possible variations.
	//    - ArrayList<JSONObject> modelComponent.getVariations(JSONObject in) creates a list of JSON objects that encode variations.
	//       - If this returns null when the component is fully optimized. 
	// 5) Add optional fitting strategy:
	//    - Default is to add sequentially/after unlocking. All experiments are added initially.
	//    - In the model specification, add option:
	//           [ {"class":"bindingMode",             "id":0, "initial": true},
	//             {"class":"bindingMode",             "id":1, "initial": false},
	//             {"class":"bindingMode",             "id":2, "initial": false},
	//             {"class":"bindingModeInteractions", "id":0, "initial": false},
	//             {"class":"experiment",              "id":0, "initial": false} ]

	
	//Variables for setting up the numerical optimization object
	////////////////////////////////////////////////////////////
	modelOptimization.Optimizer o;
	String minimizerType = null;
	
	//SGD Settings
	String sgdMethod;
	
	// Settings for combined likelihood
	///////////////////////////////////
	double lambdaL2, expBound, likelihoodThreshold;
	int nThreads, nRetries;
	boolean fixedLibrarySize;
	
	//Settings for output
	/////////////////////
	String outputPath, baseName;
	boolean verbose, storeHessian, printPSAM, printTrajectory;
	
	public LikelihoodOptimizer(JSONObject config) {
	
		JSONObject oOptSet = config.getJSONObject("optimizerSetting");
		

		//1. Determines where to write? Makes the CombinedLikelihood  
		if(oOptSet.has("output")) {
			JSONObject oOut = oOptSet.getJSONObject("output");
			verbose         = oOut.getBoolean("verbose");	
			outputPath      = oOut.has("outputPath") ? oOut.getString("outputPath") : null;
			baseName        = oOut.has("baseName")   ? oOut.getString("baseName")   : null;
			storeHessian    = oOut.getBoolean("storeHessian");
			printPSAM       = oOut.getBoolean("printPSAM");
			printTrajectory = oOut.getBoolean("printTrajectory");			
		} else {
			verbose         = false;	
			outputPath      = null;
			baseName        = null;
			storeHessian    = false;
			printPSAM       = false;
			printTrajectory = false;			
		}

		
		//2. Based on config, builds optimizer object. 
		
		minimizerType       = oOptSet.getString("minimizerType");
		lambdaL2            = oOptSet.getDouble("lambdaL2");
		expBound            = oOptSet.getDouble("expBound");
		nThreads            = oOptSet.getInt("nThreads");
		nRetries            = oOptSet.getInt("nRetries");
		likelihoodThreshold = oOptSet.getDouble("likelihoodThreshold");
		
		fixedLibrarySize = oOptSet.getBoolean("fixedLibrarySize"); 

	}
	
	//Writes the current model as a single-line JSON object. 
	public static void writeCompactJSONModel(JSONObject model, String outPath, boolean append) {
		PrintStream original	= System.out;
		
		//Changes output stream.
		if(outPath != null){
			try {
				System.setOut(new PrintStream(new FileOutputStream(outPath, append)));
				
			} catch (FileNotFoundException e) {
				System.out.println("Cannot create trajectory file at this "
						+ "location: "+outPath);
				e.printStackTrace();
			}
		
		}
		
		System.out.println(model.toString());
		
		//Returns to original stream.
		System.setOut(original);
	}
	
}