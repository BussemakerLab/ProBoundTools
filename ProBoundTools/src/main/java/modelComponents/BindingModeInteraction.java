package modelComponents;

import java.util.ArrayList;

import org.json.*;

import base.Array;
import proBoundTools.Misc;

public class BindingModeInteraction extends ModelComponent {
	
	public int i0, i1;
	BindingMode b0, b1;
	public int wPositionBiasBin0, wPositionBiasBin1;
	int maxIntOverlap, maxIntSpacing;
	boolean usePositionBias;
	boolean roundSpecificActivity, experimentSpecificInteraction, experimentSpecificActivity;
	boolean fitLogActivity; 
	
	//Experiment-specific parameters
	public ArrayList<CountTable> countTables;
	public ArrayList<ArrayList<ArrayList<double[][]>>> interactionBetas, interactionAlphas;
	public ArrayList<double[]> activityBetas, activityAlphas;
	
	
	public BindingModeInteraction(JSONObject config, int iIntIn, ArrayList<BindingMode> bindingModes) {
		super("bindingModeInteractions");
		
		iComp         = iIntIn;
		componentName = "Binding mode interaction "+iComp;
		if(verbose)
			System.out.println(">> Creating "+componentName+".");
		
		//Initially we have no experiments.
		countTables = new  ArrayList<CountTable>();
		maxFreezeLevel = 0;
		readFromJSON_settings(config);
		readFromJSON_constraints(config);

		//Reads some parameters from the binding modes, cannot be changed after being set.
		b0                    = bindingModes.get(i0);
		b1                    = bindingModes.get(i1);
		wPositionBiasBin0     = b0.wPositionBiasBin;
		wPositionBiasBin1     = b1.wPositionBiasBin;

		seed_component(config);
		
		setFreezingLevel(0);
		
	}

	public void allocateParameters() {
		
		activityBetas     = new ArrayList<double[]>();
		interactionBetas  = new ArrayList<ArrayList<ArrayList<double[][]>>>();
		
		for(int iExp=0; iExp<countTables.size(); iExp++) {
			//Creates activity variables
			activityBetas.add( fitLogActivity ? zero_d(countTables.get(iExp).nColumns) : null );

			//Adds interaction matrices.
			interactionBetas.add( zero_AAdd(2, 2, b0.maxFrames.get(iExp)/2,  b1.maxFrames.get(iExp)/2) ); 
		}
		
		updateAlphas();
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
		if(fitLogActivity)
			activityAlphas = exp_Ad(activityBetas);

		interactionAlphas = new ArrayList<ArrayList<ArrayList<double[][]>>>();
		for(int iExp=0; iExp<countTables.size(); iExp++) {
			ArrayList<ArrayList<double[][]>> temp = exp_AAdd(interactionBetas.get(iExp));
			maskInteractionMatrix(temp);
			interactionAlphas.add(temp);
		}
		
	}
	
	
	//Takes a alphaInteraction matrix, and prohibited entries.
	public void maskInteractionMatrix(ArrayList<ArrayList<double[][]>> interactionAlphas) {

		int nF0=interactionAlphas.get(0).get(0).length;
		int nF1=interactionAlphas.get(0).get(0)[0].length;
		
		for(int s0=0; s0<2; s0++) 
			for(int s1=0; s1<2; s1++) 
				for(int x0=0; x0<nF0; x0++) 
					for(int x1=0; x1<nF1; x1++) {
						//Location of binding-mode frame on positive strand
						int x0PositiveStart = (s0==0 ? x0 : nF0-1 - x0) - b0.flankLength;
						int x1PositiveStart = (s1==0 ? x1 : nF1-1 - x1) - b1.flankLength;
						int x0PositiveEnd   = x0PositiveStart + b0.k;
						int x1PositiveEnd   = x1PositiveStart + b1.k;
						
						//Computes the number of overlapping basepairs
						int spacing = Math.max(x1PositiveStart - x0PositiveEnd, 
								               x0PositiveStart - x1PositiveEnd);
						int overlap = -spacing;	
						if( ((  spacing > maxIntSpacing && !(maxIntSpacing==-1))) ||
								overlap > maxIntOverlap) 
							interactionAlphas.get(s0).get(s1)[x0][x1] = 0;
					}
		
		return;
	}
	
	//Takes a interactionPacking matrix, and sets prohibited entries to -1.
	public void maskInteractionPackingMatrix(ArrayList<ArrayList<int[][]>> interactionPacking) {

		int nF0=interactionPacking.get(0).get(0).length;
		int nF1=interactionPacking.get(0).get(0)[0].length;
		
		for(int s0=0; s0<2; s0++) 
			for(int s1=0; s1<2; s1++) 
				for(int x0=0; x0<nF0; x0++) 
					for(int x1=0; x1<nF1; x1++) {
						//Location of binding-mode frame on positive strand
						int x0PositiveStart = (s0==0 ? x0 : nF0-1 - x0) - b0.flankLength;
						int x1PositiveStart = (s1==0 ? x1 : nF1-1 - x1) - b1.flankLength;
						int x0PositiveEnd   = x0PositiveStart + b0.k;
						int x1PositiveEnd   = x1PositiveStart + b1.k;
						
						//Computes the number of overlapping basepairs
						int spacing = Math.max(x1PositiveStart - x0PositiveEnd, 
								               x0PositiveStart - x1PositiveEnd);
						int overlap = -spacing;	
						if( ((  spacing > maxIntSpacing && !(maxIntSpacing==-1))) ||
								overlap > maxIntOverlap) 
							interactionPacking.get(s0).get(s1)[x0][x1] = -1;
					}
		
		return;
	}
	
	//Takes a interactionPacking matrix, and sets prohibited entries to -1.
	public static void makeInteractionPackingHomogenous(ArrayList<ArrayList<int[][]>> interactionPacking) {

		int nF0=interactionPacking.get(0).get(0).length;
		int nF1=interactionPacking.get(0).get(0)[0].length;
		
		for(int s0=0; s0<2; s0++) {
			for(int s1=0; s1<2; s1++) { 
				for(int x0=0; x0<nF0; x0++) { 
					for(int x1=0; x1<nF1; x1++) {
						
						if(s0==s1) { //On same strand
							if(x1>=x0)
								interactionPacking.get(s0).get(s1)[x0][x1] = interactionPacking.get(0).get(0)[0     ][x1-x0           ];
							else
								interactionPacking.get(s0).get(s1)[x0][x1] = interactionPacking.get(0).get(0)[x0-x1 ][0               ];
							
						} else { //On opposite strand
							int x1p = (nF1-1) - x1;
							if(x1p > x0)
								interactionPacking.get(s0).get(s1)[x0][x1] = interactionPacking.get(0).get(1)[0     ][(nF1-1)-(x1p-x0)];
							else
								interactionPacking.get(s0).get(s1)[x0][x1] = interactionPacking.get(0).get(1)[x0-x1p][nF1-1           ];
						}
					}
				}
			}
		}
		
		return;
	}

	public void addCountTable(CountTable exp) {
		countTables.add(exp);
		allocateParameters();
	}
		
	@Override
	public void seed_component(JSONObject config) {

		String coefficientKey = "modelSeeding";
		if(config.has(coefficientKey)) {
			JSONObject oSeed = config.getJSONObject(coefficientKey); 
			if( oSeed.has(componentKey) ) {
				JSONObject oInt = oSeed.getJSONArray(componentKey).getJSONObject(iComp);

				//TODO: If the activity seed is too short, then repeat.
				
				if(oInt.has("activity")) {
					if(fitLogActivity)
						activityBetas = readFromJSON_Ad(   oInt.getJSONArray("activity"));
					else
						activityAlphas= readFromJSON_Ad(   oInt.getJSONArray("activity"));
				}
				
				if(oInt.has("positionMatrix"))
					interactionBetas  = readFromJSON_AAAdd(oInt.getJSONArray("positionMatrix"));

				if(oInt.has("spacingVector")); 
					//TODO: Implement conversion from spacing vector to position matrix.

				updateAlphas();
			}
		}
	}
	
	@Override
	int packModel_component(JSONObject packing, int iFirst) {
		
		String coefficientKey = "packing";
		addEmptyJSON_component_O(packing, coefficientKey, componentKey, iComp);
		
		int iCurr = iFirst;

		if(fitComponent) {
			Integer iLast;

			ArrayList<ArrayList<ArrayList<int[][]>>> interactionRange = range_AAAdd(interactionBetas, iCurr+1);
			iLast = last_AAAii(interactionRange);
			iCurr = iLast==null ? iCurr : iLast;
			
			//Removes prohibited interactions.
			for(int iExp=0; iExp<interactionRange.size(); iExp++) {
				//Masks interactions that are prohibited.
				maskInteractionPackingMatrix(interactionRange.get(iExp));
			}

			if(!usePositionBias) 
				for(int iExp=0; iExp<interactionRange.size(); iExp++)
					makeInteractionPackingHomogenous(interactionRange.get(iExp));
			
			//Imposes experimentSpecificInteraction=False
			if(!experimentSpecificInteraction) {
				for(int iExp=0; iExp<interactionRange.size(); iExp++) {
					for(int iS0=0; iS0<2; iS0++) {
						for(int iS1=0; iS1<2; iS1++) {
							int[][] v0 = interactionRange.get(0).get(iS0).get(iS1);
							int[][] v1 = interactionRange.get(iExp).get(iS0).get(iS1);
							if(v0.length != v1.length || v0[0].length != v1[0].length)
								throw new java.lang.RuntimeException("ERROR: All experiments must have the same number of binding frames if BindingModeInteractions.experimentSpecificInteraction=False.");

							for(int iX0=0; iX0<v0.length; iX0++)
								for(int iX1=0; iX1<v0[0].length; iX1++)
									v1[iX0][iX1] = v0[iX0][iX1];
						}
					}
				}
			}
			
			saveToJSON_interaction_AAAii(packing, coefficientKey, interactionRange);
			
			ArrayList<int[]> activityRange = range_Ad(activityBetas, iCurr+1);
			iLast = last_Ai(activityRange);
			iCurr = iLast==null ? iCurr : iLast;
			
			//Removes round-specific activities if appropriate.
			if(!roundSpecificActivity)
				for(int iExp=0; iExp<activityRange.size(); iExp++)
					for(int iRound=0; iRound<activityRange.get(iExp).length; iRound++)
						activityRange.get(iExp)[iRound] = activityRange.get(iExp)[0]; 

			
			//Removes experiment-specific activities if appropriate.
			if(!experimentSpecificActivity && activityRange.size()>1) {
				for(int iExp=1; iExp<activityRange.size(); iExp++) {
					if(activityRange.get(0).length != activityRange.get(iExp).length)
						throw new java.lang.RuntimeException("ERROR: All experiments must have the rounds if BindingModeInteracion.experimentSpecificActivity=False.");

					for(int iRound=0; iRound<activityRange.get(iExp).length; iRound++)
						activityRange.get(0)[iRound] = activityRange.get(iExp)[iRound];
				}
			}
					
			saveToJSON_activity_Ai(packing, coefficientKey, activityRange);	
		}


		
		return iCurr;
	}
	
	@Override
	public void addZeroJSON_component(JSONObject in, String coefficientKey) {
		
		addEmptyJSON_component_O(in,      coefficientKey, componentKey, iComp);
		saveToJSON_activity_Ad(in,        coefficientKey, zero_Ad(activityBetas));
		saveToJSON_interaction_AAAdd(in,  coefficientKey, zero_AAAdd(interactionBetas));
		
		return;
		
	}

	@Override
	void saveToJSON_settings(JSONObject out) {
		
		String coefficientKey = "modelSettings";
		addEmptyJSON_component_O(out, coefficientKey, componentKey, iComp);
		JSONObject oInt = out.getJSONObject("modelSettings").getJSONArray(componentKey).getJSONObject(iComp);
		oInt.put("bindingModes", new JSONArray());
		oInt.getJSONArray("bindingModes").put(i0);
		oInt.getJSONArray("bindingModes").put(i1);
		oInt.put("maxOverlap",              maxIntOverlap);
		oInt.put("maxSpacing",              maxIntSpacing);
		oInt.put("positionBias",           usePositionBias);
		oInt.put("fitLogActivity",         fitLogActivity);
		
	}

	@Override
	void readFromJSON_settings(JSONObject in) {

		String coefficientKey = "modelSettings";
		JSONObject oInt = in.getJSONObject(coefficientKey).getJSONArray(componentKey).getJSONObject(iComp);
		i0              = oInt.getJSONArray("bindingModes").getInt(0);
		i1              = oInt.getJSONArray("bindingModes").getInt(1);
		maxIntOverlap   = oInt.getInt("maxOverlap");
		maxIntSpacing   = oInt.getInt("maxSpacing");
		usePositionBias = oInt.getBoolean("positionBias");
		fitLogActivity  = oInt.getBoolean("fitLogActivity");
		
		allocateParameters();
		
	}
	
	@Override
	void saveToJSON_constraints(JSONObject out) {
		
		String coefficientKey = "modelFittingConstraints";
		addEmptyJSON_component_O(out, coefficientKey, componentKey, iComp);
		JSONObject oExp       = out.getJSONObject(coefficientKey).getJSONArray(componentKey).getJSONObject(iComp);
		oExp.put("roundSpecificActivity",         roundSpecificActivity);
		oExp.put("experimentSpecificInteraction", experimentSpecificInteraction);
		oExp.put("experimentSpecificActivity",    experimentSpecificActivity);
		
	}

	@Override
	void readFromJSON_constraints(JSONObject in) {
		
		String coefficientKey = "modelFittingConstraints";
		if(in.has(coefficientKey)) {
			JSONObject oOpt       = in.getJSONObject(coefficientKey).getJSONArray(componentKey).getJSONObject(iComp);
			roundSpecificActivity = oOpt.getBoolean("roundSpecificActivity");
			experimentSpecificInteraction
			                      = oOpt.getBoolean("experimentSpecificInteraction");
			experimentSpecificActivity
			                      = oOpt.getBoolean("experimentSpecificActivity");
		}
		
	}
	
	@Override
	void saveToJSON_parameters(JSONObject out, String coefficientKey) {

		addEmptyJSON_component_O(out,     coefficientKey, componentKey, iComp);
		saveToJSON_activity_Ad(out,       coefficientKey, activityBetas);
		saveToJSON_interaction_AAAdd(out, coefficientKey, interactionBetas);

		return;
		
	}
	
	@Override
	void readFromJSON_parameters(JSONObject in, String coefficientKey) {
		JSONObject oInt       = in.getJSONObject(coefficientKey).getJSONArray(componentKey).getJSONObject(iComp);
		activityBetas         = readFromJSON_Ad(   oInt.getJSONArray("activity"));
		interactionBetas      = readFromJSON_AAAdd(oInt.getJSONArray("positionMatrix"));
		//TODO: Implement functions to read "spacingVector"?
		updateAlphas();
	}
	
	//Methods for saving double objects
	void saveToJSON_activity_Ad(JSONObject out, String coefficientKey, ArrayList<double[]> activityBetas) {
		JSONObject oBM = out.getJSONObject(coefficientKey).getJSONArray(componentKey).getJSONObject(iComp);
		oBM.put("activity", JSONArrayConvert_Ad(activityBetas));
		return;
	}
	
	void saveToJSON_activity_d(JSONObject out, String coefficientKey, double[] activityBetas, int iExp) {
		
		JSONObject oInt = out.getJSONObject(coefficientKey).getJSONArray(componentKey).getJSONObject(iComp);
		
		JSONArray aA;
		if(!oInt.has("activity"))
			oInt.put("activity", new JSONArray());
		aA = oInt.getJSONArray("activity");
		
		//Makes sure all he
		while(aA.length() < countTables.size())
			aA.put(new JSONArray());
		
		//Saves the specific betas.
		aA.put(iExp, JSONArrayConvert_d(activityBetas));

		return;
	}
	
	void saveToJSON_interaction_AAAdd(JSONObject out, String coefficientKey, ArrayList<ArrayList<ArrayList<double[][]>>> interactionBetas) {
		JSONObject oBM = out.getJSONObject(coefficientKey).getJSONArray(componentKey).getJSONObject(iComp);
		oBM.put("positionMatrix", JSONArrayConvert_AAAdd(interactionBetas));
		return;
	}
	
	void saveToJSON_interaction_AAdd(JSONObject out, String coefficientKey, ArrayList<ArrayList<double[][]>> interactionBetas, int iExp) {
		
		JSONObject oInt = out.getJSONObject(coefficientKey).getJSONArray(componentKey).getJSONObject(iComp);
		
		JSONArray aPb;
		if(!oInt.has("positionMatrix"))
			oInt.put("positionMatrix", new JSONArray());
		aPb = oInt.getJSONArray("positionMatrix");
		
		//Makes sure all he
		while(aPb.length() < countTables.size())
			aPb.put(new JSONArray());
		
		//Saves the specific betas.
		aPb.put(iExp, JSONArrayConvert_AAdd(interactionBetas));

		return;
	}

	//Methods for saving long objects
	void saveToJSON_activity_Ai(JSONObject out, String coefficientKey, ArrayList<int[]> activityBetas) {
		JSONObject oBM = out.getJSONObject(coefficientKey).getJSONArray(componentKey).getJSONObject(iComp);
		oBM.put("activity", JSONArrayConvert_Ai(activityBetas));
		return;
	}

	
	void saveToJSON_interaction_AAAii(JSONObject out, String coefficientKey, ArrayList<ArrayList<ArrayList<int[][]>>> interactionBetas) {
		JSONObject oBM = out.getJSONObject(coefficientKey).getJSONArray(componentKey).getJSONObject(iComp);
		oBM.put("positionMatrix", JSONArrayConvert_AAAii(interactionBetas));
		return;
	}
	
	
	static public void printJSONObjectCoefficients(JSONObject model, String coefficientKey, int iBM) {
		
		JSONObject oInt = model.getJSONObject(coefficientKey).getJSONArray("bindingModeInteractions").getJSONObject(iBM);
				
		if(oInt.has("activity")) {
			JSONArray aAct = oInt.getJSONArray("activity");
			for(int iExp=0; iExp<aAct.length(); iExp++) {
				System.out.println("Activity(exp="+iExp+"):   "
						+Misc.formatVector_d(ModelComponent.readFromJSON_d(aAct.getJSONArray(iExp))));
			}
		}

		if(oInt.has("positionMatrix")) {
			JSONArray aPM = oInt.getJSONArray("positionMatrix");
			for(int iExp=0; iExp<aPM.length(); iExp++) {
				for(int s1=0; s1<2; s1++) {
					for(int s2=0; s2<2; s2++) {
						System.out.println("Position Matrix(exp="+iExp+", strand 1="+s1+", strand 2="+s2+"):");
						JSONArray aMatrix = aPM.getJSONArray(iExp).getJSONArray(s1).getJSONArray(s2);
						for(int x1=0; x1<aMatrix.length(); x1++)
							System.out.println(
									(x1==0 ? "{" : "") +
									Misc.formatVector_d(ModelComponent.readFromJSON_d(aMatrix.getJSONArray(x1)))
									+ (x1<aMatrix.length()-1 ? "," : "}"));
					}
				}
			}
		}
		
		if(oInt.has("spacingVector")) {
			JSONArray aSV = oInt.getJSONArray("spacingVector");
			for(int iExp=0; iExp<aSV.length(); iExp++)
				for(int deltaS=0; deltaS<2; deltaS++)
					System.out.println("Spacing vector(exp="+iExp+", delta strand="+deltaS+"): "+ Misc.formatVector_d(ModelComponent.readFromJSON_d(aSV.getJSONArray(iExp).getJSONArray(deltaS))));
		}
	}
	
	public void modifyBindingMode(JSONObject model, int iBM, int deltaLeft, int deltaRight, int deltaFlankLength, int newDIntMax, ArrayList<ModelComponent> componentList) {
		
		//Identifies the JSON output object. 
		JSONObject oIntCoeff    = model.getJSONObject("coefficients" ).getJSONArray("bindingModeInteractions").getJSONObject(iComp);

		//Determines how much should be padded for each binding mode.
		int dBefore1 = 0, dAfter1 = 0, dBefore2 = 0, dAfter2 = 0; 
		if(iBM==i0) {
			dBefore1 = -deltaLeft  + deltaFlankLength;
			dAfter1  = -deltaRight + deltaFlankLength;
		}
		
		if(iBM==i1) {
			dBefore2 = -deltaLeft  + deltaFlankLength;
			dAfter2  = -deltaRight + deltaFlankLength;
		}
		
		//TODO: POTENTIAL PROBLEM: If a homogeneous interaction is modified, the first row/column may become zero after modifyBindingMode. The 
		// packing may then put the entire matrix to zero. However, this should never happen since the interactions aren't activated and fitted until after
		// the binding modes are optimized.
		
		if(dBefore1!=0 || dAfter1!=0 || dBefore2!=0 || dAfter2!=0) {
			
			for(int iExp=0; iExp<interactionBetas.size(); iExp++) {
				for(int s1=0; s1<2; s1++)
					for(int s2=0; s2<2; s2++)
						oIntCoeff.getJSONArray("positionMatrix").getJSONArray(iExp).getJSONArray(s1).put(s2,JSONArrayConvert_dd(Array.padMatrix(interactionBetas.get(iExp).get(s1).get(s2), dBefore1, dAfter1,  dBefore2, dAfter2 )));

			}
		}
	}

}
