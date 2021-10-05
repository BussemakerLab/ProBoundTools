package modelComponents;

import java.util.ArrayList;

import org.json.JSONObject;

import proBoundTools.Misc;

public class ExponentialKineticsModel extends EnrichmentModel {

	// Selection = (e^delta/(1+e^delta))*e^{-alpha}+(1/(1+e^delta))*(1-e^{-alpha})
	private double[] delta;
	private boolean[] fitDelta;
	private boolean roundSpecificDelta;
	
	ExponentialKineticsModel(JSONObject config, int iExpIn, 
			CountTable inTable,
			ArrayList<BindingMode> allBindingModes, 
			ArrayList<BindingModeInteraction> allInteractions) {
		

		super(config, iExpIn, inTable, allBindingModes, allInteractions);

		
		componentName = "Exponential Kinetics enrichment model "+iComp;
		if(verbose)
			System.out.println(">> Creating "+componentName+".");
		
		readFromJSON_settings(config);
		readFromJSON_constraints(config);

		maxFreezeLevel = 0;
		for(boolean b: fitDelta)
			if(b) maxFreezeLevel=1;
		setFreezingLevel(maxFreezeLevel);
		allocateParameters();
	}
	
	public void allocateParameters() {
		
		delta = new double[nColumns];
		return;
		
	}
	
	@Override
	public void seed_component(JSONObject config) {
		
		delta                 = null;
		String coefficientKey = "modelSeeding";
		if(config.has(coefficientKey)) {
			JSONObject oSeed = config.getJSONObject(coefficientKey); 
			if( oSeed.has(componentKey) ) {
				JSONObject oEnr = oSeed.getJSONArray(componentKey).getJSONObject(iComp);
				if(oEnr.has("delta"))
					delta = readFromJSON_d(oEnr.getJSONArray("delta"));

			}
		}
		
		if(delta==null) {
			delta   = new double[nColumns];
			for(int i=0; i<nColumns; i++) {
				delta[i] = 0;
			}
		}
		
		return;
	}
	
	@Override
	int packModel_component(JSONObject packing, int iFirst) {
		
		String coefficientKey = "packing";
		addEmptyJSON_component_O(packing, coefficientKey, componentKey, iComp);
		
		int iCurr = iFirst;
		if(fitComponent&&freezingLevel==0) {
			int[] deltaRange = new int[nColumns];
			boolean valueAdded = false;
			int nFitted = 0;
			for(int iCol=0; iCol<nColumns; iCol++) {
				if(fitDelta[iCol]) {
					nFitted+=1;
					if(roundSpecificDelta||(!valueAdded)) {
						iCurr++;
						valueAdded = true;
					}
					deltaRange[iCol] = iCurr;
				} else {
					deltaRange[iCol] = -1;
				}
			}
			
			//If roundSpecificDelta=False, then either all or no detlas must be fitted
			if(!roundSpecificDelta && (nFitted>0&&nFitted<nColumns)) {
				throw new IllegalArgumentException("Either all or no detlas must be fitted when roundSpecificDelta=False for enrichment model "+iComp+", nFitted="+nFitted);
			}
			
			saveToJSON_delta_i(packing, coefficientKey, deltaRange);
		
		}
		
		
		return iCurr;
	}
	

	@Override
	public void addZeroJSON_component(JSONObject in, String coefficientKey) {
		
		addEmptyJSON_component_O(in, coefficientKey, componentKey, iComp);
		saveToJSON_delta_d(in,   coefficientKey, zero_d(delta));
		
		return;
	}
	
	@Override
	void saveToJSON_settings(JSONObject out) {
		
		super.saveToJSON_settings(out);
				
	}
	
	@Override
	void readFromJSON_settings(JSONObject in) {
		
		super.readFromJSON_settings(in);

	}
	
	
	@Override
	void saveToJSON_constraints(JSONObject out) {
		
		String coefficientKey = "modelFittingConstraints";
		addEmptyJSON_component_O(out, coefficientKey, componentKey, iComp);
		
		JSONObject oEnr        = out.getJSONObject(coefficientKey).getJSONArray(componentKey).getJSONObject(iComp);

		oEnr.put("fitDelta",           JSONArrayConvert_b(fitDelta));
		oEnr.put("roundSpecificDelta", roundSpecificDelta);
		
	}
	
	@Override
	void readFromJSON_constraints(JSONObject in) {
		
		String coefficientKey = "modelFittingConstraints";
		if(in.has(coefficientKey)) {
			JSONObject oConsEnr   = in.getJSONObject(coefficientKey).getJSONArray(componentKey).getJSONObject(iComp);
			
			if(oConsEnr.has("fitDelta")) {
				fitDelta = ModelComponent.readFromJSON_b(oConsEnr.getJSONArray("fitDelta"));
				if(fitDelta.length==1) {
					boolean val = fitDelta[0];
					fitDelta = new boolean[nColumns];
					for(int iCol=0; iCol<nColumns; iCol++)
						fitDelta[iCol] = val;
				}
			}
			else
				fitDelta = new boolean[nColumns];

			roundSpecificDelta = oConsEnr.has("roundSpecificDelta") ? oConsEnr.getBoolean("roundSpecificDelta") : true;			
		}

	}
	
	@Override
	void saveToJSON_parameters(JSONObject in, String coefficientKey) {
		
		addEmptyJSON_component_O(in, coefficientKey, componentKey, iComp);
		saveToJSON_delta_d(in,   coefficientKey, delta);
		
	}

	@Override
	void readFromJSON_parameters(JSONObject in, String coefficientKey) {
		
		JSONObject oEnr = in.getJSONObject(coefficientKey).getJSONArray(componentKey).getJSONObject(iComp);
		delta           = readFromJSON_d(oEnr.getJSONArray("delta"));
	}
	
	void saveToJSON_delta_d(JSONObject out, String coefficientKey, double[] delta) {
		
		JSONObject oBM = out.getJSONObject(coefficientKey).getJSONArray(componentKey).getJSONObject(iComp);
		oBM.put("delta", JSONArrayConvert_d(delta));

		return;
	}
	
	void saveToJSON_delta_i(JSONObject out, String coefficientKey, int[] delta) {
		
		JSONObject oBM = out.getJSONObject(coefficientKey).getJSONArray(componentKey).getJSONObject(iComp);
		oBM.put("delta", JSONArrayConvert_i(delta));

		return;
	}
	
	public void updateDeltaKappaRI(double[] deltaKappaRI, double[] alphaRI, int nRounds) {
		
		// Selection = (e^delta/(1+e^delta))*e^{-alpha}+(1/(1+e^delta))*(1-e^{-alpha})
		
		for(int r=0; r<nRounds; r++) {
			double expD  = Math.exp(delta[r]);
			double w1   = 1.0/(1.0+1.0/expD);
			double w2   = 1.0/(1.0+    expD);
			deltaKappaRI[r] = w1 * Math.exp(-alphaRI[r]) + w2 *  (1 - Math.exp(-alphaRI[r]));
		}
	}
	/// dSelection / Selection
	public void updateNablaF1(double[] nablaF1, double[] nablaN, double[] alphaRI, int nColumns) {
		for(int r=0; r<nColumns; r++) {
			double expD = Math.exp(delta[r]);
			double expA = Math.exp(alphaRI[r]);
			nablaF1[r]          = nablaN[r] * (1-expD) / (expD  + (expA-1.0));
		}
	}
	
	static public void printJSONObjectCoefficients(JSONObject model, String coefficientKey, int iExp) {

		JSONObject oEnrichment = model.getJSONObject(coefficientKey).getJSONArray("enrichmentModel").getJSONObject(iExp);
		System.out.println("delta:             "+Misc.formatVector_d(ModelComponent.readFromJSON_d(oEnrichment.getJSONArray("delta"))));
		
	}

}
