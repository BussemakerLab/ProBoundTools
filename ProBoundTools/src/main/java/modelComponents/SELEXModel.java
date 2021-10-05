package modelComponents;

import java.util.ArrayList;
import org.json.JSONObject;

public class SELEXModel extends EnrichmentModel{
	
	public boolean bindingSaturation;
	public boolean trySaturation;
	SELEXModel(JSONObject config, int iExpIn, 
			CountTable inTable,
			ArrayList<BindingMode> allBindingModes, 
			ArrayList<BindingModeInteraction> allInteractions) {
		
		super(config, iExpIn, inTable, allBindingModes, allInteractions);
		
		componentName = "SELEX enrichment model "+iComp;
		if(verbose)
			System.out.println(">> Creating "+componentName+".");
		
		bindingSaturation = false;
		trySaturation     = false;
		
		readFromJSON_settings(config);
		readFromJSON_constraints(config);

		allocateParameters();
		
		setFreezingLevel(0);

	}
	
	public void allocateParameters() {
		
		return;
		
	}
	
	@Override
	public void seed_component(JSONObject config) {

		return;
	}
	
	@Override
	int packModel_component(JSONObject packing, int iFirst) {
		
		String coefficientKey = "packing";
		addEmptyJSON_component_O(packing, coefficientKey, componentKey, iComp);
		int iCurr = iFirst;
		return iCurr;
	}
	

	@Override
	public void addZeroJSON_component(JSONObject in, String coefficientKey) {
		
		addEmptyJSON_component_O(in, coefficientKey, componentKey, iComp);
		
		return;
	}
	
	@Override
	void saveToJSON_settings(JSONObject out) {
		
		super.saveToJSON_settings(out);
		
		String coefficientKey = "modelSettings";
		JSONObject oEnr   =                  out.getJSONObject(coefficientKey).getJSONArray(componentKey).getJSONObject(iComp);
		oEnr.put("bindingSaturation",  bindingSaturation);

				
	}
	
	@Override
	void readFromJSON_settings(JSONObject in) {
		
		String coefficientKey = "modelSettings";
		super.readFromJSON_settings(in);
		JSONObject oSettEnr = in.getJSONObject(coefficientKey).getJSONArray(componentKey).getJSONObject(iComp);
		bindingSaturation     = oSettEnr.getBoolean("bindingSaturation");

	}
	
	
	@Override
	void saveToJSON_constraints(JSONObject out) {
		
		String coefficientKey = "modelFittingConstraints";
		addEmptyJSON_component_O(out, coefficientKey, componentKey, iComp);
		
		JSONObject oEnr       = out.getJSONObject(coefficientKey).getJSONArray(componentKey).getJSONObject(iComp);

		oEnr.put("trySaturation",             trySaturation);
		
	}
	
	@Override
	void readFromJSON_constraints(JSONObject in) {
		
		String coefficientKey = "modelFittingConstraints";
		if(in.has(coefficientKey)) {
			JSONObject oConsEnr   = in.getJSONObject(coefficientKey).getJSONArray(componentKey).getJSONObject(iComp);
			trySaturation          = oConsEnr.has("trySaturation") ? oConsEnr.getBoolean("trySaturation") : false;		
		}
	}
	
	@Override
	void saveToJSON_parameters(JSONObject in, String coefficientKey) {
		
		addEmptyJSON_component_O(in, coefficientKey, componentKey, iComp);
		
	}

	@Override
	void readFromJSON_parameters(JSONObject in, String coefficientKey) {
		
	}
	
	public void updateDeltaKappaRI(double[] deltaKappaRI, double[] alphaRI, int nRounds) {
		for(int r=0; r<nRounds; r++) {					
			deltaKappaRI[r] = alphaRI[r];
			if(bindingSaturation)
				deltaKappaRI[r] /= (1+alphaRI[r]);
		}
	}
	
	public void updateNablaF1(double[] nablaF1, double[] nablaN, double[] alphaRI, int nColumns) {
		for(int r=0; r<nColumns; r++) {
			nablaF1[r]      = nablaN[r] * (  1.0/alphaRI[r]     );
			if(bindingSaturation)
				nablaF1[r] += nablaN[r] * ( -1.0/(1+alphaRI[r]) );		
		}
	}
	
}
