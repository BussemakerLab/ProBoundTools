package modelComponents;

import java.util.ArrayList;

import org.json.JSONObject;

import proBoundTools.Misc;


public class RhoGammaModel extends EnrichmentModel {
	
	private double[] rho, gamma;
	boolean  fitRho, fitGamma, roundSpecificRho, roundSpecificGamma;

	RhoGammaModel(JSONObject config, int iExpIn, 
			CountTable inTable,
			ArrayList<BindingMode> allBindingModes, 
			ArrayList<BindingModeInteraction> allInteractions) {
		
		super(config, iExpIn, inTable, allBindingModes, allInteractions);
		
		componentName = "Rho-gamma enrichment model "+iComp;
		if(verbose)
			System.out.println(">> Creating "+componentName+".");

		modelType = "RhoGamma";
		
		readFromJSON_settings(config);
		readFromJSON_constraints(config);
		
		allocateParameters();
		
		setFreezingLevel(0);
		
	}
	
	public void allocateParameters() {
		
		rho   = new double[nColumns];
		gamma = new double[nColumns];
		return;
		
	}
	
	@Override
	public void seed_component(JSONObject config) {

		rho   = null;
		gamma = null;
		String coefficientKey = "modelSeeding";
		if(config.has(coefficientKey)) {
			JSONObject oSeed = config.getJSONObject(coefficientKey); 
			if( oSeed.has(componentKey) ) {
				JSONObject oEnr = oSeed.getJSONArray(componentKey).getJSONObject(iComp);
				if(oEnr.has("rho"))
					rho   = readFromJSON_d(oEnr.getJSONArray("rho"));
				if(oEnr.has("gamma"))
					gamma = readFromJSON_d(oEnr.getJSONArray("gamma"));

			}
		}
		if(rho==null) {
			rho   = new double[nColumns];
			for(int i=0; i<nColumns; i++)
				rho[i]   = roundSpecificRho ? i : 1;
		}
		
		if(gamma==null) {
			gamma   = new double[nColumns];
			for(int i=0; i<nColumns; i++) {
				gamma[i] = 0;
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

			if(fitRho) {
				int[] rhoRange                     = range_d(rho,          iCurr+1);
				iLast = last_i(rhoRange);
				iCurr = iLast==null ? iCurr : iLast;
				if(!roundSpecificRho)
					for(int iExp=0; iExp<rhoRange.length; iExp++)
						rhoRange[iExp] = rhoRange[0];
				saveToJSON_rho_i(packing, coefficientKey, rhoRange);
			}
			
			if(fitGamma) {
				int[] gammaRange                   = range_d(gamma,          iCurr+1);
				iLast = last_i(gammaRange);
				iCurr = iLast==null ? iCurr : iLast;
				if(!roundSpecificGamma)
					for(int iExp=0; iExp<gammaRange.length; iExp++)
						gammaRange[iExp] = gammaRange[0];
				saveToJSON_gamma_i(packing, coefficientKey, gammaRange);
			}
		}

		return iCurr;
	}
	

	@Override
	public void addZeroJSON_component(JSONObject in, String coefficientKey) {
		
		addEmptyJSON_component_O(in, coefficientKey, componentKey, iComp);
		saveToJSON_rho_d(in,   coefficientKey, zero_d(rho));
		saveToJSON_gamma_d(in, coefficientKey, zero_d(gamma));
		
		return;
	}
	
	@Override
	void saveToJSON_settings(JSONObject out) {
		
		super.saveToJSON_settings(out);
		
		//TODO: Decide if we should simply have cumulative enrichment.
		//oEnr.put("cumulativeEnrichment",  cumulativeEnrichment);
		//oEnr.put("exponentialEnrichment", exponentialEnrichment);
		
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

		oEnr.put("fitRho",             fitRho);
		oEnr.put("fitGamma",           fitGamma);
		oEnr.put("roundSpecificRho",   roundSpecificRho);
		oEnr.put("roundSpecificGamma", roundSpecificGamma);
		
	}
	
	@Override
	void readFromJSON_constraints(JSONObject in) {
		
		String coefficientKey = "modelFittingConstraints";
		if(in.has(coefficientKey)) {
			JSONObject oConsEnr   = in.getJSONObject(coefficientKey).getJSONArray(componentKey).getJSONObject(iComp);
			fitRho                = oConsEnr.has("fitRho")             ? oConsEnr.getBoolean("fitRho")            : false;
			fitGamma              = oConsEnr.has("fitGamma")           ? oConsEnr.getBoolean("fitGamma")          : false;
			roundSpecificRho      = oConsEnr.has("roundSpecificRho")   ? oConsEnr.getBoolean("roundSpecificRho")  : true;
			roundSpecificGamma    = oConsEnr.has("roundSpecificGamma") ? oConsEnr.getBoolean("roundSpecificGamma"): true;			
		}
	}
	
	@Override
	void saveToJSON_parameters(JSONObject in, String coefficientKey) {
		
		addEmptyJSON_component_O(in, coefficientKey, componentKey, iComp);
		saveToJSON_rho_d(in,   coefficientKey, rho);
		saveToJSON_gamma_d(in, coefficientKey, gamma);
		
	}

	@Override
	void readFromJSON_parameters(JSONObject in, String coefficientKey) {
		
		JSONObject oEnr = in.getJSONObject(coefficientKey).getJSONArray(componentKey).getJSONObject(iComp);
		rho             = readFromJSON_d(oEnr.getJSONArray("rho"));
		gamma           = readFromJSON_d(oEnr.getJSONArray("gamma"));
		
	
	}
	
	void saveToJSON_rho_d(JSONObject out, String coefficientKey, double[] rho) {
		
		JSONObject oBM = out.getJSONObject(coefficientKey).getJSONArray(componentKey).getJSONObject(iComp);
		oBM.put("rho", JSONArrayConvert_d(rho));

		return;
	}
	
	void saveToJSON_gamma_d(JSONObject out, String coefficientKey, double[] gamma) {
		
		JSONObject oBM = out.getJSONObject(coefficientKey).getJSONArray(componentKey).getJSONObject(iComp);
		oBM.put("gamma", JSONArrayConvert_d(gamma));

		return;
	}
	
	void saveToJSON_rho_i(JSONObject out, String coefficientKey, int[] rho) {
		
		JSONObject oBM = out.getJSONObject(coefficientKey).getJSONArray(componentKey).getJSONObject(iComp);
		oBM.put("rho", JSONArrayConvert_i(rho));

		return;
	}
	
	void saveToJSON_gamma_i(JSONObject out, String coefficientKey, int[] gamma) {
		
		JSONObject oBM = out.getJSONObject(coefficientKey).getJSONArray(componentKey).getJSONObject(iComp);
		oBM.put("gamma", JSONArrayConvert_i(gamma));

		return;
	}

	public void updateDeltaKappaRI(double[] deltaKappaRI, double[] alphaRI, int nRounds) {
		
		for(int r=0; r<nRounds; r++) {					
			deltaKappaRI[r] = 1.0;
			deltaKappaRI[r] *= Math.pow(alphaRI[r],     rho[r]);
			deltaKappaRI[r] *= Math.pow(1.0+alphaRI[r], gamma[r]);
		}
	}
	
	public void updateNablaF1(double[] nablaF1, double[] nablaN, double[] alphaRI, int nRounds) {
		for(int r=0; r<nColumns; r++) {
			nablaF1[r] = nablaN[r] * ( rho[r]/alphaRI[r] + gamma[r]/(1+alphaRI[r]) );
		}
	}
	
	static public void printJSONObjectCoefficients(JSONObject model, String coefficientKey, int iExp) {

		JSONObject oEnrichment = model.getJSONObject(coefficientKey).getJSONArray("enrichmentModel").getJSONObject(iExp);
		System.out.println("rho:               "+Misc.formatVector_d(ModelComponent.readFromJSON_d(oEnrichment.getJSONArray("rho"))));
		System.out.println("gamma:             "+Misc.formatVector_d(ModelComponent.readFromJSON_d(oEnrichment.getJSONArray("gamma"))));
		
	}
}
