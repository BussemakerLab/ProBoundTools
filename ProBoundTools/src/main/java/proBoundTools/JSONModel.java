package proBoundTools;

import java.io.*;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Set;

import org.everit.json.schema.Schema;
import org.everit.json.schema.ValidationException;
import org.everit.json.schema.loader.SchemaLoader;

import org.json.*;

import modelComponents.BindingMode;
import modelComponents.BindingModeInteraction;
import modelComponents.CountTable;
import modelComponents.ExponentialKineticsModel;
import modelComponents.RhoGammaModel;


public class JSONModel {

	public JSONObject oModel;
	
	String jsonSchemaPath;
	
	//Constructor reading whole file
	public JSONModel(String jsonPath, String jsonSchemaPathIn) {
		
		// Loads JSON model.
    	oModel = loadJSONFile(jsonPath);
    	
    	// Validates the basic JSON Schema
    	jsonSchemaPath = jsonSchemaPathIn;
    	
    	validateModelSchema();
    	
    	//TODO: Note that we need to validate the model twice to fully fill in empty objects.
    	validateModelSchema();

    	//Validation of binding modes.
    	validateModelSettings(oModel);

		
	}
	
	//Constructor reading a single line.
	public JSONModel(String jsonPath, String jsonSchemaPathIn, int iLine) {
		
		// Loads JSON model.
    	oModel = loadJSONLine(jsonPath, iLine);
    	
    	// Validates the basic JSON Schema
    	jsonSchemaPath = jsonSchemaPathIn;
    	
    	validateModelSchema();
    	
    	//TODO: Note that we need to validate the model twice to fully fill in empty objects.
    	validateModelSchema();

    	//Validation of binding modes.
    	validateModelSettings(oModel);

	}
	
	void validateModelSchema() {
    	validateSchemaFile_O(jsonSchemaPath, oModel);
	}
	
	/*public void validateConfigruationFile(String settingsSchemaPath) {
		validateSchemaFile(settingsSchemaPath, oModel);
	}*/
	
	public void writeModel(String outPath) {

		try {
    		FileWriter file = new FileWriter(outPath);
    		oModel.write(file);
    		file.close();
    	} catch (java.lang.Exception e) {
        	throw new java.lang.RuntimeException("Can't write to JSON vile '"+outPath+"'.");
    	}
    	    	
	}
		
	/////////////////////////////
	// CHECKING modelSettings //  
	////////////////////////////
	static void validateModelSettings(JSONObject oModel) {
		
		JSONObject oModelSettings       = oModel.getJSONObject("modelSettings");
		JSONArray aBindingModeSettings  = oModelSettings.getJSONArray("bindingModes");
		JSONArray aInteractionsSettings = oModelSettings.getJSONArray("bindingModeInteractions");

		int nBindingModes = aBindingModeSettings.length();
		int nInteractions = aInteractionsSettings.length();

		// Checking binding modes 
		////////////////////////////
		
		// Should be OK by default.

		// Checking bindignModeInteractions
		////////////////////////////////////
		for(int iInt=0; iInt<nInteractions; iInt++){
			///Checking so the interacting binding modes exist. 
			JSONObject oInt = aInteractionsSettings.getJSONObject(iInt);
			JSONArray aModes = oInt.getJSONArray("bindingModes");
			int iBM0=aModes.getInt(0), iBM1=aModes.getInt(1);
			if( Math.max(iBM0, iBM1) >= nBindingModes ) 
	        	throw new java.lang.RuntimeException("Interactions defined for non-existent binding modes ("+iBM0+","+iBM1+").");
		}

	}
	
	public static void validateConfig(String jsonSchemaPath, String jsonSettingsSchemaPath, JSONObject oModel) {
		
		
		//Checks so all entries are of correct type.
		validateSchemaFile_O(jsonSchemaPath, oModel);
		validateSchemaFile_O(jsonSettingsSchemaPath, oModel);
		
		if(!oModel.has("modelSettings"))
			throw new java.lang.IllegalArgumentException("ERROR: The config file does not have a 'modelSettings' object.");
		JSONObject oModelSettings       = oModel.getJSONObject("modelSettings");
		
		//Makes sure the arrays 'bindingModes', 'countTable', 'bindingModeInteractions', and 'enrichmentModel' exist.
		if(! (oModelSettings.has("bindingModes")))
			throw new java.lang.IllegalArgumentException("ERROR: The config file does not have an 'modelSettings.bindingModes' array.");
		if(!oModelSettings.has("countTable"))
			throw new java.lang.IllegalArgumentException("ERROR: The config file does not have an 'modelSettings.countTable' array.");
		if(!oModelSettings.has("enrichmentModel"))
			oModelSettings.put("enrichmentModel", new JSONArray());
		if(!oModelSettings.has("bindingModeInteractions"))
			oModelSettings.put("bindingModeInteractions", new JSONArray());
		
		JSONArray aBindingModeSettings  = oModelSettings.getJSONArray("bindingModes");
		JSONArray aInteractionsSettings = oModelSettings.getJSONArray("bindingModeInteractions");
		JSONArray aEnrichmentSettings   = oModelSettings.getJSONArray("enrichmentModel");
		JSONArray aCountTableSettings   = oModelSettings.getJSONArray("countTable");
		
		int nExperiments  = aCountTableSettings.length();
		int nBindingModes = aBindingModeSettings.length();
		int nInteractions = aInteractionsSettings.length();

		//Makes sure we have the correct number of enrichment models.
		while(aEnrichmentSettings.length() < nExperiments)
			aEnrichmentSettings.put(new JSONObject());
		if(aEnrichmentSettings.length()>nExperiments)
			throw new java.lang.IllegalArgumentException("ERROR: The array modelSettings.enrichmentModel has more entries than modelSettings.countTable.");
		
		//Makes sure we have the correct number of modelFittingConstraints
		if(!oModel.has("modelFittingConstraints")) oModel.put("modelFittingConstraints", new JSONObject());
		JSONObject oModelConstr = oModel.getJSONObject("modelFittingConstraints");
		System.out.println(oModelConstr.toString());
		//Checks so bindingModes, bindingModeInteractions, countTables, and enrichmentModel have the correct number of entries.
		String[] entries = new String[4];
		entries[0] = "bindingModes"; entries[1] = "bindingModeInteractions"; entries[2] = "countTable"; entries[3] = "enrichmentModel";
		for(int i=0; i<entries.length; i++) {
			
			//Adds new entries if necessary.
			if(!oModelConstr.has(entries[i])) oModelConstr.put(entries[i], new JSONArray());
			JSONArray aEntry   = oModelConstr.getJSONArray(entries[i]);
			JSONArray aSetting = oModelSettings.getJSONArray(entries[i]);
			System.out.println("Entry="+entries[i]+", aEntry="+aEntry.toString());
			while(aEntry.length() < aSetting.length())
				aEntry.put(new JSONObject());
			
			//Checks so there are not more entries in modelFittingConstraints than in modelSettings
			if(aEntry.length() > aSetting.length())
				throw new java.lang.IllegalArgumentException("ERROR: modelFittingConstraints."+entries[i]+" has more entries ("+aEntry.length()+") than modelSettings."+entries[i]+" ("+aSetting.length()+").");
			
		}
				
		//Runs schema validation to fill in default values.  
		validateSchemaFile_O(jsonSchemaPath, oModel);
		validateSchemaFile_O(jsonSchemaPath, oModel);
		
		// Checks countTable
		////////////////////

		//Reads the number of columns
		ArrayList<Integer> nColumns = new ArrayList<Integer>();
		for(int iCnt=0; iCnt<aCountTableSettings.length(); iCnt++) {
			int nColTemp = aCountTableSettings.getJSONObject(iCnt).getInt("nColumns");
			nColumns.add(nColTemp);
			//Checks so all modeled columns exist
			if(aCountTableSettings.getJSONObject(iCnt).has("modeledColumns")) {
				JSONArray modCols = aCountTableSettings.getJSONObject(iCnt).getJSONArray("modeledColumns");
				if(modCols.length()==1) {
					if(modCols.getInt(0)!=-1)
						throw new java.lang.IllegalArgumentException("ERROR: Count table "+iCnt+" only has one modeled columns specified  \"modeledColumns\": "+modCols.getInt(0)+".");
				} else {
					for(int i=0; i<modCols.length(); i++)
						if(modCols.getInt(i)<0 || modCols.getInt(i)>nColTemp)
							throw new java.lang.IllegalArgumentException("ERROR: Count table "+iCnt+" has "+nColTemp+" columns, but \"modeledColumns\" indicates that column i="+modCols.getInt(i)+" should be modeled.");
				}
			}
		}

		// Checks enrichmentModel
		/////////////////////////
		for(int iExp=0; iExp<aEnrichmentSettings.length(); iExp++) {
			JSONObject oEnr = aEnrichmentSettings.getJSONObject(iExp);
			JSONArray oBM = oEnr.getJSONArray("bindingModes");
			JSONObject oCnt = aCountTableSettings.getJSONObject(iExp);
			 
			//If bindingModes=[-1], then build new array containing all binding modes.
			if(oBM.length()==1 && oBM.getInt(0)==-1) {
				oBM = new JSONArray();
				for(int iBM=0; iBM<nBindingModes; iBM++)
					oBM.put(iBM);
			}
				
			//Checks so all binding modes are valid
			for(int iBM=0; iBM<oBM.length(); iBM++) {

				//Checks so the binding mode exists.
				if(oBM.getInt(iBM) >= aBindingModeSettings.length())
					throw new java.lang.IllegalArgumentException("ERROR: For enrichment models "+iExp+", "+iBM+" is not a valid binding mode.");

				//Checks so the binding mode flank lengths used in this enrichment model not larger than flanks sequences in the count table.
				JSONObject oBMSet = aBindingModeSettings.getJSONObject(oBM.getInt(iBM));
				//Computes the maximum allowed length.
				int maxFlankLength         = Math.min(
						oCnt.getString("rightFlank").length(),  
						oCnt.getString("leftFlank").length());
				if(oBMSet.getInt("flankLength") > maxFlankLength)
					throw new java.lang.RuntimeException("Flank length of binding mode "+ iBM +" exceeds the length of the specified flanking sequence in dataset "+iExp+".");


			}
			
			JSONArray oInt = oEnr.getJSONArray("bindingModeInteractions");
			
			//If bindingModeInteractions=[-1], then build new array containing all interactions.
			if(oInt.length()==1 && oInt.getInt(0)==-1) {
				oInt = new JSONArray();
				for(int iInt=0; iInt<nInteractions; iInt++)
					oInt.put(iInt);
			}
			

			//Checks so all binding mode interactions are valid
			for(int iInt=0; iInt<oInt.length(); iInt++) {
				
				//Checks so the interaction exists
				if(oInt.getInt(iInt) >= aInteractionsSettings.length())
					throw new java.lang.IllegalArgumentException("ERROR: For enrichment models "+iExp+", "+iInt+" is not a valid binding mode interaction.");
				JSONObject oIntSet = aInteractionsSettings.getJSONObject(oInt.getInt(iInt));
				
				//Checks to the respective binding modes are included. 
				boolean hasBM0=false, hasBM1=false;
				for(int iBM=0; iBM<oBM.length(); iBM++) {
					if(oBM.getInt(iBM) == oIntSet.getJSONArray("bindingModes").getInt(0))
						hasBM0 = true;
					if(oBM.getInt(iBM) == oIntSet.getJSONArray("bindingModes").getInt(1))
						hasBM1 = true;
				}
				if(!hasBM0 || !hasBM1)
					throw new java.lang.IllegalArgumentException("ERROR: For enrichment models, binding mode interaction "+iInt+" is not a valid since the corresponding binding modes are not included in the enrichment model.");
			}

		}

		//Checks binding modes
		//////////////////////
		//Checks so the binding motif size and shift are not optimized if a symmetry string is specified.
		for(int iBM=0; iBM<oModelConstr.getJSONArray("bindingModes").length(); iBM++) {
			JSONObject oBM = oModelConstr.getJSONArray("bindingModes").getJSONObject(iBM);
			String symmetryString = oBM.getString("symmetryString");
			if(!symmetryString.equals("null")) {
				if(oBM.getBoolean("optimizeSize"))
					throw new java.lang.IllegalArgumentException("For binding mode "+iBM+": cannot have optimizeSize=true when symmetryString!=\"null\".");
				if(oBM.getBoolean("optimizeMotifShift"))
					throw new java.lang.IllegalArgumentException("For binding mode "+iBM+": cannot have optimizeMotifShift=true when symmetryString!=\"null\".");
				if(oBM.getBoolean("optimizeMotifShiftHeuristic"))
					throw new java.lang.IllegalArgumentException("For binding mode "+iBM+": cannot have optimizeMotifShiftHeuristic=true when symmetryString!=\"null\".");
			}
		}
		
		//Checks binding mode interactions
		//////////////////////////////////
		
		//Checks so the interactions are between valid binding modes.
		for(int iInt=0; iInt<aInteractionsSettings.length(); iInt++) {
			JSONArray aBM = aInteractionsSettings.getJSONObject(iInt).getJSONArray("bindingModes");
			int iBM0 = aBM.getInt(0), iBM1 = aBM.getInt(1);
			if( iBM0>=nBindingModes || iBM1>=nBindingModes )
				throw new java.lang.IllegalArgumentException("For binding mode interaction "+iInt+", one of the binding modes ("+iBM0+","+iBM1+") does not exist.");
		}
				
		
	}

    static void validateModelCoefficients(JSONObject oModel, String coefficientKey) {
    	
    	JSONObject coefficients         = oModel.getJSONObject(coefficientKey);
		JSONObject oModelSettings       = oModel.getJSONObject("modelSettings");
    
		//Validates the basics of the settings object, fills in default values.
		validateModelSettings(oModel);
    	
		////////////////////////////////
		// CHECKING modelCoefficients //  
		////////////////////////////////
    	
    	Integer nExperiments        = null;  //Number of experiments
    	boolean hasActivity         = false; //True if there are activities
    	ArrayList<Integer> nColumns = null;  //Number of requred ativities for each binding modes

    	
    	
		// CHECKING countTable   
		/////////////////////////
    	if(oModelSettings.has("countTable")) {
    		if(!coefficients.has("countTable"))
    			throw new java.lang.RuntimeException("ERROR: modelSettings countTable but "+coefficientKey+" does not.");
    		
    		JSONArray aCTSett = oModelSettings.getJSONArray("countTable");
    		JSONArray aCTCoef = coefficients.getJSONArray("countTable");
    	    		
    		if( aCTSett.length() != aCTCoef.length())
    			throw new java.lang.RuntimeException("ERROR: The number of experiments do not match in modelSettings and in "+coefficientKey+".");
    		nExperiments=aCTSett.length();
    		
    		//Saves the number of columns
    		nColumns = new ArrayList<Integer>();
    		for(int iExp=0; iExp<nExperiments; iExp++ ) {
    			JSONArray aEta = aCTCoef.getJSONObject(iExp).getJSONArray("h");
    			nColumns.add(aEta.length());
    		}
    	}
    	
		// CHECKING enrichmentModel   
		/////////////////////////
		JSONArray aBindingModeSettings  = oModelSettings.getJSONArray("bindingModes");
		JSONArray aBindingModesCoef     = coefficients.getJSONArray("bindingModes");
		int nBindingModes               = aBindingModeSettings.length();		
		
    	JSONArray aInteractionsSettings = oModelSettings.getJSONArray("bindingModeInteractions");
    	JSONArray aInteractionCoef      = coefficients.getJSONArray("bindingModeInteractions");
		int nInteractions               = aInteractionsSettings.length();
		
    	ArrayList<Set<String>> bmModificationSets  = null;
    	
    	if(oModelSettings.has("enrichmentModel")) {
    		
    		JSONArray aEnrichment = oModelSettings.getJSONArray("enrichmentModel");
    		
			if(nExperiments!=null && aEnrichment.length() != nExperiments)
	        	throw new java.lang.RuntimeException(
	        			"The number of experiments in the enrichmentModel ("+aEnrichment.length()+") does not "
	        					+ "match then number of experiments ("+nExperiments+").");

			bmModificationSets  = new ArrayList<Set<String>>();
        	for(int iBM=0; iBM<nBindingModes; iBM++)
        		bmModificationSets.add(new HashSet<String>());

    		for(int iEnr=0; iEnr<aEnrichment.length(); iEnr++) {
    			
    			JSONObject oEnr = aEnrichment.getJSONObject(iEnr);
    			JSONArray aModifications = oEnr.getJSONArray("modifications");
    			JSONArray aBM            = oEnr.getJSONArray("bindingModes");
    			JSONArray aInt           = oEnr.getJSONArray("bindingModeInteractions");
    			
    			//Checks if the binding modes in the enrichment models is valid, and records the modifications used for each binding mode.  
				for(int iBM=0; iBM<aBM.length(); iBM++) {
					int iBMset = aBM.getInt(iBM);
					for(int iMod=0; iMod<aModifications.length(); iMod++) {
    					if(iBMset >= nBindingModes)
    						throw new java.lang.RuntimeException("ERROR: A binding mode in enrichment model "+iEnr+" does not exist.");
    					bmModificationSets.get(iBMset).add(aModifications.getString(iMod));
    				}
    			}
				
				//Checks so all interactions in the enrichment model exist.
				for(int iInt=0; iInt<aInt.length(); iInt++)
					if(aInt.getInt(iInt) >= nInteractions)
						throw new java.lang.RuntimeException("ERROR: A binding mode interactions in enrichment model "+iEnr+" does not exist.");
				
				//If the enrichment model is of SELEX type, checks the length of rho and gamma:
				if(oEnr.getString("modelType").equals("RhoGamma") && oModelSettings.has("countTable")) {
					if(oEnr.has("rho")) {
						if( oEnr.getJSONArray("rho").length() != nExperiments )
							throw new java.lang.RuntimeException("ERROR: The rho vector has the incorrect number of components.");
					} else {
						throw new java.lang.RuntimeException("ERROR: There are no rho vector specified.");
					}
					
					if(oEnr.has("gamma")) {
						if( oEnr.getJSONArray("gamma").length() != nExperiments )
							throw new java.lang.RuntimeException("ERROR: The gamma vector has the incorrect number of components.");
					} else {
						throw new java.lang.RuntimeException("ERROR: There are no gamma vector specified.");
					}
				}
    		}
    	}


		// CHECKING bindingModes   
		/////////////////////////
		
		//TODO: We learn the number of experiments from the first binding mode and make sure that everything else is consistent with this.
    	
		if( aBindingModesCoef.length() != nBindingModes )
        	throw new java.lang.RuntimeException("The number of binding modes ("+nBindingModes+") does not match the number binding mode coefficients ("+aBindingModesCoef.length()+").");

		
		//Checks the individual binding modes.
		for(int iBM=0; iBM<nBindingModes; iBM++){
			JSONObject oBM = aBindingModesCoef.getJSONObject(iBM);
			int size       = aBindingModeSettings.getJSONObject(iBM).getInt("size");
			int diIntMax   = aBindingModeSettings.getJSONObject(iBM).getInt("dinucleotideDistance");

			// CHECKING mononucleotide
			///////////////////////////
			JSONArray aMonoBeta = oBM.getJSONArray("mononucleotide");
			if( aMonoBeta.length() != size*4 )
	        	throw new java.lang.RuntimeException("Binding mode "+iBM+" has " + aMonoBeta.length() +" but "+size*4+" coefficeints were expected.");

			// CHECKING dinucleotide
			/////////////////////////
			JSONArray aDiBetaAll = oBM.getJSONArray("dinucleotide");
			if( aDiBetaAll.length() != diIntMax ) 
	        	throw new java.lang.RuntimeException("Binding mode "+iBM+" has " + aDiBetaAll.length() +" dinucleotide spacings but "+diIntMax+" was/were expected.");

			for(int iSpacing=0; iSpacing<aDiBetaAll.length(); iSpacing++) {
				JSONArray aDiBeta = aDiBetaAll.getJSONArray(iSpacing);
				int expectedBetas = (size-1-iSpacing)*16;
				if( aDiBeta.length() != expectedBetas )
		        	throw new java.lang.RuntimeException("Binding mode "+iBM+" has " + aDiBeta.length() +" dinculeotide coefficients for spacing " + (iSpacing+1) + " but "+ expectedBetas +" coefficeints were expected.");
			}

			// CHECKING activity
			/////////////////////
			if(iBM==0 && oBM.has("activity"))
				hasActivity = true;
			if(hasActivity && !oBM.has("activity"))
				throw new java.lang.RuntimeException("When parsing coeffficients, either all binding modes must have activities, or none.");
				
			if(hasActivity) {
				JSONArray aActivity = oBM.getJSONArray("activity");
				if( nExperiments == null) {
					nExperiments = aActivity.length();
				} else {
					if( nExperiments != aActivity.length())
			        	throw new java.lang.RuntimeException("There are "+nExperiments+" experiments, but only "+aActivity.length()+" binding mode activity vectors.");
				}
				
				if( nColumns == null) {
					nColumns = new ArrayList<Integer>();
					for(int iExp=0; iExp<nExperiments; iExp++)
						nColumns.add(aActivity.getJSONArray(iExp).length());
					
				} else {
					for(int iExp=0; iExp<nExperiments; iExp++) {
						if( aActivity.getJSONArray(iExp).length() != nColumns.get(iExp) )
				        	throw new java.lang.RuntimeException("Binding mode "+iBM+" has " + aActivity.getJSONArray(iExp).length() + 
				        			" activities in experiment "+iExp+" but should have "+nColumns.get(iExp)+".");
					}
					
				}
			}

				
			// CHECKING modifications (if these are defined in the enrichmentModel
			//////////////////////////////////////////////////////////////////////
			if( bmModificationSets != null ){
				JSONArray aModifications          = oBM.getJSONArray("modifications");
				JSONArray aModificationSetting    = aBindingModeSettings.getJSONObject(iBM).getJSONArray("modifications");
				int nMods                         = aModifications.length();
				Set<String> requiredModifications = new HashSet<String>(bmModificationSets.get(iBM));

				for(int iMod=0; iMod<nMods; iMod++) {
					JSONObject oMod = aModifications.getJSONObject(iMod);
					String modName = oMod.getString("name");

					//Identifying the corresponding modification in settings.
					int iBmMatch = -1;
					for(int iModSet=0; iModSet<nMods; iModSet++) {
						String modNameSetting = aModificationSetting.getJSONObject(iModSet).getString("name");
						if(!modNameSetting.equals(modName)) {
							continue;
						} else {
							iBmMatch = iModSet;
							requiredModifications.remove(modNameSetting);
							break;
						}
					}
					if(iBmMatch == -1)
						throw new java.lang.RuntimeException("Binding mode "+iBM+" has coefficients for the modification '" 
								+ modName + "' but this is not defined in the modelSettings.");
				}

				if(!requiredModifications.isEmpty())
					throw new java.lang.RuntimeException("Binding mode "+iBM+" should have coefficients for the modification '" 
							+ requiredModifications + "' but these are not defined.");
			}


			// CHECKING positionBias
			/////////////////////////
			boolean bmUsePositionBias = aBindingModeSettings.getJSONObject(iBM).getBoolean("positionBias");
			JSONArray aPositionBetas = oBM.getJSONArray("positionBias");
			if( bmUsePositionBias && aPositionBetas.length()==0 )
				throw new java.lang.RuntimeException("Binding mode "+iBM+" should have position bias coefficients but does not.");
			if( !bmUsePositionBias && aPositionBetas.length()>0 )
	        	throw new java.lang.RuntimeException("Binding mode "+iBM+" should not have position bias but it has coefficients.");
		}

		// CHECKING bindingModeInteractions   
		///////////////////////////////////
		if( aInteractionCoef.length() != nInteractions )
        	throw new java.lang.RuntimeException("The number of binding mode interactions ("+nInteractions+") does not match the number binding mode interaction coefficients ("+aInteractionCoef.length()+").");
		
		//Checks the individual interactions
		for(int iInt=0; iInt<nInteractions; iInt++){
			
			JSONObject oInt = aInteractionCoef.getJSONObject(iInt);
			// CHECKING interaction activity
			////////////////////////////////
			if(hasActivity) {
				
				if(!oInt.has("activity"))
					throw new java.lang.RuntimeException("Interaction "+iInt+" does not have an activity defined.");
							
				JSONArray aActivity = oInt.getJSONArray("activity");

				if( nExperiments != aActivity.length())
					throw new java.lang.RuntimeException("There are "+nExperiments+" experiments, but only "+aActivity.length()+" binding mode activity vectors.");

				for(int iExp=0; iExp<nExperiments; iExp++) {
					if( aActivity.getJSONArray(iExp).length() != nColumns.get(iExp) && 
							aActivity.getJSONArray(iExp).length() != 1)
						throw new java.lang.RuntimeException("Binding mode interaction has freeActivity=true "+iInt+" and " + aActivity.getJSONArray(iExp).length() + 
								" activities in experiment "+iExp+" but should have either 1 or "+nColumns.get(iExp)+".");

				} 
				
			}

			// CHECKING positionBias/interactions
			/////////////////////////////////////

			boolean usePositionBias = aInteractionsSettings.getJSONObject(iInt).getBoolean("positionBias");
			if(usePositionBias) {
				JSONArray aPositionMatrix = aInteractionCoef.getJSONObject(iInt).getJSONArray("positionMatrix");
				if( aPositionMatrix.length() == 0)
					throw new java.lang.RuntimeException("Binding mode interaction "+iInt+" has positionBias=true the matrix 'positionMatrix' is empty.");
			} else {
				JSONArray aSpacingVector = aInteractionCoef.getJSONObject(iInt).getJSONArray("spacingVector");
				if( aSpacingVector.length() == 0)
					throw new java.lang.RuntimeException("Binding mode interaction "+iInt+" has positionBias=false the matrix 'spacingVector' is empty.");
			}
		}
		
    }
    
    public static void validateSchemaFile_O(String schemaPath, JSONObject jsonObject) {

        try {
        	//Loads general schema.
        	JSONObject rawSchema = new JSONObject(new JSONTokener(new FileInputStream(new File(schemaPath))));
        	Schema schema        = SchemaLoader.builder().useDefaults(true).schemaJson(rawSchema).build().load().build();
        	schema.validate(jsonObject);
        	
        } catch (ValidationException e) {
        	        	
        	List<String> ev= e.getAllMessages();
        	for(int i=0; i<ev.size(); i++)
        		System.err.println(ev.get(i));
        	
        	throw new java.lang.RuntimeException("Error when validating JSON Schema in file  '" + schemaPath + "':\n"+e.getMessage() );
        } catch (org.json.JSONException e) {
        	System.err.println("printStackTrace:");
        	e.printStackTrace();
        	throw new java.lang.RuntimeException("JSON Exception: '"+schemaPath+"'.");
        } catch (java.lang.Exception e) {
        	System.err.println(e);
        	e.printStackTrace();
        	throw new java.lang.RuntimeException("Couldn't open JSON Schema file '"+schemaPath+"'.");
		}
    }
    
    
    public static void validateSchemaFile_A(String schemaPath, JSONArray jsonArray) {

        try {
        	
        	//Loads general schema.
        	JSONObject rawSchema = new JSONObject(new JSONTokener(new FileInputStream(new File(schemaPath))));
        	Schema schema        = SchemaLoader.builder().useDefaults(true).schemaJson(rawSchema).build().load().build();
        	schema.validate(jsonArray);
        	
        } catch (ValidationException e) {
        	        	
        	List<String> ev= e.getAllMessages();
        	for(int i=0; i<ev.size(); i++)
        		System.err.println(ev.get(i));
        	
        	throw new java.lang.RuntimeException("Error when validating JSON Schema in file  '" + schemaPath + "':\n"+e.getMessage() );
        } catch (org.json.JSONException e) {
        	System.err.println("printStackTrace:");
        	e.printStackTrace();
        	throw new java.lang.RuntimeException("JSON Exception: '"+schemaPath+"'.");
        } catch (java.lang.Exception e) {
        	System.err.println(e);
        	e.printStackTrace();
        	throw new java.lang.RuntimeException("Couldn't open JSON Schema file '"+schemaPath+"'.");
		}
    }
    
    public static JSONObject loadJSONFile(String jsonPath) {

    	JSONObject parsedJSONObject = null;

    	try {

    		File jsonFile = new File(jsonPath);
    		InputStream is = new FileInputStream(jsonFile);
    		BufferedReader streamReader = new BufferedReader(new InputStreamReader(is, "UTF-8"));
    		StringBuilder responseStrBuilder = new StringBuilder();

    		String inputStr;
    		while ((inputStr = streamReader.readLine()) != null)
    			responseStrBuilder.append(inputStr);

    		parsedJSONObject = new JSONObject(responseStrBuilder.toString());

    		streamReader.close();

    	} catch (JSONException e) {
    		throw new java.lang.RuntimeException("JSONException reading input JSON file '"+jsonPath+"': "+e.getMessage(), e);

    	} catch (java.lang.Exception e) {
    		throw new java.lang.RuntimeException("Problem opening and reading the file '"+jsonPath+"'.", e);
    	}

    	return parsedJSONObject;
    }
    
    public static JSONObject loadJSONLine(String jsonPath) {
    	return loadJSONLine(jsonPath, -1);
    }
    
    public static JSONObject loadJSONLine(String jsonPath, int iLine) {
    	
    	JSONObject parsedJSONObject = null;
    	
    	try {
    		
        	File jsonFile = new File(jsonPath);
        	InputStream is = new FileInputStream(jsonFile);
            BufferedReader streamReader = new BufferedReader(new InputStreamReader(is, "UTF-8"));
            //StringBuilder responseStrBuilder = new StringBuilder();
            
            String line = null;
            int iCurrent = 0;
            Double lBest = null;
            while ((line = streamReader.readLine()) != null) {
            	
            	if(iLine<0) {
            		//Keeps the line with the lowest -log(likelihood)
            		JSONObject tempObject = new JSONObject(line);
            		if(tempObject.has("metadata")&&tempObject.getJSONObject("metadata").has("logLikelihood")) {
            			//Only keep the best model if  metadata.logLikelihood is defined.
            			double l = tempObject.getJSONObject("metadata").getDouble("logLikelihood");
            			if(lBest==null || l<lBest) {
            				lBest = l;
            				parsedJSONObject = tempObject;
            			} 
            		} else {
            			//Keeps the last line if metadata.logLikelihood is defined.
            			parsedJSONObject = tempObject;
            		}
            	} else {
                	if(iCurrent==iLine) {
                		parsedJSONObject = new JSONObject(line);
                		break;
                	} else {
                		iCurrent++;
                	}            		
            	}
            	

            }      
            streamReader.close();
            
            if(parsedJSONObject==null) {
            	throw new java.lang.RuntimeException("End of file was reached before line "+iLine);
            }
            
            
            
    	} catch (JSONException e) {
    		throw new java.lang.RuntimeException("JSONException reading input JSON file '"+jsonPath+"': "+e.getMessage(), e);
    		
    	} catch (java.lang.Exception e) {
        	throw new java.lang.RuntimeException("Problem opening and reading the file '"+jsonPath+"'.", e);
		}
        
        return parsedJSONObject;
    }
    
    static public void printJSONState(JSONObject fullModel, String coefficients) {
    	JSONObject oCoeff = fullModel.getJSONObject(coefficients);
    	
    	if(oCoeff.has("countTable")) {
    		JSONArray aCount = oCoeff.getJSONArray("countTable");
        	
        	for(int iExp=0; iExp<aCount.length(); iExp++) {
        		System.out.println("");
        		System.out.println("Count Table "+iExp+":");
        		System.out.println("---------------");
        		CountTable.printJSONObjectCoefficients(fullModel, coefficients, iExp);

        	}
    	}
    	
    	if(oCoeff.has("enrichmentModel")) {
    		JSONArray aEnr = oCoeff.getJSONArray("enrichmentModel");
        	
        	for(int iExp=0; iExp<aEnr.length(); iExp++) {
        		String modelType = fullModel.getJSONObject("modelSettings").getJSONArray("enrichmentModel").getJSONObject(iExp).getString("modelType");
        		if(modelType.equals("RhoGamma") || modelType.equals("ExponentialKinetics")) {
        			System.out.println("");
        			System.out.println("Enrichment model"+iExp+":");
        			System.out.println("--------------------");
        			if(modelType.equals("RhoGamma"))
        				RhoGammaModel.printJSONObjectCoefficients(fullModel, coefficients, iExp);
        			else
            			ExponentialKineticsModel.printJSONObjectCoefficients(fullModel, coefficients, iExp);

        		}
        	}
    	}
    	
    	
    	if(oCoeff.has("bindingModes")) {
    		JSONArray aBM = oCoeff.getJSONArray("bindingModes");
        	
        	for(int iBM=0; iBM<aBM.length(); iBM++) {
        		System.out.println("");
        		System.out.println("Binding mode "+iBM+":");
        		System.out.println("---------------");
        		BindingMode.printJSONObjectCoefficients(fullModel, coefficients, iBM);
        	}
    	}
    	
    	if(oCoeff.has("bindingModeInteractions")) {
    		JSONArray aInt = oCoeff.getJSONArray("bindingModeInteractions");
        	
        	for(int iInt=0; iInt<aInt.length(); iInt++) {
        		System.out.println("");
        		System.out.println("Binding mode interaction "+iInt+":");
        		System.out.println("---------------------------");
        		BindingModeInteraction.printJSONObjectCoefficients(fullModel, coefficients, iInt);
        	}
    	}
    }
    
    public static void addEmptyModelFittingConstraints(String jsonSchemaPath, JSONObject oModel) {

    	// Gets the model settings
    	JSONObject oSet = oModel.getJSONObject("modelSettings");
    	
    	// Creates model fitting constraints.
    	if(!oModel.has("modelFittingConstraints"))
    		oModel.put("modelFittingConstraints", new JSONObject());
    	JSONObject oMFC = oModel.getJSONObject("modelFittingConstraints");

    	//Loops over "bindingModes", "bindingModeInteractions", "countTable", and "enrichmentModel"
    	Iterator<String> comp = oSet.keys();
    	while(comp.hasNext()) {
    		//If the next entry in 'modelSettings' is a JSONArray...
    		String k  = comp.next();
			if ( oSet.get(k) instanceof JSONArray )  {
	    		//... retrieves the relevant object...
	    		JSONArray aSet = oSet.getJSONArray(k);
	    		//...and creates the corresponding array of model fitting constraints.
	    		if(!oMFC.has(k)) 
	    			oMFC.put(k, new JSONArray());
	    		JSONArray aMFC = oMFC.getJSONArray(k);
	    		
	    		//Creates empty model fitting constraints..
	    		while(aMFC.length()<aSet.length())
	    			aMFC.put(new JSONObject());				
			}


    		
    	}
    	
    	//Validates the config to fill in missing values.
    	validateSchemaFile_O(jsonSchemaPath, oModel);
		
    }
}
