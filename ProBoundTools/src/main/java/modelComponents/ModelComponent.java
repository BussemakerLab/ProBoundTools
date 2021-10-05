package modelComponents;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Random;
import java.util.Set;

import org.json.*;

abstract public class ModelComponent {
	public int iComp;             // Index of model component (iBM, iInt, iTable, iEnr).
	int iFirst, nParameters;      // Index of first parameter in the combined vector.
	int nRawParameters;           // Number of raw parameters used before collapsing to free parameters.
	Math M0, M1;                  // Matrices for compressing the current block
	public String componentName;  // Name of the object.
	public boolean verbose;
	
	//Parameters determining scoring and fitting. 
	public boolean fitComponent;         // Determines if these parameters should be optimized. Must be false if includeComponent=false
	public boolean includeComponent;     // Indicates if this component should be included. 
	abstract public void setComponentFiting(boolean fit); //Functions activating/deactivating the fitting, accoring to includeComponent and freezingLevel
	
	//Parameters/functions for freezing some of the parameters:
	public int freezingLevel;             // The freezing level indicates how many components are frozen, with zero being totally free.
	protected int maxFreezeLevel;            // Maximum freezing level.
	//parameters/functions for generating variations of the settings:
	
	//public boolean componentOptimized;   // Indicates if this component should be included.
	String componentKey;                   //Key of the component in the JSON object
	abstract public void seed_component(JSONObject config);       //Seeds binding model component

	//Functions for reading and writing JSON objects, and for packing.
	abstract void saveToJSON_settings(JSONObject out                          );        // Writes the settings of a ModelComponent to a JSON object
	abstract void readFromJSON_settings(JSONObject in                         );        // Reads the settings of a ModelComponent from a JSON object
	abstract void saveToJSON_constraints(JSONObject out                          );     // Writes the parameter constraints of a ModelComponent to a JSON object
	abstract void readFromJSON_constraints(JSONObject in                         );     // Reads the parameter constraints of a ModelComponent from a JSON object
	abstract void saveToJSON_parameters(JSONObject out ,        String coefficientKey); // Writes the parameters of a ModelComponent to a JSON object
	abstract void readFromJSON_parameters(JSONObject in,        String coefficientKey); // Reads the parameters of a ModelComponent to a JSON object
	abstract int packModel_component(JSONObject packing,        int iFirst);            // Adds one ModelComponet to a JSON packing object.
	public abstract void addZeroJSON_component(JSONObject in,   String coefficientKey); // Adds components containing zero.
	
	ModelComponent(String componentKeyIn) {
		componentKey         = componentKeyIn;
		fitComponent         = false;
		includeComponent     = false;
		verbose              = false;
	}
	
	//Sets the internal parameters associated to with the new freezing level.
	public void setFreezingLevel(int level) {
		freezingLevel = level;
		setComponentFiting(fitComponent);
	}
	
	// General function for writing the component metadata.
	public void saveToJSON_settings_metadata(JSONObject out) {
		String coefficientKey = "modelSettings";
		JSONObject oComp = out.getJSONObject(coefficientKey).getJSONArray(componentKey).getJSONObject(iComp);
		oComp.put("componentName",       componentName);
		oComp.put("includeComponent",    includeComponent);
		oComp.put("freezingLevel",       freezingLevel);

	}
	
	// General function for reading the component metadata.
	public void readFromJSON_settings_metadata(JSONObject in) {

		String coefficientKey = "modelSettings";
		JSONObject oComp = in.getJSONObject(coefficientKey).getJSONArray(componentKey).getJSONObject(iComp);
		
		if(oComp.has("componentName"))
			componentName       = oComp.getString("componentName");
		
		if(oComp.has("includeComponent"))
			includeComponent    = oComp.getBoolean("includeComponent");
		
		if(oComp.has("freezingLevel"))
			freezingLevel       = oComp.getInt("freezingLevel");

	}


	// This functions sdds an empty placeholder OBJECT so that in["coefficients" or "modelSettings"][componentKey][iLast] exists
	void addEmptyJSON_component_O(JSONObject in, String coefficientKey, String componentKey, int iLast) {
		
		JSONArray oComponent;
		
		//Creates coefficients if required. 
		if(!in.has(coefficientKey))
			in.put(coefficientKey, new JSONObject());
		
		JSONObject oCoeff = in.getJSONObject(coefficientKey);
		
		//Creates binding modes if required. 
		if(!oCoeff.has(componentKey))
			oCoeff.put(componentKey, new JSONArray());
		oComponent = oCoeff.getJSONArray(componentKey);
		
		//Makes sure that iBM points to a (possible empty) object.
		while(oComponent.length() <= iLast)
			oComponent.put(new JSONObject());

	}
	
	//Saves a list of model components to a JSON object.
	public static JSONObject saveToJSON(ArrayList<ModelComponent> component, String coefficientKey) {
		
		//Creates output object
		JSONObject out = createEmptyJSONModel(coefficientKey);

		for(int i=0; i<component.size(); i++)
			component.get(i).saveToJSON_settings(out);
		
		for(int i=0; i<component.size(); i++)
			component.get(i).saveToJSON_parameters(out, coefficientKey);

		for(int i=0; i<component.size(); i++)
			component.get(i).saveToJSON_constraints(out);
		

		return out;
	}
	
	//Reads a settings and parameters in a JSON object and updates the model components.   
	public static void readFromJSON(JSONObject in, String coefficientKey, ArrayList<ModelComponent> component) {
		
		for(int i=0; i<component.size(); i++) {
			component.get(i).readFromJSON_settings(in                  );
			component.get(i).readFromJSON_parameters(in, coefficientKey);
			component.get(i).readFromJSON_constraints(in               );

		}
		
	}
	
	public static JSONObject createEmptyJSONModel(String coefficientKey) {
		JSONObject out    = new JSONObject();
		addEmptyJSONObject(out, coefficientKey);
		addEmptyJSONObject(out, "modelSettings");
		return out;
	}
	
	public static void addEmptyJSONObject(JSONObject baseObject, String key) {
		JSONObject newObject = new JSONObject();
		baseObject.put(key, newObject);
		newObject.put("bindingModes",                   new JSONArray());
		newObject.put("bindingModeInteractions",        new JSONArray());
		newObject.put("enrichmentModel",                new JSONArray());
		newObject.put("countTable",                     new JSONArray());
	}


	// Creates a packing object by calling void packModel_component(JSONObject packing, int iFirst)
	// for each model component sequentially.
	public static JSONObject packModel(ArrayList<ModelComponent> component) {
		
		String coefficientKey = "packing";
		
		JSONObject packing = createEmptyJSONModel(coefficientKey);
		
		int iLast = -1;
		for(int iComp=0; iComp < component.size(); iComp++) {
			component.get(iComp).saveToJSON_settings(packing);
			iLast = component.get(iComp).packModel_component(packing, iLast);
		}
		
		removePackingGaps(packing);

		return packing;
	};  

	

	//By default, an included component is not fit and the 
	public void setComponentInclusion(boolean include) {
		includeComponent = include;
		freezingLevel    = maxFreezeLevel;
		setComponentFiting(false);
	}
	
	//Saves a list of model regularization toa JSON object.
	public static JSONObject regularizationToJSON(ArrayList<ModelComponent> component, String coefficientKey) {
		
		//Creates output object
		JSONObject out = createEmptyJSONModel(coefficientKey);

		for(int i=0; i<component.size(); i++)
			component.get(i).saveToJSON_settings(out);
		
		return out;
	}
	
	//Updates a packing object so that there are no gaps in the packing
	///////////////////////////////////////////////////////////////////
	public static int removePackingGaps(JSONObject packing) {
		Set<Integer> indexSet = new HashSet<Integer>();

		addIndicesToSet_O(packing.getJSONObject("packing"), indexSet);
		
		//Sorts all indices (except for -1, which indicates 'no parameter')
		indexSet.remove((int)-1);
		List<Integer> sortedIndexSet = new ArrayList<Integer>(indexSet);
		Collections.sort(sortedIndexSet);
		
		//Creates dictionary to update indices
		HashMap<Integer,Integer> indexUpdate = new HashMap<Integer,Integer>();
		indexUpdate.put((int)-1, (int)-1);
		for(int i=0; i<sortedIndexSet.size(); i++)
			indexUpdate.put(sortedIndexSet.get(i), (int) i);
		
		//Updates the indices in the packing
		updateIndices_O(packing.getJSONObject("packing"), indexUpdate);
		
		return sortedIndexSet.size();
		
	}

	//Adds all indices in a JSONObject to the index set
	private static void addIndicesToSet_O(JSONObject in, Set<Integer> indexSet) {
		//Iterator<?> keys = in.keys();
		Iterator<String> keys = in.keys();

		while( keys.hasNext() ) {
			String key = (String)keys.next();
			if ( in.get(key) instanceof JSONObject )      addIndicesToSet_O(in.getJSONObject(key), indexSet);
			else if ( in.get(key) instanceof JSONArray )  addIndicesToSet_A(in.getJSONArray(key),  indexSet);
			else                                          indexSet.add(in.getInt(key));
		}
	}

	//Array version of addIndicesToSet_O
	private static void addIndicesToSet_A(JSONArray in, Set<Integer> indexSet) {

		for(int i=0; i<in.length(); i++) {
			if ( in.get(i) instanceof JSONObject )      addIndicesToSet_O(in.getJSONObject(i), indexSet);
			else if ( in.get(i) instanceof JSONArray )  addIndicesToSet_A(in.getJSONArray(i),  indexSet);
			else                                        indexSet.add(in.getInt(i));
		}
	}
	
	//Updates all indices in a JSONObject using a hashmap
	private static void updateIndices_O(JSONObject in, HashMap<Integer,Integer> indexMap) {
		Iterator<String> keys = in.keys();

		while( keys.hasNext() ) {
			String key = (String)keys.next();
			if ( in.get(key) instanceof JSONObject )      updateIndices_O(in.getJSONObject(key), indexMap);
			else if ( in.get(key) instanceof JSONArray )  updateIndices_A(in.getJSONArray(key),  indexMap);
			else                                          in.put(key, indexMap.get(in.getInt(key)));
		}
	}
	
	//Array version of updateIndices_O
	private static void updateIndices_A(JSONArray in, HashMap<Integer,Integer> indexMap) {
		
		for(int i=0; i<in.length(); i++) {
			if ( in.get(i) instanceof JSONObject )      updateIndices_O(in.getJSONObject(i), indexMap);
			else if ( in.get(i) instanceof JSONArray )  updateIndices_A(in.getJSONArray(i),  indexMap);
			else                                        in.put(i, indexMap.get(in.getInt(i)));
		}
	}

	// Gets the maximum value of all indices
	////////////////////////////////////////
	public static Integer maxIndex_O(JSONObject in) {

		Integer maxIndex = null;
		
		Iterator<String> keys = in.keys();

		while( keys.hasNext() ) {
			Integer newIndex = null;
			String key = (String)keys.next();
			if ( in.get(key) instanceof JSONObject )      newIndex = maxIndex_O(in.getJSONObject(key));
			else if ( in.get(key) instanceof JSONArray )  newIndex = maxIndex_A(in.getJSONArray(key));
			else                                          newIndex =            in.getInt(key);
			
			if(maxIndex==null || (newIndex!=null && newIndex>maxIndex)) maxIndex = newIndex;
		}
		
		return maxIndex;
	}
	
	//Array version of updateIndices_O
	public static Integer maxIndex_A(JSONArray in) {
		
		Integer maxIndex = null;

		for(int i=0; i<in.length(); i++) {
			Integer newIndex = null;
			if ( in.get(i) instanceof JSONObject )      newIndex = maxIndex_O(in.getJSONObject(i));
			else if ( in.get(i) instanceof JSONArray )  newIndex = maxIndex_A(in.getJSONArray(i));
			else                                        newIndex =            in.getInt(i);
			
			if(maxIndex==null || (newIndex!=null && newIndex>maxIndex)) maxIndex = newIndex;
		}
		
		return maxIndex;
	}

	// FUNCTIONS FOR PACKING PARAMETERS
	///////////////////////////////////

	public static void packParameters_O(JSONObject model, JSONObject packing, double[] param) {
		
		Iterator<String> keys = packing.keys();

		while( keys.hasNext() ) {
			String key = (String)keys.next();
			if (      packing.get(key) instanceof JSONObject )  packParameters_O(model.getJSONObject(key), packing.getJSONObject(key), param);
			else if ( packing.get(key) instanceof JSONArray  )  packParameters_A(model.getJSONArray(key),  packing.getJSONArray(key),  param);
			else if ( packing.getInt(key) != -1 )               param[packing.getInt(key)] = model.getDouble(key);
			
		}
		
		return ;
	}
	
	public static void packParameters_A(JSONArray model, JSONArray packing, double[] param) {
		for(int i=0; i<packing.length(); i++) {
			if (      packing.get(i) instanceof JSONObject )  packParameters_O(model.getJSONObject(i), packing.getJSONObject(i), param);
			else if ( packing.get(i) instanceof JSONArray  )  packParameters_A(model.getJSONArray(i),  packing.getJSONArray(i),  param);
			else if ( packing.getInt(i) != -1 )               param[packing.getInt(i)] = model.getDouble(i);
		}
		
		return ;
	}

	// FUNCTIONS FOR UNPACKING PARAMETERS
	/////////////////////////////////////
	public static void unpackParameters_O(double[] param, JSONObject packing, JSONObject model) {
		
		Iterator<String> keys = packing.keys();

		while( keys.hasNext() ) {
			String key = (String)keys.next();
			if (      packing.get(key) instanceof JSONObject ) {
				if(!model.has(key))
					model.put(key, new JSONObject());				
				unpackParameters_O(param, packing.getJSONObject(key), model.getJSONObject(key));
				
			} else if ( packing.get(key) instanceof JSONArray  ) {
				
				if(!model.has(key))
					model.put(key, new JSONArray());
				unpackParameters_A(param, packing.getJSONArray(key),  model.getJSONArray(key) );
				
			} else if ( packing.getInt(key) != -1 ) {
				
				model.put(key, param[packing.getInt(key)]);
				
			}
		}
	}
	
	public static void unpackParameters_A(double[] param, JSONArray packing, JSONArray model) {
		
		for(int i=0; i<packing.length(); i++) {
			if (      packing.get(i) instanceof JSONObject ) {

				while(model.length() <= i)
					model.put(new JSONObject());
				unpackParameters_O(param, packing.getJSONObject(i), model.getJSONObject(i));
				
			} else if ( packing.get(i) instanceof JSONArray  ) {
				
				while(model.length() <= i)
					model.put(new JSONArray());
				unpackParameters_A(param, packing.getJSONArray(i),  model.getJSONArray(i));
				
			} else if ( packing.getInt(i) != -1 ) {
				
				while(model.length() <= i)
					model.put(0.0);
				model.put(i, param[packing.getInt(i)]);
			}
		}
	}

	// FUNCTIONS FOR ADDING (POSSIBLY WITH A WEIGHT) A NEW TERM TO A SUM-OBJECT
	////////////////////////////////////////////////////////////////////////////
	public static void addJSONObject_O(JSONObject sum, JSONObject newTerm) { addJSONObject_O(sum, newTerm, 1);}
	public static void addJSONObject_A(JSONArray sum,  JSONArray newTerm)  { addJSONObject_A(sum, newTerm, 1);}
	
	public static void addJSONObject_O(JSONObject sum, JSONObject newTerm, double weight) {
		
		Iterator<String> keys = newTerm.keys();

		while( keys.hasNext() ) {
			String key = (String)keys.next();
			
			//Creates a new vale in sum if the current term is not there.
			if(!sum.has(key)) {
				if (      newTerm.get(key) instanceof JSONObject ) sum.put(key, new JSONObject());
				else if ( newTerm.get(key) instanceof JSONArray  ) sum.put(key, new JSONArray());
				else                                               sum.put(key, (double) 0);
			}
			
			//Recursively adds to the sum.
			if (      newTerm.get(key) instanceof JSONObject ) addJSONObject_O(sum.getJSONObject(key), newTerm.getJSONObject(key), weight);
			else if ( newTerm.get(key) instanceof JSONArray  ) addJSONObject_A(sum.getJSONArray(key),  newTerm.getJSONArray(key),  weight);
			else                                               sum.put(key,    sum.getDouble(key) +    newTerm.getDouble(key) *    weight);
		}
	}
	
	public static void addJSONObject_A(JSONArray sum, JSONArray newTerm, double weight) {
		
		for(int i=0; i<newTerm.length(); i++) {
			//Creates a new vale in sum if the current term is not there.
			if(i>=sum.length()) {
				if (      newTerm.get(i) instanceof JSONObject ) sum.put(new JSONObject());
				else if ( newTerm.get(i) instanceof JSONArray  ) sum.put(new JSONArray() );
				else                                             sum.put((double) 0      );
			}
			
			//Recursively adds to the sum.
			if (      newTerm.get(i) instanceof JSONObject ) addJSONObject_O(sum.getJSONObject(i), newTerm.getJSONObject(i), weight);
			else if ( newTerm.get(i) instanceof JSONArray  ) addJSONObject_A(sum.getJSONArray(i),  newTerm.getJSONArray(i),  weight);
			else                                             sum.put(i,      sum.getDouble(i) +    newTerm.getDouble(i) *    weight);
		}
	}
	
	// FUNCTIONS FOR MULTIPLYING JSON OBJECTS (Only common factors are used)
	//////////////////////////////////////////
		
	public static JSONObject multiplyJSONObject_O(JSONObject a, JSONObject b) {
		
		JSONObject out = new JSONObject(); 
		Set<String> commonKeys = new HashSet<String>();
		commonKeys.addAll(a.keySet());
		commonKeys.retainAll(b.keySet());
		

		for(String key : commonKeys ) {

			if(a.get(key).getClass() != b.get(key).getClass())
				throw new java.lang.RuntimeException("The fractors "+a.get(key).toString()+ " and "+b.get(key).toString()+" must be of the same type.");

			//Recursively adds to the sum.
			if (      a.get(key) instanceof JSONObject ) out.put(key, multiplyJSONObject_O(a.getJSONObject(key), b.getJSONObject(key)));
			else if ( a.get(key) instanceof JSONArray  ) out.put(key, multiplyJSONObject_A(a.getJSONArray( key), b.getJSONArray( key)));
			else                                         out.put(key,                      a.getDouble(    key)* b.getDouble(    key));
		}
		return out;
	}
	
	public static JSONArray multiplyJSONObject_A(JSONArray a, JSONArray b) {
		
		JSONArray out = new JSONArray(); 
		
		for(int i=0; i<Math.min(a.length(), b.length()); i++) {
			//Creates a new vale in sum if the current term is not there.
			if(a.get(i).getClass() != b.get(i).getClass())
				throw new java.lang.RuntimeException("The fractors "+a.get(i).toString()+ " and "+b.get(i).toString()+" must be of the same type.");

			//Recursively adds to the sum.
			if (      a.get(i) instanceof JSONObject ) out.put(multiplyJSONObject_O(a.getJSONObject(i), b.getJSONObject(i)));
			else if ( a.get(i) instanceof JSONArray  ) out.put(multiplyJSONObject_A(a.getJSONArray( i), b.getJSONArray( i)));
			else                                       out.put(                     a.getDouble(    i)* b.getDouble(    i));

		}
		return out;
	}
	
	// FUNCTIONS FOR SUMMING OVER A JSON OBJECT
	///////////////////////////////////////////
		
	public static double trJSONObject_O(JSONObject in) {
		
		double sum = 0; 

		for(String key : in.keySet() ) {

			//Recursively adds to the sum.
			if (      in.get(key) instanceof JSONObject ) sum += trJSONObject_O(in.getJSONObject(key));
			else if ( in.get(key) instanceof JSONArray  ) sum += trJSONObject_A(in.getJSONArray( key));
			else                                          sum +=                in.getDouble(    key);
		}
		return sum;
	}
	
	public static double trJSONObject_A(JSONArray in) {
		
		double sum = 0; 

		for(int i=0; i<in.length(); i++) {

			//Recursively adds to the sum.
			if (      in.get(i) instanceof JSONObject ) sum += trJSONObject_O(in.getJSONObject(i));
			else if ( in.get(i) instanceof JSONArray  ) sum += trJSONObject_A(in.getJSONArray( i));
			else                                        sum +=                in.getDouble(    i);
		}
		return sum;
	}
		
	
	// FUNCTIONS FOR EXPONENTIATING A JSON OBJECT
	///////////////////////////////////////////
		
	public static JSONObject expJSONObject_O(JSONObject a) {
		
		JSONObject out = new JSONObject(); 
		for(String key : a.keySet() ) {
			//Recursively exponentiates
			if (      a.get(key) instanceof JSONObject ) out.put(key, expJSONObject_O(a.getJSONObject(key)));
			else if ( a.get(key) instanceof JSONArray  ) out.put(key, expJSONObject_A(a.getJSONArray( key)));
			else                                         out.put(key, Math.exp(a.getDouble(key))           );
		}
		return out;
	}
	
	public static JSONArray expJSONObject_A(JSONArray a) {
		JSONArray out = new JSONArray(); 
		for(int i=0; i<a.length(); i++) {

			//Recursively exponentiates
			if (      a.get(i) instanceof JSONObject ) out.put(expJSONObject_O(a.getJSONObject(i)));
			else if ( a.get(i) instanceof JSONArray  ) out.put(expJSONObject_A(a.getJSONArray( i)));
			else                                       out.put(Math.exp(a.getDouble(i))           );

		}
		return out;
	}
	
	

	// FUNCTIONS FOR SCALAR MULTIPLYING A JSON OBJECT
	///////////////////////////////////////////

	public static JSONObject scalarMultiplyJSONObject_O(JSONObject a, double b) {

		JSONObject out = new JSONObject(); 
		for(String key : a.keySet() ) {
			//Recursively exponentiates
			if (      a.get(key) instanceof JSONObject ) out.put(key, scalarMultiplyJSONObject_O(a.getJSONObject(key), b));
			else if ( a.get(key) instanceof JSONArray  ) out.put(key, scalarMultiplyJSONObject_A(a.getJSONArray( key), b));
			else                                         out.put(key, a.getDouble(key) * b                               );
		}
		return out;
	}

	public static JSONArray scalarMultiplyJSONObject_A(JSONArray a, double b) {
		JSONArray out = new JSONArray(); 
		for(int i=0; i<a.length(); i++) {

			//Recursively exponentiates
			if (      a.get(i) instanceof JSONObject ) out.put(scalarMultiplyJSONObject_O(a.getJSONObject(i), b));
			else if ( a.get(i) instanceof JSONArray  ) out.put(scalarMultiplyJSONObject_A(a.getJSONArray( i), b));
			else                                       out.put(a.getDouble(i) * b                                );

		}
		return out;
	}

	//FUNCTIONS FOR ADDING A SCALAR TO A JSON OBJECT
	///////////////////////////////////////////

	public static JSONObject scalarAddJSONObject_O(JSONObject a, double b) {

		JSONObject out = new JSONObject(); 
		for(String key : a.keySet() ) {
			//Recursively exponentiates
			if (      a.get(key) instanceof JSONObject ) out.put(key, scalarAddJSONObject_O(a.getJSONObject(key), b));
			else if ( a.get(key) instanceof JSONArray  ) out.put(key, scalarAddJSONObject_A(a.getJSONArray( key), b));
			else                                         out.put(key, a.getDouble(key) + b                          );
		}
		return out;
	}

	public static JSONArray scalarAddJSONObject_A(JSONArray a, double b) {
		JSONArray out = new JSONArray(); 
		for(int i=0; i<a.length(); i++) {

			//Recursively exponentiates
			if (      a.get(i) instanceof JSONObject ) out.put(scalarAddJSONObject_O(a.getJSONObject(i), b));
			else if ( a.get(i) instanceof JSONArray  ) out.put(scalarAddJSONObject_A(a.getJSONArray( i), b));
			else                                       out.put(a.getDouble(i) + b                          );

		}
		return out;
	}
	

	//Clones a nested JSON object
	/////////////////////////////
	
	public static JSONObject clone_JSON_O(JSONObject inObject) {
		JSONObject outObject = new JSONObject();
		clone_JSON_O(outObject, inObject);
		return outObject;
	}
	
	public static void clone_JSON_O(JSONObject outObject, JSONObject inObject) {
		
		Iterator<String> keys = inObject.keys();

		while( keys.hasNext() ) {
			String key = (String)keys.next();
			
			if ( inObject.get(key) instanceof JSONObject ) {
				
				outObject.put(key, new JSONObject());				
				clone_JSON_O(outObject.getJSONObject(key), inObject.getJSONObject(key));
				
			} else if ( inObject.get(key) instanceof JSONArray  ) {
				
				outObject.put(key, new JSONArray());
				clone_JSON_A(outObject.getJSONArray(key), inObject.getJSONArray(key));
	
			} else if ( inObject.get(key) instanceof String  ) {
				outObject.put(key, inObject.getString(key));
			} else if ( inObject.get(key) instanceof Integer  ) {
				outObject.put(key, inObject.getInt(key));
			} else if ( inObject.get(key) instanceof Boolean  ) {
				outObject.put(key, inObject.getBoolean(key));
			} else if ( inObject.get(key) instanceof Number  ) {
				outObject.put(key, inObject.getNumber(key));
			} else 
				throw new java.lang.RuntimeException("Can't clone JSONObject: "+inObject.toString());
		}
	}
	
	public static JSONArray clone_JSON_A(JSONArray inArray) {
		JSONArray outArray = new JSONArray();
		clone_JSON_A(outArray, inArray);
		return outArray;
	}
	
	
	public static void clone_JSON_A(JSONArray outArray, JSONArray inArray) {
		
		for(int i=0; i<inArray.length(); i++) {
			//Creates a new vale in sum if the current term is not there.
			if ( inArray.get(i) instanceof JSONObject ) {
				
				outArray.put(new JSONObject());				
				clone_JSON_O(outArray.getJSONObject(i), inArray.getJSONObject(i));
				
			} else if ( inArray.get(i) instanceof JSONArray  ) {
				
				outArray.put(new JSONArray());
				clone_JSON_A(outArray.getJSONArray(i), inArray.getJSONArray(i));
	
			} else if ( inArray.get(i) instanceof String  ) {
				outArray.put(i, inArray.getString(i));
			} else if ( inArray.get(i) instanceof Integer  ) {
				outArray.put(i, inArray.getInt(i));
			} else if ( inArray.get(i) instanceof Boolean  ) {
				outArray.put(i, inArray.getBoolean(i));
			} else if ( inArray.get(i) instanceof Number  ) {
				outArray.put(i, inArray.getNumber(i));
			} else 
				throw new java.lang.RuntimeException("Can't clone JSONArray: "+inArray.toString());
		}
	}
	

	
	
	
	/////////////////////////////////////////////
	//                                         //
	//   CODE FOR MANIPULATING NESTED ARRAYS   //
	//   ===================================   //
	/////////////////////////////////////////////
	// Code for converting from nested ArrayList/bool[] to nested JSONArray.
	/////////////////////////////////////////////////////////////////////////

	static public JSONArray JSONArrayConvert_b(boolean[] values) {
		JSONArray newEntry = new JSONArray();
		for(int i=0; i<values.length; i++)
			newEntry.put(values[i]);
		return newEntry;
	}
	
	
	// Code for converting from nested ArrayList/double[] to nested JSONArray.
	/////////////////////////////////////////////////////////////////////////

	static public JSONArray JSONArrayConvert_d(double[] values) {
		JSONArray newEntry = new JSONArray();
		for(int i=0; i<values.length; i++)
			newEntry.put(values[i]);
		return newEntry;
	}
	
	static public JSONArray JSONArrayConvert_dd(double[][] values) {
		JSONArray newEntry = new JSONArray();
		for(int i=0; i<values.length; i++) 
			newEntry.put(JSONArrayConvert_d(values[i]));
		return newEntry;
	}
	
	static public JSONArray JSONArrayConvert_Ad(ArrayList<double[]> values) {
		JSONArray newEntry = new JSONArray();
		for(int i=0; i<values.size(); i++) 
			newEntry.put(JSONArrayConvert_d(values.get(i)));
		return newEntry;
	}
	
	static public JSONArray JSONArrayConvert_AAd(ArrayList<ArrayList<double[]>> values) {
		JSONArray newEntry = new JSONArray();
		for(int i=0; i<values.size(); i++) 
			newEntry.put(JSONArrayConvert_Ad(values.get(i)));
		return newEntry;
	}
	
	static public JSONArray JSONArrayConvert_Add(ArrayList<double[][]> values) {
		JSONArray newEntry = new JSONArray();
		for(int i=0; i<values.size(); i++) 
			newEntry.put(JSONArrayConvert_dd(values.get(i)));
		return newEntry;
	}
	
	static public JSONArray JSONArrayConvert_AAdd(ArrayList<ArrayList<double[][]>> values) {
		JSONArray newEntry = new JSONArray();
		for(int i=0; i<values.size(); i++) 
			newEntry.put(JSONArrayConvert_Add(values.get(i)));
		return newEntry;
	}
	
	static public JSONArray JSONArrayConvert_AAAdd(ArrayList<ArrayList<ArrayList<double[][]>>> values) {
		JSONArray newEntry = new JSONArray();
		for(int i=0; i<values.size(); i++) 
			newEntry.put(JSONArrayConvert_AAdd(values.get(i)));
		return newEntry;
	}
	
	// Code for converting from nested ArrayList/int[] to nested JSONArray.
	/////////////////////////////////////////////////////////////////////////
	static public JSONArray JSONArrayConvert_i(int[] values) {
		JSONArray newEntry = new JSONArray();
		for(int i=0; i<values.length; i++)
			newEntry.put(values[i]);
		return newEntry;
	}
	
	static public JSONArray JSONArrayConvert_ii(int[][] values) {
		JSONArray newEntry = new JSONArray();
		for(int i=0; i<values.length; i++) 
			newEntry.put(JSONArrayConvert_i(values[i]));
		return newEntry;
	}
	
	static public JSONArray JSONArrayConvert_Ai(ArrayList<int[]> values) {
		JSONArray newEntry = new JSONArray();
		for(int i=0; i<values.size(); i++) 
			newEntry.put(JSONArrayConvert_i(values.get(i)));
		return newEntry;
	}
	
	static public JSONArray JSONArrayConvert_AAi(ArrayList<ArrayList<int[]>> values) {
		JSONArray newEntry = new JSONArray();
		for(int i=0; i<values.size(); i++) 
			newEntry.put(JSONArrayConvert_Ai(values.get(i)));
		return newEntry;
	}
	
	static public JSONArray JSONArrayConvert_Aii(ArrayList<int[][]> values) {
		JSONArray newEntry = new JSONArray();
		for(int i=0; i<values.size(); i++) 
			newEntry.put(JSONArrayConvert_ii(values.get(i)));
		return newEntry;
	}
	
	static public JSONArray JSONArrayConvert_AAii(ArrayList<ArrayList<int[][]>> values) {
		JSONArray newEntry = new JSONArray();
		for(int i=0; i<values.size(); i++) 
			newEntry.put(JSONArrayConvert_Aii(values.get(i)));
		return newEntry;
	}
	
	static public JSONArray JSONArrayConvert_AAAii(ArrayList<ArrayList<ArrayList<int[][]>>> values) {
		JSONArray newEntry = new JSONArray();
		for(int i=0; i<values.size(); i++) 
			newEntry.put(JSONArrayConvert_AAii(values.get(i)));
		return newEntry;
	}
	
	// Code for converting from nested JSONArray to nested ArrayList/double[].
	//////////////////////////////////////////////////////////////////////////	
	
	static public double[] readFromJSON_d(JSONArray array) {
		double[] out = new double[array.length()];
		for(int i=0; i<array.length(); i++)
			out[i] = array.getDouble(i);
		return out;
	}
	
	static public double[][] readFromJSON_dd(JSONArray array) {
		double[][] out = new double[array.length()][];
		for(int i=0; i<array.length(); i++)
			out[i] = readFromJSON_d(array.getJSONArray(i));
		return out;
	}

	static public ArrayList<double[]> readFromJSON_Ad(JSONArray array) {
		ArrayList<double[]> out = new ArrayList<double[]>();
		for(int i=0; i<array.length(); i++)
			out.add(readFromJSON_d(array.getJSONArray(i)));
		return out;
	}
	
	static public ArrayList<ArrayList<double[]>> readFromJSON_AAd(JSONArray array) {
		ArrayList<ArrayList<double[]>> out = new ArrayList<ArrayList<double[]>>();
		for(int i=0; i<array.length(); i++)
			out.add(readFromJSON_Ad(array.getJSONArray(i)));
		return out;
	}
	
	static public ArrayList<double[][]> readFromJSON_Add(JSONArray array) {
		ArrayList<double[][]> out = new ArrayList<double[][]>();
		for(int i=0; i<array.length(); i++)
			out.add(readFromJSON_dd(array.getJSONArray(i)));
		return out;
	}
	
	static public ArrayList<ArrayList<double[][]>> readFromJSON_AAdd(JSONArray array) {
		ArrayList<ArrayList<double[][]>> out = new ArrayList<ArrayList<double[][]>>();
		for(int i=0; i<array.length(); i++)
			out.add(readFromJSON_Add(array.getJSONArray(i)));
		return out;
	}
	
	static public ArrayList<ArrayList<ArrayList<double[][]>>> readFromJSON_AAAdd(JSONArray array) {
		ArrayList<ArrayList<ArrayList<double[][]>>> out = new ArrayList<ArrayList<ArrayList<double[][]>>>();
		for(int i=0; i<array.length(); i++)
			out.add(readFromJSON_AAdd(array.getJSONArray(i)));
		return out;
	}
	
	// Code for converting from nested JSONArray to nested ArrayList/int[].
	//////////////////////////////////////////////////////////////////////////	
	
	static public int[] readFromJSON_i(JSONArray array) {
		int[] out = new int[array.length()];
		for(int i=0; i<array.length(); i++)
			out[i] = array.getInt(i);
		return out;
	}
	
	// Code for converting from nested JSONArray to nested ArrayList/boolean[].
	//////////////////////////////////////////////////////////////////////////	
	
	static public boolean[] readFromJSON_b(JSONArray array) {
		boolean[] out = new boolean[array.length()];
		for(int i=0; i<array.length(); i++)
			out[i] = array.getBoolean(i);
		return out;
	}
	
	// Code for creating a copy of a nested ArrayList/int[] with constant entries
	////////////////////////////////////////////////////////////////////////////
	static public int[] constant_i(double[] values, int cons) {
		int[] out = new int[values.length];
		for(int i=0; i<out.length; i++)out[i] = cons;
		return out;
	}	
	
	static public int[][] constant_ii(double[][] values, int cons) {
		int[][] out = new int[values.length][];
		for(int i=0; i<values.length; i++)
			out[i] = constant_i(values[i], cons);
		return out;
	}	
	
	static public ArrayList<int[]> constant_Ai(ArrayList<double[]> values, int cons) {
		ArrayList<int[]> out = new ArrayList<int[]>();
		for(int i=0; i<values.size(); i++)
			out.add(constant_i(values.get(i), cons));
		return out;
	}	
	
	static public ArrayList<ArrayList<int[]>> constant_AAi(ArrayList<ArrayList<double[]>> values, int cons) {
		ArrayList<ArrayList<int[]>> out = new ArrayList<ArrayList<int[]>>();
		for(int i=0; i<values.size(); i++)
			out.add(constant_Ai(values.get(i), cons));
		return out;
	}	

	static public ArrayList<int[][]> constant_Aii(ArrayList<double[][]> values, int cons) {
		ArrayList<int[][]> out = new ArrayList<int[][]>();
		for(int i=0; i<values.size(); i++)
			out.add(constant_ii(values.get(i), cons));
		return out;
	}	
	
	static public ArrayList<ArrayList<int[][]>> constant_AAii(ArrayList<ArrayList<double[][]>> values, int cons) {
		ArrayList<ArrayList<int[][]>> out = new ArrayList<ArrayList<int[][]>>();
		for(int i=0; i<values.size(); i++)
			out.add(constant_Aii(values.get(i), cons));
		return out;
	}	
	
	static public ArrayList<ArrayList<ArrayList<int[][]>>> constant_AAAii(ArrayList<ArrayList<ArrayList<double[][]>>> values, int cons) {
		ArrayList<ArrayList<ArrayList<int[][]>>> out = new ArrayList<ArrayList<ArrayList<int[][]>>>();
		for(int i=0; i<values.size(); i++)
			out.add(constant_AAii(values.get(i), cons));
		return out;
	}	
	
	// Code for creating a copy of a nested ArrayList/double[] with constant entries
	////////////////////////////////////////////////////////////////////////////
	static public double[] constant_d(double[] values, double cons) {
		double[] out = new double[values.length];
		for(int i=0; i<out.length; i++)out[i] = cons;
		return out;
	}	
	
	static public double[][] constant_dd(double[][] values, double cons) {
		double[][] out = new double[values.length][];
		for(int i=0; i<values.length; i++)
			out[i] = constant_d(values[i], cons);
		return out;
	}	
	
	static public ArrayList<double[]> constant_Ad(ArrayList<double[]> values, double cons) {
		ArrayList<double[]> out = new ArrayList<double[]>();
		for(int i=0; i<values.size(); i++)
			out.add(constant_d(values.get(i), cons));
		return out;
	}	
	
	static public ArrayList<ArrayList<double[]>> constant_AAd(ArrayList<ArrayList<double[]>> values, double cons) {
		ArrayList<ArrayList<double[]>> out = new ArrayList<ArrayList<double[]>>();
		for(int i=0; i<values.size(); i++)
			out.add(constant_Ad(values.get(i), cons));
		return out;
	}	

	static public ArrayList<double[][]> constant_Add(ArrayList<double[][]> values, double cons) {
		ArrayList<double[][]> out = new ArrayList<double[][]>();
		for(int i=0; i<values.size(); i++)
			out.add(constant_dd(values.get(i), cons));
		return out;
	}	
	
	static public ArrayList<ArrayList<double[][]>> constant_AAdd(ArrayList<ArrayList<double[][]>> values, double cons) {
		ArrayList<ArrayList<double[][]>> out = new ArrayList<ArrayList<double[][]>>();
		for(int i=0; i<values.size(); i++)
			out.add(constant_Add(values.get(i), cons));
		return out;
	}	
	
	static public ArrayList<ArrayList<ArrayList<double[][]>>> constant_AAAdd(ArrayList<ArrayList<ArrayList<double[][]>>> values, double cons) {
		ArrayList<ArrayList<ArrayList<double[][]>>> out = new ArrayList<ArrayList<ArrayList<double[][]>>>();
		for(int i=0; i<values.size(); i++)
			out.add(constant_AAdd(values.get(i), cons));
		return out;
	}	
	
	// Code for creating a copy of a nested ArrayList/double[] with zero entries
	////////////////////////////////////////////////////////////////////////////
	static public double[] zero_d(double[] values)
	{ return constant_d(values, 0); }	
	
	static public double[][] zero_dd(double[][] values)
	{ return constant_dd(values, 0); }	
	
	static public ArrayList<double[]> zero_Ad(ArrayList<double[]> values)
	{ return constant_Ad(values, 0); }	
	
	static public ArrayList<ArrayList<double[]>> zero_AAd(ArrayList<ArrayList<double[]>> values) 
	{ return constant_AAd(values, 0); }	

	static public ArrayList<double[][]> zero_Add(ArrayList<double[][]> values) 
	{ return constant_Add(values, 0); }
	
	static public ArrayList<ArrayList<double[][]>> zero_AAdd(ArrayList<ArrayList<double[][]>> values) 
	{ return constant_AAdd(values, 0); }	
	
	static public ArrayList<ArrayList<ArrayList<double[][]>>> zero_AAAdd(ArrayList<ArrayList<ArrayList<double[][]>>> values) 
	{ return constant_AAAdd(values, 0); }	

	// Code for creating nested arrays of different dimensions
	//////////////////////////////////////////////////////////
	static public double[] zero_d(int d1) {
		return new double[d1];
	}	
	
	static public double[][] zero_dd(int d1, int d2) {
		double[][] out = new double[d1][];
		for(int i=0; i<d1; i++)
			out[i] = zero_d(d2);
		return out;
	}	
	
	static public ArrayList<double[]> zero_Ad(int d1, int d2) {
		ArrayList<double[]> out = new ArrayList<double[]>();
		for(int i=0; i<d1; i++)
			out.add(zero_d(d2));
		return out;
	}	
	
	static public ArrayList<ArrayList<double[]>> zero_AAd(int d1, int d2, int d3) {
		ArrayList<ArrayList<double[]>> out = new ArrayList<ArrayList<double[]>>();
		for(int i=0; i<d1; i++)
			out.add(zero_Ad(d2,d3));
		return out;
	}	

	static public ArrayList<double[][]> zero_Add(int d1, int d2, int d3) {
		ArrayList<double[][]> out = new ArrayList<double[][]>();
		for(int i=0; i<d1; i++)
			out.add(zero_dd(d2, d3));
		return out;
	}	
	
	static public ArrayList<ArrayList<double[][]>> zero_AAdd(int d1, int d2, int d3, int d4) {
		ArrayList<ArrayList<double[][]>> out = new ArrayList<ArrayList<double[][]>>();
		for(int i=0; i<d1; i++)
			out.add(zero_Add(d2,d3,d4));
		return out;
	}	
	
	static public ArrayList<ArrayList<ArrayList<double[][]>>> zeroCopy_AAAdd(int d1, int d2, int d3, int d4, int d5) {
		ArrayList<ArrayList<ArrayList<double[][]>>> out = new ArrayList<ArrayList<ArrayList<double[][]>>>();
		for(int i=0; i<d1; i++)
			out.add(zero_AAdd(d2, d3, d4, d5));
		return out;
	}	

	// Code for creating a exponentiated copy of a nested ArrayList/double[]
	////////////////////////////////////////////////////////////////////////
	
	static public ArrayList<Double> exp_A(ArrayList<Double> values) {
		ArrayList<Double> out = new ArrayList<Double>();
		for(int i=0; i<values.size(); i++)
			out.add(Math.exp(values.get(i)));

		return out;	
	}	
	
	static public double[] exp_d(double[] values) {
		double[] out = new double[values.length];
		for(int i=0; i<values.length; i++)
			out[i] = Math.exp(values[i]);
		return out;	
	}	
	
	static public double[][] exp_dd(double[][] values) {
		double[][] out = new double[values.length][];
		for(int i=0; i<values.length; i++)
			out[i] = exp_d(values[i]);
		return out;
	}	
	
	static public ArrayList<double[]> exp_Ad(ArrayList<double[]> values) {
		ArrayList<double[]> out = new ArrayList<double[]>();
		for(int i=0; i<values.size(); i++)
			out.add(exp_d(values.get(i)));
		return out;
	}	
	
	static public ArrayList<ArrayList<double[]>> exp_AAd(ArrayList<ArrayList<double[]>> values) {
		ArrayList<ArrayList<double[]>> out = new ArrayList<ArrayList<double[]>>();
		for(int i=0; i<values.size(); i++)
			out.add(exp_Ad(values.get(i)));
		return out;
	}	

	static public ArrayList<double[][]> exp_Add(ArrayList<double[][]> values) {
		ArrayList<double[][]> out = new ArrayList<double[][]>();
		for(int i=0; i<values.size(); i++)
			out.add(exp_dd(values.get(i)));
		return out;
	}	
	
	static public ArrayList<ArrayList<double[][]>> exp_AAdd(ArrayList<ArrayList<double[][]>> values) {
		ArrayList<ArrayList<double[][]>> out = new ArrayList<ArrayList<double[][]>>();
		for(int i=0; i<values.size(); i++)
			out.add(exp_Add(values.get(i)));
		return out;
	}	
	
	static public ArrayList<ArrayList<ArrayList<double[][]>>> exp_AAAdd(ArrayList<ArrayList<ArrayList<double[][]>>> values) {
		ArrayList<ArrayList<ArrayList<double[][]>>> out = new ArrayList<ArrayList<ArrayList<double[][]>>>();
		for(int i=0; i<values.size(); i++)
			out.add(exp_AAdd(values.get(i)));
		return out;
	}	
	
	//Code for creating a exponentiated copy of a nested ArrayList/double[]
	////////////////////////////////////////////////////////////////////////
	static public double[] clone_d(double[] values) {
		double[] out = new double[values.length];
		for(int i=0; i<values.length; i++)
			out[i] = values[i];
		return out;	
	}	

	static public double[][] clone_dd(double[][] values) {
		double[][] out = new double[values.length][];
		for(int i=0; i<values.length; i++) 
			out[i] =  values[i]==null ? null : clone_d(values[i]);

		return out;
	}	

	static public ArrayList<double[]> clone_Ad(ArrayList<double[]> values) {
		ArrayList<double[]> out = new ArrayList<double[]>();
		
		for(int i=0; i<values.size(); i++) 
			out.add( values.get(i)==null ? null : clone_d(values.get(i)) );

		return out;
	}	

	static public ArrayList<ArrayList<double[]>> clone_AAd(ArrayList<ArrayList<double[]>> values) {
		ArrayList<ArrayList<double[]>> out = new ArrayList<ArrayList<double[]>>();
		for(int i=0; i<values.size(); i++)
			out.add( values.get(i)==null ? null : clone_Ad(values.get(i)) );
//			out.add(clone_Ad(values.get(i)));
		return out;
	}	
	
	static public ArrayList<ArrayList<ArrayList<double[]>>> clone_AAAd(ArrayList<ArrayList<ArrayList<double[]>>> values) {
		ArrayList<ArrayList<ArrayList<double[]>>> out = new ArrayList<ArrayList<ArrayList<double[]>>>();
		for(int i=0; i<values.size(); i++)
			out.add( values.get(i)==null ? null : clone_AAd(values.get(i)) );
//			out.add(clone_AAd(values.get(i)));
		return out;
	}	

	static public ArrayList<double[][]> clone_Add(ArrayList<double[][]> values) {
		ArrayList<double[][]> out = new ArrayList<double[][]>();
		for(int i=0; i<values.size(); i++)
			out.add( values.get(i)==null ? null : clone_dd(values.get(i)) );
//			out.add(clone_dd(values.get(i)));
		return out;
	}	

	static public ArrayList<ArrayList<double[][]>> clone_AAdd(ArrayList<ArrayList<double[][]>> values) {
		ArrayList<ArrayList<double[][]>> out = new ArrayList<ArrayList<double[][]>>();
		for(int i=0; i<values.size(); i++)
			out.add( values.get(i)==null ? null : clone_Add(values.get(i)) );
//			out.add(clone_Add(values.get(i)));
		return out;
	}	

	static public ArrayList<ArrayList<ArrayList<double[][]>>> clone_AAAdd(ArrayList<ArrayList<ArrayList<double[][]>>> values) {
		ArrayList<ArrayList<ArrayList<double[][]>>> out = new ArrayList<ArrayList<ArrayList<double[][]>>>();
		for(int i=0; i<values.size(); i++)
			out.add( values.get(i)==null ? null : clone_AAdd(values.get(i)) );
//			out.add(clone_AAdd(values.get(i)));
		return out;
	}	


	//Code for creating a exponentiated copy of a nested ArrayList/double[]
	////////////////////////////////////////////////////////////////////////
	static public int[] range_d(double[] values, int iFirst) {
		int[] out = new int[values.length];
		for(int i=0; i<values.length; i++)
			out[i] = iFirst+i;
		return out;	
	}	

	static public int[][] range_dd(double[][] values, int iFirst) {
		int[][] out = new int[values.length][];
		int iCurr   = iFirst;
		for(int i=0; i<values.length; i++) {
			out[i] = range_d(values[i], iCurr);
			Integer temp = last_i(out[i]);
			iCurr = temp==null ? iCurr : temp+1;

		}
		return out;
	}	

	static public ArrayList<int[]> range_Ad(ArrayList<double[]> values, int iFirst) {
		ArrayList<int[]> out = new ArrayList<int[]>();
		int iCurr = iFirst;
		for(int i=0; i<values.size(); i++) {
			out.add(range_d(values.get(i), iCurr));
			Integer temp = last_i(out.get(i));
			iCurr = temp==null ? iCurr : temp+1;

		}
		return out;
	}	

	static public ArrayList<ArrayList<int[]>> range_AAd(ArrayList<ArrayList<double[]>> values, int iFirst) {
		ArrayList<ArrayList<int[]>> out = new ArrayList<ArrayList<int[]>>();
		int iCurr = iFirst;
		for(int i=0; i<values.size(); i++) {
			out.add(range_Ad(values.get(i), iCurr));
			Integer temp = last_Ai(out.get(i));
			iCurr = temp==null ? iCurr : temp+1;
		}
		return out;
	}	

	static public ArrayList<int[][]> range_Add(ArrayList<double[][]> values, int iFirst) {
		ArrayList<int[][]> out = new ArrayList<int[][]>();
		int iCurr = iFirst;
		for(int i=0; i<values.size(); i++) {
			out.add(range_dd(values.get(i), iCurr));
			Integer temp = last_ii(out.get(i));
			iCurr = temp==null ? iCurr : temp+1;
		}
		return out;
	}	

	static public ArrayList<ArrayList<int[][]>> range_AAdd(ArrayList<ArrayList<double[][]>> values, int iFirst) {
		ArrayList<ArrayList<int[][]>> out = new ArrayList<ArrayList<int[][]>>();
		int iCurr = iFirst;
		for(int i=0; i<values.size(); i++) {
			out.add(range_Add(values.get(i), iCurr));
			Integer temp = last_Aii(out.get(i));
			iCurr = temp==null ? iCurr : temp+1;
		}
		return out;
	}	

	static public ArrayList<ArrayList<ArrayList<int[][]>>> range_AAAdd(ArrayList<ArrayList<ArrayList<double[][]>>> values, int iFirst) {
		ArrayList<ArrayList<ArrayList<int[][]>>> out = new ArrayList<ArrayList<ArrayList<int[][]>>>();
		int iCurr = iFirst;
		for(int i=0; i<values.size(); i++) {
			out.add(range_AAdd(values.get(i), iCurr));
			Integer temp = last_AAii(out.get(i));
			iCurr = temp==null ? iCurr : temp+1;
		}
		return out;
	}	


	//Code for getting the last element, returns null if there are no elements.
	///////////////////////////////////////////////////////////////////////////
	static public Integer last_i(int[] values) {
		return values.length == 0 ? null : values[values.length-1];
	}
	
	static public Integer last_ii(int[][] values) {
		Integer out = null; 
		for(int i=values.length-1; i>=0&&out==null; i--)
			out = last_i(values[i]);
		return out;
	}
	
	static public Integer last_Ai(ArrayList<int[]> values) {
		Integer out = null; 
		for(int i=values.size()-1; i>=0&&out==null; i--)
			out = last_i(values.get(i));
		return out;
	}
	
	static public Integer last_AAi(ArrayList<ArrayList<int[]>> values) {
		Integer out = null; 
		for(int i=values.size()-1; i>=0&&out==null; i--)
			out = last_Ai(values.get(i));
		return out;
	}
	
	static public Integer last_Aii(ArrayList<int[][]> values) {
		Integer out = null; 
		for(int i=values.size()-1; i>=0&&out==null; i--)
			out = last_ii(values.get(i));
		return out;
	}
	
	static public Integer last_AAii(ArrayList<ArrayList<int[][]>> values) {
		Integer out = null; 
		for(int i=values.size()-1; i>=0&&out==null; i--)
			out = last_Aii(values.get(i));
		return out;
	}
	
	static public Integer last_AAAii(ArrayList<ArrayList<ArrayList<int[][]>>> values) {
		Integer out = null; 
		for(int i=values.size()-1; i>=0&&out==null; i--)
			out = last_AAii(values.get(i));
		return out;
	}
	
	

	//Sums over nested array.
	///////////////////////////////////////////////////////////////////////////
	static public double tr_d(double[] values) {
		double sum = 0;
		for(int i=0; i<values.length; i++)
			sum += values[i];
		return sum;	
	}	

	static public double tr_dd(double[][] values) {
		double sum = 0;
		for(int i=0; i<values.length; i++) 
			sum += tr_d(values[i]);
		return sum;
	}	

	static public double tr_Ad(ArrayList<double[]> values) {
		double sum = 0;
		for(int i=0; i<values.size(); i++) 
			sum += tr_d(values.get(i));
		return sum;
	}	
	
	static public double tr_AAd(ArrayList<ArrayList<double[]>> values) {
		double sum = 0;
		for(int i=0; i<values.size(); i++) 
			sum += tr_Ad(values.get(i));
		return sum;
	}	

	static public double tr_Add(ArrayList<double[][]> values) {
		double sum = 0;
		for(int i=0; i<values.size(); i++) 
			sum += tr_dd(values.get(i));
		return sum;
	}		

	static public double tr_AAdd(ArrayList<ArrayList<double[][]>> values) {
		double sum = 0;
		for(int i=0; i<values.size(); i++) 
			sum += tr_Add(values.get(i));
		return sum;
	}			

	static public double tr_AAAdd(ArrayList<ArrayList<ArrayList<double[][]>>> values) {
		double sum = 0;
		for(int i=0; i<values.size(); i++) 
			sum += tr_AAdd(values.get(i));
		return sum;
	}			
	
	

	//Code for creating a random copy nested arary with dimensions matching the input.
	///////////////////////////////////////////////////////////////////////////////////
	static public double[] random_d(double[] values, double epsilon) {
		double[] out = new double[values.length];
		Random generator	= new Random();
		for(int i=0; i<values.length; i++)
			out[i] = (generator.nextDouble()-0.5)*epsilon;
		return out;	
	}	

	static public double[][] random_dd(double[][] values, double epsilon) {
		double[][] out = new double[values.length][];
		for(int i=0; i<values.length; i++)
			out[i] = random_d(values[i], epsilon);
		return out;
	}	

	static public ArrayList<double[]> random_Ad(ArrayList<double[]> values, double epsilon) {
		ArrayList<double[]> out = new ArrayList<double[]>();
		for(int i=0; i<values.size(); i++)
			out.add(random_d(values.get(i), epsilon));
		return out;
	}	

	static public ArrayList<ArrayList<double[]>> random_AAd(ArrayList<ArrayList<double[]>> values, double epsilon) {
		ArrayList<ArrayList<double[]>> out = new ArrayList<ArrayList<double[]>>();
		for(int i=0; i<values.size(); i++)
			out.add(random_Ad(values.get(i), epsilon));
		return out;
	}	

	static public ArrayList<double[][]> random_Add(ArrayList<double[][]> values, double epsilon) {
		ArrayList<double[][]> out = new ArrayList<double[][]>();
		for(int i=0; i<values.size(); i++)
			out.add(random_dd(values.get(i), epsilon));
		return out;
	}	

	static public ArrayList<ArrayList<double[][]>> random_AAdd(ArrayList<ArrayList<double[][]>> values, double epsilon) {
		ArrayList<ArrayList<double[][]>> out = new ArrayList<ArrayList<double[][]>>();
		for(int i=0; i<values.size(); i++)
			out.add(random_Add(values.get(i), epsilon));
		return out;
	}	

	static public ArrayList<ArrayList<ArrayList<double[][]>>> random_AAAdd(ArrayList<ArrayList<ArrayList<double[][]>>> values, double epsilon) {
		ArrayList<ArrayList<ArrayList<double[][]>>> out = new ArrayList<ArrayList<ArrayList<double[][]>>>();
		for(int i=0; i<values.size(); i++)
			out.add(random_AAdd(values.get(i), epsilon));
		return out;
	}	
	

	//Generates integer arrays with values between 0 and (n-1)
	///////////////////////////////////////////////////////////////////////////////////
	static public int[] random_i(double[] values, int n) {
		int[] out = new int[values.length];
		Random generator	= new Random();
		for(int i=0; i<values.length; i++)
			out[i] = generator.nextInt(n);
		return out;	
	}	

	static public int[][] random_ii(double[][] values, int n) {
		int[][] out = new int[values.length][];
		for(int i=0; i<values.length; i++)
			out[i] = random_i(values[i], n);
		return out;
	}	

	static public ArrayList<int[]> random_Ai(ArrayList<double[]> values, int n) {
		ArrayList<int[]> out = new ArrayList<int[]>();
		for(int i=0; i<values.size(); i++)
			out.add(random_i(values.get(i), n));
		return out;
	}	

	static public ArrayList<ArrayList<int[]>> random_AAi(ArrayList<ArrayList<double[]>> values, int n) {
		ArrayList<ArrayList<int[]>> out = new ArrayList<ArrayList<int[]>>();
		for(int i=0; i<values.size(); i++)
			out.add(random_Ai(values.get(i), n));
		return out;
	}	

	static public ArrayList<int[][]> random_Aii(ArrayList<double[][]> values, int n) {
		ArrayList<int[][]> out = new ArrayList<int[][]>();
		for(int i=0; i<values.size(); i++)
			out.add(random_ii(values.get(i), n));
		return out;
	}	

	static public ArrayList<ArrayList<int[][]>> random_AAii(ArrayList<ArrayList<double[][]>> values, int n) {
		ArrayList<ArrayList<int[][]>> out = new ArrayList<ArrayList<int[][]>>();
		for(int i=0; i<values.size(); i++)
			out.add(random_Aii(values.get(i), n));
		return out;
	}	

	static public ArrayList<ArrayList<ArrayList<int[][]>>> random_AAAii(ArrayList<ArrayList<ArrayList<double[][]>>> values, int n) {
		ArrayList<ArrayList<ArrayList<int[][]>>> out = new ArrayList<ArrayList<ArrayList<int[][]>>>();
		for(int i=0; i<values.size(); i++)
			out.add(random_AAii(values.get(i), n));
		return out;
	}	
	
}



