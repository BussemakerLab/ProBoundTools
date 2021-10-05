package modelComponents;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;

public class BindingModeSymmetry {

	public Map<Integer, Boolean> symmetryCMap = new HashMap<Integer, Boolean>(); //Is the base-feature complement-symmetric?
	public ArrayList<Integer> symmetryIDList = new ArrayList<Integer>(); //IDs of bases along binding site
	public ArrayList<Boolean> symmetryBarrierList = new ArrayList<Boolean>(); //Indicators of barriers across binding site. 
	int k;

	String orderedLetters;
	int nLetters;
	String letterComplement;
	//int[] rcSymDiRowShift;
	HashMap<Integer,Integer> iIntToExtMono;
	HashMap<Integer,Integer> iIntToExtDi;
	HashMap<Integer,Integer> iIntRC;
	
	/*  FORMAT OF SYMMETRY STRING:
	 *  ==========================
	 *  
	 *  FORMAT 1: Binding block definition
	 *  	- Symmetry string = "Block_1|Block_2;...|Block_N", where for each Block i we have: 
	 *  
	 *              Block_i = "ID_i:L_i:RC_i", where
	 *              
	 *               ID_i : Integer ID of block. Negative integer denotes reverse complement. 
	 *                      ID=0 is gives block with vanishing readout coefficients.
	 *               L_i  : Length of block.
	 *               RC_i : 1=(Reverse-complement symmetric block), 0=(Unconstrainted)
	 *              
	 *      - Example: 1:5:0|0:2:1|-1:5:0  gives a binding model with two identical separated blocks:
	 *           >>>>>NN<<<<<        
	 *                 
	 *  FORMAT 2: Base definition
	 *  	- The symmetry string consists of a series of bases and barriers:
	 *  	- Bases are indicated using letters, numbers, and "."
	 *  	- Each letter defines a set of equivalent bases.
	 *  	- Lower case indicates the complement of upper case
	 *  	- Numbers indicate bases that complement map to themselves (e.g. complement-symmetric).
	 *  	- "." indicate a base with vanishing readout   
	 *  	- Pipe sign "|" indicate barriers that the interactions cannot span.
	 *  	- Example 1: "ABCDE|..|edcba" is the same as above
	 *  	- Example 2: "AB1ba|..|AB1ba" is the same as above but each block is reverse-complement symmetric.
	 *  
	 */		

	//NOTE:
	// The constructor defines the binding mode symmetry by:
	//
	// 1) Saving a base-id integer for each position in the binding mode (negative being RC of positive)
	// 2) Indicating if each base ID is complement-symmetric.
	// 3) Saving all barriers.
	//
	// (The constructor does not construct the symmetry matrices, this is done by constructMMatrices() (OLd) and constructPackingVectors)
	public BindingModeSymmetry(String symString, int k, String orderedLetters, String letterComplement) {
		
		// 1. Saves the alphabet
		this.orderedLetters   = orderedLetters;
		nLetters              = orderedLetters.length();
		this.letterComplement = letterComplement;
		
		// 2. Creates a vectors that indicates how rows should be shifted to avoid empty rows
		//    in the reverse-compliment symmetric dinucleotide matrices.
		/*rcSymDiRowShift = new int[nLetters*nLetters];
		int tempShift   = 0;
		for(int i1=0; i1<nLetters; i1++) {
			for(int i2=0; i2<nLetters; i2++) {
				int iDiF = nLetters *              i1 +              i2; //Index of forward 2mer
				int iDiR = nLetters * (nLetters-i2-1) + (nLetters-i1-1); //Index of reverse 2mer
				if(iDiR<iDiF) // In this case the indicator drops of the diagonal and creates an empty row which fixed using the shift.
					tempShift+=1;
				rcSymDiRowShift[iDiF] = tempShift;
			}
		}*/
		
		// 3. Determines the order of the mononucleotides and dinucleotides on the 'external' representation (i.e., used by BindingMode):
		HashMap<String,Integer> extMonoToIndex = new HashMap<String,Integer>();
		HashMap<String,Integer> extDiToIndex   = new HashMap<String,Integer>();
		for(int i1=0; i1<nLetters; i1++) { 
			extMonoToIndex.put(  ""+orderedLetters.charAt(i1),                                       i1);
			for(int i2=0; i2<nLetters; i2++)
				extDiToIndex.put(""+orderedLetters.charAt(i1)+orderedLetters.charAt(i2), i1*nLetters+i2);
		}
		

		// 4. Determines the order of the mono- and di-nucleotides in the 'internal' representation.
		iIntRC                                 = new HashMap<Integer,Integer>();
		HashMap<Integer,String> intIndexToMono = new HashMap<Integer,String>();
		HashMap<Integer,String> intIndexToDi   = new HashMap<Integer,String>();
		String[] splitComplement = letterComplement.split(",");

//		if(splitComplement.length*2 != nLetters)
//				throw new IllegalArgumentException("The number of letters in the alphabet specification '"+orderedLetters+"' must match the number of letters in the complement specification '"+letterComplement+"'.");

		boolean selfMap = false;
		for(int i=0; i<splitComplement.length; i++) {
			//Checks the format of the complement transformation
			String[] compChar = splitComplement[i].split("-");
			if( compChar.length!=2 || compChar[0].length()!=1 || compChar[1].length()!=1 ) 
				throw new IllegalArgumentException("Each complement pair contain exctly two characters, but '"+splitComplement[i]+"' is given.");
			
			//Checks so all characters are either self-mapping or paired.
			boolean currentSelfMap = compChar[0].equals(compChar[1]);
			if(i==0) {
				selfMap = currentSelfMap;

			} else {
				if(currentSelfMap!=selfMap)
					throw new IllegalArgumentException("Either all or none of the letters must be self-complementary:"+letterComplement+"'.");
			}
			
			//Saves letters
			if(currentSelfMap) {
				intIndexToMono.put(i,            compChar[0] );
				iIntRC.put(        i,            i           );
			} else {
				intIndexToMono.put(i,            compChar[0] );
				intIndexToMono.put(nLetters-i-1, compChar[1] );
				iIntRC.put(        i,            nLetters-i-1);
				iIntRC.put(        nLetters-i-1, i           );
			}
		}
		
		for(int i1=0; i1<nLetters; i1++)
			for(int i2=0; i2<nLetters; i2++)
				intIndexToDi.put(nLetters*i1+i2, intIndexToMono.get(i1)+intIndexToMono.get(i2));
				
		//5. Create a map from the internal to external representation.
		iIntToExtMono = new HashMap<Integer,Integer>();
		iIntToExtDi   = new HashMap<Integer,Integer>();
		for(int i=0; i<nLetters; i++)
			iIntToExtMono.put(i, extMonoToIndex.get(intIndexToMono.get(i)));
		for(int i=0; i<nLetters*nLetters; i++)
			iIntToExtDi.put(i,   extDiToIndex.get(intIndexToDi.get(i)));
		
		//DEFAULT: No symmetry if symmetryStirng=(null or "null").
		if(symString == null || symString.toLowerCase().equals("null")){//Binding mode unconstrained if symString is "null"
			this.k=k;
			for(int ID=1; ID<=k; ID++){
				symmetryCMap.put(ID, false);
				symmetryIDList.add(ID);
				symmetryBarrierList.add(false);
			}
			return;
		}
		
		
		
		
		//CAEE 1: USING FORMAT 1 ABOVE.
		if(symString.split(":").length>1){
			String[] blockStrings = symString.split("\\|");

			Map<Integer, Integer[]> blockVectors = new HashMap<Integer, Integer[]>(); //Is the base-feature complement-symmetric?
			Map<Integer, Boolean> rcSym = new HashMap<Integer, Boolean>();

			int nID=0;//Current number of unique base IDs

			for(int i=0; i<blockStrings.length; i++){
				//Parses block string
				String[] entries = blockStrings[i].split(":");
				if(entries.length != 3)
					throw new IllegalArgumentException("Each block must have exactly 3 integer values.");
				int id = Integer.parseInt(entries[0]);
				int L  = Integer.parseInt(entries[1]);
				if(L<=0)
					throw new IllegalArgumentException("Length of block must larger than zero.");
				boolean rc = Integer.parseInt(entries[2]) == 1 ? true : false;

				if(rc)id = Math.abs(id); //Can't have reverse of reverse-symmetric block. 

				if(id<0 && rcSym.get(-id)!= null && rcSym.get(-id) == true )
					throw new IllegalArgumentException("Blocks are inconsistent.");
				Integer[] currentBlockIDVector;
				if(id == 0) {
					//Uses null-filled block vector if ID=0;
					currentBlockIDVector = new Integer[L];
				} else {
					//Makes sure that the block stored
					if(blockVectors.get(id) == null){

						//Adds new block to list of blocks
						if(rc){
							//Adds a reverse-symmetric block.
							Integer[] newBlockVector = new Integer[L];
							blockVectors.put( Math.abs(id), newBlockVector);
							rcSym.put( Math.abs(id), true);

							for(int j=0; j<L;j++) {
								//L=4: 0,1<>2,3
								//L=3" 0,<1>2 
								if(j<(L+1)/2){//If we are in the 1st half of the block (including center)
									nID++;
									if(2*j+1 == L) //Checks if the new ID is in the center.
										symmetryCMap.put( nID,  true);
									else {
										symmetryCMap.put( nID, false);
										symmetryCMap.put(-nID, false);
									}
									newBlockVector[j]       =  nID;
								} else {
									newBlockVector[j]       = -(nID-(j-L/2));
								}
							}
						} else {
							//Adds a non-symmetric block.
							Integer[] newBlockVectorRC = new Integer[L];
							Integer[] newBlockVector = new Integer[L];
							blockVectors.put( Math.abs(id), newBlockVector);
							blockVectors.put(-Math.abs(id), newBlockVectorRC);
							rcSym.put( Math.abs(id), false);
							rcSym.put(-Math.abs(id), false);

							for(int j=0; j<L;j++) {
								nID++;
								newBlockVector[j]       =  nID;
								newBlockVectorRC[L-1-j] = -nID;
								symmetryCMap.put(nID,  false);
								symmetryCMap.put(-nID, false);
							}
						}
					} else {
						//If the block ID already exists, check so the new and old block are consistent.
						if(rc != rcSym.get(id) || L != blockVectors.get(id).length)
							throw new IllegalArgumentException("Blocks are inconsistent.");

					}
					currentBlockIDVector = blockVectors.get(id);
				}

				//Adds block to binding mode.
				for(int j=0; j<currentBlockIDVector.length; j++){
					//Adds base IDs to list.
					symmetryIDList.add(currentBlockIDVector[j]);
					//Adds barriers to list (true for last element):
					symmetryBarrierList.add(j==currentBlockIDVector.length-1);
				}
			}

		}
		//CASE 2: USING INPUT FORMAT 2 
		else {
			Map<String, Integer> charIDRep = new HashMap<String, Integer>();//Numeric representation of character

			int nID=0;//Current number of unique IDs 
			for(int i=0; i<symString.length(); i++){
				String current = symString.substring(i,i+1);
				if(current.equals("|")){
					//Records Barrier
					if(nID==0) //Ignores barrier at first position.
						continue;
					if(symmetryBarrierList.size() < symmetryIDList.size())  //Record barrier
						symmetryBarrierList.add(true);
					else
						throw new IllegalArgumentException("Invalid barriers.");
				} else {
					if(symmetryBarrierList.size() < symmetryIDList.size()) //Record lack of barrier.
						symmetryBarrierList.add(false);
					//Parses base data
					if(current.equals(".")){ //Adds empty base.
						symmetryIDList.add(null);
					} else { 
						//Adds new non-trivial base
						if(charIDRep.get(current) == null){
							nID++;
							charIDRep.put(current.toUpperCase(),  nID);
							if(current.toUpperCase().equals(current.toLowerCase()) ) {//Self reverse-complement symmetric
								symmetryCMap.put(charIDRep.get(current.toUpperCase()),true);
							} else {
								charIDRep.put(current.toLowerCase(), -nID);
								symmetryCMap.put(charIDRep.get(current.toUpperCase()),false);
								symmetryCMap.put(charIDRep.get(current.toLowerCase()),false);

							}
						}	
						//Records base.
						symmetryIDList.add(charIDRep.get(current));
					}
				}
			}
		}

		if (symmetryIDList.size() == k)
			this.k = k;
		else
			throw new IllegalArgumentException("Symmetry string does not match k. L(String)="+symString+", k="+k);
	}


	public ArrayList<int[]> constructPackingVectors(int iFirst) {
		ArrayList<int[]> packing = new ArrayList<int[]>();

		// Looping over spacings
		int iCurrent = iFirst;
		for(int dInt=0; dInt<k;dInt++){
			int nDof = dInt>0 ? nLetters*nLetters : nLetters;
			//Adds a new packing vector
			int[] currentPacking = new int[nDof*(k-dInt)];

			//Creating a list of valid interactions
			//=================================================			
			Map<String, Integer> dinuclRep = new HashMap<String, Integer>();    // Maps a dinucleotide string (e.g. 1,-2) to an integer 
			Map<Integer, Boolean> rcSymRep = new HashMap<Integer, Boolean>();   //Indicates if a dinucleotide (integer) is reverse-complement symmetric. 
			ArrayList<Integer> intIDList   = new ArrayList<Integer>();          //List of interactions in integer notation. 
			int nInt                       = buildFeatureMap(dinuclRep, rcSymRep, intIDList, dInt);
			
			//Packing list
			//================	
			//Creating a list with the offset for each intID.
			ArrayList<Integer> intIDofToPack0 = new ArrayList<Integer>();
			for(int id=1; id<=nInt; id++) {
				intIDofToPack0.add(iCurrent);
				
				//Creating submatrices.
				iCurrent += nDof; 
				
			}

			//Filling in matrix.
			for(int x=0; x<k-dInt; x++) { 
				int i0 = intIDList.get(x)==null ? 0 : intIDofToPack0.get(Math.abs(intIDList.get(x))-1); //Fist possible possible packing index.

				if(dInt == 0) { 
					//Mononucleotide
					for(int iLet=0; iLet<nLetters; iLet++) {
						
						//Determines forward and reverse indices.
						int iMonoF =            iLet;
						int iMonoR = iIntRC.get(iLet);
						int iCol   = nDof*x+iIntToExtMono.get(iMonoF);
						
						//Determines the packing
						int packedIndex;
						if(intIDList.get(x) != null) {
							if(rcSymRep.get(intIDList.get(x)))  //self-complement base.
								packedIndex = i0 + Math.min(             iMonoF,  iMonoR);
							else
								packedIndex = i0 + (intIDList.get(x)>0 ? iMonoF : iMonoR);
						} else {
							packedIndex     = -1;
						}
						
						//Saves the packing
						currentPacking[iCol] = packedIndex;
					}
				} else {        
					//Dinucleotide interactions
					for(int iLet1=0; iLet1<nLetters; iLet1++) {
						for(int iLet2=0; iLet2<nLetters; iLet2++) {
							
							//Determines forward and reverse indices.
							int iDiF        = nLetters *            iLet1 +             iLet2;  //Index of forward 2mer
							int iDiR        = nLetters * iIntRC.get(iLet2) + iIntRC.get(iLet1); //Index of reverse 2mer
							int iCol        = nDof*x + iIntToExtDi.get(iDiF);
							//Determines the packing
							int packedIndex;
							if(intIDList.get(x) != null) {
								if(rcSymRep.get(intIDList.get(x))) { //Interaction self-complement symmetric.
									packedIndex = i0 +  Math.min(iDiF, iDiR);
								} else {
									//Forward feature = Diagonal matrix, Reverse feature = reverse identity matrix
									packedIndex = i0 + (intIDList.get(x)>0 ? iDiF : iDiR);
								}
							} else {
								packedIndex     = -1;
							}
							
							//Saves the packing
							currentPacking[iCol] = packedIndex;
						}
					}
				}
			}
			packing.add(currentPacking);
		}

		return packing;
	}


	//Creating a list of valid interactions
	int buildFeatureMap(Map<String, Integer>  dinuclRep, Map<Integer, Boolean> rcSymRep, ArrayList<Integer> intIDList, int dInt) {
		//=================================================			
		//dinuclRep -  Maps a dinucleotide string (e.g. 1,-2) to an integer 
		//rcSymRep  -  Indicates if a dinucleotide (integer) is reverse-complement symmetric. 
		//intIDList -  List of interactions in integer notation. 

		int nInt = 0; //Number of dinucleotide independent interactions.
		for(int i=0; i<k-dInt; i++){ 
			//Makes sure both bases non-zero.
			if(symmetryIDList.get(i)==null || symmetryIDList.get(i+dInt)==null){
				intIDList.add(null);
				continue;
			}
			//Makes sure the interaction doesn't span a barrier.
			boolean barrierSpanned = false;
			for(int x=i;x<i+dInt;x++)
				if(symmetryBarrierList.get(x)){
					barrierSpanned = true;
					break;					
				}
			if(barrierSpanned){
				intIDList.add(null);
				continue;
			}

			//Creates key strings for the interactions.. 
			int n1          = symmetryIDList.get(i);
			int n2          = symmetryIDList.get(i+dInt);
			int n1C         = symmetryCMap.get(n1) ? n1 : -n1;  //Complement base.
			int n2C         = symmetryCMap.get(n2) ? n2 : -n2;
			String intKey   = ""+ n1 +","+ n2;
			String intKeyRC = ""+n2C+"," +n1C;

			//Storing information about interaction if it is new. 
			if(dinuclRep.get(intKey) == null){
				nInt++;

				//Adding information about the new binding mode 				
				dinuclRep.put(intKey,nInt);

				//Determining symmetry of new interaction
				Boolean isRCSym = intKey.equals(intKeyRC);
				rcSymRep.put(nInt, isRCSym);

				//Adding information to RC interaction. 
				if (!isRCSym){
					dinuclRep.put(intKeyRC,-nInt);
					rcSymRep.put(-nInt,    isRCSym);
				}			
			}

			//Adding information about the current interaction to the binding mode. 
			intIDList.add(dinuclRep.get(intKey));				
		}
		return nInt;
	}
}
