package sequenceTools;

import java.util.*;
import java.util.regex.*;

import proBoundTools.Misc;

import java.util.Arrays;
import java.math.BigInteger;

public class LongSequence implements Comparable<LongSequence>
{	
	private HashMap<Character,Integer> ALPHABET_LOOKUP; //translates a ACGTacgt character to a number from 0-3
	private HashMap<Integer,Character> ALPHABET_LOOKUP_REVERSE; //translates a number to a ACGT character
	private long[] ALPHABET_COMPLEMENT; //Lookup table for complementing internal alphabet representation;
	private long ALPHABET_SIZE;
	//private static long mask = ( 1 << 31) -1;
	private static long bitMask = 0;
	//private static long comMask = 0;
	private String input;
	private long value[];
	private int bBits;
	private int basesPerLong; //Number of letters stored in each long
	private int length;

	/**
	 * creates a sequence from part of a string.
	 * @param base
	 * @param start
	 * @param length
	 */
	private LongSequence(String in, int bBits,HashMap<Character,Integer> t1,HashMap<Integer,Character> t2, long[] ac ,long as,int bPL,long bm) {
		input = in;
		this.length = in.length();
		this.bBits = bBits;
		basesPerLong = bPL;
		bitMask = bm;
		ALPHABET_LOOKUP = t1;
		ALPHABET_LOOKUP_REVERSE = t2;
		ALPHABET_COMPLEMENT = ac;
		ALPHABET_SIZE = as;
		setCompMask(length,bBits);
		value = new long[(length + basesPerLong - 1) / basesPerLong]; //int ceiling a/b = (a +b - 1)/b
		if (!in.equals(""))
			encode(in,0,in.length());
	}

	private LongSequence(String base, int start, int length,int bBits, HashMap<Character,Integer> t1,HashMap<Integer,Character> t2, long[] ac,long as,int bPL,long bm)
	{
		this.length = length;
		this.bBits = bBits;
		basesPerLong = bPL;
		bitMask = bm;
		if (start + length > base.length())
			throw new RuntimeException("start= "+start + ", length= "+length + " are not valid parameters!");
		ALPHABET_LOOKUP = t1;
		ALPHABET_LOOKUP_REVERSE = t2;
		ALPHABET_COMPLEMENT = ac;
		ALPHABET_SIZE = as;
		setCompMask(length,bBits);
		value = new long[(length + basesPerLong - 1) / basesPerLong]; //int ceiling a/b = (a +b - 1)/b
		if(!base.equals(""))
			encode(base,start,length);//string to binary
	}

	private LongSequence(long[] in, int length,int bBits,HashMap<Character,Integer> t1,HashMap<Integer,Character> t2, long[] ac,long as,int bPL,long bm)//the number representation of the sequence
	{
		value = Arrays.copyOf(in,in.length);
		this.length=length;
		this.bBits = bBits;
		basesPerLong = bPL;
		bitMask = bm;
		ALPHABET_LOOKUP = t1;
		ALPHABET_LOOKUP_REVERSE = t2;
		ALPHABET_COMPLEMENT = ac;
		ALPHABET_SIZE = as;
		setCompMask(length,bBits);
	}

	private void setCompMask(int length,int bBits) {		
		StringBuilder sb = new StringBuilder(Long.toBinaryString(ALPHABET_SIZE - 1));

		if (sb.length() < bBits) 
			for(int i= bBits - sb.length();i >0;i--) 
				sb.insert(0,'0');
		for(int i=basesPerLong;i > 1; i /= 2)
			sb.append(sb);
		//comMask = (new BigInteger(sb.toString(), 2)).longValue();
	}

	private void encode(String input, int start, int length)
	{
		int idx = 0;

		for(int i=start;i<length+start;i++)
		{
			char c=input.charAt(i);
			if(!ALPHABET_LOOKUP.containsKey(c))
				throw new RuntimeException("Invalid char:"+c);

			int tmp = ALPHABET_LOOKUP.get(c);
			if (idx % basesPerLong != 0)
				value[idx/basesPerLong] <<= bBits;

			value[idx++/basesPerLong] |=tmp;
		}
		value[(length-1)/basesPerLong] <<= (basesPerLong - length % basesPerLong) *bBits;
	}

	public long getAlphabetCharCode(char c)
	{
		return ALPHABET_LOOKUP.get(c);
	}
	
	public char getAlphabetChar(int i)
	{
		return ALPHABET_LOOKUP_REVERSE.get(i);
	}
	
	public boolean isValidAlphabetSymbol(char c)
	{
		return ALPHABET_LOOKUP.get(c)!= null;
	}

	private long[] encode()
	{
		if(value == null) {
			encode(input, 0, length);
		}
		return Arrays.copyOf(value,value.length);
	}
	
	public long[] getValue()
	{
		if(value == null) {
			encode();
		}
		return Arrays.copyOf(value,value.length);
	}
	
	/**
	 * Constructs the ACGT string from this sequence.
	 * @return
	 */
	public String getString()
	{
		if(input ==null)
		{
			int l=this.length;
			long tmp;
			int shift;
			StringBuffer sb=new StringBuffer();

			for(int i = l-1; i >= 0;i--) {
				shift =  bBits * (basesPerLong - 1 - i % basesPerLong);
				tmp= (value[i / basesPerLong] >>> shift ) & bitMask;
				sb.append(ALPHABET_LOOKUP_REVERSE.get((int)tmp));
			}

			int numOfZeors = this.length - sb.length();
			if(numOfZeors!=0)
			{
				StringBuffer sb2=new StringBuffer();
				while(numOfZeors>0)
				{
					sb2.append(ALPHABET_LOOKUP_REVERSE.get(0));
					numOfZeors--;
				}
				sb.reverse();
				sb2.append(sb);
				input=sb2.toString();
			}
			else
			{
				sb.reverse();
				input=sb.toString();
			}
		}
		return input;
	}
	
	public int getLength()
	{
		return this.length;
	}

	public int getbBits()
	{
		return this.bBits;
	}

	public HashMap<Character,Integer> getALPHABET_LOOKUP() {
		return new HashMap<Character,Integer>(ALPHABET_LOOKUP);
	}

	public HashMap<Integer,Character> getALPHABET_LOOKUP_REVERSE() {
		return new HashMap<Integer,Character>(ALPHABET_LOOKUP_REVERSE);
	}

	public long[] getALPHABET_COMPLEMENT() {
		long[] out = new long[ALPHABET_COMPLEMENT.length];
		for(int i=0; i<ALPHABET_COMPLEMENT.length; i++)
			out[i] = ALPHABET_COMPLEMENT[i];
		return out;
	}
	
	public long getALPHABET_SIZE() {
		return ALPHABET_SIZE;
	}

	public void setInput(String s) {
		input = s;
		length = s.length();
		setCompMask(length,bBits);
		value = new long[(length + basesPerLong - 1) / basesPerLong]; //int ceiling a/b = (a +b - 1)/b
		encode(input,0,length);//string to binary
	}

	//@Override
	public int compareTo(LongSequence o)
	{
		int val = this.getString().compareTo(o.getString());
		if (val > 0)
			return 1;
		if (val < 0)
			return -1;

		return 0;

	}
	
	
	@Override
	public int hashCode()
	{
		return Arrays.hashCode(value);
	}

	@Override
	public boolean equals(Object obj)
	{
		if (this == obj)
			return true;

		LongSequence other = (LongSequence) obj;
		if (length != other.length)
			return false;
	
		long[] v1 = value;
		long[] v2 = other.getValue();

		for(int i=0;i < v1.length;i++)
			if (v1[i] != v2[i])
				return false;
		
		return true;
	}

	public String toString()
	{
		return this.getString();
	}
	
	/**
	 * Returns the reverse complement of the current sequence
	 * 
	 * For example: AAAACGT ---> ACGTTTT
	 * 
	 * @return
	 */
	public LongSequence getComplement()
	{
		return new LongSequence(complement(),length,bBits,ALPHABET_LOOKUP,ALPHABET_LOOKUP_REVERSE, ALPHABET_COMPLEMENT,ALPHABET_SIZE,basesPerLong,bitMask);
	}

	public LongSequence getReverse()
	{
		return new LongSequence(reverse(),length,bBits,ALPHABET_LOOKUP,ALPHABET_LOOKUP_REVERSE, ALPHABET_COMPLEMENT,ALPHABET_SIZE,basesPerLong,bitMask);
	}

	public LongSequence getReverseComplement()
	{
		return new LongSequence(reverseComplement(),length,bBits,ALPHABET_LOOKUP,ALPHABET_LOOKUP_REVERSE, ALPHABET_COMPLEMENT, ALPHABET_SIZE,basesPerLong,bitMask);
	}
    /**
     * Returns an array of features for the current sequence given an order
     * @param order  inclusive order of the features (0th order - single bases, 1st order - dimers, etc.)
     * @return map of features with their counts
     */
    private static Map<String, List<LongSequence>> getFeatures(LongSequence seq, int order, Map<String, List<LongSequence>> features) {
        if(order<0) return null;

        if(features==null) features = new HashMap<String, List<LongSequence>>();
        String str = seq.getString();

        for(int j=order+1; j>0; j--) {
            List<LongSequence> sl = seq.slidingWindow(j,0);
			String key = str + ":" + j;
			if (!features.containsKey(key))
				features.put(key,sl);
        }

        return features;
    }

    //OLD, DEFUNCT, AND UNUSED FUNCTION FOR GENERATING AN FEATURE-COUNT LOOKUP TABLE FOR ALL k-mers UP TO A LENGTH.
    /*public static Map<String, Integer> generateFeaturesMatrix(int k) {
        if(k<0)return null;
        char[] bases=new char[256];

        Map<String,Integer> featuresMatrix = new HashMap<String, Integer>();
        List<String> keys = new ArrayList<String>();
        // add initial values
        for(char b: bases) {
            keys.add(b + "");
        }

        for(int i=1; i<=k; i++) {
           int size = keys.size();
           for(int j=0; j<size; j++) {
               for(char b: bases) {
                   keys.add(keys.get(j) + b);
               }
           }
        }

        for(String key: keys) featuresMatrix.put(key,0);

        return featuresMatrix;
    }*/

    public List<LongSequence> slidingWindow(int kmerLength, int leftOffset)
    {
        int ws = kmerLength;
		long lMask = bitMask << (basesPerLong -1)* bBits;
		long loShift = (basesPerLong - 1 - (ws - 1) % basesPerLong) *bBits;
		long lastCarry, newCarry;
        // size of output items
        int maxWin = length - leftOffset - ws + 1; //of sliding window

        // DebugLog.log("slidingWindow:"+base.substring(start, end));
        if (maxWin <= 0)
            throw new RuntimeException("Window size= " + ws + ", left offset= " + leftOffset+ " are not valid inputs!");

        LinkedList<LongSequence> result = new LinkedList<LongSequence>();
        long[] win = new long[(ws + basesPerLong - 1)/basesPerLong];

		int outIndex = 0;
		
       	for(int bi= leftOffset; bi < leftOffset + ws; bi++) {
            int shift = (basesPerLong - 1 - bi % basesPerLong) *bBits;
			long tmp = (value[bi/basesPerLong] >>> shift) & bitMask;
			if (outIndex % basesPerLong != 0)
				win[outIndex/basesPerLong] <<= bBits;
			win[outIndex++ / basesPerLong] |= tmp;
		}
		win[(ws-1)/basesPerLong] <<= (basesPerLong -  ws % basesPerLong) *bBits; 
		result.add(new LongSequence(win,ws,bBits,ALPHABET_LOOKUP,ALPHABET_LOOKUP_REVERSE,ALPHABET_COMPLEMENT,ALPHABET_SIZE,basesPerLong,bitMask));
        for (int i = 1; i < maxWin; i++) {
			lastCarry = (win[win.length-1] & lMask) >>> ((basesPerLong- 1)*bBits);
            for(int j= 0;j < win.length - 1;j++) {
  			  	newCarry = (win[j] & lMask) >>> ((basesPerLong- 1)*bBits);
				win[j] <<= bBits;
				win[j] = win[j] | lastCarry;
				lastCarry = newCarry;
            }

            int shift = (basesPerLong - 1 - (leftOffset + ws - 1 + i) % basesPerLong) *bBits;
			long tmp = ((value[(leftOffset + ws - 1 + i)/basesPerLong] >>> shift) & bitMask);
			tmp <<= loShift;
			win[win.length-1] <<= bBits;
			win[(ws-1)/basesPerLong] |= tmp; // fix the left over in last array elements

			result.add(new LongSequence(win,ws,bBits,ALPHABET_LOOKUP,ALPHABET_LOOKUP_REVERSE,ALPHABET_COMPLEMENT,ALPHABET_SIZE,basesPerLong,bitMask));
		}

        return result;
    }

    public List<long[]> getSlidingWindowCodes(int kmerLength, int leftOffset)
    {
        int ws = kmerLength;
		long lMask = bitMask << (basesPerLong -1)* bBits;
		long loShift = (basesPerLong - 1 - (ws - 1) % basesPerLong) *bBits;
		long lastCarry, newCarry;
        // size of output items
        int maxWin = length - leftOffset - ws + 1; //of sliding window

        // DebugLog.log("slidingWindow:"+base.substring(start, end));
        if (maxWin <= 0)
            throw new RuntimeException("Window size= " + ws + ", left offset= " + leftOffset+ " are not valid inputs!");

        LinkedList<long[]> result = new LinkedList<long[]>();
        long[] win = new long[(ws + basesPerLong - 1)/basesPerLong];

		int outIndex = 0;
		
       	for(int bi= leftOffset; bi < leftOffset + ws; bi++) {
            int shift = (basesPerLong - 1 - bi % basesPerLong) *bBits;
			long tmp = (value[bi/basesPerLong] >>> shift) & bitMask;
			if (outIndex % basesPerLong != 0)
				win[outIndex/basesPerLong] <<= bBits;
			win[outIndex++ / basesPerLong] |= tmp;
		}
		win[(ws-1)/basesPerLong] <<= (basesPerLong -  ws % basesPerLong) *bBits; 
		result.add(Arrays.copyOf(win,win.length));
        for (int i = 1; i < maxWin; i++) {
			lastCarry = (win[win.length-1] & lMask) >>> ((basesPerLong- 1)*bBits);
            for(int j= 0;j < win.length - 1;j++) {
  			  	newCarry = (win[j] & lMask) >>> ((basesPerLong- 1)*bBits);
				win[j] <<= bBits;
				win[j] = win[j] | lastCarry;
				lastCarry = newCarry;
            }

            int shift = (basesPerLong - 1 - (leftOffset + ws - 1 + i) % basesPerLong) *bBits;
			long tmp = ((value[(leftOffset + ws - 1 + i)/basesPerLong] >>> shift) & bitMask);
			tmp <<= loShift;
			win[win.length-1] <<= bBits;
			win[(ws-1)/basesPerLong] |= tmp; // fix the left over in last array elements

			result.add(Arrays.copyOf(win,win.length));
		}

        return result;
    }

    private static Map<String, List<long[]>> getFeaturesCodes(LongSequence seq, int order, Map<String, List<long[]>> features) {
        if(order<0) return null;

        if(features==null) features = new HashMap<String, List<long[]>>();
        String str = seq.getString();

        for(int j=order+1; j>0; j--) {
            List<long[]> sl = seq.getSlidingWindowCodes(j,0);
			String key = str + ":" + j;
			if (!features.containsKey(key))
				features.put(key,sl);
        }

        return features;
    }
 
	public long[] complement() {
		
		long[] output = new long[value.length];

		for(int i=0; i <length; i++) {
			if (i % basesPerLong != 0)
				output[i/basesPerLong] <<= bBits;
			int shift = (basesPerLong - 1 - i % basesPerLong) *bBits;
			long tmp = (value[i/basesPerLong] >>> shift) & bitMask;
			output[i/basesPerLong] |= ALPHABET_COMPLEMENT[(int)tmp];
		} 
		output[(length-1)/basesPerLong] <<= (basesPerLong - length % basesPerLong) *bBits;
		return output;
		
	}

	public long[] reverse() {
		long[] output = new long[value.length];
		for(int i=0,j=length-1; i <length; i++,j--) {
			if (i % basesPerLong != 0)
				output[i/basesPerLong] <<= bBits;
			int shift = (basesPerLong - 1 - j % basesPerLong) *bBits;
			long tmp = (value[j/basesPerLong] >>> shift) & bitMask;
			output[i/basesPerLong] |= tmp;
		} 
		output[(length-1)/basesPerLong] <<= (basesPerLong - length % basesPerLong) *bBits;

		return output;
	}	

	public long[] reverseComplement() {
		long[] output = new long[value.length];

		for(int i=0,j=length-1; i <length; i++,j--) {
			if (i % basesPerLong != 0)
				output[i/basesPerLong] <<= bBits;
			int shift = (basesPerLong - 1 - j % basesPerLong) *bBits;
			long tmp = (value[j/basesPerLong] >>> shift) & bitMask;
			output[i/basesPerLong] |= ALPHABET_COMPLEMENT[(int)tmp];
		} 
		output[(length-1)/basesPerLong] <<= (basesPerLong - length % basesPerLong) *bBits;
		return output;
	}

	public static class SequenceClass {
		private int bBits; // Bits per letter
		private HashMap<Character,Integer> ALPHABET_LOOKUP; 
		private HashMap<Integer,Character> ALPHABET_LOOKUP_REVERSE; 
		private long[] ALPHABET_COMPLEMENT;
		private long ALPHABET_SIZE;
		private String repScheme;
		private long bitMask;
		private int basesPerLong;

		public SequenceClass(String s) {
			repScheme = s;
			setUpAlphaTab(repScheme);
			setbBits();
		}

		public SequenceClass() {
			repScheme = "ATCG";
			setUpAlphaTab(repScheme);
			setbBits();
		}

		public SequenceClass(HashMap<Character,Integer> t1,HashMap<Integer,Character> t2, long[] ac,long as) {
			if (t1 == null || t2 == null || ac==null) {
				repScheme = "A-T,C-G";
				setUpAlphaTab(repScheme);
			}
			else {
				ALPHABET_LOOKUP = new HashMap<Character,Integer>();
				ALPHABET_LOOKUP_REVERSE = new HashMap<Integer,Character>();
				ALPHABET_LOOKUP.clear();
				ALPHABET_LOOKUP_REVERSE.clear();
				ALPHABET_LOOKUP.putAll(t1);
				ALPHABET_LOOKUP_REVERSE.putAll(t2);
				ALPHABET_COMPLEMENT = ac;
				ALPHABET_SIZE = as;
			}
			setbBits();
		}

		public int getbBits() {
			return bBits;
		}
	
		public long getBitMask() {
			return bitMask;
		}

		public HashMap<Character,Integer> getALPHABET_LOOKUP() {
			return new HashMap<Character,Integer>(ALPHABET_LOOKUP);
		}

		public HashMap<Integer,Character> getALPHABET_LOOKUP_REVERSE() {
			return new HashMap<Integer,Character>(ALPHABET_LOOKUP_REVERSE);
		}
		
		public long[] getALPHABET_COMPLEMENT() {
			long[] out = new long[ALPHABET_COMPLEMENT.length];
			for(int i=0; i<ALPHABET_COMPLEMENT.length; i++)
				out[i] = ALPHABET_COMPLEMENT[i];
			return out;
		}

		public long getALPHABET_SIZE() {
			return ALPHABET_SIZE;
		}

		public int getBasesPerLong() {
			return basesPerLong;
		}
		
		
		// Nine characters (each 8bits) are currently stored in the maximally left-shifted position.
		//| Long 1 | Long 2 |
		//|12345678|9.......|
		//Note: When basesPerLong*bBits!=128, there will be empty bits on the 'left' in all longs. In the case of the last
		//      long there could also be zeros on the 'right' if it is only partially filled.
		private void setbBits() {
			bBits = (int) Math.ceil(Math.log(ALPHABET_SIZE) / Math.log(2));
			bBits = (bBits < 2)? 2 : bBits;
//			while (Long.SIZE % bBits != 0)
//				bBits++; // bBits must be divisor of Long.SIZE
			basesPerLong = Long.SIZE / bBits;
			bitMask = (long)(Math.pow(2,bBits) - 1);
		}

		public LongSequence build(String in) {
			return new LongSequence(in,bBits,ALPHABET_LOOKUP,ALPHABET_LOOKUP_REVERSE,ALPHABET_COMPLEMENT,ALPHABET_SIZE,basesPerLong,bitMask);
		}

		public LongSequence build(String base, int start, int length) {
			return new LongSequence(base,start,length,bBits,ALPHABET_LOOKUP,ALPHABET_LOOKUP_REVERSE,ALPHABET_COMPLEMENT,ALPHABET_SIZE,basesPerLong,bitMask);
		}

		public LongSequence build(long[] input,int length) {
			return new LongSequence(input,length, bBits,ALPHABET_LOOKUP,ALPHABET_LOOKUP_REVERSE,ALPHABET_COMPLEMENT,ALPHABET_SIZE,basesPerLong,bitMask);
		}

		public Map<String, List<LongSequence>> getFeatures(String in, int order, Map<String, List<LongSequence>> features) {
			LongSequence sq = new LongSequence(in, bBits,ALPHABET_LOOKUP,ALPHABET_LOOKUP_REVERSE,ALPHABET_COMPLEMENT,ALPHABET_SIZE,basesPerLong,bitMask);
			return LongSequence.getFeatures(sq,order,features);
		}

		public Map<String, List<long[]>> getFeaturesCodes(String in, int order, Map<String, List<long[]>> features) {
			LongSequence sq = new LongSequence(in, bBits, ALPHABET_LOOKUP,ALPHABET_LOOKUP_REVERSE,ALPHABET_COMPLEMENT,ALPHABET_SIZE,basesPerLong,bitMask);
			return LongSequence.getFeaturesCodes(sq,order,features);
		}

		//Saves converts the string 'input' into forward and reverse bit vector representations.
		public void doubleStrand(String input, long[] fv, long[] rv) 
		{
			int length = input.length();
			if (input.equals(""))
				return;

			for(int i=0;i<length;i++) {
				char fc = input.charAt(i);
				char rc = input.charAt(length - i - 1);

				if(!ALPHABET_LOOKUP.containsKey(fc))
					throw new RuntimeException("Invalid char: "+fc);

				if (i % basesPerLong != 0) {
					fv[i/basesPerLong] <<= bBits;
					rv[i/basesPerLong] <<= bBits;
				}

				fv[i++/basesPerLong] |=                     ALPHABET_LOOKUP.get(fc);;
				rv[i++/basesPerLong] |= ALPHABET_COMPLEMENT[ALPHABET_LOOKUP.get(rc)];;
			}
			fv[(length-1)/basesPerLong] <<= (basesPerLong - length % basesPerLong) *bBits;
			rv[(length-1)/basesPerLong] <<= (basesPerLong - length % basesPerLong) *bBits;
		}

		//Given a sequence and the left and right flanking sequences, all stored as long vectors, writes forward and reverse long vectors
		//for the concatenated sequence.
		public void sequenceWithFlank(long[] input, int inLen, int offset, long[] lf, String lfStr, long[] rf, String rfStr, long[] fv, long[] rv) 
		{
			int length = inLen + lfStr.length() + rfStr.length();
			int i = 0, j= length - 1, shift = 0, fShift = 0, rShift = 0;
			long tmp;

			//Left Flank
			for(int k=0; k < lfStr.length();k++) {
				shift  = (basesPerLong - 1 - k % basesPerLong) *bBits;
				fShift = (basesPerLong - 1 - i % basesPerLong) *bBits;
				rShift = (basesPerLong - 1 - j % basesPerLong) *bBits;
				tmp = (lf[k/basesPerLong] >>> shift) & bitMask;                    //Extracts the current letter from the left flank
				fv[i++/basesPerLong] |= (                         tmp  << fShift); //Writes it in the correct position in the forward out vector
				rv[j--/basesPerLong] |= (ALPHABET_COMPLEMENT[(int)tmp] << rShift); //Writes it in the correct position in the reverse out vector

			}
			//Variable region
			for(int k=offset; k < inLen; k++) {
				shift = (basesPerLong - 1 - k % basesPerLong) *bBits;
				fShift = (basesPerLong - 1 - i % basesPerLong) *bBits;
				rShift = (basesPerLong - 1 - j % basesPerLong) *bBits;
				tmp = (input[k/basesPerLong] >>> shift) & bitMask;
				fv[i++/basesPerLong] |= (                         tmp  << fShift);
				rv[j--/basesPerLong] |= (ALPHABET_COMPLEMENT[(int)tmp] << rShift);
			}

			//Right flank
			for(int k=0; k < rfStr.length(); k++) {
				shift = (basesPerLong - 1 - k % basesPerLong) *bBits;
				fShift = (basesPerLong - 1 - i % basesPerLong) *bBits;
				rShift = (basesPerLong - 1 - j % basesPerLong) *bBits;
				tmp = (rf[k/basesPerLong] >>> shift) & bitMask;
				fv[i++/basesPerLong] |= (                         tmp  << fShift) ;
				rv[j--/basesPerLong] |= (ALPHABET_COMPLEMENT[(int)tmp] << rShift);
			}
		}

		private void setUpAlphaTab(String s) {

			if (! s.matches("([A-Za-z]\\-[A-Za-z],?)+") ) {
				throw new RuntimeException("Error alphabet must be in the format X-Y.");
			}

			ALPHABET_LOOKUP                = new HashMap<Character,Integer>();
			ALPHABET_LOOKUP_REVERSE        = new HashMap<Integer,Character>();
			HashMap<Integer,Integer> compl = new HashMap<Integer,Integer>();

			s = s.replace(",","");
			Pattern p = Pattern.compile("([A-Za-z]\\-[A-Za-z],?)+?");
			Matcher m = p.matcher(s);
			int iLetter = 0;
			while (m.find()) {		
				//Gets the new letters in the new pair
				Character c1=m.group(0).charAt(0), c2=m.group(0).charAt(2);
				
				//Adds the characters to the lookup tables if they have not yet been encountered.
				if(!ALPHABET_LOOKUP.containsKey(c1)) {
					ALPHABET_LOOKUP.put(c1, iLetter);
					ALPHABET_LOOKUP_REVERSE.put(iLetter, c1);
					iLetter++;
				}
				
				if(!ALPHABET_LOOKUP.containsKey(c2)) {
					ALPHABET_LOOKUP.put(c2, iLetter);
					ALPHABET_LOOKUP_REVERSE.put(iLetter, c2);
					iLetter++;
				}
				
				//Stores information about the complement in a temporary hash map
				int i1=ALPHABET_LOOKUP.get(c1), i2=ALPHABET_LOOKUP.get(c2);
				if(!compl.containsKey(i1))
					compl.put(i1, i2);
				if(!compl.containsKey(i2))
					compl.put(i2, i1);
			}
			
			ALPHABET_SIZE       = iLetter;
			
			//Converts the hash-map complement lookup into an array.
			ALPHABET_COMPLEMENT = new long[(int)iLetter];
			for(int i=0; i<iLetter; i++)
				ALPHABET_COMPLEMENT[i] = compl.get(i);
		
		}
	}	
}
