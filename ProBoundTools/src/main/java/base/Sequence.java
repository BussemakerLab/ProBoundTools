package base;

import java.util.*;

//THIS CODE IS FROM THE SELEX BIOCONDUCTOR PACKAGE

public class Sequence implements Comparable<Sequence>
{	
	
	private static int ATCG_LOOKUP[]; //translates a ACGTacgt character to a number from 0-3
	private static char ATCG_LOOKUP_REVERSE[]; //translates a number to a ACGT character
	private static long ATCG_COMPLEMENT[]; //complement mapping
	private static long mask = ( 1 << 31) -1;

	private String input;
	private long value = -1;
	private int length;
	
	static
	{
		//
		int size=256;
		ATCG_LOOKUP=new int[size];
		for(int i=0;i<size;i++)
		{
			ATCG_LOOKUP[i]=-1;
		}
		ATCG_LOOKUP['a']=0;
		ATCG_LOOKUP['A']=0;
		ATCG_LOOKUP['c']=1;
		ATCG_LOOKUP['C']=1;
		ATCG_LOOKUP['g']=2;
		ATCG_LOOKUP['G']=2;
		ATCG_LOOKUP['t']=3;
		ATCG_LOOKUP['T']=3;
		
		ATCG_LOOKUP_REVERSE=new char[]{'A','C','G','T'};
		
		ATCG_COMPLEMENT = new long[4];
		ATCG_COMPLEMENT[0]=3;//A --> T
		ATCG_COMPLEMENT[1]=2;//C --> G
		ATCG_COMPLEMENT[2]=1;//G --> C
		ATCG_COMPLEMENT[3]=0;//T --> A
	}
	
	/**
	 * creates a sequence from a string.
	 * @param in
	 */
	public Sequence(String in)
	{
//		input = in.toLowerCase();
		input = in;
		this.length = in.length();
		if(in.length()>32)
		{
			throw new RuntimeException("The length of input string["+in+"] is longer than 31.");
		}
	}

	/**
	 * creates a sequence from part of a string.
	 * @param base
	 * @param start
	 * @param length
	 */
	public Sequence(String base, int start, int length)
	{
		this.length = length;
		encode(base,start,length);//string to binary
	}
	
	/**
	 * Creates a sequence from a long value.
	 * @param in
	 * @param length
	 */
	public Sequence(long in, int length)//the number representation of the sequence
	{
		this.value=in;
		this.length=length;
	}
	
	private void encode(String input, int start, int length)
	{
		long v = 0;
		for(int i=start;i<length+start;i++)
		{
			char c=input.charAt(i);
			int tmp = ATCG_LOOKUP[c];
			if(tmp==-1)
			{
				throw new RuntimeException("Invalid char:"+c);
			}
			v <<= 2;
			v |=tmp;
		}
		
		//output
		this.value=v;
	}
	
	public static long getCharCode(char c)
	{
		return ATCG_LOOKUP[c];
	}
	
	public static char getChar(int i)
	{
		return ATCG_LOOKUP_REVERSE[i];
	}
	
	public static boolean isValidSymbol(char c)
	{
		return ATCG_LOOKUP[c]!=-1;
	}
	
	public static boolean isValidSymbolEnhanced(char c)
	{
		return (c<255 && Sequence.isValidSymbol(c));
	}
	
	private long encode()
	{
		if(value==-1)
		{
			encode(input, 0, length);
		}
		return value;
	}
	
	public long getValue()
	{
		if(value==-1)
		{
			encode();
		}
		return value;
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
			long v=this.value;
			StringBuffer sb=new StringBuffer();
			while(l!=0)
			{
				long tmp= v&3;
				sb.append(ATCG_LOOKUP_REVERSE[(int)tmp]);
				v >>>= 2;
				l--;
			}
			int numOfZeors = this.length - sb.length();
			if(numOfZeors!=0)
			{
				StringBuffer sb2=new StringBuffer();
				while(numOfZeors>0)
				{
					sb2.append(ATCG_LOOKUP_REVERSE[0]);
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
	
	//@Override
	public int compareTo(Sequence o)
	{
//		if(this.length != o.length)
//		{
//			return this.length-o.length;
//		}
//		else
		{
			long v1=getValue();
			long v2=o.getValue();
			
			int v1_u= (int) ((v1) >> 30); 
			int v2_u= (int) ((v2) >> 30);
			if(v1_u != v2_u)
			{
				return (v1_u-v2_u);
			}
			
			int v1_l= (int)(v1& mask); 
			int v2_l= (int)(v2& mask);
			
    		return (v1_l - v2_l);
		}
	}
	
	
	@Override
	public int hashCode()
	{
		Long v=this.getValue();
		if(v> Integer.MAX_VALUE)
		{
			return v.hashCode();
		}
		else
		{
			return (int)(long)(v);
		}
		//return v.hashCode();
	}

	@Override
	public boolean equals(Object obj)
	{
		if (this == obj)
			return true;
		Sequence other = (Sequence) obj;
		if (length != other.length)
		{
			System.out.println(length + " =?= "+ other.length);
			System.out.println(this.getString() + " =?= "+ other.getString());
			return false;
		}
		return this.getValue() == other.getValue();
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
	public Sequence getReverseComplement()
	{
		long val =  this.getValue();
		int len =  this.getLength();
		
		long val2 = 0;
		for(int i=0;i<len;i++)
		{
			long tmp= val & 3;
			val  >>= 2;
			val2 <<= 2;
			tmp = ATCG_COMPLEMENT[(int)tmp];
			val2 = tmp|val2; 
		}
		
		Sequence seq = new Sequence(val2,len);
		return seq;
	}

    /**
     * Returns an array of features for the current sequence given an order
     * @param order  inclusive order of the features (0th order - single bases, 1st order - dimers, etc.)
     * @return map of features with their counts
     */
    public static Map<String, Integer> getFeatures(Sequence seq, byte order, Map<String, Integer> features) {
        if(order<0) return null;

        if(features==null) features = new HashMap<String, Integer>();
        String str = seq.getString();
        int length = seq.getLength();

        //List<Sequence> subSequenceList = Sequence.slidingWindow(str,0,length,0, true, order+1);
        for(int j=order+1; j>0; j--) {
            //subSequenceList.addAll(Sequence.slidingWindow(str,0,length, 0, true, j));
            List<Sequence> subSequenceList = Sequence.slidingWindow(str,0,length, 0, true, j);
            for(Sequence s : subSequenceList) {
                String key = s.getString();
                Integer currVal =  features.get(key);
                int value = currVal==null? 1 : currVal+1;
                features.put(key,value);
            }
        }
        return features;

    }

    public static Map<String, Integer> generateFeaturesMatrix(byte order) {
        if(order<0)return null;

        char[] bases = ATCG_LOOKUP_REVERSE;

        Map<String,Integer> featuresMatrix = new HashMap<String, Integer>();
        List<String> keys = new ArrayList<String>();
        // add initial values
        for(char b: bases) {
            keys.add(b + "");
        }

        for(int i=1; i<=order; i++) {
           int size = keys.size();
           for(int j=0; j<size; j++) {
               for(char b: bases) {
                   keys.add(keys.get(j) + b);
               }
           }
        }

        for(String key: keys) featuresMatrix.put(key,0);

        return featuresMatrix;
    }

    public static List<Sequence> slidingWindow(String base, int start, int end,int leftOffset,boolean useSlidingWindow, int kmerLength)
    {
        int len = end - start;
        int c = kmerLength;

        int kmerMask = 0;
        for (int i = 0; i < kmerLength; i++)
        {
            kmerMask <<= 2;
            kmerMask |= 3;
        }

        // size of output items
        int size = len - c + 1  - leftOffset; //of sliding window

        // DebugLog.log("slidingWindow:"+base.substring(start, end));
        if (size <= 0)
        {
            String subString = base.substring(start, end);
            // DebugLog.log(subString);
            throw new FatalException("Raw input [" + subString + "](length:"
                    + len + ") not long enough for required " + c
                    + "-mer counting. Offset=["+leftOffset+"] UseSlidingWindow=["+ useSlidingWindow+ "]");
        }

        if(!useSlidingWindow)
        {
            size=1;
        }

        StringBuffer sb=new StringBuffer();
        LinkedList<Sequence> result = new LinkedList<Sequence>();
        long currentChar=0;

        for (int i = 0; i < size; i++)
        {
            // result[i]=new Sequence(seq, i, c);

            if (i == 0)
            {
                result.add( new Sequence(base, start + i + leftOffset, c) );
            } else
            {
                long v = result.peekLast().getValue(); // result[i - 1].getValue();
                currentChar = Sequence.getCharCode(base.charAt(start + i + c - 1 + leftOffset));
                long newValue = (v << 2 | currentChar) & kmerMask;
                result.add( new Sequence(newValue, c) );
            }

        }
        return result;
    }
    
	public long reverseComplement() {
		long curr = this.getValue();
		long mask = 3;
		long output = 0;
		curr = ~curr;
		output = output | (curr & mask);
		for (int i=1; i<this.length; i++) {
			curr = curr >> 2;
			output = output << 2;
			output = output | (curr & mask);
		}
		return output;
	}
	
}
