package proBoundTools;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.Reader;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Date;
import java.util.HashSet;
import java.util.Iterator;
import java.util.zip.GZIPInputStream;

import org.json.JSONObject;

import base.ArraysMergerHeap;
import base.CountObject;
import base.MersenneTwisterFast;
import sequenceTools.*;
import modelComponents.MultiRoundData;

public class Misc {
	
	public static String getCurrentTimeStamp() {
		SimpleDateFormat format = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss");	    
	    return format.format(new Date());
	}
	
	public static void msg(String str, boolean verbose) {
		if(verbose)
			System.out.println(str);
			
	}
	

	public static MultiRoundData readDataFile(JSONObject config, LongSequence.SequenceClass sc, int iExp) {
		
		//Dataset info
		JSONObject oData        = config.getJSONObject("modelSettings").getJSONArray("countTable").getJSONObject(iExp);
		int l					= oData.getInt("variableRegionLength");
		String lFlank			= oData.getString("leftFlank");
		String rFlank			= oData.getString("rightFlank");
		String inputFileType    = oData.getString("inputFileType");
		int nColumns            = oData.getInt("nColumns");
		String[] samplePaths    = oData.getString("countTableFile").split(",");

		boolean verbose         = false;

		MultiRoundData out      = null;
		
		if (inputFileType.toLowerCase().equals("dat")){
			
			//Building count table from "dat" files.
			////////////////////////////////////////
			ArrayList<MultiRoundData> datasets = new ArrayList<MultiRoundData>();
			Misc.msg(">> Starting to read 'dat' files ("+Misc.getCurrentTimeStamp()+").", verbose);
			if(samplePaths.length != nColumns)
				throw new java.lang.RuntimeException("Incorrect number of files. " + 
						samplePaths.length + " were given but " + nColumns + " were expected.");

			for(int r=0; r<samplePaths.length; r++) {
				Misc.msg(">>> Starting to read file (column ="+r+"):"+samplePaths[r], verbose);
				datasets.add(readSeqFile(samplePaths[r], l, lFlank, rFlank, sc));
			}
			Misc.msg(">> Starting to combine tables. ("+Misc.getCurrentTimeStamp()+")", verbose);
			out = joinSortedTables(datasets, lFlank, rFlank, sc);	

		} else if (inputFileType.toLowerCase().equals("tsv") || inputFileType.toLowerCase().equals("tsv.gz") ){

			boolean isGzip = inputFileType.toLowerCase().equals("tsv.gz");
			
			//Building count table from "tsv" file.
			///////////////////////////////////////
			Misc.msg(">> Reading tsv file: "+samplePaths[0]+" ("+Misc.getCurrentTimeStamp()+")", verbose);
			out = readProbeCountTable(samplePaths[0], l, lFlank, rFlank, sc, isGzip);
			
			if(out.countPerRound.length != nColumns)
				throw new java.lang.RuntimeException("Incorrect number of columns in count table. The table has " + 
						out.countPerRound.length + " columns but " + nColumns + " were expected.");

		} else if (inputFileType.toLowerCase().equals("fastq")){
			
			//Building count table from "fastq" files.
			//////////////////////////////////////////
			Misc.msg(">> Reading fastq files ("+Misc.getCurrentTimeStamp()+"):", verbose);
			if(samplePaths.length != nColumns)
				throw new java.lang.RuntimeException("Incorrect number of files. " + 
						samplePaths.length + " were given but " + nColumns + " were expected.");
			ArrayList<MultiRoundData> datasets = new ArrayList<MultiRoundData>();

			for(int r=0; r<samplePaths.length; r++){
				Misc.msg(">>> Starting to read "+r, verbose);
				FileParser fp = new FileParser(samplePaths[r], lFlank, rFlank, sc);
				fp.createCountingTable();
				int nReads = 0;
				for (int i : fp.getCounts())nReads += i;
				datasets.add(new MultiRoundData(nReads, lFlank, rFlank, fp.getCounts(), fp.getSequences(), sc));
			} 
			Misc.msg(">> Starting to combine tables ("+Misc.getCurrentTimeStamp()+").", verbose);
			out = joinSortedTables(datasets, lFlank, rFlank, sc);
		} else 
			throw new java.lang.RuntimeException("No valid input format specified.");
	
		Misc.msg(">> Done ("+Misc.getCurrentTimeStamp()+").", verbose);

		return out;

	}
	

	//This function joins multiple sequence-sorted count lists into one combined count table.
	public static MultiRoundData joinSortedTables(ArrayList<MultiRoundData> datasets, String lFlank, 
			String rFlank, LongSequence.SequenceClass sc){

		int nColumns                          = datasets.size();
		int[] countPerColumn                  = new int[nColumns];

		//List containing the current position in all sequence lists. 
		int [] iList                          = new int[nColumns];

		//Length of all lists
		long [] NList = new long[nColumns];
		for(int r=0; r<nColumns; r++)NList[r] = datasets.get(r).counts.size();

		//Creates the long-probes sequences if necessary.
		for(int r=0; r<nColumns; r++)
			if(datasets.get(r).longProbes==null)
				datasets.get(r).encodeLongToLongSequence(sc);
		
		//Creates output lists. (Could be moved to the constructor)
		ArrayList<LongSequence> longProbes    = new ArrayList<LongSequence>();
		ArrayList<int[]> counts               = new ArrayList<int[]>();
		counts.add(                             new int[nColumns]); //Add first sequence with zero count 

		//Finding the first/smallest sequence across all input lists 
		LongSequence currentLongSeq           = datasets.get(0).longProbes.get(0);
		for(int r=0; r<nColumns; r++)
			if(currentLongSeq.compareTo(datasets.get(r).longProbes.get(0)) > 0)
				currentLongSeq=datasets.get(r).longProbes.get(0);
		longProbes.add(currentLongSeq);



		//Crawls along the lists, moving incrementing the index in the round with the smallest sequence.
		while(true){
			//Identifying which round has the sequence with the lowest value
			Integer rMin            = null;
			LongSequence minLongSeq = null;
			for(int r=0;r<nColumns;r++){
				//Only considers the rounds that have not reached the end			
				if(iList[r]<NList[r]){
					if(minLongSeq == null || minLongSeq.compareTo(datasets.get(r).longProbes.get(iList[r])) > 0) { 
						rMin       = r;
						minLongSeq = datasets.get(r).longProbes.get(iList[r]);
					}			
				}
			}
			//Exits loop if all sequences have been added to the combined table
			if(rMin == null) break;

			//Checks so the input list actually was sorted
			if(minLongSeq.compareTo(currentLongSeq) < 0  )
				throw new java.lang.RuntimeException("Input sequence files are not sorted and cannot be joined into a single table.");

			//Creates new entries in t output lists if the new sequence is larger than the last.
			if(minLongSeq.compareTo(currentLongSeq) > 0 ){
				currentLongSeq = minLongSeq;
				longProbes.add(currentLongSeq);
				counts.add(new int[nColumns]);
			}

			//Adds the count to the output table
			int currentCount                   = (int) datasets.get(rMin).counts.get(iList[rMin]);
			counts.get(counts.size()-1)[rMin] += currentCount;
			//Adds the running read-per-round count.
			countPerColumn[rMin]              += currentCount;

			//Moves increments index in the list we just processed
			iList[rMin]+=1;

		}
		
		return new MultiRoundData(lFlank, rFlank, counts, longProbes, countPerColumn, sc);
	}
	

	//Read TSV count table.
	public static MultiRoundData readProbeCountTable(String path, int l, String lFlank, String rFlank, LongSequence.SequenceClass sc) {
		return readProbeCountTable(path, l, lFlank, rFlank, sc, false);
	}
	
	public static MultiRoundData readProbeCountTable(String path, int l, String lFlank, String rFlank, LongSequence.SequenceClass sc, boolean isGZIP){

		ArrayList<int[]> counts            = new ArrayList<int[]>();
		ArrayList<LongSequence> longProbes = new ArrayList<LongSequence>();

		int nColumns        = 0;
		boolean firstLine   = true;
		int[] countPerRound = null;


		String[] fields;
		int [] tempCounts;

		try {
			BufferedReader br;
			if( isGZIP ) {
				InputStream fileStream = new FileInputStream(path);
				InputStream gzipStream = new GZIPInputStream(fileStream);
				Reader decoder = new InputStreamReader(gzipStream, "UTF-8");
				br = new BufferedReader(decoder);
			} else {
				br = new BufferedReader(new FileReader(path));
			}
			
			
			for(String line; (line = br.readLine()) != null; ) {
				//Splits line
				fields            = line.replaceAll("\"", "").replaceAll("\'", "").split("\t");

				//Checks table formating.
				if(firstLine){
					l             = fields[0].length();
					nColumns      = fields.length - 1;
					countPerRound = new int[nColumns];
					firstLine     = false;					
				} else {
					if(fields[0].length() != l || fields.length != nColumns+1){
						System.out.println("ERROR: The sequences do not have uniform length or the number of table columns vary.");
						throw new IllegalArgumentException();
					}
				}

				//Processes counts
				tempCounts = new int[nColumns];
				for(int r=0; r<nColumns; r++){
					int c = Integer.parseInt(fields[r+1]);
					tempCounts[r] = c;
					countPerRound[r] += c;
				}			    	
				counts.add(tempCounts);			
				longProbes.add(sc.build(fields[0]));
				
			}
			br.close();
		} catch (Exception e) {
			e.printStackTrace();
		}

		MultiRoundData out = new MultiRoundData(lFlank, rFlank, counts, longProbes, countPerRound, sc);

		return out;

	}
	
	

	//Reads a single binary file, returns 32bp MultiRoundData object.
	public static MultiRoundData readSeqFile(String samplePath, int l, String lFlank, String rFlank, LongSequence.SequenceClass sc) {

		int currCount, totCount = 0;
		long currSeq;
		
		ArrayList<Long> probes              = new ArrayList<Long>(10*1000*1000);
		ArrayList<Integer> counts	        = new ArrayList<Integer>(10*1000*1000);

		CountObject obj;
		ArraysMergerHeap.MinHeapNode node = new ArraysMergerHeap.MinHeapNode(samplePath);
		
		//First extract all probes from file and store
		{
			while ((obj = node.peek()) != null) {
				currSeq 	= obj.getKey().getValue();
				currCount	= obj.getCount(); 
				totCount	+= currCount;
				probes.add(currSeq);
				counts.add(currCount); 
				node.pop();
			}
		}

		MultiRoundData out = new MultiRoundData(l, totCount, lFlank, rFlank, counts, probes);
		
		return out;
	}
	

	//Formats double vectors
	/////////////////////////
	public static String formatVector_d(ArrayList<Double> in, String separator){
		return formatVector_d(in,separator, "","", 4);
	}

	public static String formatVector_d(ArrayList<Double> in){
		return formatVector_d(in,",", "{","}", 4);
	}
	
	public static String formatVector_d(ArrayList<Double> in, String separator, String left, String right){
		return formatVector_d(in, separator, left, right, 4);
	}
	
	public static String formatVector_d(ArrayList<Double> in, String separator, String left, String right, int nDigits){
		String seq = left; 
		for(int i=0;i<in.size();i++){
			if(i>0)seq+=separator;
			seq+= String.format("%."+nDigits+"f", in.get(i));
		}
		return seq+right;
	}
	
	/// Formating uwing the format 1.2e-4
	public static String formatVectorE_d(ArrayList<Double> in, String separator){
		return formatVectorE_d(in,separator, "","", 4);
	}

	public static String formatVectorE_d(ArrayList<Double> in){
		return formatVectorE_d(in,",", "{","}", 4);
	}
	
	public static String formatVectorE_d(ArrayList<Double> in, String separator, String left, String right){
		return formatVectorE_d(in, separator, left, right, 4);
	}
	
	public static String formatVectorE_d(ArrayList<Double> in, String separator, String left, String right, int nDigits){
		String seq = left; 
		for(int i=0;i<in.size();i++){
			if(i>0)seq+=separator;
			seq+= String.format("%."+nDigits+"e", in.get(i));
		}
		return seq+right;
	}
	
	///
	
	public static String formatVectorE_d(double[] in, String separator){
		return formatVectorE_d(in,separator, "","", 4);
	}

	public static String formatVectorE_d(double[] in){
		return formatVectorE_d(in,",", "{","}", 4);
	}
	
	public static String formatVectorE_d(double[] in, String separator, String left, String right){
		return formatVectorE_d(in, separator, left, right, 4);
	}
	
	public static String formatVectorE_d(double[] in, String separator, String left, String right, int nDigits){
		String seq = left; 
		for(int i=0;i<in.length;i++){
			if(i>0)seq+=separator;
			seq+= String.format("%."+nDigits+"e", in[i]);
		}
		return seq+right;
	}
	
	////
	
	public static String formatVector_d(double[] in, String separator){
		return formatVector_d(in,separator, "","", 4);
	}
	
	public static String formatVector_d(double[] in){
		return formatVector_d(in,",", "{","}", 4);
	}
	
	public static String formatVector_d(double[] in, String separator, String left, String right){
		return formatVector_d(in, separator, left, right, 4);
	}
	
	
	public static String formatVector_d(double[] in, String separator, String left, String right, int nDigits){
		String seq = left; 
		if(in==null) {
			seq = "null";
		} else {
			for(int i=0;i<in.length;i++){
				if(i>0)seq+=separator;
				seq+= String.format("%."+nDigits+"f", in[i]);
			}
			seq += right;
		}

		return seq;
	}
	
	//Formats integer vectors
	/////////////////////////
	public static String formatVector_i(int[] in, String separator){
		return formatVector_i(in,separator, "","");
	}
	
	public static String formatVector_i(int[] in){
		return formatVector_i(in,",", "{","}");
	}
	
	public static String formatVector_i(int[] in, String separator, String left, String right){
		String seq = left; 
		for(int i=0;i<in.length;i++){
			if(i>0)seq+=separator;
			seq+= String.format("%d", in[i]);
		}
		return seq+right;
	}
	
	
	//Formats long vectors
	/////////////////////////
	public static String formatVector_l(long[] in, String separator){
		return formatVector_l(in,separator, "","");
	}
	
	public static String formatVector_l(long[] in){
		return formatVector_l(in,",", "{","}");
	}
	
	public static String formatVector_l(long[] in, String separator, String left, String right){
		String seq = left; 
		for(int i=0;i<in.length;i++){
			if(i>0)seq+=separator;
			seq+= String.format("%d", in[i]);
		}
		return seq+right;
	}
		
	public static int[] randomSample(int N, int k, MersenneTwisterFast generator) {
		int idx;
		int[] output = new int[k];
		
		// Change sampling method based on number of samples required
		if (((double) k)/((double) N) > .33) {
			// Use Reservoir Sampling; start by filling array
	        for (int i=0; i<k; i++) {
	                output[i] = i;
	        }
	
	        // Replace elements with decreasing probability
	        for (int i=k; i<N; i++) {
	            idx = generator.nextInt(i+1);
	            if (idx<k) {
	            	output[idx] = i;
	            }
	        }
		} else {
			HashSet<Integer> h = new HashSet<Integer>();
	
	        for (int i=0; i<k; i++) {
	                h.add(generator.nextInt(N));
	        }
	        while(h.size()<k-1) {
	                h.add(generator.nextInt(N));
	        }
	        Iterator<Integer> i = h.iterator();
	        idx = 0;
			while (i.hasNext()) {
			        output[idx] = i.next();
			        idx++;
			}
		}
		return output;
	}

	public static double median(ArrayList<Double> l)
	{
		int iCenter = l.size()/2; //Index of center element for odd length, index after center element for even length.
		Collections.sort(l);
		if (l.size()%2 == 0)
			return (l.get(iCenter-1)+l.get(iCenter))/2;
		else
			return l.get(iCenter);
	}
	
	public static String joinStrings(String[] in) {
		StringBuilder s = new StringBuilder("");
		for(String si: in)
			s.append(si);
		return s.toString();
	}
}
