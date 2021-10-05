package modelComponents;

import java.util.ArrayList;
import base.*;
import java.io.PrintWriter;
import sequenceTools.*;

public class MultiRoundData{

	public int l;
	public String leftFlank;
	public String rightFlank;

	public ArrayList<Long> probes;
	public ArrayList<LongSequence> longProbes;
	public LongSequence.SequenceClass sc;
	
	// Single-round variables
	public ArrayList<Integer> counts;
	public int nCount; 
	
	//Multi-round variables
	public ArrayList<int[]> countTable;
	public int[] countPerRound;	
	
	//Constructor for single multi-round object
	public MultiRoundData(int l, String leftFlank, String rightFlank, ArrayList<int[]> countTable, 
			ArrayList<Long> probes, int[] countPerRound) {
		
		this.l			= l;
		this.leftFlank	= leftFlank;
		this.rightFlank = rightFlank;
		this.probes		= probes;

		this.countTable = countTable;
		this.countPerRound = countPerRound;
		this.longProbes = null;
	}	
	
	//Constructor for single multi-round object
	public MultiRoundData(String leftFlank, String rightFlank, ArrayList<int[]> countTable, 
			ArrayList<LongSequence> longProbes, int[] countPerRound, LongSequence.SequenceClass sc) {
		
		this.l			= longProbes.get(0).getLength();

		this.leftFlank	= leftFlank;
		this.rightFlank = rightFlank;
		this.longProbes = longProbes;

		this.countTable = countTable;
		this.countPerRound = countPerRound;
		this.sc         = sc;

		this.probes     = null;
		
		for(LongSequence p : longProbes)
			if(p.getLength() != this.l){
				System.out.println("WARNING: Probes differ in length.");
				break;
			}
		
	}	
		
	//Constructor for single round data.
	public MultiRoundData(int l, int nCount, String leftFlank, String rightFlank, ArrayList<Integer> counts, 
			ArrayList<Long> probes) {
		this.l			= l;
		this.nCount		= nCount;
		this.leftFlank	= leftFlank;
		this.rightFlank = rightFlank;
		this.counts		= counts;
		this.probes		= probes;
		this.longProbes = null;
	}
	
	//Constructor for single round data.
	public MultiRoundData(int nCount, String leftFlank, String rightFlank, ArrayList<Integer> counts, 
			ArrayList<LongSequence> longProbes, LongSequence.SequenceClass sc) {
		this.l			= longProbes.get(0).getLength(); //OBS: This should be changed if we generalize to varying length.
		this.nCount		= nCount;
		this.leftFlank	= leftFlank;
		this.rightFlank = rightFlank;
		this.counts		= counts;
		this.longProbes = longProbes;
		this.sc         = sc;
		
		this.probes     = new ArrayList<Long>(); //OBS: This will fail for long sequences and should be removed.
		for(int i=0; i<longProbes.size(); i++)
			this.probes.add((new Sequence(longProbes.get(i).toString())).getValue());
		
		for(LongSequence p : longProbes)
			if(p.getLength() != this.l){
				System.out.println("ERROR: Probes differ in length.");
				break;
			}
	}

	
	public void encodeLongToLongSequence(LongSequence.SequenceClass sc) {
		this.sc = sc;
		this.longProbes = new ArrayList<LongSequence>();
		for(int i=0; i<probes.size(); i++)			
			longProbes.add(sc.build((new Sequence(probes.get(i), l)).getString()));
	}
	
	public void writeMultiRoundTable(String path) {
		try {
			PrintWriter writer = new PrintWriter(path, "UTF-8");

			for(int i=0;i<longProbes.size();i++){
				writer.println(longProbes.get(i).getString()+","+formatVector(countTable.get(i)));
			}
			writer.close();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	public String formatVector(int[] in){
		String seq = ""; 
		for(int i=0;i<in.length;i++){
			if(i>0)seq+=",";
			seq+= String.format("%d", in[i]);
		}
		return seq;
	}

	public static void transliterate(MultiRoundData table, LongSequence.SequenceClass sc, ArrayList<String> trIn, ArrayList<String> trOut) {

		for(int iSeq = 0; iSeq < table.longProbes.size(); iSeq++) {
			String trStr   = table.longProbes.get(iSeq).getString();
			for(int iStr=0; iStr<trIn.size(); iStr++) 
				trStr = trStr.replaceAll(trIn.get(iStr), trOut.get(iStr));
			table.longProbes.set(iSeq, sc.build(trStr));
			if(table.probes!=null)
				table.probes.set(iSeq, (new Sequence(trStr)).getValue());
		}
	}
}