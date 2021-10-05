package sequenceTools;
import java.util.regex.*;
import java.util.*;
//import sequenceTools.*;
//import base.DebugLog;
/*import base.SELEXConfigReader;

import config.ExperimentReference;
import config.InputDataSetStats;
import config.SELEXSequencingConfig;
import config.Sample;
import config.SequencingRunInfo;*/

import java.io.BufferedReader;
import java.io.InputStreamReader;
import java.io.FileInputStream;
import java.util.zip.GZIPInputStream;

public class FileParser {
	private String path;
	private String leftFlank;
	private String rightFlank;
	//private int fileLen;
	private LongSequence.SequenceClass sc;
	private ArrayList<Integer> seqCnt = new ArrayList<Integer>();
	private ArrayList<LongSequence> seqArr = new ArrayList<LongSequence>();

	public FileParser(String path, String leftFlank, String rightFlank, LongSequence.SequenceClass sc) {
		this.path = path;
		this.leftFlank = leftFlank;
		this.rightFlank = rightFlank;
		//this.fileLen = fileLen;
		this.sc = sc;
	}

	public ArrayList<Integer> getCounts() {
		return seqCnt;
	}

	public ArrayList<LongSequence> getSequences() {
		return seqArr;
	}

	public void createCountingTable() {
		BufferedReader br = null;
		seqCnt.clear();
		seqArr.clear();
		//HashMap<Character,Integer> table = sc.getALPHABET_LOOKUP();

		ArrayList<LongSequence> sqArr = new ArrayList<LongSequence>();

		try {
			if (path.endsWith(".gz")) {
				br = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(path))));
			}
			else {
				br = new BufferedReader(new InputStreamReader(new FileInputStream(path)));
			}

			int c;
			StringBuffer sb = new StringBuffer();
			String ms = leftFlank + "([a-zA-Z]+)" + rightFlank;
			Pattern p = Pattern.compile(ms);
			Matcher m; 

			while((c = br.read()) != -1)	{
				if (c == '\n') {
					String line = sb.toString();
					sb.setLength(0);

					m= p.matcher(line);
					if (m.find()) {
						String seqStr =  m.group(1);
						sqArr.add(sc.build(seqStr));
					}
				}
				else {
					sb.append((char)c);
				}
			}

			if(br!= null) {
				br.close();
			}

			Collections.sort(sqArr);
			LongSequence cur = sqArr.get(0);
			LongSequence tmp;
			int count = 1;

			for(int i=1; i < sqArr.size();i++) {
				tmp = sqArr.get(i);

				if(tmp.compareTo(cur) !=0) {
					seqArr.add(cur);
					seqCnt.add(count);

					cur = tmp;

					count = 1;
				}
				else {
					count++;
				}
			}

			//Handle the last group of items
			seqArr.add(cur);
			seqCnt.add(count);

		}
		catch(Exception e) {
			System.out.println(e.getMessage());
		}
	}

}