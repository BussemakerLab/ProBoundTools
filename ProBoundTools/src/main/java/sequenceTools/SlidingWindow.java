package sequenceTools;
import java.util.*;
import java.util.Arrays;

import proBoundTools.Misc;


public class SlidingWindow {
	private String leftFlank = "";
	private String rightFlank = "";
	private long[] lfCode;
	private long[] rfCode;
	//private LongSequence sq;
	//private double[] beta;
	private int bBits;
	private long bitMask;
	private int bPL;
	private int ws;
	private int nn_feat = 1;
	private int nuc_feat = 1;
	private LongSequence.SequenceClass sc;
	private HashMap<Character, Integer> ALPHABET_LOOKUP;
	//private HashMap<Integer, Character> ALPHABET_LOOKUP_REVERSE_DB;
	private double[] pBetas;
	//private int[] betaFeatOff;
	private int[] betaCellLen;
	//private int[] betaNNLen;
	private int[][][] betaOffset;
	private int[][] winsize_per_feat;
	private int[] nnsize_per_nuc;

	/*private SlidingWindow() {
		// disable the default constructor
	}*/

	public SlidingWindow(LongSequence.SequenceClass sc, ArrayList<ArrayList<double[]>> betas, String betaMap) {
		this.sc = sc;
		bitMask = sc.getBitMask();
		bPL = sc.getBasesPerLong();
		bBits = sc.getbBits();
		setBetaArray(sc.getALPHABET_LOOKUP(), (int) sc.getALPHABET_SIZE(), betas, betaMap);
		ALPHABET_LOOKUP = sc.getALPHABET_LOOKUP();
		//ALPHABET_LOOKUP_REVERSE_DB = sc.getALPHABET_LOOKUP_REVERSE();
	}

	public SlidingWindow(String leftFlank, String rightFlank, LongSequence.SequenceClass sc,
			ArrayList<ArrayList<double[]>> betas, String betaMap) {
		this.sc = sc;
		this.leftFlank = leftFlank;
		this.rightFlank = rightFlank;
		lfCode = sc.build(leftFlank).getValue();
		rfCode = sc.build(rightFlank).getValue();
		bitMask = sc.getBitMask();
		bPL = sc.getBasesPerLong();
		
		//Number of bits used to represent a letter
		bBits = sc.getbBits();
		setBetaArray(sc.getALPHABET_LOOKUP(), (int) sc.getALPHABET_SIZE(), betas, betaMap);
		ALPHABET_LOOKUP = sc.getALPHABET_LOOKUP();
		//ALPHABET_LOOKUP_REVERSE_DB = sc.getALPHABET_LOOKUP_REVERSE();
	}

	private void setBetaArray(HashMap<Character, Integer> alpha, int as, ArrayList<ArrayList<double[]>> rawBetas,
			String betaMap) {
		char c;
		int betaLen;

		if (rawBetas.size() == 0)
			throw new RuntimeException("Error! Collection of beta arrays is empty !");

		if (as != betaMap.length())
			throw new RuntimeException("Error! Beta mapping does not contain all alphabet characters !");

		if (rawBetas.get(0).get(0).length % betaMap.length() != 0)
			throw new RuntimeException("Error! Beta array size is not consistent with beta mapping (mono beta length = "+rawBetas.get(0).get(0).length+", beta mapping size = "+betaMap.length()+") !");

		ws = rawBetas.get(0).get(0).length / betaMap.length();
		if (rawBetas.size() > 1)
			nn_feat = (rawBetas.get(1).size() != 0) ? rawBetas.get(1).size() : 1;
		nuc_feat = rawBetas.size();
		for (int i = 1; i < rawBetas.size(); i++) {
			if (ws - i - rawBetas.get(i).size() < 0)
				throw new RuntimeException("(" + i + ",) has too many neighbor parameters !");

			for (int j = 0; j < rawBetas.get(i).size(); j++) {
				double rawLen = rawBetas.get(i).get(j).length / Math.pow(betaMap.length(), i + 1);
				if (rawLen % 1 != 0) {
					String errMsg = "Error! (" + i + "," + j + ")-beta array length ("+ rawBetas.get(i).get(j).length+") not multiple of  "+(int)Math.pow(betaMap.length(), i + 1)+"!";
					throw new RuntimeException(errMsg);
				}
					
				if (rawLen != ws -i - j) {
					String errMsg = "Error! (" + i + "," + j + ")-beta has incorrect number of k-mer cells ! " +
									"Expected  " + (ws -i - j) + " cells, got "+(int)rawLen;
					throw new RuntimeException(errMsg);
				}
			}
		}

		int rawBase = betaMap.length();
 
		//Number of spots occupied by each letter in the non-packed representation (base for long representation)
		int cvtBase = (int) Math.pow(2, bBits);
		if (bBits * rawBetas.size() > Long.SIZE)
			throw new RuntimeException("Error ! " + rawBetas.size() + "-nucleotide is too big to support");

		// Initialize motif beta array related index
		
		//Position of features in the flattened array
		//INDEXING: [0/1/2=mono/di/tri][di/tri-spacing][window index]
		betaOffset = new int[nuc_feat][nn_feat][ws]; 
		
		//Length of 1-mers features, 2-mer reatures in flattened fector
		//For ACGT:  1-mer = 4, 2-mer = 16
		//In general: 2^{bBits * k} where bBits=(bits/letter) and k=length of k=mer. 
		betaCellLen = new int[nuc_feat];//length of a specific nuc_feat array  
		for (int i = 0; i < nuc_feat; i++) {
			betaCellLen[i] = (int) Math.pow(cvtBase, i + 1);
		}

		//Initialize window size per feat array
		winsize_per_feat = new int[nuc_feat][nn_feat];
		nnsize_per_nuc = new int[nuc_feat];
		for(int i =0; i < nuc_feat; i ++) { //Loops over mono/di/tri
			nnsize_per_nuc[i] = rawBetas.get(i).size(); //Gap length for mono/di/tri
		}

		winsize_per_feat[0][0] = ws;
		for(int i=1; i < nuc_feat; i++) {
			for(int j=0; j < nn_feat; j++) {
				winsize_per_feat[i][j] = (ws - i - j > 0)? ws - i - j : 0; 
			}
		}

		//Initialize the 3D offset array
		//For nuc_feat = 0, no neighbors
		for(int k=0; k < ws; k++) { //Loops over 1-mer windows and places 1-mer "cells"
			betaOffset[0][0][k] = betaCellLen[0] * k;
		}

		//When nuc_feat > 1, set up first element of nuc_feat = 1
		if(nuc_feat > 1) {
			for(int k=0; k < winsize_per_feat[1][0]; k++) { //Loops 2-mer windows and places 2-mer "cells".
				betaOffset[1][0][k] = betaOffset[0][0][ws-1] + betaCellLen[0] + betaCellLen[1] * k;
			}
		}

		for(int i=1; i < nuc_feat; i++) {                       //Loops over length k of k-mers
			for(int j=1; j < nnsize_per_nuc[i]; j++) {          //Loops over gaps in k-mers (Starting with 1-letter gap)
				for(int k=0; k < winsize_per_feat[i][j]; k++) { //Loops over location of gapped k-mer within winndow.
					int last_win_len = winsize_per_feat[i][j-1] - 1;
					betaOffset[i][j][k] = betaOffset[i][j-1][last_win_len] + (k + 1) * betaCellLen[i]; 
				}
			}

			if (i + 1 < nuc_feat) {
				for(int k=0; k < winsize_per_feat[i+1][0]; k++) {
					int last_win_len = winsize_per_feat[i][nnsize_per_nuc[i] - 1] - 1;
					int last_nn_len = nnsize_per_nuc[i] - 1;
					betaOffset[i+1][0][k] = betaOffset[i][last_nn_len][last_win_len] 
										  + betaCellLen[i] + k * betaCellLen[i + 1];
				}
			}
		}

		int last_nn_len = nnsize_per_nuc[nuc_feat - 1] - 1;
		int last_win_len = winsize_per_feat[nuc_feat - 1][last_nn_len] - 1;

		betaLen = betaOffset[nuc_feat - 1][last_nn_len][last_win_len] + betaCellLen[nuc_feat - 1];
		pBetas = new double[betaLen];

		int[] stbMap = new int[as];

		// Create map of external alphabet to internal alphabets
		for (int i = 0; i < betaMap.length(); i++) {
			c = betaMap.charAt(i);
			if (!alpha.containsKey(c))
				throw new RuntimeException("Beta mapping string " + betaMap + " cotains invalid character " + c + "!");

			stbMap[i] = alpha.get(c);
		}

		for (int i = 0; i < nuc_feat; i++) {
			int nn = (i == 0) ? 1 : nn_feat;

			for (int m = 0; m < nn; m++) {
				int rLen = (int) Math.pow(rawBase, i + 1);

				// place elements of rawBeta array into correct position in
				// internal array
				for (int j = 0; j < winsize_per_feat[i][m]; j++) {
					for (int k = 0; k < rLen; k++) {
						int rawIndex = k;
						int cvtIndex = 0;
						int cvtMult = 1;
						int pow = i;

						while (pow >= 0) {
							cvtIndex += cvtMult * (stbMap[rawIndex % rawBase]);
							cvtMult = cvtMult * cvtBase;
							rawIndex = rawIndex / rawBase;
							pow--;
						}
						pBetas[betaOffset[i][m][j] + cvtIndex] = rawBetas.get(i).get(m)[j * rLen + k];
					}
				}
			}
		}

	}

	public double[] getBeta() {
		return Arrays.copyOf(pBetas, pBetas.length);
	}

	public ArrayList<ArrayList<Double> > slideSN(long[] seqCode, int seqLen, int offset) {
		ArrayList<ArrayList<Double>> ret = new ArrayList<ArrayList<Double>>();
		int length = 0;
		double fwSum = 0, efwSum = 0;
		double rcSum = 0, ercSum = 0;
//		double ts = 0;
		long[] fwCode, rcCode;

		for(int i=0; i < 4; i++)
			ret.add(new ArrayList<Double>());

		length = seqLen + leftFlank.length() + rightFlank.length();
		fwCode = new long[(length + bPL - 1) / bPL];
		rcCode = new long[fwCode.length];

		sc.sequenceWithFlank(seqCode, seqLen, offset, lfCode, leftFlank, rfCode, rightFlank, fwCode, rcCode);

		// loop through all possible windows include Flanks
		for (int i = 0; i < length - ws + 1; i++) {
			fwSum = 0;
			rcSum = 0;

			for (int j = 0; j < ws; j++) {
				int shift = (bPL - 1 - (j + i) % bPL) * bBits;
				fwSum += pBetas[j * betaCellLen[0] + (int) ((fwCode[(j + i) / bPL] >>> shift) & bitMask)];
				rcSum += pBetas[j * betaCellLen[0] + (int) ((rcCode[(j + i) / bPL] >>> shift) & bitMask)];
			}

			efwSum = Math.exp(fwSum);
			ercSum = Math.exp(rcSum);

			ret.get(0).add(efwSum);
			ret.get(1).add(ercSum);
			ret.get(2).add(efwSum + ercSum);
		}

		return ret;
	}

	public ArrayList<ArrayList<Double> > slidePN(long[] seqCode, int seqLen, int offset, int p, int nn) {
		ArrayList<ArrayList<Double> > ret = new ArrayList<ArrayList<Double> >(0);
		int length = 0;
		double fwSum = 0, efwSum = 0;
		double rcSum = 0, ercSum = 0;
		double totSum = 0;
		long[] rcCode, fwCode;
		//double ts = 0;
		int[] fbidx = new int[p];
		int[] rbidx = new int[p];
		int fbpos = 0, rbpos = 0;

		if (p > nuc_feat)
			throw new RuntimeException("Nuc feature p is too big!");

		if (nn > nn_feat)
			throw new RuntimeException("Nearest neighbor feature n is too big!");

		for(int i=0; i < 4; i++) {
			ret.add(new ArrayList<Double>());
		}

		length = seqLen + leftFlank.length() + rightFlank.length();
		fwCode = new long[(length + bPL - 1) / bPL];
		rcCode = new long[fwCode.length];

		sc.sequenceWithFlank(seqCode, seqLen, offset, lfCode, leftFlank, rfCode, rightFlank, fwCode, rcCode);
		// loop through all possible windows include Flanks
		//System.out.println("nn_feat = " + nn);
		for (int i = 0; i < length - ws + 1; i++) {
			fwSum = 0;
			rcSum = 0;

			for (int j = 0; j < ws; j++) {
				for (int m = p - 1; m >= 1; m--) {
					if (j >= m) {
						for (int k = nn - 1; k >= 0; k--) {
							if (j + k >= ws)
								continue; //skip the nn that is out of bound

							int shift = (bPL - 1 - (j + i + k) % bPL) * bBits;
							fbpos = (int) ((fwCode[(j + i + k) / bPL] >>> shift) & bitMask);
							rbpos = (int) ((rcCode[(j + i + k) / bPL] >>> shift) & bitMask);

							fbidx[m] = (int) ((fbidx[m - 1] << bBits) | fbpos);
							rbidx[m] = (int) ((rbidx[m - 1] << bBits) | rbpos);
							//System.out.println("j=" + (j - m));
							//System.out.println("Slide " + dbCvtIdx2NucStr(fbidx[m],m+1));
							//System.out.println("Index=" + (betaOffset[m][k][j - m] + fbidx[m]));
							//System.out.println("pBetas[betaOffset[m][k][j]=" + pBetas[betaOffset[m][k][j-m] + fbidx[m]] );
							assert(pBetas[betaOffset[m][k][j - m] + fbidx[m]] == fbidx[m]);
							assert(pBetas[betaOffset[m][k][j - m] + rbidx[m]] == rbidx[m]);
							//since we consider nucleotide from right to left
							//the actual position of di/tri/etc is j - m
							fwSum += pBetas[betaOffset[m][k][j - m] + fbidx[m]];
							rcSum += pBetas[betaOffset[m][k][j - m] + rbidx[m]];
						}
					}
				}
				int shift = (bPL - 1 - (j + i) % bPL) * bBits;
				fbpos = (int) ((fwCode[(j + i) / bPL] >>> shift) & bitMask);
				rbpos = (int) ((rcCode[(j + i) / bPL] >>> shift) & bitMask);
				fbidx[0] = fbpos;
				rbidx[0] = rbpos;

				fwSum += pBetas[betaOffset[0][0][j] + fbidx[0]];
				rcSum += pBetas[betaOffset[0][0][j] + rbidx[0]];
				assert(betaOffset[0][0][j] + fbidx[0] == fbidx[0]);
				assert(betaOffset[0][0][j] + rbidx[0] == rbidx[0]);
				//System.out.println("j=" + j);
				//System.out.println("Slide " + dbCvtIdx2NucStr(fbidx[0],1));
				//System.out.println("Index=" + (betaOffset[0][0][j] + fbidx[0]));
				//System.out.println("pBetas[betaOffset[0][0][j]=" + pBetas[betaOffset[0][0][j] + fbidx[0]]);
			}

			efwSum = Math.exp(fwSum);
			ercSum = Math.exp(rcSum);
			totSum += efwSum + ercSum;

			ret.get(0).add(efwSum);
			ret.get(1).add(ercSum);
			ret.get(2).add(efwSum + ercSum);
		}

		ret.get(3).add(totSum);

		return ret;
	}

	public void swGradient(long[] seqCode, int seqLen, double gradW, ArrayList<ArrayList<Double> > alpha, double[] gradient) {
		
		int length = 0;
		long[] rcCode, fwCode;
//		double ts = 0;
		int[] fbidx = new int[nuc_feat];
		int[] rbidx = new int[nuc_feat];
		int fbpos = 0, rbpos = 0;

		if(gradient.length != pBetas.length) {
			throw new RuntimeException("Error! Gradient length is not correct ! Expected =" + pBetas.length);
		}

		length = seqLen + leftFlank.length() + rightFlank.length();
		fwCode = new long[(length + bPL - 1) / bPL];
		rcCode = new long[fwCode.length];

		sc.sequenceWithFlank(seqCode, seqLen, 0, lfCode, leftFlank, rfCode, rightFlank, fwCode, rcCode);
		// loop through all possible windows include Flanks

		//Loop over windows;
		for (int i = 0; i < length - ws + 1; i++) {
			// Note that (j-m) is the (j-m)-th window
			// betNNLen[x] is the first index of x+1 neighbors
			//Loop first-base-position in window
			for (int j = 0; j < ws; j++) {
				for (int m = nuc_feat - 1; m >= 1; m--) {
					if (j >= m) {
						for (int k = nn_feat - 1; k >= 0; k--) {
							if (j + k >= ws)
								continue;

							int shift = (bPL - 1 - (j + i + k) % bPL) * bBits;
							fbpos = (int) ((fwCode[(j + i + k) / bPL] >>> shift) & bitMask);
							rbpos = (int) ((rcCode[(j + i + k) / bPL] >>> shift) & bitMask);

							fbidx[m] = (int) ((fbidx[m - 1] << bBits) | fbpos);
							rbidx[m] = (int) ((rbidx[m - 1] << bBits) | rbpos);
							//System.out.println("Slide " + dbCvtIdx2NucStr(fbidx[m],m+1));
							//System.out.println("Index=" + (betaOffset[m][k][j - m] + fbidx[m]));
							//System.out.println("alpha=" + alpha.get(0).get(i) * gradW);
							gradient[betaOffset[m][k][j-m] + fbidx[m]] += alpha.get(0).get(i) * gradW;
							gradient[betaOffset[m][k][j-m] + rbidx[m]] += alpha.get(1).get(i) * gradW;
						}
					}
				}
				//Determines the shift in the long.
				int shift = (bPL - 1 - (j + i) % bPL) * bBits;
				fbpos = (int) ((fwCode[(j + i) / bPL] >>> shift) & bitMask);
				rbpos = (int) ((rcCode[(j + i) / bPL] >>> shift) & bitMask);
				fbidx[0] = fbpos;
				rbidx[0] = rbpos;
				//System.out.println("Slide " + dbCvtIdx2NucStr(fbidx[0],1));
				//System.out.println("Index=" + (betaOffset[0][0][j] + fbidx[0]));
				//System.out.println("alpha=" + alpha.get(0).get(i) * gradW);
				gradient[betaOffset[0][0][j] + fbidx[0]] += alpha.get(0).get(i) * gradW;
				gradient[betaOffset[0][0][j] + rbidx[0]] += alpha.get(1).get(i) * gradW;
			}
		}

	}

	public void swHessian(long[] seqCode, int seqLen, double gradW, ArrayList<ArrayList<Double> > alpha, double[][] hessian) {
		int length = 0;
		long[] rcCode, fwCode;
//		double ts = 0;
		int[] fbidx = new int[nuc_feat];
		int[] rbidx = new int[nuc_feat];
		int fbpos = 0, rbpos = 0;
		int totCount = 0;

		length = seqLen + leftFlank.length() + rightFlank.length();
		fwCode = new long[(length + bPL - 1) / bPL];
		rcCode = new long[fwCode.length];

		sc.sequenceWithFlank(seqCode, seqLen, 0, lfCode, leftFlank, rfCode, rightFlank, fwCode, rcCode);

		// count how many elements we have to add
		for (int j = 0; j < ws; j++) {
			for (int m = nuc_feat - 1; m >= 1; m--) {
				if (j >= m) {
					for (int k = nn_feat - 1; k >= 0; k--) {
						if (j + k >= ws)
							continue;

						++totCount;
					}
				}
			}

			++totCount;
		}

		// loop through all possible windows include Flanks
		int[] fTmp = new int[totCount];
		int[] rTmp = new int[totCount];
		
		for (int i = 0; i < length - ws + 1; i++) {

			// Note that (j-m) is the (j-m)-th window
			// betNNLen[x] is the first index of x+1 neighbors
			int count = 0;

			for (int j = 0; j < ws; j++) {
				for (int m = nuc_feat - 1; m >= 1; m--) {
					if (j >= m) {
						for (int k = nn_feat - 1; k >= 0; k--) {
							if (j + k >= ws)
								continue;

							int shift = (bPL - 1 - (j + i + k) % bPL) * bBits;
							fbpos = (int) ((fwCode[(j + i + k) / bPL] >>> shift) & bitMask);
							rbpos = (int) ((rcCode[(j + i + k) / bPL] >>> shift) & bitMask);

							fbidx[m] = (int) ((fbidx[m - 1] << bBits) | fbpos);
							rbidx[m] = (int) ((rbidx[m - 1] << bBits) | rbpos);

							fTmp[count] = betaOffset[m][k][j-m] + fbidx[m];
							rTmp[count++] = betaOffset[m][k][j-m] + rbidx[m];
						}
					}
				}
				int shift = (bPL - 1 - (j + i) % bPL) * bBits;
				fbpos = (int) ((fwCode[(j + i) / bPL] >>> shift) & bitMask);
				rbpos = (int) ((rcCode[(j + i) / bPL] >>> shift) & bitMask);
				fbidx[0] = fbpos;
				rbidx[0] = rbpos;

				fTmp[count] = betaOffset[0][0][j] + fbidx[0];
				rTmp[count++] = betaOffset[0][0][j] + rbidx[0];
			}

			//Score all combination with repition of fTmp and rTmp to Hessian 2D matrix
			for(int a = 0; a < totCount; a++) {
				for(int b=0; b < totCount; b++) {
					hessian[fTmp[a]][fTmp[b]] += alpha.get(0).get(i) * gradW;
					hessian[rTmp[a]][rTmp[b]] += alpha.get(1).get(i) * gradW;
				}
			}
		}
	}
	 
	//Creates a list of gradient indices with non-zero design matrix X for each offset j in the sequence 
	public ArrayList<ArrayList<ArrayList<Integer>>> swNonzeroXList(long[] seqCode, int seqLen) {
		
		ArrayList<ArrayList<Integer>> forwardIndices = new ArrayList<ArrayList<Integer>>();
		ArrayList<ArrayList<Integer>> reverseIndices = new ArrayList<ArrayList<Integer>>();
		int length = 0;
		long[] rcCode, fwCode;
		int[] fbidx = new int[nuc_feat];
		int[] rbidx = new int[nuc_feat];
		int fbpos = 0, rbpos = 0;


		length = seqLen + leftFlank.length() + rightFlank.length();
		fwCode = new long[(length + bPL - 1) / bPL];
		rcCode = new long[fwCode.length];

		sc.sequenceWithFlank(seqCode, seqLen, 0, lfCode, leftFlank, rfCode, rightFlank, fwCode, rcCode);
		// loop through all possible windows include Flanks

		for (int i = 0; i < length - ws + 1; i++) {
			forwardIndices.add(new ArrayList<Integer>());
			reverseIndices.add(new ArrayList<Integer>());

			// Note that (j-m) is the (j-m)-th window
			// betNNLen[x] is the first index of x+1 neighbors
			for (int j = 0; j < ws; j++) {
				
				for (int m = nuc_feat - 1; m >= 1; m--) {
					if (j >= m) {
						for (int k = nn_feat - 1; k >= 0; k--) {
							if (j + k >= ws)
								continue;

							int shift = (bPL - 1 - (j + i + k) % bPL) * bBits;
							fbpos = (int) ((fwCode[(j + i + k) / bPL] >>> shift) & bitMask);
							rbpos = (int) ((rcCode[(j + i + k) / bPL] >>> shift) & bitMask);

							fbidx[m] = (int) ((fbidx[m - 1] << bBits) | fbpos);
							rbidx[m] = (int) ((rbidx[m - 1] << bBits) | rbpos);
							
							forwardIndices.get(i).add(betaOffset[m][k][j-m] + fbidx[m]);
							reverseIndices.get(i).add(betaOffset[m][k][j-m] + rbidx[m]);
						}
					}
				}
				int shift = (bPL - 1 - (j + i) % bPL) * bBits;
				fbpos = (int) ((fwCode[(j + i) / bPL] >>> shift) & bitMask);
				rbpos = (int) ((rcCode[(j + i) / bPL] >>> shift) & bitMask);
				fbidx[0] = fbpos;
				rbidx[0] = rbpos;
				
				forwardIndices.get(i).add(betaOffset[0][0][j] + fbidx[0]);
				reverseIndices.get(i).add(betaOffset[0][0][j] + rbidx[0]);

			}
		}
		ArrayList<ArrayList<ArrayList<Integer>>> out = new ArrayList<ArrayList<ArrayList<Integer>>>();
		out.add(forwardIndices);
		out.add(reverseIndices);
		return out;
	}

	public void setLeftFlank(String lf) {
		leftFlank = lf;
		lfCode = sc.build(lf).getValue();
	}

	public void setRightFlank(String rf) {
		rightFlank = rf;
		rfCode = sc.build(rf).getValue();
	}

	public void setFlanks(String lf, String rf) {
		leftFlank = lf;
		rightFlank = rf;
		lfCode = sc.build(lf).getValue();
		rfCode = sc.build(rf).getValue();
	}

	public ArrayList<ArrayList<Double> > slideSN(long[] seqCode, int seqLen) {
		return slideSN(seqCode, seqLen, 0);
	}

	public ArrayList<ArrayList<Double> > slidePN(long[] seqCode, int seqLen, int p) {
		return slidePN(seqCode, seqLen, 0, p, nn_feat);
	}

	public ArrayList<ArrayList<Double> > slideSN(LongSequence sq, int offset) {
		return slideSN(sq.getValue(), sq.getLength(), offset);
	}

	public ArrayList<ArrayList<Double> > slideSN(LongSequence sq) {
		return slideSN(sq.getValue(), sq.getLength(), 0);
	}

	public ArrayList<ArrayList<Double> > slidePN(LongSequence sq, int offset, int p) {
		return slidePN(sq.getValue(), sq.getLength(), offset, p, nn_feat);
	}

	public ArrayList<ArrayList<Double> > slidePN(LongSequence sq, int p) {
		return slidePN(sq.getValue(), sq.getLength(), 0, p, nn_feat);
	}

	public ArrayList<ArrayList<Double> > slideSN(String s, int offset) {
		LongSequence sq = sc.build(s);

		return slideSN(sq.getValue(), sq.getLength(), offset);
	}

	public ArrayList<ArrayList<Double> > slideSN(String s) {
		LongSequence sq = sc.build(s);

		return slideSN(sq.getValue(), sq.getLength(), 0);
	}

	public ArrayList<ArrayList<Double> > slidePN(String s, int offset, int p) {
		LongSequence sq = sc.build(s);

		return slidePN(sq.getValue(), sq.getLength(), offset, p, nuc_feat);
	}

	public ArrayList<ArrayList<Double> > slidePN(String s, int p) {
		LongSequence sq = sc.build(s);

		return slidePN(sq.getValue(), sq.getLength(), 0, p, nuc_feat);
	}

	/*private String dbCvtIdx2NucStr(long index, int len) {
		StringBuilder sb = new StringBuilder();
		int cvtBase = (int) Math.pow(2, bBits);

		for(int i=0; i < len; i++) {
			sb.insert(0, ALPHABET_LOOKUP_REVERSE_DB.get((int)(index % cvtBase)));
			index /= cvtBase;
		}

		return sb.toString();
	}*/

	//Converts a gradient vector from the internal representation 
	//(which users 2^bBits positions per letter and an internal alphabet ordering) 
	//to the packed and re-ordered external representation.
	public double[] convertGradient(double[] grad, String newAlphaMap) {
		int[] stbMap    = new int[ALPHABET_LOOKUP.size()];
		double[] output = new double[grad.length];

		//Creates map index in internal alphabet to location in external alphabet 
		for (int i = 0; i < newAlphaMap.length(); i++) {
			char c = newAlphaMap.charAt(i);
			if (!ALPHABET_LOOKUP.containsKey(c))
				throw new RuntimeException("Beta mapping string " + newAlphaMap + " cotains invalid character " + c + "!");

			stbMap[ALPHABET_LOOKUP.get(c)] = i;
		}

		int nAlpha    = newAlphaMap.length();                      // Number of letters in the (external) alphabet
		int nInt      = (int) Math.pow(2, bBits);                  // Number of possible letters in the internal long representation.
		int iOutBlock = 0;                                         // Index of the current block in the output vector.
		for (int i = 0; i < nuc_feat; i++) {                       // Loops over k-mer = (monomer, dimer, trimer)
			int rLen = (int) Math.pow(nAlpha, i + 1);              // Number of k-mers 
			int nn   = (i == 0) ? 1 : nn_feat;                     // Maximum gap length (??)
			for (int m = 0; m < nn; m++) {                         // Loops over k-mer gap length.

				// Place the elements of the internal gradient vector (which uses the same ordering as beta)
				// into a new vector that uses the same ordering as the input vector that was used to build beta.
				for (int j = 0; j < winsize_per_feat[i][m]; j++) { // Loops over location of gapped k-mer in window.
					for (int k = 0; k < rLen; k++) {               // Loops over k-mers.
						int kMer     = k;                          // Index of k-mer using packed internal alphabet ordering.
						int intIndex = 0;                          // Index of k-mer using non-packed internal ordering. 
						int cvtIndex = 0;                          // Index of k-mer using external alphabet ordering.
						int intMult  = 1;                          // Current base multiplier in the internal k-mer representation.
						int cvtMult  = 1;                          // Current base multiplier in the external k-mer representation.
						for(int pow=0; pow<=i; pow++) {            // Loops across base power factors
							int rightLetter = kMer % nAlpha;        // Extracts the rightmost letter;
							intIndex += intMult *        rightLetter;  // Adds to the internal index with the appropriate base multiplier
							cvtIndex += cvtMult * stbMap[rightLetter]; // Adds to the converted index with the appropriate base multiplier 
							cvtMult  *= nAlpha;                    // Goes to the next base multiplier
							intMult  *= nInt;                      // Goes to the next base multiplier in the internal representation.
							kMer     /= nAlpha;                    // Removes the rightmost letter in the k-mer, shifts to the right.
						}

						output[iOutBlock + cvtIndex] = grad[betaOffset[i][m][j] + intIndex]; //Copies the current entry to the relevant position in the output vector.
					}
					iOutBlock += rLen;                             //Moves forward in the output vector
				}
			}
		}
		
		return Arrays.copyOfRange(output, 0, iOutBlock);
	}
	
	/*

	//This function has a bug. The corresponding bug was fixed in convertGradient. Fix before using
	public double[][] convertHessian(double[][] h, String newAlphaMap) {
		int[] stbMap = new int[ALPHABET_LOOKUP.size()];
		double[][] output = new double[h.length][h[0].length];

		for (int i = 0; i < newAlphaMap.length(); i++) {
			char c = newAlphaMap.charAt(i);
			if (!ALPHABET_LOOKUP.containsKey(c))
				throw new RuntimeException("Beta mapping string " + newAlphaMap + " cotains invalid character " + c + "!");

			stbMap[ALPHABET_LOOKUP.get(c)] = i;
		}

		int rawBase = newAlphaMap.length();
		int cvtBase = rawBase;

		for (int i1 = 0; i1 < nuc_feat; i1++) {
			int nn1 = (i1 == 0) ? 1 : nn_feat;

			for (int m1 = 0; m1 < nn1; m1++) {
				int rLen1 = (int) Math.pow(rawBase, i1 + 1);

				// place elements of rawBeta array into correct position in
				// internal array
				for (int j1 = 0; j1 < winsize_per_feat[i1][m1]; j1++) {
					for (int k1 = 0; k1 < rLen1; k1++) {
						int rawIndex1 = k1;
						int cvtIndex1 = 0;
						int cvtMult1 = 1;
						int pow1 = i1;

						while (pow1 >= 0) {
							cvtIndex1 += cvtMult1 * (stbMap[rawIndex1 % rawBase]);
							cvtMult1 = cvtMult1 * cvtBase;
							rawIndex1 = rawIndex1 / rawBase;
							pow1--;
						}

						for (int i2 = 0; i2 < nuc_feat; i2++) {
							int nn2 = (i2 == 0) ? 1 : nn_feat;

							for (int m2 = 0; m2 < nn2; m2++) {
								int rLen2 = (int) Math.pow(rawBase, i2 + 1);

								// place elements of rawBeta array into correct position in
								// internal array
								for (int j2 = 0; j2 < winsize_per_feat[i2][m2]; j1++) {
									for (int k2 = 0; k2 < rLen2; k2++) {
										int rawIndex2 = k2;
										int cvtIndex2 = 0;
										int cvtMult2 = 1;
										int pow2 = i2;

										while (pow2 >= 0) {
											cvtIndex2 += cvtMult2 * (stbMap[rawIndex2 % rawBase]);
											cvtMult2 = cvtMult2 * cvtBase;
											rawIndex2 = rawIndex2 / rawBase;
											pow2--;
										}

										output[betaOffset[i1][m1][j1] + cvtIndex1][betaOffset[i2][m2][j2] + cvtIndex2] 
											= h[betaOffset[i1][m1][j1] + cvtIndex1][betaOffset[i2][m2][j2] + cvtIndex2];
									}
								}
							}
						}
					}
				}
			}
		}

		return output;
	}*/
}