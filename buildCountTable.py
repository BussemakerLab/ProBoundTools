#!/usr/bin/env python

import os.path
import sys
import argparse
import gzip
import random as rn
import numpy as np

############ VARIOUS FUNCTION ########

def err(msg):
	sys.stderr.write("ERROR: %s\n"%msg)
	sys.exit(1)

def disp(msg, verbose):
	if verbose:
		sys.stdout.write("%s\n"%msg)
	return

############## MAIN CODE #############

def main():

	#Parsing arguments
	##################
	parser = argparse.ArgumentParser(description='Combines multiple lists of sequences into a count table')
	parser.add_argument('-i', metavar='seq.txt',   help='Input sequences.',                  required=True, nargs="+" )
	parser.add_argument('-m', metavar='max_reads', help='Maximum number of sequences',       required=False, type=int, default=None)
	parser.add_argument('-o', metavar='output',    help='Output file (instead of pipeline)', required=False )
	parser.add_argument("--gzip",    help="Input and output is gzipped", action="store_true")
	parser.add_argument("--verbose", help="Increase output verbosity", action="store_true")
	args  = parser.parse_args()

	if args.m is None:
		nMax = np.inf
	else:
		nMax = args.m

	nColumns = len(args.i)
	countTable = {}           #Count table: [probe][round]

	for r in range(nColumns):

		#Streams all reads and selects nMax probes
		seqFile   = args.i[r]
		probeList = []
		iLine     = 0

		#Opens the input file
		if args.gzip:
			try:
				fIn = gzip.open(seqFile)
			except FileNotFoundError:
				err("The gzip file %s cannot be found"%seqFile)
		else:
			try:
				fIn = open(seqFile)
			except FileNotFoundError:
				err("The file %s cannot be found."%seqFile)

		for fl in fIn:
			#Randomly draws nMax reads from the full dataset
			if iLine<nMax:
				#Records all of the first nMax probes
				if args.gzip:
					probeList += [ fl.decode().rstrip() ]	
				else:
					probeList += [ fl.rstrip() ]	
			else:
				#After that, the new probe is kept with probability nMax/(iLine+1) and one random probe is removed from probeList
				if rn.random() < nMax / (iLine+1):
					probeList[rn.randint(0,nMax-1)] = fl.rstrip()
			iLine+=1

		#Closes the input file
		fIn.close()

		#Builds count table
		for probe in probeList:
			if probe not in countTable:
				countTable[probe] =  [0] * nColumns
			countTable[probe][r] += 1

	#Determines the output file:			
	if args.o is None:
		for probe in sorted(countTable.keys()):
			sys.stdout.write(("%s\t%s\n"%(probe, "\t".join([ "%d"%c for c in countTable[probe]]))))

	else:
		#Opens the output file
		if args.gzip:
			try:
				fOut = gzip.open(args.o, 'w')
			except FileNotFoundError:
				err("Cannot open the gzip output file %s"%args.o)
		else:
			try:
				fOut =    open(args.o, 'w')
			except FileNotFoundError:
				err("Cannot open the output file: %s"%args.o)

		#Writes table.
		disp('Writes to file %s'%args.o, args.verbose)
		for probe in sorted(countTable.keys()):
			if args.gzip:
				fOut.write(("%s\t%s\n"%(probe, "\t".join([ "%d"%c for c in countTable[probe]]))).encode())
			else:
				fOut.write(("%s\t%s\n"%(probe, "\t".join([ "%d"%c for c in countTable[probe]]))))

		fOut.close()

##################  END ################

main()
