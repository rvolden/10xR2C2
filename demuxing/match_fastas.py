#!/usr/bin/env python3
# Roger Volden

'''
Used to match up the fasta entries between the demux fasta and the
full length consensi. This is because when demuxing, sometimes reads
get dropped, which throws off the script that sorts reads into cells.

Assumes the demux fasta file has fewer lines than the flc file.

Usage:
python3 match_fastas.py \
    kmer_demuxed.fasta \
    R2C2_postprocessed.fasta \
    >R2C2_matched.fasta
'''

import sys

def readFasta(inFile):
    readDict = {}
    for line in open(inFile):
        line = line.rstrip()
        if not line:
            continue
        if line.startswith('>'):
            lastHead = line[1:].split('_')[0]
            readDict[lastHead] = ''
        else:
            readDict[lastHead] += line
    return readDict

def readDemux(inFile, flcDict):
    for line in open(inFile):
        line = line.rstrip()
        if line.startswith('>'):
            header = line[1:].split('_')[0]
            print('>' + line[1:].split('|')[0])
            print(flcDict[header])

def main():
    demux, flc = sys.argv[1], sys.argv[2]

    flcDict = readFasta(flc)

    readDemux(demux, flcDict)

if __name__ == '__main__':
    main()
