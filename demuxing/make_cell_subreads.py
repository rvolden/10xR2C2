#!/usr/bin/env python3
# Roger Volden

'''
Take the demuxed fasta files and make dictionary like this:
{header:cell_#_barcode.fasta}
Then read through the massive subreads file to be able to
query each header in the dictionary to get the right cell.

Usage:
    python3 make_cell_subreads.py path/to/demuxed/fastas all_subreads.fastq outDir
'''

import os
import sys
from progress.bar import Bar

def readFasta(inFile, headerDict):
    '''Add to a dictionary of {header: cell ID}'''
    with open(inFile) as f:
        for line in f:
            line = line.rstrip()
            if not line:
                continue
            if line[0] == '>':
                # {header:cell_#_barcode}
                headerDict[line[1:].split('_')[0]] = inFile.split('.')[0].split('/')[-1]
    return headerDict

def filterSubreads(headerDict, subreads, outDir):
    lineNum = 0
    header, seq, quality = '', '', ''
    print('Starting to read in the subreads')
    for line in open(subreads):
        line = line.rstrip()
        if not line: continue
        if lineNum % 4 == 0 and line[0] == '@':
            if header:
                if header in headerDict:
                    outSub = open(outDir + '/' + headerDict[header] + '_subs.fastq', 'a+')
                    outSub.write('@' + header + '\n')
                    outSub.write(seq + '\n+\n')
                    outSub.write(quality + '\n')
                    outSub.close()
                    print(headerDict[header])
            header = line[1:].split('_')[0]
        if lineNum % 4 == 1:
            seq = line
        if lineNum % 4 == 3:
            quality = line
        lineNum += 1
    if header:
        if header in headerDict:
            outSub = open(outDir + '/' + headerDict[header] + '_subs.fastq', 'a+')
            outSub.write('@' + header + '\n')
            outSub.write(seq + '\n+\n')
            outSub.write(quality + '\n')
            outSub.close()
            print(headerDict[header])

def main():
    fastaDir, subreads = sys.argv[1], sys.argv[2]
    outDir = sys.argv[3]
    # fileList has all of the demuxed reads by cell
    fileList = [fastaDir + '/' + x for x in os.listdir(fastaDir) if x.endswith('.fasta') and x.startswith('cell_')]
    headerDict = {}
    b = Bar('Reading in fasta files', max=len(fileList))
    for file in fileList:
        headerDict = readFasta(file, headerDict)
        b.next()
    b.finish()
    filterSubreads(headerDict, subreads, outDir)

if __name__ == '__main__':
    main()
