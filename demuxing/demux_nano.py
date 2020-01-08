#!/usr/bin/env python3
# Roger Volden

'''
Use to demultiplex R2C2 reads based on the 1500 most
common barcodes.

This will need to read in the 1500 barcodes and the
R2C2 reads. R2C2 reads will be divided into two files:
labeled_10x_sequences.fasta and restOfTheRead.fasta.
Labeled 10x sequences will only contain barcode info,
while the rest of the read is in the other file. I can
read these in at the same time and demux them accordingly.

General usage:
    python3 demux_nano.py bc.fasta labeled_10x_sequences.fasta restOfTheRead.fasta
'''

import sys, os

def readBC(inFile):
    '''
    Read in the barcode fasta file and return a dict:
    readDict = {"cell_#" : "barcode_sequence"}
    '''
    readDict = {}
    for line in open(inFile):
        line = line.rstrip()
        if line[0] == '>':
            cellNum = line[1:].split('_')[1]
            c = 'cell_' + cellNum
            readDict[c] = ''
        else:
            readDict[c] += line
    return readDict

def readFasta(fq1, fq2, bcDict):
    '''
    Reads through two fasta files at the same time.
    I can do this because the first fasta contains the 10X
    sequences, while I have the entire consensus sequence
    in the other file.

    Doesn't return anything, just writes the appropriate
    sequences to each cell's file.
    '''
    readDict, lastHead = {}, ''
    for lineA, lineB in zip(open(fq1), open(fq2)):
        lineA, lineB = lineA.rstrip(), lineB.rstrip()
        if not lineA or not lineB:
            continue
        if lineA[0] == '>' and lineB[0] == '>':
            if lineA.split('|')[0] != lineB:
                print('A' + lineA.split('|')[0] + '\n' + lineB)
            barcode = lineA[1:].split('|')[-1]
            lastHead = lineA[1:].split('|')[0]
            if not barcode:
                readDict[lastHead] = ['del']
            if lastHead not in readDict:
                readDict[lastHead] = []
        elif len(readDict[lastHead]) == 1:
            del readDict[lastHead]
        else:
            cell = 'cell_' + barcode.split('_')[1]
            out = os.getcwd() + '/demuxed/' + cell + '_' + bcDict[cell] + '.fasta'
            writing = open(out, 'a+')
            writing.write(">" + lastHead + "\n" + lineB + "\n")
            writing.close()

def main():
    bcDict = readBC(sys.argv[1])

    # write out a barcode guide to make sure I can easily
    # go back and forth between cell # and barcode sequence
    # so I can match cells between datasets.
    out = os.getcwd() + '/demuxed/bcGuide'
    guide = open(out, 'w+')
    for k, v in bcDict.items():
        guide.write(k + '\t' + v + '\n')
    guide.close()
    sys.stderr.write("Finished reading BC file\n")

    readFasta(sys.argv[2], sys.argv[3], bcDict)
    sys.stderr.write("Finished reading FASTA\n")

if __name__ == '__main__':
    main()
