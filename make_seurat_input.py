#!/usr/bin/env python3
# Roger Volden

'''
Copy the cellranger output format for Seurat.
cellranger output:
    hg38/
        barcodes.tsv
        genes.tsv
        matrix.mtx

File formats:
    barcodes.tsv
        one barcode per line + '-1'
    genes.tsv
        gene id \t gene name
    matrix.mtx
        %%MatrixMarket matrix coordinate real general
        %
        # genes  # cells  # of nonzero entries
        gene index  cell index  count

Usage:
    python3 make_seurat_input.py
        -a gencode.v29.annotation.gtf
        -b bcGuide
        -e featureCountsOutput
        -o path/to/output
'''

import os
import sys
import argparse

def argParser():
    '''Parses arguments'''
    parser = argparse.ArgumentParser(description = 'Copies the cellranger output for use with Seurat.',
                                     add_help = True,
                                     prefix_chars = '-')
    parser.add_argument('-a', type=str, action='store',
            help='Gencode annotation file (gtf). eg. gencode.v29.annotation.gtf')
    parser.add_argument('-b', type=str, action='store',
            help='bcGuide file produced from demuxing.')
    parser.add_argument('-e', type=str, action='store',
            help='Output file from featureCounts.')
    parser.add_argument('-o', type=str, action='store',
            help='Output path. eg. hg38')
    return vars(parser.parse_args())

args = argParser()
anno = args['a']
bcGuide = args['b']
exp = args['e']
outPath = args['o']
if not outPath:
    outPath = ""
if not outPath.endswith('/'):
    outPath += '/'

def makeGenes(inFile):
    '''
    Takes a gtf file and returns a list of tuples (gene ID, gene name, index)
    genes = [(gene ID, name, index), ...]
    This list of genes is getting returned because the indeces are needed to make
    the matrix.mtx file.

    Also writes out the ID and name information to genes.tsv
    '''
    idDict, geneDict = {}, {}
    for line in open(inFile):
        if line[0] == '#':
            continue
        line = line.rstrip().split('\t')
        line = [x.strip() for x in line[8].split(';')]
        id, name = '', ''
        for section in line:
            if section[:7] == 'gene_id':
                id = section.split('"')[1]
            if section[:9] == 'gene_name':
                name = section.split('"')[1]
            if id and name:
                break
        idDict[id] = name

    genesOut = open(outPath + 'genes.tsv', 'w+')
    index = 1
    sys.stderr.write('Writing genes to {}genes.tsv\n'.format(outPath))
    for gID, gName in idDict.items():
        genesOut.write(gID + '\t' + gName + '\n')
        geneDict[gID] = index
        index += 1
    genesOut.close()
    return geneDict

def makeBC(inFile):
    '''
    Makes the barcodes.tsv file
    Literally just adds a '-1' to the end of the barcode
    '''
    numCells = 0
    bcOut = open(outPath + 'barcodes.tsv', 'w+')
    sys.stderr.write('Writing barcodes to {}barcodes.tsv\n'.format(outPath))
    for line in open(inFile):
        line = line.rstrip().split('\t')
        bcOut.write(line[1] + '-1\n')
        numCells += 1
    bcOut.close()
    return numCells

def modifyFC(fcIn):
    '''
    Reads the featureCounts output and throws it into a dictionary
    countDict = {geneID: [0, 1, 0, 0, ...], ...}
    The index of the list +1 is the cell #
    '''
    first = True
    countDict = {}
    sys.stderr.write('Reading in the featureCounts output\n')
    for line in open(fcIn):
        if line.startswith('#'): # featureCounts command info
            continue
        if first: # featureCounts column names
            first = False
            continue
        line = line.rstrip().split('\t')
        gene = line[0]
        counts = line[6:] # only need the count columns
        countDict[gene] = counts
    return countDict

def makeMtx(countDict, geneDict, numCells):
    '''
    Makes the matrix.mtx file
    countDict = {geneID: [0, 1, 0, 0, ...]} where the index of the list +1 is the cell #
    geneDict = {geneID: index, ...}
    '''
    nonZero = 0
    for countList in countDict.values():
        nonZero += len([x for x in countList if int(x) > 0])

    mtxOut = open(outPath + 'matrix.mtx', 'w+')
    sys.stderr.write('Writing matrix to {}matrix.mtx\n'.format(outPath))
    mtxOut.write("%%MatrixMarket matrix coordinate real general\n%\n")
    mtxOut.write("{0} {1} {2}\n".format(str(len(geneDict)), str(numCells), str(nonZero)))

    for gene, counts in countDict.items():
        for i in range(len(counts)):
            count = counts[i]
            cell = str(i + 1)
            if count == '0':
                continue
            else:
                mtxOut.write(str(geneDict[gene]) + ' ' + cell + ' ' + count + '\n')

def checkArgs():
    '''
    Prepares environment for writing the output
    and makes sure everything is working
    '''
    if not anno:
        sys.stderr.write('Missing annotation input (-a).\n')
        sys.exit(1)
    if not bcGuide:
        sys.stderr.write('Missing bcGuide input (-b).\n')
        sys.exit(1)
    if not exp:
        sys.stderr.write('Missing featureCounts expression input (-e).\n')
        sys.exit(1)
    if len(outPath) < 2:
        sys.stderr.write('Missing path to output (-o).\n')
        sys.exit(1)
    if not os.path.exists(outPath):
        sys.stderr.write('Output dir does not exist, making it.\n')
        sys.stderr.write('mkdir {}\n'.format(outPath))
        os.mkdir(outPath)

def main():
    checkArgs()
    # first make the genes.tsv
    geneDict = makeGenes(anno)

    # then make the barcodes.tsv
    numCells = makeBC(bcGuide)

    # change the featureCounts output to something easier to parse
    countDict = modifyFC(exp)

    # make the matrix.mtx file
    makeMtx(countDict, geneDict, numCells)

if __name__ == '__main__':
    main()
