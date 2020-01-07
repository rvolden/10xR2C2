#!/usr/bin/env python3
# Roger Volden

'''
This will generate a cell to celltype file with the following format:
    cell number \t cell type \t cell barcode

The cell number is taken from the bcGuide, which is produced during
the demuxing process.
The cell type will need to be determined manually based on expression.
The cell barcode is also taken from the bcGuide file.
Coordinates will need to be determined manually for the cell type
boundaries.

Usage:
    python3 cell_to_celltype.py bcGuide 42_*coords
'''

import sys

def readBCGuide(bcGuide):
    '''{barcode:cell#}'''
    bcDict = {}
    for line in open(bcGuide):
        line = line.rstrip().split()
        bcDict[line[1]] = line[0]
    return bcDict

def readCoordsAndSort(bcDict, embeddings):
    print('# cell\ttype\tbarcode')
    for line in open(embeddings):
        line = line.rstrip().split()
        bc, x, y = line[0], float(line[1]), float(line[2])
        cellType = ''

        if x > 30:
            cellType = 'bCell'
        elif x > 7 and y > 14:
            cellType = 'monocyte'
        else:
            cellType = 'tCell'

        print(bcDict[bc] + '\t' + cellType + '\t' + bc)

def main():
    bcDict = readBCGuide(sys.argv[1])
    readCoordsAndSort(bcDict, sys.argv[2])

main()
