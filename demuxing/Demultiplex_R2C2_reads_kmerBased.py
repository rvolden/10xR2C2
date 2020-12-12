#!/usr/bin/env python3
# Roger Volden

'''
Adds a pipe to the end of the sequence header and puts the header
of the barcode that it matched. This will later be parsed to demux
the full length R2C2 read.

python3 Demultiplex_R2C2_reads_kmerBased.py \
    --input_fasta_file R2C2_10x_postprocessed.fasta \
    --output_path . \
    --tenx_index_file 1500_most_frequent_bcs.fasta
'''

import sys
import argparse
import editdistance as ld
from tqdm import tqdm

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input_fasta_file', type=str)
parser.add_argument('-o', '--output_path', type=str)
parser.add_argument('-n', '--tenx_index_file', type=str)

args = parser.parse_args()
output_path = args.output_path + '/'
input_file = args.input_fasta_file
bc_file = args.tenx_index_file

def read_fasta(inFile, bc):
    readDict = {}
    bcKmerDict = {}
    for line in open(inFile):
        line = line.rstrip()
        if not line:
            continue
        if line.startswith('>'):
            readDict[line[1:]] = ''
            lastHead = line[1:]
        else:
            readDict[lastHead] += line
            if bc:
                kmers = get_kmer_list(line)
                for kmer in kmers:
                    if not bcKmerDict.get(kmer):
                        bcKmerDict[kmer] = [(lastHead, readDict[lastHead])]
                    else:
                        bcKmerDict[kmer].append((lastHead, readDict[lastHead]))
    if bc:
        return bcKmerDict
    else:
        return readDict

def reverse_complement(sequence):
    bases = {'A':'T', 'C':'G', 'G':'C', 'T':'A', 'N':'N', '-':'-'}
    return ''.join([bases[x] for x in list(sequence)])[::-1]

def get_kmer_list(bc):
    kmer_list = [bc[:4] + '___',
                 '_' + bc[4:8] + '__',
                 '__' + bc[8:12] + '_',
                 '___' + bc[12:16]]
    return kmer_list

def demultiplex(reads, bcKmerDict):
    indexed_reads, counter = {}, 0
    for read, complete_sequence in tqdm(reads.items()):
        # if counter % 10000 == 0:
        #     print(str(counter) + ' of ' + str(len(reads)))
        sequence = complete_sequence[:16]

        bc_set = set()
        readKmers = get_kmer_list(sequence)
        for kmer in readKmers:
            if bcKmerDict.get(kmer):
                barcodes = bcKmerDict[kmer] # list of tuples
                for b in barcodes:
                    bc_set.add(b) # b = (header, bc sequence)
        dists, bc_match, maxDist = [], '', 3
        for barcode in bc_set:
            dist = ld.eval(barcode[1], sequence)
            dists.append((barcode, dist))

        if dists:
            sorted_bc_list = sorted(dists, key=lambda x: x[1])

            if sorted_bc_list[0][1] < maxDist and \
               sorted_bc_list[0][1] < sorted_bc_list[1][1]-1:
                bc_match = sorted_bc_list[0][0][0]

            read += '|' + bc_match
            indexed_reads[read] = complete_sequence
        counter += 1
    return indexed_reads

def write_fasta_file(path, reads):
    out = open(path + '/kmer_demuxed.fasta', 'w+')
    for name, sequence in reads.items():
        out.write('>%s\n%s\n' %(name, sequence))

def main():
    bcKmerDict = read_fasta(bc_file, True)
    reads = read_fasta(input_file, False)
    indexed_reads = demultiplex(reads, bcKmerDict)
    write_fasta_file(output_path, indexed_reads)

main()
