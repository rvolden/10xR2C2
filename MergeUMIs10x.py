import argparse
import editdistance
import numpy as np
import os
import sys
import multiprocessing as mp
import subprocess

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input_path', type=str, help='This should point to the folder that contains your demultiplexed R2C2 fasta and subread (ending on _subs.fastq) files for each cell.')
parser.add_argument('-o', '--output_path', type=str, help='Merged fasta and subread files will be written to this folder')
parser.add_argument('-c', '--config_file', type=str, help='Same config file used for C3POa')
parser.add_argument('-m', '--score_matrix', type=str, help='Same matrix file used for C3POa')
parser.add_argument('-t', '--threads', type=str, help='defines the number of threads the multiprocessing will use')

args = parser.parse_args()
path = args.output_path+'/'
temp = path+'/temp'
os.system('mkdir '+temp)
input_path = args.input_path+'/'

config_file = args.config_file
score_matrix = args.score_matrix
threads = int(args.threads)
subsample = 200

def reverse_complement(sequence):
  Seq = ''
  complement = {'A':'T', 'C':'G', 'G':'C', 'T':'A', 'N':'N', '-':'-'}
  for item in sequence[::-1]:
    Seq = Seq + complement[item]
  return Seq

def configReader(configIn):
    '''Parses the config file.'''
    progs = {}
    for line in open(configIn):
        if line.startswith('#') or not line.rstrip().split():
            continue
        line = line.rstrip().split('\t')
        progs[line[0]] = line[1]
    # should have minimap, poa, racon, water, consensus
    # check for extra programs that shouldn't be there
    possible = set(['poa', 'minimap2', 'gonk', 'consensus', 'racon', 'blat','emtrey', 'psl2pslx'])
    inConfig = set()
    for key in progs.keys():
        inConfig.add(key)
        if key not in possible:
            raise Exception('Check config file')
    # check for missing programs
    # if missing, default to path
    for missing in possible-inConfig:
        if missing == 'consensus':
            path = 'consensus.py'
        else:
            path = missing
        progs[missing] = path
        sys.stderr.write('Using ' + str(missing)
                         + ' from your path, not the config file.\n')
    return progs

progs = configReader(config_file)
poa = progs['poa']
minimap2 = progs['minimap2']
racon = progs['racon']
consensus = progs['consensus']


def determine_consensus(fasta, fastq, temp_folder, process_count):
    '''Aligns and returns the consensus'''
    corrected_consensus = ''
    pc = process_count
    out_F = fasta
    fastq_reads = read_fastq_file(fastq)
    out_Fq = temp_folder + '/subsampled.' + pc + '.fastq'
    out = open(out_Fq,'w')
    indexes = np.random.choice(np.arange(0, len(fastq_reads), 1), min(len(fastq_reads), subsample), replace=False)
    subsample_fastq_reads = []
    for index in indexes:
        subsample_fastq_reads.append(fastq_reads[index])

    for read in subsample_fastq_reads:
        out.write('@' + read[0] + '\n' + read[2] + '\n+\n' + read[3] + '\n')
    out.close()

    poa_cons = temp_folder + '/consensus.' + pc + '.fasta'
    final = temp_folder + '/corrected_consensus.' + pc + '.fasta'
    overlap = temp_folder +'/overlaps.' + pc + '.sam'
    pairwise = temp_folder + '/prelim_consensus.' + pc + '.fasta'

    max_coverage = 0
    reads = read_fasta(out_F)
    repeats = 0
    qual = []
    raw = []
    before = []

    after = []
    for read in reads:
        best = read
    out_cons_file = open(poa_cons, 'w')
    out_cons_file.write('>' + best + '\n' + reads[best].replace('-', '') + '\n')
    out_cons_file.close()

    final = poa_cons
    for i in np.arange(1, 2, 1):
        if i == 1:
            input_cons = poa_cons
            output_cons = poa_cons.replace('.fasta', '_' + str(i) + '.fasta')
        else:
            input_cons = poa_cons.replace('.fasta', '_' + str(i-1) + '.fasta')
            output_cons = poa_cons.replace('.fasta', '_' + str(i) + '.fasta')

        minimap2_command = '%s -t 1 --secondary=no -ax map-ont \
                            %s %s >%s 2>./minimap2_messages' \
                            %(minimap2, input_cons, out_Fq,overlap)
        minimap2_process = subprocess.run(minimap2_command, shell=True)

        racon_command = '%s -q 5 -t 1 %s %s %s >%s 2>./racon_messages.txt' \
                         %(racon, out_Fq, overlap, input_cons, output_cons)
        racon_process = subprocess.run(racon_command, shell=True)
        final = output_cons

    reads = read_fasta(final)
    if len(reads) == 0:
        reads = read_fasta(poa_cons)
    for read in reads:
        corrected_consensus = reads[read]

    return corrected_consensus

def read_fasta(infile):
    reads = {}
    sequence = ''

    for line in open(infile):
      if line:
        a = line.strip()
        if len(a) > 0:
            if a[0] == '>':
                if sequence != '':
                    reads[name] = sequence
                name = ('-').join(a[1:].split('_')[0].split('-')[:5])
                sequence = ''
            else:
                sequence += a
    if sequence != '':
        reads[name] = sequence
    return reads

def read_subreads(seq_file):
    lineNum = 0
    lastPlus  =  False
    read_dict = {}
    for line in open(seq_file):
        line = line.rstrip()
        # make an entry as a list and append the header to that list
        if lineNum % 4 == 0:
            if line[0] == '@':
                if lastPlus:
                    if root_name not in read_dict:
                        read_dict[root_name] = []  # chrom_reads needs to contain root_names
                    read_dict[root_name].append((root_name,seq,qual))
                name = line[1:]
                root_name = ('-').join(name.split('_')[0].split('-')[:5])

        if lineNum % 4 == 1:
            seq = line
        if lineNum % 4 == 2:
            lastPlus = True
        if lineNum % 4 == 3 and lastPlus:
            qual = line

        lineNum += 1
    read_dict[root_name].append((root_name, seq, qual))
    return read_dict

def read_fastq_file(seq_file):
    '''
    Takes a FASTQ file and returns a list of tuples
    In each tuple:
        name : str, read ID
        seed : int, first occurrence of the splint
        seq : str, sequence
        qual : str, quality line
        average_quals : float, average quality of that line
        seq_length : int, length of the sequence
    '''

    lineNum = 0
    lastPlus = False
    read_list = []
    for line in open(seq_file):
        line = line.rstrip()
        if not line:
            continue
        # make an entry as a list and append the header to that list
        if lineNum % 4 == 0 and line[0] == '@':
            if lastPlus == True:
                read_list.append((name, '', seq, qual, average_quals, seq_length))
            name = line[1:]
        if lineNum % 4 == 1:
            seq = line
            seq_length = len(seq)
        if lineNum % 4 == 2:
            lastPlus = True
        if lineNum % 4 == 3 and lastPlus:
            qual = line
            quals = []
            for character in qual:
                number = ord(character) - 33
                quals.append(number)
            average_quals = np.average(quals)
        lineNum += 1
    return read_list

def make_consensus(Molecule, subreads, process_count):
    subread_file = temp + '/' + process_count + '.temp_subreads.fastq'
    fastaread_file = temp + '/' + process_count + '.temp_consensusreads.fasta'
    subs = open(subread_file, 'w')
    fasta = open(fastaread_file, 'w')
    raw_count = 0
    combined_root_name = list(Molecule)[0][1:].split('\n')[0]
    for read in Molecule:
        fasta.write(read)
        root_name = read[1:].split('\n')[0]
        if root_name in subreads:
            raw = subreads[root_name]
            for entry in raw:
                if entry[1]:
                    raw_count += 1
                    subs.write('@' + combined_root_name + '_' + str(raw_count) + '\n' + entry[1] + '\n+\n' + entry[2] + '\n')
    subs.close()
    fasta.close()
    if len(read_fastq_file(subread_file)) > 0:
        corrected_consensus = determine_consensus(fastaread_file, subread_file, temp, process_count)
        return '>%s\n%s\n' %(combined_root_name, corrected_consensus)
    else:
        return list(Molecule)[0].split('\n')[0] + '~N\n' + list(Molecule)[0].split('\n')[1] + '\n'

def group_reads(reads, subreads, final, final_UMI_only, final_subreads, matched_reads, process_count):
    UMI_dict = {}
    unique_number = 0
    for read,sequence in reads.items():
        umi = reverse_complement(sequence)[56:66]
        if umi[4:10] == 'TTTTTT':
            unique_number += 1
            umi = str(unique_number)
        if umi not in UMI_dict:
            UMI_dict[umi] = []
        UMI_dict[umi].append('>' + read + '\n' + sequence + '\n')

    for umi, Molecule in UMI_dict.items():
        if len(Molecule) == 1:
            final.write(list(Molecule)[0])
            if list(Molecule)[0][1:].split('\n')[0] in subreads:
                for subread in subreads[list(Molecule)[0][1:].split('\n')[0]]:
                    root_name, sequence, qual = subread[0], subread[1], subread[2]
                    final_subreads.write('@' + root_name + '\n' + sequence + '\n+\n' + qual + '\n')
        elif len(Molecule) > 1:
            matched_reads.write(umi + '\t' + str(umi.count('TT')) + '\t')
            for entry in list(Molecule):
                matched_reads.write(entry[1:].split('\n')[0] + '_' + str(len(entry.split('\n')[1])) + ',')
            matched_reads.write('\n')
            new_read = make_consensus(list(Molecule), subreads, process_count)
            if new_read != 'nope':
                final.write(new_read)
                final_UMI_only.write(new_read)
                combined_root_name = new_read[1:].split('\n')[0]
                for molecule in list(Molecule):
                    name = molecule.split('\n')[0][1:]
                    if name in subreads:
                        for subread in subreads[name]:
                            root_name, sequence, qual = subread[0], subread[1], subread[2]
                            final_subreads.write('@' + combined_root_name + '\n' + sequence + '\n+\n' + qual + '\n')

def processing(fasta_file, subreads_file, process_count):
    final = open(path + '/' + process_count + '.merged.fasta', 'w')
    final_UMI_only = open(path + '/' + process_count + '.UMI_only.fasta', 'w')
    final_subreads = open(path + '/' + process_count + '.merged.subreads.fastq', 'w')
    matched_reads = open(path + '/' + process_count + '.matched_reads.txt', 'w')
    print('reading reads required for subprocess ' + process_count)
    reads = read_fasta(fasta_file)
    print(len(reads))
    print('reading subreads ' + process_count)
    subreads = read_subreads(subreads_file)
    print(len(subreads))
    print('grouping and merging consensus reads ' + process_count)
    group_reads(reads, subreads, final, final_UMI_only, final_subreads, matched_reads, process_count)

def main():
    pool = mp.Pool(processes=threads)
    print('kmer-matching UMIs')
    print('reading consensus reads')
    count = 0
    big_count = 0
    sub_reads = {}
    fileList = []
    for file in os.listdir(input_path):
        if '.fasta' in file and 'cell' in file and 'merged' not in file:
            fileList.append(file)
    for file in sorted(fileList, key=lambda x: int(x.split('_')[1])):
        fasta_file = input_path + '/' + file
        sub_reads = input_path + '/' + file.split('.')[0] + '_subs.fastq'
        pool.apply_async(processing, [fasta_file, sub_reads, file.split('.')[0]])
        count = 0
        sub_reads = {}
    pool.close()
    pool.join()

main()
