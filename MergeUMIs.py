import argparse
import editdistance
import numpy as np
import os
import sys
import mappy as mm
from tqdm import tqdm

parser = argparse.ArgumentParser()
parser.add_argument('-f', '--fasta_file', type=str)
parser.add_argument('-s', '--subreads_file', type=str)
parser.add_argument('-o', '--output_path', type=str)
parser.add_argument('-u', '--umi_file', type=str)
parser.add_argument('-c', '--config_file', type=str)
parser.add_argument('-m', '--score_matrix', type=str)

args = parser.parse_args()
path = args.output_path + '/'
fasta_file = args.fasta_file
subreads_file = args.subreads_file

umi_file = args.umi_file
config_file= args.config_file
score_matrix = args.score_matrix
subsample = 200

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

def determine_consensus(name, fasta, fastq, temp_folder):
    '''Aligns and returns the consensus'''
    corrected_consensus = ''
    out_F = fasta
    fastq_reads = read_fastq_file(fastq, False)
    out_Fq = temp_folder + '/subsampled.fastq'
    out = open(out_Fq, 'w')
    indexes = np.random.choice(np.arange(0, len(fastq_reads), 1), min(len(fastq_reads), subsample), replace=False)

    subsample_fastq_reads = []
    for index in indexes:
        subsample_fastq_reads.append(fastq_reads[index])

    for read in subsample_fastq_reads:
        out.write('@' + read[0] + '\n' + read[1] + '\n+\n' + read[2] + '\n')
    out.close()

    poa_cons = temp_folder + '/consensus.fasta'
    final = temp_folder + '/corrected_consensus.fasta'
    overlap = temp_folder + '/overlaps.sam'
    pairwise = temp_folder + '/prelim_consensus.fasta'

    max_coverage, repeats = 0, 0
    reads = read_fasta(out_F)
    qual, raw, before, after = [], [], [], []

    header_log = path + 'header_associations.tsv'
    header_fh = open(header_log, 'a+')
    headers = []
    for read in reads:
        info = read.split('_')
        # print(info)
        coverage = int(info[3])
        headers.append(info[0])
        qual.append(float(info[1]))
        raw.append(int(info[2]))
        repeats += int(info[3])
        before.append(int(info[4]))
        after.append(int(info[5].split('|')[0]))

        if coverage >= max_coverage:
            best = read
            max_coverage = coverage

    print('\t'.join(headers), file=header_fh)
    header_fh.close()

    out_cons_file = open(poa_cons, 'w')
    out_cons_file.write('>' + best + '\n' + reads[best].replace('-', '') + '\n')
    out_cons_file.close()

    final = poa_cons
    input_cons = poa_cons
    output_cons = poa_cons.replace('.fasta', '_1.fasta')

    os.system('%s --secondary=no -ax map-ont \
              %s %s > %s 2> ./minimap2_messages.txt'
              % (minimap2, input_cons, out_Fq, overlap))
    os.system('%s -q 5 -t 1 \
               %s %s %s >%s 2>./racon_messages.txt'
               %(racon, out_Fq, overlap, input_cons, output_cons))
    final = output_cons

    # print(final)
    reads = read_fasta(final)
    for read in reads:
        corrected_consensus = reads[read]

    return corrected_consensus, repeats, headers[0], round(np.average(qual), 2), int(np.average(raw)), int(np.average(before)), int(np.average(after))

def read_subreads(seq_file, chrom_reads):
    for read in mm.fastx_read(seq_file, read_comment=False):
        root_name = read[0].split('_')[0]
        if root_name in chrom_reads:
            # root_name : [(header, seq, qual), ...]
            chrom_reads[root_name].append(read) # read = (header, seq, qual)
    return chrom_reads

def read_fasta(infile):
    reads = {}
    for read in mm.fastx_read(infile, read_comment=False):
        reads[read[0]] = read[1]
    return reads

def read_fastq_file(seq_file, check):
    '''
    Takes a FASTQ file and returns a list of tuples
    In each tuple:
        name : str, read ID
        seed : int, first occurrence of the splint
        seq : str, sequence
        qual : str, quality line
        average_quals : float, average quality of that line
        seq_length : int, length of the sequence
    Has a check mode where if it sees one read, it'll return True
    '''
    read_list = []
    for read in mm.fastx_read(seq_file, read_comment=False):
        split_name = read[0].split('_')
        name, seed = split_name[0], 0
        seq, qual = read[1], read[2]
        if check:
            return True
        avg_q = np.average([ord(x)-33 for x in qual])
        s_len = len(seq)
        read_list.append((read[0], seq, qual, avg_q, s_len))
    return read_list

def make_consensus(Molecule, UMI_number, subreads):
    subread_file = path + '/temp_subreads.fastq'
    fastaread_file = path + '/temp_consensus_reads.fasta'
    subs = open(subread_file, 'w')
    fasta = open(fastaread_file, 'w')
    for read in Molecule:
        fasta.write(read)
        # print(read)
        root_name = read[1:].split('_')[0]
        raw = subreads[root_name] # [(header, seq, qual), ...], includes rootname_subread_n
        for entry in raw:
            subs.write('@' + entry[0] + '\n' + entry[1] + '\n+\n' + entry[2] + '\n')
    subs.close()
    fasta.close()
    if read_fastq_file(subread_file, True):
        corrected_consensus, repeats, name, qual, raw, before, after = determine_consensus(str(UMI_number), fastaread_file, subread_file, path)
        return '>%s_%s_%s_%s_%s_%s|%s\n%s\n' %(name, str(qual), str(raw), str(repeats), str(before), str(after), str(UMI_number), corrected_consensus)
    else:
        return ''

def parse_reads(reads, sub_reads, UMIs):
    group_dict, chrom_reads= {}, {}
    groups = []
    UMI_group, previous_start, previous_end = 0, 0, 0
    for name, group_number in sub_reads.items():
        root_name = name.split('_')[0]
        UMI5 = UMIs[name][0]
        UMI3 = UMIs[name][1]
        if root_name not in chrom_reads:
            chrom_reads[root_name] = []
            if group_number not in group_dict:
                group_dict[group_number] = []
            group_dict[group_number].append((name, UMI5, UMI3, reads[name]))
    for group_number in sorted(group_dict):
        group = group_dict[group_number]
        groups.append(list(set(group)))
    return groups, chrom_reads

def group_reads(groups, reads, subreads, UMIs, final, final_UMI_only, matched_reads):
    UMI_group = 0
    # print(len(groups))
    for group in groups:
        group = list(set(group))
        UMI_dict = {}
        set_dict = {}
        group_counter = 0
        if len(group) > 1:
            group_counter += 1
            # print('test', group_counter, len(group))
            UMI_counter = 0
            for i in range(0, len(group), 1):
                UMI_dict[group[i][0]] = set()
            # if len(group) == 2:
            #     print(group[0][1], group[0][2])
            #     print(group[1][1], group[1][2])

            for i in range(len(group)):
                UMI_counter += 1
                UMI_dict[group[i][0]].add(UMI_counter)
                for j in range(i+1, len(group)):
                    if np.abs(len(group[i][3])-len(group[j][3]))/len(group[i][3]) < 0.1:
                        status = 'both'
                        if len(group[i][1]) > 0 and len(group[j][1]) > 0:
                            dist5 = editdistance.eval(group[i][1], group[j][1])

                        else:
                            dist5 = 15
                            status = 'single'
                        if len(group[i][2]) > 0 and len(group[j][2]) > 0:
                            dist3 = editdistance.eval(group[i][2], group[j][2])
                        else:
                            dist3 = 15
                            status = 'single'

                        match = 0
                        if status == 'both':
                            if dist5 + dist3 <= 2:
                                match = 1
                        if match == 1:
                            UMI_dict[group[j][0]] = UMI_dict[group[j][0]]|UMI_dict[group[i][0]]
                            UMI_dict[group[i][0]] = UMI_dict[group[j][0]]|UMI_dict[group[i][0]]

            for i in range(0, len(group), 1):
                for j in range(i+1, len(group), 1):
                    if np.abs(len(group[i][3])-len(group[j][3]))/len(group[i][3]) < 0.1:
                        status = 'both'
                        if len(group[i][1]) > 0 and len(group[j][1]) > 0:
                            dist5 = editdistance.eval(group[i][1], group[j][1])
                        else:
                            dist5 = 15
                            status = 'single'
                        if len(group[i][2]) > 0 and len(group[j][2]) > 0:
                            dist3 = editdistance.eval(group[i][2], group[j][2])
                        else:
                            dist3 = 15
                            status = 'single'

                        match = 0
                        if status == 'both':
                            if dist5 + dist3 <= 2:
                                match = 1
                        if match == 1:
                            UMI_dict[group[j][0]] = UMI_dict[group[j][0]]|UMI_dict[group[i][0]]
                            UMI_dict[group[i][0]] = UMI_dict[group[j][0]]|UMI_dict[group[i][0]]

            for entry in UMI_dict:
                counter_set = UMI_dict[entry]
                if not set_dict.get(tuple(counter_set)):
                    UMI_group += 1
                    set_dict[tuple(counter_set)] = UMI_group

            read_list = []
            for i in range(0, len(group), 1):
                UMI_number = set_dict[tuple(UMI_dict[group[i][0]])]
                read_list.append(('>%s|%s\n%s\n' % (group[i][0], str(UMI_number), group[i][3]), UMI_number, group[i][1], group[i][2]))

            previous_UMI = ''
            Molecule = set()
            for read, UMI_number, umi5, umi3 in sorted(read_list, key=lambda x:int(x[1])):
                matched_reads.write(str(UMI_number) + '\t' + read.split('|')[0] + '\t' + umi5 + '\t' + umi3 + '\n')
                if UMI_number != previous_UMI:
                    if len(Molecule) == 1:
                        final.write(list(Molecule)[0])
                    elif len(Molecule) > 1:
                        new_read = make_consensus(list(Molecule), previous_UMI, subreads)
                        if not new_read:
                            continue
                        final.write(new_read)
                        final_UMI_only.write(new_read)
                    Molecule = set()
                    Molecule.add(read)
                    previous_UMI = UMI_number

                elif UMI_number == previous_UMI:
                    Molecule.add(read)

            if len(Molecule) == 1:
                final.write(list(Molecule)[0])
            elif len(Molecule) > 1:
                new_read = make_consensus(list(Molecule), previous_UMI, subreads)
                if not new_read:
                    continue
                # print('new_read', new_read)
                # print('written')
                final.write(new_read)
                # print('wrote')
                final_UMI_only.write(new_read)
                # print('wrote')
        elif len(group) > 0:
            UMI_group += 1
            final.write('>%s|%s\n%s\n' % (group[0][0], str(UMI_group), group[0][3]))
    # print(group_counter)

def read_UMIs(UMI_file):
    UMI_dict, group_dict, kmer_dict = {}, {}, {}
    group_number = 0
    for line in open(UMI_file):
        a = line.strip('\n').split('\t')
        name, UMI5, UMI3 = a[0], a[1], a[2]
        kmer_list = ['', '', '', '']

        UMI_dict[name] = (UMI5, UMI3)
        if UMI5[5:10] == 'TATAT':
            kmer_list[0] = UMI5[:5]
            kmer_list[1] = UMI5[10:]
        if UMI3[5:10] == 'ATATA':
            kmer_list[2] = UMI3[:5]
            kmer_list[3] = UMI3[10:]

        combination_list = []
        for x in range(4):
            for y in range(x+1, 4):
                first_kmer = kmer_list[x]
                second_kmer = kmer_list[y]
                if first_kmer != '' and second_kmer != '':
                    combination = '_' * x+first_kmer + '_' * (y-x)+second_kmer + '_' * (3-y)
                    combination_list.append(combination)
        match = ''
        for combination in combination_list:
            if combination in kmer_dict:
                match = kmer_dict[combination]
                for combination1 in combination_list:
                    kmer_dict[combination1] = match
                break
        if match == '':
            group_number += 1
            group_dict[group_number] = []
            match = group_number
            for combination in combination_list:
                kmer_dict[combination] = match
        group_dict[match].append(name)
    return group_dict, UMI_dict

def processing(reads, sub_reads, UMIs, groups, final, final_UMI_only, matched_reads):
    annotated_groups, chrom_reads = parse_reads(reads, sub_reads, UMIs)
    # print('reading subreads')
    subreads = read_subreads(subreads_file, chrom_reads)
    # print('grouping and merging consensus reads')
    group_reads(annotated_groups, reads, subreads, UMIs, final, final_UMI_only, matched_reads)

def main():
    final = open(path + '/R2C2_full_length_consensus_reads_UMI_merged.fasta', 'w')
    final_UMI_only = open(path + '/R2C2_full_length_consensus_reads_UMI_only.fasta', 'w')
    matched_reads = open(path + '/matched_reads.txt', 'w')
    print('kmer-matching UMIs')
    groups, UMIs = read_UMIs(umi_file)
    print('reading consensus reads')
    reads = read_fasta(fasta_file)
    count = 0
    sub_reads = {}
    for group in tqdm(sorted(groups)):
        count += len(groups[group])
        for name in groups[group]:
            sub_reads[name] = group
        if count > 500000:
            # print('processing')
            processing(reads, sub_reads, UMIs, groups, final, final_UMI_only, matched_reads)
            count = 0
            sub_reads = {}
    processing(reads, sub_reads, UMIs, groups, final, final_UMI_only, matched_reads)

if __name__ == '__main__':
    main()
