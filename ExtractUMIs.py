import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-i5', '--input5_fasta', type=str)
parser.add_argument('-i3', '--input3_fasta', type=str)
parser.add_argument('-i', '--input_fasta', type=str)
parser.add_argument('-o', '--output_path', type=str)
parser.add_argument('-x', '--splint_number', type=str)
args = parser.parse_args()
output_path = args.output_path + '/'
input_reads = args.input_fasta
input5 = args.input5_fasta
input3 = args.input3_fasta
splint = args.splint_number

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
            name = a[1:].split()[0]
            sequence = ''
          else:
            sequence += a
    if sequence != '':
        reads[name] = sequence
    return reads

def getFlanking(splint):
    flankingDict = {}
    flankingDict['1'] = ('CCATA', 'ATCAC', 'ATCGC', 'TTAGT')
    flankingDict['2'] = ('AATAA', 'ATTGC', 'TCACG', 'CTGCC')
    flankingDict['3.5'] = ('TGGGT', 'TAAAA', 'CAGCT', 'ATATT')
    flankingDict['4.5'] = ('TCCGT', 'TACGA', 'AGGCG', 'ACCTG')
    flankingDict['5'] = ('GATAG', 'TTCTG', 'AAGAG', 'GAACC')
    flankingDict['6'] = ('CTGGT', 'ACGTC', 'ATTAG', 'TAGCC')
    flankingDict['7'] = ('TTATA', 'TGGAC', 'AGAGG', 'CAGAT')
    flankingDict['8'] = ('AGAAT', 'TTCTC', 'AGGCT', 'AGGGA')
    flankingDict['9'] = ('CGGGG', 'AAGAT', 'GTAAC', 'ACCTA')
    flankingDict['10'] = ('TTCAG', 'ATTAA', 'CTGGT', 'AGCCT')
    flankingDict['11'] = ('CGATA', 'ATTCT', 'TGGTG', 'TCATA')
    flankingDict['12'] = ('CAAGT', 'CATAC', 'AAGTC', 'CACAA')
    flankingDict['13'] = ('ATCTG', 'ATGCA', 'ATGAC', 'CGTAG')
    flankingDict['14'] = ('TGGAC', 'ACTTC', 'TGAGG', 'TTAAT')
    flankingDict['15'] = ('CCTGT', 'CATGG', 'TCGTA', 'AATTT')
    flankingDict['16'] = ('GTACA', 'TCTAT', 'AGAAT', 'TTTCT')

    return flankingDict[splint]

def reverse_complement(sequence):
  Seq = ''
  complement = {'A':'T', 'C':'G', 'G':'C', 'T':'A', 'N':'N', '-':'-'}
  for item in sequence[::-1]:
    Seq = Seq + complement[item]
  return Seq

def find_UMI(sequence, frame1, frame2):
    UMI = ''
    for i in range(0, max(len(sequence)-25, 1), 1):
        match1 = sequence[i:i+5]
        if match1 == frame1:

            for j in range(i+12, i+25, 1):
                match2 = sequence[j:j+5]
                if match2 == frame2:
                    UMI = sequence[i+5:j]
    return UMI

print('reading left adapter')
reads5 = read_fasta(input5)
print('reading right adapter')
reads3 = read_fasta(input3)
print('reading full_length_reads')
reads = read_fasta(input_reads)

out = open(output_path + '/R2C2_full_length_consensus_reads.UMI', 'w')

for read in reads:
    name = read
    if read in reads5:
        sequence5 = reads5[read]
    else:
        sequence5 = ''

    if read in reads3:
        sequence3 = reads3[read]
    else:
        sequence3 = ''

    length5 = len(sequence5)
    length3 = len(sequence3)
    UMI5 = ''
    UMI3 = ''

    flanking1, flanking2, flanking3, flanking4 = getFlanking(splint)

    UMI5 = find_UMI(sequence3, flanking1, flanking2) #Splint 1
    UMI3 = find_UMI(sequence5, flanking3, flanking4) #Splint 1
    if UMI5=='' and UMI3=='':
        UMI3 = find_UMI(sequence3, flanking3, flanking4) #Splint 1
        UMI5 = find_UMI(sequence5, flanking1, flanking2) #Splint 1
    out.write('%s\t%s\t%s\n' %(name, UMI5, UMI3))
