# 10xR2C2
[![Published in Genome Biology](https://img.shields.io/badge/Published%20in-Genome-Biology-blue.svg)](https://doi.org/10.1186/s13059-022-02615-z)
[![DOI](https://zenodo.org/badge/232209750.svg)](https://zenodo.org/badge/latestdoi/232209750)  

Scripts for analyzing 10x R2C2 data

Contents
- [Demultiplexing](#Demultiplexing)
- [Formatting for Seurat](#seurat)
- [Merging UMIs](#umis)

## Demultiplexing ##
To preface, I'm well aware that this isn't the most efficient way of demuxing these reads.
However, it works and isn't worth rewriting.
For this, you need R2C2 reads and R2C2 10x reads.
R2C2 reads are the full length consensus sequences, whereas the R2C2 10x reads are the beginning/end of the sequence that contains the cell barcode and UMI information.
There are a few steps for demuxing 10x R2C2 reads:

#### Determine which barcodes are the most frequent ####
```bash
python3 detBarcodes.py R2C2_10x_postprocessed.fasta 737K-august-2016.txt >1500_most_frequent_bcs.fasta
```
This takes postprocessed reads that contain the 10x barcode information, tallies the occurrences of each barcode, makes a rank ordered histogram of barcode counts, and prints a fasta of the 1500 most frequent barcodes.
The figure making part will need to be uncommented if this gets used with new data, as the number of barcodes is determined manually.
The most frequent barcodes are then given to `Demultiplex_R2C2_reads_kmerBased.py` to determine which reads belong to which barcodes.

#### Add barcoding info to 10x reads ####
```bash
python3 Demultiplex_R2C2_reads_kmerBased.py -i R2C2_10x_postprocessed.fasta -o . -n 1500_most_frequent_bcs.fasta
```
This script outputs `kmer_demuxed.fasta`, which is `R2C2_10x_postprocessed.fasta` with the corresponding barcode added to the header.
The reads are going to be demuxed into cells using the 10x postprocessed fasta (contains cell barcodes).
I do this by reading the file with the barcode data and the full length read at the same time.
To enable that, I need to match the order between the 10x sequences and the full length consensus sequences.

#### Match fasta files ####
```bash
python3 match_fastas.py kmer_demuxed.fasta R2C2_postprocessed.fasta >R2C2_matched.fasta
```

After matching the fasta files, you can actually separate the reads into individual cells.

#### Separate reads into cells ####
```bash
python3 demux_nano.py 1500_most_frequent_bcs.fasta kmer_demuxed.fasta R2C2_matched.fasta
```
This will write cell fasta files to an output directory called `demuxed`.
It will also create a bcGuide, which is needed for downstream analysis.
To merge and analyze specific loci, we need to also demultiplex the subreads associated with each consensus read.
This can be done using `make_cell_subreads.py`.

#### Separating subreads into cells ####
```bash
python3 make_cell_subreads.py /path/to/demuxed/fastas all_subreads.fastq /path/to/output
```
This will collect the base fasta headers for each cell and extract the associated subreads.
It will output to the desired directory and add `_subs.fastq` to the original cell fasta file name.

## Formatting for Seurat <a name="seurat"></a>
[Seurat](https://satijalab.org/seurat/) has a Read10X function that takes a directory that contains three files: genes.tsv, barcodes.tsv, and matrix.mtx.
Usually [Cell Ranger](https://github.com/10XGenomics/cellranger) (10X Genomics) outputs this file structure for you.
The `make_seurat_input.py` script will replicate the output given by Cell Ranger when given an annotation file, barcode guide, and [featureCounts](http://bioinf.wehi.edu.au/featureCounts/) (Subread) output.
The script assumes that featureCounts was run on all cells in one command.

#### Usage ####
```bash
python3 make_seurat_input.py -a gencode.v29.annotation.gtf -b bcGuide -e featureCounts.out -o /path/to/output
```

#### Options ####
```
-a  path to the annotation gtf
-b  path to the bcGuide file (made during demuxing)
-e  path to the featureCounts output
-o  desired output directory (will be made if it doesn't exist)
```

## Cell to cell type ##
The `cell_to_celltype.py` script is used by `merge_psls.py`.

#### Usage ####
```bash
python3 cell_to_celltype.py bcGuide embeddings
```

The bcGuide should come from the demultiplexing script.
The embeddings file should come from Seurat. It can be extracted by doing:

```R
sink('embeddings')
Embeddings(object=pbmc[["tsne"]])
sink()
```

The script outputs a file with the following format:
```
cell number \t cell type \t cell barcode
```

**This script requires heavy modification to fit your dataset. View source code for details.**

## Merging UMIs <a name="umis"></a>
The point of merging UMIs for R2C2 data is to reduce redundancy as well as increase the accuracy for that molecule.
UMI merging identifies all of the subreads that belong to that specific molecule to redo the consensus calling.
We have two different UMI merging scripts, one for the R2C2 splint UMI, one for the 10X barcode UMI.
Before the R2C2 UMI merging script, UMIs need to be extracted from the reads.
This is done by using the `ExtractUMIs.py` script.

```bash
python3 ExtractUMIs.py -i5 5primeAdapter.fasta -i3 3primeAdapter.fasta -i consensus_reads.fasta -o /path/to/output -x splint_number >umi_file.txt
```

```bash
python3 MergeUMIs.py -f consensus_reads.fasta -s subreads.fastq -o /path/to/output -u umi_file.txt -c config -m scores.mat
```

```bash
python3 MergeUMIs10x.py -i /path/to/input -o /path/to/output -c config_file -m scores.mat -t nThreads
```
