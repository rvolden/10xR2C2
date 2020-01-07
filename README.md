# 10xR2C2
Scripts for analyzing 10x R2C2 data

## Demultiplexing ##

## Formatting for Seurat ##
Seurat has a Read10X function that takes a directory that contains three files: genes.tsv, barcodes.tsv, and matrix.mtx.
Usually Cell Ranger (10X Genomics) outputs this file structure for you.
The ```make_seurat_input.py``` script will replicate the output given by Cell Ranger when given an annotation file, barcode guide, and featureCounts (Subread) output.
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
The ```cell_to_celltype.py``` script is used by ```merge_psls.py```.

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

This script requires heavy modification to fit your dataset. View source code for details.
