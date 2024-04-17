# Input data
This directory contains input data for the analysis.

## Basic information about sequences and alignments

These files are used for the basic processing of the deep sequencing data to call variants by barcode and count barcodes:
  
   - `wildtype_sequence_x.fasta`: for each RBD background x, the sequence of the unmutated RBDs used for initializing CodonVariantTables.


## Files related to processing NGS sequencing data

  - [`codon_variant_table_x.csv`: for each RBD background x, barcode-variant lookup tables for all backgrounds generated from the primary DMS study of each background

  - `barcode_runs_x.csv`: for each RBD background x, list of the Illumina runs used to count the barcodes for different samples.

## For visualizing escape data:

These files are used for visualizing the antibody- or serum-escape data:

  - [site_color_schemes.csv](site_color_schemes.csv): Schemes for how to color sites (can be used in escape profiles). Here are details on these schemes.

  - `escape_profiles_config_x.yaml`: for each RBD backgorund x, information on how to plot the escape profiles; manually edit this to alter their plotting.

  - [dms-view_metadata.md](dms-view_metadata.md): Used for rendering the dms-view page to visualize data.


## Data from GISAID surveillance

  - [210801_mutation_counts.csv](210801_mutation_counts.csv): The counts of all RBD mutations in the SARS-CoV-2 spike alignment downloaded from [GISAID](https://www.gisaid.org/) on Aug. 1, 2021. These data are used to help set filters on the ACE2 binding and RBD expression scores, to ideally retain mutations that are observed an appreciable number of times in GISAID sequences.

## PDB files in [pdbs](pdbs/) subdirectory

  - [6M0J](pdbs/6M0J.pdb): Wuhan-Hu-1 RBD bound to huACE2.

  - [7LYQ](pdbs/7LYQ.pdb): B.1.351 spike trimer

  - [7LYQ_RBD_ACE2](pdbs/7LYQ_RBD_ACE2.pdb): RBD from 7LYQ aligned with 6M0J and shown bound to ACE2

  - [7LYQ_RBD](pdbs/7LYQ_RBD.pdb): RBD from 7LYQ
