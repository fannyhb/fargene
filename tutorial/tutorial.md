# fARGene tutorial

This tutorial gives examples on how to use `fARGene` for identification of resistance genes in fragmented metagenomic data and in longer sequences. It also gives examples on how to to create and optimize new models using `fargene_model_creation`.

## Getting started

Install fARGene by following the instructions given [here](../README.md).

The tutorial data can be found in `/path/to/downloaded/fargene/tutorial/tutorialdata/`

Consider the two FASTQ files `reads_1.fastq` and `reads_2.fastq`.
They are simulated paired-end data consisting of fragments of length 100 bases randomly sampled from a *klebsiella pneumoniae* plasmid carrying the metallo-beta-lactamase blaNDM-1 (Yong et al 2009).
The data furthermore consist of the above mentioned plasmid as full-length in the FASTA file `klebisella_plasmid.fasta`. 

For the model creation and optimization, the tutorial data consists of the two files `class_b1_b2.fasta` and `mbl_superfamily.fasta`.
Here `class_b1_b2.fasta` includes 35 genes belonging to MBL subclass B1 and B2, while `mbl_superfamily.fasta` includes 38 sequences from gene homologs without the resistance phenotye.


## Analysis

### Fragmented data usage

To identify and assemble the NDM gene in the fragmented data, run the following command

```
fargene -i /path/to/tutorialdata/*.fastq --meta --hmm-model class_b_1_2 -o tutorial_output_1
```

The file `results_summary.txt` will then display

```
Summary Results fARGene      
The used HMM-model was: class_B_1_2.hmm
Number of input files:  2              
                                       
Metagenomic Input                      
------------------------               
Number of retrieved contigs:    1      
Number of predicted ORFs:   1
```

and the directory `tutorial_output_1/predictedGenes` now consists of four files

  - predicted-orfs-amino.fasta
  - predicted-orfs.fasta
  - retrieved-contigs.fasta
  - retrieved-contigs-peptides.fasta
  
Where `predicted-orfs*` contain the assembled full-length NDM gene as amino and nucleotide sequences. If you are interested in the surroundings of the gene, the file `retrieved-contigs.fasta` contains the from the retrieved paired-end fragments assembled contig. For this data, the contig is 935 nucleotides long while the NDM gene is 816 nucleotides long.

### Longer sequences as input

To analyze the *klebsiella pneumoniae* plasmid, run 

```
fargene -i /path/to/tutorialdata/klebisella_plasmid.fasta --hmm-model class_b_1_2 -o tutorial_output_2
```

The `results_summary.txt` now displays

```
Summary fARGene      
The used HMM-model was: class_B_1_2.hmm
Number of input files:  1              
                                       
Genomic Input                          
------------------------               
Number of predicted genes:  1          
Number of predicted ORFs:   1
```

and the directory `tutorial_output_2/predictedGenes` now consists of four files

  - klebisella_plasmid-class_B_1_2-predicted-orfs-amino.fasta
  - klebisella_plasmid-class_B_1_2-predicted-orfs.fasta
  - klebisella_plasmid-class_B_1_2-filtered-peptides.fasta
  - klebisella_plasmid-class_B_1_2-filtered.fasta

Where `*predicted-orfs*` contain the assembled full-length NDM gene as amino and nucleotide sequences.

## Model creation and optimization

To create a new model and estimate its sensitivity and specificity for full-length genes and fragments of length 100 nucleotides

```
fargene_model_creation -rin class_b1_b2.fasta -nin mbl_superfamily -o model_creation_b_1_2
```

Then the following output is produced

```
Estimating sensitivity for full length genes...                                         
Estimating specificity...                                                               
                                                                                        
For a preferred sensitivity of 1 and max false positive rate of 0.0,                    
the suggested minimal score and corresponding sensitivity and specificity are:          
score=91.2      sensitivity=1.0000      specificity=1.0000                              
                                                                                        
Estimating sensitivity for fragments of lengths 33 AA...                                
Estimating specificity...                                                               
                                                                                        
For a max false positive rate of 0.05,                                                  
the suggested minimal score, score/AA and corresponding sensitivity and specificity are:
score=10.10   score/AA=0.3061 sensitivity=0.8401      specificity=0.9501
```

To decide which threshold scores to set it's helpful to look at the files `results/resulting_sensitivity_specificity_b1_b2_full_length.txt` and `results/resulting_sensitivity_specificity_b1_b2_33.txt` where the sensitivity and specificity is reported for each score.

For full-length genes, we see from `resulting_sensitivity_specificity_b1_b2_full_length.txt` that we can choose any score between 92 and 161 and get a sensitivity and specificity of 1. 

For fragments of 33aa, we see from `resulting_sensitivity_specificity_b1_b2_33.txt` that there is a trade off between the sensitivity and the specificity. If we want a minimal specificity of 0.96, the corresponding sensitivity would be 0.81. The domain score would then be chosen to be 12. The estimation was done on fragments of length 33 aa and therefore the meta-score would be set to 0.3636 (12/33).

To get a better overview of the estimation of the sensitivity and specificity, it is useful to look at the figure `results/resulting_sensitivity_specificity_b1_b2_full_length.png`

![Full-length](misc/resulting_sensitivity_specificity_b1_b2_full_length.png?raw=true "Sensitivity and 1-specificity full-length.")

and figure `results/resulting_sensitivity_specificity_b1_b2_33.png`


![Alt text](misc/resulting_sensitivity_specificity_b1_b2_33.png?raw=true "Sensitivity and 1-specificity fragments.")

where the sensitivity and 1-specificity are plotted as a function of the domain score.
