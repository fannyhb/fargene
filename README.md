# fARGene

fARGene (Fragmented Antibiotic Resistance Gene iENntifiEr ) is a tool that takes either fragmented metagenomic data or longer sequences as input and predicts and delivers full-length antiobiotic resistance genes as output. The tool includes developed and optimized models for a number or resistance gene types,
and the functionality to create and optimize models of your own choice of resistance genes.

The current version of the tool includes developed and optimized models for identification of the following resistance genes

 - Class A beta-lactamases
 - Subclass B1 and B2 beta-lactamases
 - Subclass B3 beta-lactamases
 - Class C beta-lactamases
 - Class D beta-lactamases
 - *qnr*
 
## Table of  contents

* [Getting Started](#getting-started)
   * [Prerequisites](#prerequisites)
   * [Installing](#installing)
* [Data analysis](#data-analysis)
   * [Easy usage](#easy-usage)
   * [Options and usage](#options-and-usage)
   * [Output](#output)
   * [Examples](#examples)
* [Model creation and optimization](#model-creation-and-optimization)
   * [Easy usage](#easy-usage-1)
   * [Options and usage](#options-and-usage-1)
   * [Output](#output-1)
* [Tutorial](#tutorial)
* [Other included tools](#other-included-tools)
* [License](#license)

## Getting Started

These instructions will get you a copy of the most up-to-date version of fARGene. 

### Prerequisites

- Python 2.x
- [EMBOSS transeq](http://emboss.sourceforge.net/download/)
- [seqtk](https://github.com/lh3/seqtk)
- [HMMER](http://hmmer.org/)
- For short-read data:
  - [SPAdes](http://cab.spbu.ru/software/spades/) 3.7.0 or later
  - [Trim Galore!](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/) (optional)
  - [ORFfinder](https://www.ncbi.nlm.nih.gov/orffinder/) (optional)
- For long-read data
  - [prodigal](https://github.com/hyattpd/Prodigal) (optional)

Some of these requirements are optional but might affect the results. If you for example skip the installation of Trim Galore!, the option `--no-quality-filtering` must be used. If you skip the installation of ORFfinder/prodigal for short-read/long-read data, the option `--no-orf-prediction` must be used. fARGene expects these tools to be available in `$PATH`.

For the model creation package you additionally need the following packages:

- numpy
- matplotlib
- [ClustalO](http://www.clustal.org/omega/)

### Installing

#### Installing from source
```
git clone https://github.com/fannyhb/fargene.git
cd fargene
python setup.py install
```

Note:

setup.py will look for and try to install numpy and matplotlib so make sure that you either:
- have these packages installed
- run setup.py as root or with sudo
- install the program in a [conda](https://conda.io/docs/user-guide/install/download.html) environment. 

#### Installing from conda

```
conda install -c conda-forge -c bioconda fargene
```

## Data analysis

### Easy usage

**Fragmented metagenomic input**

```
fargene -i path/to/paired_end_fastqfiles/*.fastq --meta --hmm-model class_a -o output_dir -p num_of_processes
```

Note that the FASTQ files must be paired-end.

**Genomes or longer contigs as nucleotide sequences as input**

```
fargene -i path/to/fastafile(s)/*.fasta --hmm-model class_a -o output_dir
```

**Genomes or longer contigs as protein sequences as input**

```
fargene -i path/to/fastafile(s)/*.fasta --protein --hmm-model class_a -o output_dir
```

Where `hmm-model` can be any of the pre-defined models:
 - Class A beta-lactamases `--hmm-model class_a`
 - Subclass B1 and B2 beta-lactamases `--hmm-model class_b_1_2`
 - Subclass B3 beta-lactamases `--hmm-model class_b_3`
 - Class C beta-lactamases `--hmm-model class_c`
 - Class D beta-lactamases
   - `--hmm-model class_d_1`
   - `--hmm-model class_d_2`
 - *qnr* `--hmm-model qnr`
   
If you choose to use your own profile hidden Markov model you need to specify the score as follows:

```
--hmm-model /path/to/hmmfile.hmm --score full_length_threshold_score --meta-score fragmented_threshold_score_per_AA
```

Note that the meta score is given as score per aa.
If the input is not fragmented data then the `--meta-score` is not required.

It is also possible to change the scores for the predifened models, just add the option `--score new_score` and/or `--meta-score new_meta_score`.

### Options and usage

Run `fargene --help` to get all the options of fargene.

```
usage: fargene [-h] --infiles INFILES [INFILES ...] --hmm-model HMM_MODEL
               [--score LONG_SCORE] [--meta] [--meta-score META_SCORE]
               [--output OUTDIR] [--force] [--tmp-dir TMP_DIR] [--protein]
               [--processes PROCESSES] [--min-orf-length MIN_ORF_LENGTH]
               [--retrieve-whole] [--no-orf-predict] [--no-quality-filtering]
               [--no-assembly] [--orf-finder] [--store-peptides] [--rerun]
               [--amino-dir AMINO_DIR] [--fasta-dir FASTA_DIR]
               [--translation-format TRANS_FORMAT] [--loglevel {DEBUG,INFO}]
               [--logfile LOGFILE]

Searches and retrieves new and previously known genes from fragmented
metagenomic data and genomes. Copyright (c) Fanny Berglund 2018.

optional arguments:
  -h, --help            show this help message and exit
  --infiles INFILES [INFILES ...], -i INFILES [INFILES ...]
                        Input file(s) to be searched. Could either be in FASTA
                        or FASTQ format.
  --hmm-model HMM_MODEL
                        The Hidden Markov Model that should be used to analyse
                        the data. Could either be one of the pre-defined
                        models or the path to a custom HMM.
  --score LONG_SCORE, -sl LONG_SCORE
                        The threshold score for a sequence to be classified as
                        a (almost) complete gene (default: None).
  --meta                If the input data is paired end metagenomic data
                        (default: False).
  --meta-score META_SCORE, -sm META_SCORE
                        The threshold score for a fragment to be classified as
                        a positive. Expressed as score per amino acid
                        (default: None).
  --output OUTDIR, -o OUTDIR
                        The output directory for the whole run (default:
                        ./fargene_output).
  --force, -f           Overwrite output directory if it exists (default:
                        False).
  --tmp-dir TMP_DIR     Directory for (sometimes large) intermediate files.
                        (default: OUT_DIR/tmpdir)
  --protein             If the input sequence(s) is amino acids (default:
                        False).
  --processes PROCESSES, -p PROCESSES
                        Number of processes to be used when processing
                        metagenomic data (default: 1).
  --min-orf-length MIN_ORF_LENGTH
                        The minimal length for a retrieved predicted ORF (nt).
                        (default: 90% of the length of the chosen hmm.)
  --retrieve-whole      Use this flag if the whole sequence where a hit is
                        detected should be retrieved (default: False).
  --no-orf-predict      Do not perform ORF prediction.
  --no-quality-filtering
                        Use if no quality control should be performed on the
                        metagenomic data (default: False).
  --no-assembly         Use if you want to skip the assembly and retrieval of
                        contigs for metagenomic data (default: False).
  --orf-finder          Use NCBI ORFfinder instead of prodigal for ORF
                        prediction of genomes/contigs (default: False).
  --store-peptides, -sp
                        Store the translated sequences. Useful if you plan to
                        redo the analysis using a different model and want to
                        skip the preprocessing steps (default: False).
  --rerun               Use of you want to redo the analysis or do the
                        analysis using a different model and have kept either
                        the nucletide or amino acid sequences. Please note
                        that this only works if the input data is the same for
                        both runs (default: False).
  --amino-dir AMINO_DIR
                        Where the amino acid sequences generated by the method
                        are located. Only to be used in combination with
                        --rerun
  --fasta-dir FASTA_DIR
                        Where the nucleotide sequences in FASTA generated by
                        previous runs of the method are located. Only to be
                        used in combination with --rerun
  --translation-format TRANS_FORMAT
                        The translation format that transeq should use.
                        (default: pearson)
  --loglevel {DEBUG,INFO}
                        Set logging level (default: INFO).
  --logfile LOGFILE     Logfile (default: fargene_analysis.log).
```

### Output

The two most important output files are `predicted-orfs.fasta` and `predicted-orfs-amino.fasta` which are located in `output_dir/predictedGenes/`.


This directory also contains `retrieved-contigs.fasta` and `retrieved-contigs-peptides.fasta`. Where the first file contains the (complete) contigs that passed the final full-length classification, and the second file contains the parts of the contigs that passed the final classification step that aligned with the HMM, as amino acid sequences. Note that the second file is a prediction of where the genes are located on the contig and does usually not include the start and/or the stop of the genes. 


A summary of the analysis can be found in `output_dir/results_summary.txt` and the logfile is found in `output_dir/novelGeneFinder.log`.

Below is a summary of the remaining output:

|File/Directory| Description|
|--------------|------------|
|hmmsearchresults/| All the output files from `hmmsearch`.|
|retrievedFragments/| Each fragment that were classified as positive, together with its read-pair.|
|retrievedFragments/all_retrieved_[12].fastq| All quality controlled retrieved fragments gathered in two files.|
|retrievedFragments/trimmedReads/| The quality controlled retrieved fragments from each input file. |
|tmpdir/| Various files that can be deleted if you don't want to redo the analysis.|
|tmpdir/* positives.out| List of ids for sequences/genes that were classified as positives.|
|tmpdir/infile(s).fasta| The from FASTQ to FASTA converted input files.|
|tmpdir/infiles(s)-amino.fasta| The translated input sequences. Are only saved if option `--store-peptides` is used.|
|spades_assembly/| The output from the SPAdes assembly.|

#### For genomes or longer contigs as input

The output is basically the same as for the metagenomic input. The most important difference is the file
`input_file-hmm_model_name-filtered.fasta` located in `output_dir/predictedGenes`. This file contains the sequences that passed the final classification step, but only the parts that where predicted by the HMM to be part of the gene. Since this prediction is rather conservative, the start/stop of the genes are usually not included here. The file `input_file-hmm_model_name-filtered-peptides.fasta` is the above file translated in the same frame as the gene is predicted to be located. 

### Examples

#### Analyze the same dataset with several models

If you want to search the same dataset for different types of genes it is recommended to use the `--rerun` option in order to avoid repeating the time consuming preprocessing steps. This can be used on the converted FASTA files generated by the first analysis with fargene as shown below

First analysis
```
fargene -i path/to/paired_end_fastqfiles/*.fastq --meta --hmm-model class_a -o class_a_out -p num_of_processes
```

To repeat the analysis with another model

```
fargene -i path/to/paired_end_fastqfiles/*.fastq --meta --hmm-model class_c -o class_c_out -p num_of_processes
--rerun --fasta-dir class_a_out/tmpdir/
```

Depending on the size of your storage and the size of the dataset, the option `--store-peptides` can be used in the first analysis. This will significantly speed up the analysis for the following models, but be aware of the potential huge amounts of data that will be stored. The commands can then look as follows

First analysis

```
fargene -i path/to/paired_end_fastqfiles/*.fastq --meta --hmm-model class_a -o class_a_out -p num_of_processes 
--store-peptides
```

To repeat the analysis with another model

```
fargene -i path/to/paired_end_fastqfiles/*.fastq --meta --hmm-model class_c -o class_c_out -p num_of_processes
--rerun --amino-dir class_a_out/tmpdir/
```

## Model creation and optimization

### Easy usage

To create and optimize threshold score for a new model using the default values (full-length genes and 10 000 fragments of length 33 aa)

```
fargene_model_creation --reference-sequences file_with_reference_sequences.fasta
--negative-sequences file_with_negative_sequences.fasta --modelname name_of_my_new_model
--output output_dir
```

Where both the reference and negative sequences should be protein sequences.

### Options and usage

Run `fargene_model_creation --help` to get all the options of the model creation and optimization

```
usage: fargene_model_creation [-h] --reference-sequences REFERENCE_SEQUENCES        
                              [--negative-sequences NEGATIVE_SEQUENCES]             
                              [--output OUTPUT_DIR] [--modelname MODELNAME]         
                              [--fragment-lengths FRAGMENT_LENGTHS]                 
                              [--num-fragments NUM_FRAGMENTS] [--only-sens]         
                              [--only-spec] [--only-full-length] [--only-fragments] 
                                                                              
A program to create and optimize profile hidden Markov models                 
                                                                              
optional arguments:                                                           
  -h, --help            show this help message and exit                       
  --reference-sequences REFERENCE_SEQUENCES, -rin REFERENCE_SEQUENCES         
                        The sequences that the model should be built of.      
  --negative-sequences NEGATIVE_SEQUENCES, -nin NEGATIVE_SEQUENCES            
                        The sequences that should be used as the negative     
                        dataset. Should preferable be similar sequences but   
                        without the desired phenotype.                        
  --output OUTPUT_DIR, -o OUTPUT_DIR                                          
                        The directory where the output should be saved.       
  --modelname MODELNAME                                                       
                        The name of the new model                             
  --fragment-lengths FRAGMENT_LENGTHS, -l FRAGMENT_LENGTHS                    
                        The length (aa) of the fragments that should be used  
                        to determine the threshold score for metagenomic      
                        input. (default: 33 AA)                               
  --num-fragments NUM_FRAGMENTS                                               
                        The number of fragments that should be created from   
                        each gene. (default: 10 000)                          
  --only-sens           Should be used if only sensitivity of the model should
                        be estimated.                                         
  --only-spec           Should be used if only the specificity of the model   
                        should be estimated.                                  
  --only-full-length    Should be used if you only want to optimize the       
                        threshold score for full length genes.                
  --only-fragments      Should be used if you only want to optimize the       
                        threshold score for full fragmented genes.                     
```

### Output

The created HMM is located in `output_dir/hmmerModel` and is called `modelname.hmm`.

The two most important files from the estimation of sensitivity and specificity are called `resulting_sensitivity_specificity_modelname_full_length.txt` and `resulting_sensitivity_specificity_modelname_fragmentlength.txt` and are located in `output_dir/results/`. These contain the estimated senstivity and specificity as a function of the domain score from the `hmmsearch`.

This results is also visualized in the figures `resulting_sensitivity_specificity_modelname_full_length.png` `resulting_sensitivity_specificity_modelname_fragmentlength.png`.

The combined results from the `hmmsearch` of full-length sequences are called `modelname-hmmsearch-refrence-sequences-full-length.txt`and `modelname-hmmsearch-negative-sequences-full-length.txt` and are useful to detect genes responsible for potenial outline scores.

## Tutorial
 
For a tutorial of how to use fargene click [here](tutorial/tutorial.md).

## Other included tools

fARGene comes with a third tool called `pick_long_reads` which can extract sequences longer than a certain length.

Usage for nucleotide sequences:

```
pick_long_reads -i all_contigs.fasta --length 800 -o contigs_longer_than_800nt.fasta
```

Usage for protein sequences:

```
pick_long_reads -i all_contigs-amino.fasta --length 250 --cut-stars -o contigs_longer_than_250aa.fasta
```

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details.
