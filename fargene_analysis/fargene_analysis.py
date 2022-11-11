#!/usr/bin/env python2.7
from os import path, makedirs, getcwd
from sys import argv
from collections import defaultdict
from multiprocessing import Pool, cpu_count
#from distutils.spawn import find_executable
from shutil import which
import argparse
import logging
import itertools
import importlib

from Transformer import Transformer
from HmmModel import HmmModel
from predict_orfs import predict_orfs_orfFinder, predict_orfs_prodigal
from ResultsSummary import ResultsSummary
import utils

def parse_args(argv):
    desc = 'Searches and retrieves new and previously known genes from fragmented metagenomic data and genomes'
    copyright = 'Copyright (c) Fanny Berglund 2018.'
    parser = argparse.ArgumentParser(description=desc+'. '+copyright)
    parser.add_argument('--infiles', '-i', nargs='+', required=True,
                        help='Input file(s) to be searched. Could either be in FASTA or FASTQ format.')
    parser.add_argument('--hmm-model', dest='hmm_model', required=True,
                        help='The Hidden Markov Model that should be used to analyse the data.'\
                        ' Could either be one of the pre-defined models or the path to a custom HMM.')
    parser.add_argument('--score','-sl', dest='long_score', required=False,
                        help = 'The threshold score for a sequence to be classified as a (almost) complete gene (default: %(default)s).')

    parser.add_argument('--meta', action='store_true',
                        help='If the input data is paired end metagenomic data (default: %(default)s).')
    parser.add_argument('--meta-score','-sm', dest='meta_score', type=float,
                        help = 'The threshold score for a fragment to be classified as a positive. '\
                        'Expressed as score per amino acid (default: %(default)s).')
    
    parser.add_argument('--output','-o', dest='out_dir', metavar='OUTDIR',
                        help='The output directory for the whole run (default: %(default)s).')
    parser.add_argument('--force','-f',action='store_true',
                        help='Overwrite output directory if it exists (default: %(default)s).')
    
    parser.add_argument('--tmp-dir', dest='tmp_dir',
                        help='Directory for (sometimes large) intermediate files. '\
                                '(default: OUT_DIR/tmpdir)')

    parser.add_argument('--protein', action='store_true', dest='protein',
                        help= 'If the input sequence(s) is amino acids (default: %(default)s).')


    parser.add_argument('--processes','-p', type=int, default=1, dest='processes',
                        help = 'Number of processes to be used when processing metagenomic data (default: %(default)s).')

    parser.add_argument('--min-orf-length', type=int, dest='min_orf_length' ,
                        help='The minimal length for a retrieved predicted ORF (nt). '\
                                '(default: 90%% of the length of the chosen hmm.)')
   
    parser.add_argument('--retrieve-whole', action='store_true', dest='retrieve_whole',
                        help='Use this flag if the whole sequence where a hit is detected should be retrieved (default: %(default)s).')

    parser.add_argument('--no-orf-predict', action='store_false', dest='orf_predict' ,
                        help='Do not perform ORF prediction.')
    parser.add_argument('--no-quality-filtering', default=False, action='store_true', dest='no_quality_filtering',
                        help = 'Use if no quality control should be performed on the metagenomic data (default: %(default)s).')
    parser.add_argument('--no-assembly', action='store_true', dest='no_assembly',
                        help = 'Use if you want to skip the assembly and retrieval of contigs for metagenomic data (default: %(default)s).')
    parser.add_argument('--orf-finder', action='store_true', dest='orf_finder',
                        help = 'Use NCBI ORFfinder instead of prodigal for ORF prediction of genomes/contigs (default: %(default)s).')

    parser.add_argument('--store-peptides','-sp', default=False, action='store_true', dest='store_peptides',
                        help = 'Store the translated sequences. Useful if you plan to redo '\
                               'the analysis using a different model and want to skip the preprocessing steps '\
                               '(default: %(default)s).')
    parser.add_argument('--rerun', action='store_true',
                        help = 'Use of you want to redo the analysis or do the analysis using a different model '\
                                'and have kept either the nucletide or amino acid sequences. '\
                                'Please note that this only works if the input data is the same for both runs '\
                                '(default: %(default)s).')
    parser.add_argument('--amino-dir', dest='amino_dir',
                        help = 'Where the amino acid sequences generated by the method are located.'\
                                ' Only to be used in combination with --rerun')
    parser.add_argument('--fasta-dir', dest='fasta_dir',
                        help = 'Where the nucleotide sequences in FASTA generated by previous runs of the method are located. '\
                                'Only to be used in combination with --rerun')

    parser.add_argument('--translation-format', default='pearson', dest='trans_format',
            help= 'The translation format that transeq should use. (default: %(default)s)')

    parser.add_argument('--loglevel', choices=['DEBUG', 'INFO'], default='INFO', type=str,
                        help='Set logging level (default: %(default)s).')
    parser.add_argument('--logfile', type=str, default='fargene_analysis.log',
                        help='Logfile (default: %(default)s).')

    parser.set_defaults(
            meta = False,
            retrieve_whole = False,
            protein = False,
            tmp_dir = False,
            hmm_model = None,
            long_score = None,
            meta_score = None,
            sensitive = False,
            no_assembly = False,
            transformer = True,
            orf_predict = True,
            min_orf_length = None,
            rerun = False,
            amino_dir = False,
            force = False,
            orf_finder = False,
            out_dir = './fargene_output')

    options = parser.parse_args()

    logger = logging.getLogger(__name__)
    if options.loglevel == 'DEBUG':
	    logger.setLevel(logging.DEBUG)
    else:
	    logger.setLevel(logging.INFO)
    logging_format_file = logging.Formatter('%(asctime)s %(levelname)s: %(message)s')
    logging_format_console = logging.Formatter('%(levelname)s: %(message)s')
    file_handler = logging.FileHandler(options.logfile)
    console_handler = logging.StreamHandler()
    file_handler.setFormatter(logging_format_file)
    console_handler.setFormatter(logging_format_console)
    logger.addHandler(file_handler)
    logger.addHandler(console_handler)

    return options, logger

def main():
    
    options, logger = parse_args(argv)
    check_executables_in_path(options, logger)

    outdir = path.abspath(options.out_dir)
    options.hmm_out_dir = '%s/hmmsearchresults' %(outdir)
    options.res_dir = '%s/retrievedFragments' %(outdir)
    options.final_gene_dir = '%s/predictedGenes' %(outdir)
    options.assembly_dir = '%s/spades_assembly' %(outdir)
    if not options.tmp_dir:
        options.tmp_dir = '%s/tmpdir' %(outdir)
    
    if path.isdir(options.out_dir) and not options.force:
        msg = ('The directory {0} already exists. To overwrite use the'
                ' --force flag').format(options.out_dir)
        logger.error(msg)
        logger.info('Exiting pipeline')
        exit()
    elif path.isdir(options.out_dir) and options.force:
        if not utils.remove_files(options.final_gene_dir, options.res_dir):
            logger.info('Exiting pipeline')
            exit()
    else:
        utils.create_dir(options.out_dir)

    for infile in options.infiles:
        if not path.isfile(infile):
            msg ='The provided input file {0} does not exist.'.format(infile)
            logger.critical(msg)
            logger.info('Exiting pipeline')
            exit()


    check_arguments(options, logger)

    utils.create_dir(options.hmm_out_dir)
    utils.create_dir(options.tmp_dir)
    utils.create_dir(options.final_gene_dir)

    if options.meta:
        utils.create_dir(options.res_dir)
        if not options.no_quality_filtering:
            options.trimmed_dir = '%s/trimmedReads' %(path.abspath(options.res_dir))
            utils.create_dir(options.trimmed_dir)

    summaryFile = '%s/results_summary.txt' %(outdir)
    Results = ResultsSummary(summaryFile, len(options.infiles), options.hmm_model)

    logger.info('Starting fARGene')
    logger.info('Starting pipeline, planning to analyze %s files', len(options.infiles))
    logger.info('Running on %s processes' %str(options.processes))
    

    if not options.meta:
        parse_fasta_input(options, Results, logger)
        Results.write_summary(False)
        numGenes = Results.retrievedSequences
        retrieved = 'possible genes'
    else:
        options.protein = False
        parse_fastq_input(options, Results, logger)
        Results.write_summary(True)
        numGenes = Results.retrievedContigs
        retrieved = 'retrieved contigs'
    logger.info('Done with pipeline')
    
    msg = ('fARGene is done.\n'
           'Total number of {}: {}\n'
           'Total number of predicted ORFS longer than {} nt: {}\n'
           'Output can be found in {}'
           ).format(retrieved, numGenes, options.min_orf_length, Results.predictedOrfs, outdir)
    logger.info(msg)

def check_arguments(options, logger):
    predefined = False
    model_location = path.dirname(__file__)+ '/models'
    preDefinedModels = [
            HmmModel("b1", model_location + "/B1.hmm", 135.8, float(0.2424)),
            HmmModel("class_b_1_2", model_location + "/class_B_1_2.hmm", 127, float(0.3636)),
            HmmModel("class_b_3", model_location + "/class_B_3.hmm", 103, float(0.30303)),
            HmmModel("class_a", model_location + "/class_A.hmm", 105, float(0.2424)),
            HmmModel("class_c", model_location + "/class_C.hmm", 248, float(0.30303)),
            HmmModel("class_d_1", model_location + "/class_D_1.hmm", 182, float(0.3030)),
            HmmModel("class_d_2", model_location + "/class_D_2.hmm", 234, float(0.3030)),
            HmmModel("qnr", model_location + "/qnr.hmm", 150, float(0.51515)),
            HmmModel("tet_efflux", model_location + "/tet_efflux.hmm", 400, float(0.3030)),
            HmmModel("tet_rpg", model_location + "/tet_rpg.hmm", 471, float(0.4545)),
            HmmModel("tet_enzyme", model_location + "/tet_enzyme.hmm", 300, float(0.4545)),
            HmmModel("erm_type_a", model_location + "/erm_typeA.hmm", 200, float(0.4545)),
            HmmModel("erm_type_f", model_location + "/erm_typeF.hmm", 175, float(0.4242)),
            HmmModel("mph", model_location + "/mph.hmm", 100, float(0.3030)),
            HmmModel("aminoglycoside_model_a", model_location + "/aminoglycoside_model_a.hmm", 100, float(0.2727)),
            HmmModel("aminoglycoside_model_b", model_location + "/aminoglycoside_model_b.hmm", 100, float(0.5152)),
            HmmModel("aminoglycoside_model_c", model_location + "/aminoglycoside_model_c.hmm", 100, float(0.2424)),
            HmmModel("aminoglycoside_model_d", model_location + "/aminoglycoside_model_d.hmm", 100, float(0.3333)),
            HmmModel("aminoglycoside_model_e", model_location + "/aminoglycoside_model_e.hmm", 100, float(0.3333)),
            HmmModel("aminoglycoside_model_f", model_location + "/aminoglycoside_model_f.hmm", 100, float(0.3636)),
            HmmModel("aminoglycoside_model_g", model_location + "/aminoglycoside_model_g.hmm", 100, float(0.3030)),
            HmmModel("aminoglycoside_model_h", model_location + "/aminoglycoside_model_h.hmm", 100, float(0.1818)),
            HmmModel("aminoglycoside_model_i", model_location + "/aminoglycoside_model_i.hmm", 100, float(0.2121))
            ]

    for model in preDefinedModels:
        if (options.hmm_model).lower()== model.name:
            options.hmm_model = model.path
            if not options.long_score:
                options.long_score = model.long_score
            if not options.meta_score:
                options.meta_score = model.meta_score
            predefined = True

    if not path.isfile(options.hmm_model):
        names = "\n".join([str(hmmModel.name) for hmmModel in preDefinedModels])
        msg = ("\nThe HMM file {0} could not be found.\n"
                 "Either provide a valid path to a HMM or choose "
                 "one of the following pre-defined models:\n{1}").format(
                         options.hmm_model,names)
        logger.critical(msg)
        logger.info('Exiting pipeline')
        exit()

    if not predefined:
        if options.long_score is None:
            msg = "No threshold score for whole genes was given.\n"+\
            "Please provide one using the option --score"     
            logger.critical(msg)
            logger.info('Exiting pipeline')
            exit()
        if options.meta and options.meta_score is None:
            msg = "No threshold score for metagenomic fragments was given.\n"+\
            "Please provide one using the option --meta-score"     
            logger.critical(msg)
            logger.info('Exiting pipeline')
            exit()

    topFile = options.infiles[0]
    if options.meta:
        if not utils.is_fastq(topFile):
            msg = "If using the meta options, the input files must be FASTQ"
            logger.critical(msg)
            logger.info('Exiting pipeline')
            exit()
    else:
        if not utils.is_fasta(topFile):
            msg = "If not using the meta option, the input file(s) must be FASTA"
            logger.critical(msg)
            logger.info('Exiting pipeline')
            exit()              

    if not options.min_orf_length:
        options.min_orf_length = utils.decide_min_ORF_length(options.hmm_model)

    if options.rerun:
        fastqBaseName = path.splitext(path.basename(options.infiles[0]))[0]
        if options.amino_dir:
            options.amino_dir = path.abspath(options.amino_dir)
        else:
            options.amino_dir = path.abspath(options.tmp_dir)
        peptideFile = '%s/%s-amino.fasta' %(options.amino_dir, fastqBaseName)
        if path.isfile(peptideFile):
            return
        else:
            if options.fasta_dir:
                options.fasta_dir = path.abspath(options.fasta_dir)
            else:
                options.fasta_dir = path.abspath(options.tmp_dir)
            fastaFile = '%s/%s.fasta' %(path.abspath(options.fasta_dir), fastqBaseName)
            if not path.isfile(fastaFile):
                msg = 'Neither nucleotide or amino sequences exists as FASTA.\n'\
                      'Please provide path to amino or nucleotide sequences or remove flag --rerun'
                logger.critical(msg)
                logger.info('Exiting pipeline')
                exit()

def check_executables_in_path(options, logger):
    executables = ['transeq','seqtk','hmmsearch']
    if options.meta:
        if not options.no_assembly:
            executables.append('spades.py')
        if not options.no_quality_filtering:
            executables.append('trim_galore')
        if options.orf_predict:
            executables.append('ORFfinder')
    else:
        if options.orf_predict:
            executables.append('prodigal')
    for executable in executables:
        if not which(executable):
            msg = ('Did not find {} in path.\n'
                    'Exiting pipeline').format(executable)
            logger.critical(msg)
            exit()

def parse_fasta_input(options, Results, logger):
    modelName = path.splitext(path.basename(options.hmm_model))[0]
    logger.info('Parsing FASTA files')
    frame = '6'
    for fastafile in options.infiles:
        fastaBaseName = path.splitext(path.basename(fastafile))[0]
        hmmOut = '%s/%s-%s-hmmsearched.out' %(path.abspath(options.hmm_out_dir), fastaBaseName,modelName)
        fastaOut = '%s/%s-%s-filtered.fasta' %(path.abspath(options.final_gene_dir), fastaBaseName,modelName)
        aminoOut = '%s/%s-%s-filtered-peptides.fasta' %(path.abspath(options.final_gene_dir), fastaBaseName,modelName)
        orfFile = '%s/%s-%s-predicted-orfs.fasta' %(path.abspath(options.final_gene_dir), fastaBaseName,modelName)
        orfAminoFile = '%s/%s-%s-predicted-orfs-amino.fasta' %(path.abspath(options.final_gene_dir), fastaBaseName,modelName)
        hitFile = '%s/%s-positives.out' %(path.abspath(options.tmp_dir), fastaBaseName)
        elongated_fasta ='%s/%s-gene-elongated.fasta' %(path.abspath(options.tmp_dir), fastaBaseName)
        if options.protein:
            utils.perform_hmmsearch(fastafile, options.hmm_model, hmmOut, options)
            utils.classifier(hmmOut, hitFile, options)
            hitDict = utils.create_dictionary(hitFile, options)
            utils.retrieve_fasta(hitDict, fastafile, fastaOut, options)
        else: 
            if options.store_peptides:
                peptideFile ='%s/%s-amino.fasta' %(path.abspath(options.tmp_dir), fastaBaseName)
                utils.translate_sequence(fastafile, peptideFile, options, frame)
                logger.info('Performing hmmsearch')
                utils.perform_hmmsearch(peptideFile, options.hmm_model, hmmOut, options)
            else:
                utils.translate_and_search(fastafile, options.hmm_model, hmmOut, options)
            utils.classifier(hmmOut,hitFile,options)
            hitDict = utils.create_dictionary(hitFile, options)
            utils.retrieve_fasta(hitDict, fastafile, fastaOut, options)
            if not path.isfile(fastaOut):
                logger.critical('Could not find file %s', fastaOut)
#                exit()
            else:
                utils.retrieve_surroundings(hitDict, fastafile, elongated_fasta)
                if path.isfile(elongated_fasta):
                    if not options.orf_finder:
                        tmpORFfile = '%s/%s-long-orfs.fasta' %(options.tmp_dir,fastaBaseName)
                        predict_orfs_prodigal(elongated_fasta, options.tmp_dir, tmpORFfile, options.min_orf_length) 
                        orfFile = utils.retrieve_predicted_orfs(options, tmpORFfile)
                    else:
                        tmpORFfile = '%s/%s-long-orfs.fasta' %(options.tmp_dir, fastaBaseName)
                        predict_orfs_orfFinder(elongated_fasta,options.tmp_dir, tmpORFfile, options.min_orf_length)
                        orfFile = utils.retrieve_predicted_orfs(options, tmpORFfile)
                if options.store_peptides:
                    options.retrieve_whole = False
                    utils.retrieve_peptides(hitDict, peptideFile, aminoOut, options)
                else:
                    tmpFastaOut = utils.make_fasta_unique(fastaOut, options)
                    utils.retrieve_predicted_genes_as_amino(options, tmpFastaOut, aminoOut, frame='6')
        Results.count_hits(hitFile)
    if path.isfile(orfFile):                                      
        if not options.orf_finder:                                
            Results.count_orfs_genomes(orfFile)                   
        else:                                                     
            Results.predictedOrfs = Results.count_contigs(orfFile)
                                                              
    return orfFile


def parse_fastq_input(options, Results, logger):
    '''
    If the input is .fastq
    Pooled:
        1) Converts to .fasta with seqtk
        2) Translates to peptides and pipes to hmmsearch
        3) Parses the output ffrom hmmsearch and classifies the reads > hitFile
    4) Parses the hitFile and saves to dictionary,
        assumes paired end.
        If read_id_X from fileY_1 is a hit then it is saved
        as [fileY].append(read_id_X) and vice verse
        doing this after classification to save RAM
    5) Retrieves the hits in fastq using seqtk
    '''
    logger.info('Starting parse_fastq_input')
    modelName = path.splitext(path.basename(options.hmm_model))[0]
    fastqPath = path.dirname(path.abspath(options.infiles[0])) # Assuming the path is the same to every input fastqfile
    if options.processes > cpu_count():
        options.processes = cpu_count()

    if not options.rerun:
        logger.info('Converting FASTQ to FASTA')
        for fastqfile in options.infiles:
            fastqBaseName = path.splitext(path.basename(fastqfile))[0]
            fastafile = '%s/%s.fasta' %(path.abspath(options.tmp_dir), fastqBaseName)
            if not path.isfile(fastafile):
                utils.convert_fastq_to_fasta(fastqfile, fastafile)
            elif path.getsize(fastafile) == 0:
                utils.convert_fastq_to_fasta(fastqfile, fastafile)

   
    p = Pool(options.processes)
    
    logger.info('Processing and searching input files. This may take a while...')

    try:
        bases_files = p.map(pooled_processing_fastq, zip((options.infiles), itertools.repeat(options)))  
    except KeyboardInterrupt:
        logger.warning('\nCaught a KeyboardInterrupt. Terminating...')
        p.terminate()
        p.join()
        exit()
        
    fastqDict = defaultdict(list)
    transformer = Transformer()
    transformer.find_file_difference(options.infiles[0], options.infiles[1])
    transformer.find_header_endings(options.infiles[0], options.infiles[1])
    transformer.verify_transform_is_working(options.infiles[0],options.infiles[1])

    logger.info('Retrieving hits from input files.')
    
    for fastqbase_hitfile in bases_files:
        fastqDict = utils.add_hits_to_fastq_dictionary(fastqbase_hitfile[1],
                fastqDict, fastqbase_hitfile[0], options, transformer)
    logger.info('Retrieving fastqfiles')
    utils.retrieve_paired_end_fastq(fastqDict, fastqPath, options, transformer) 
    
    if not options.no_quality_filtering:
        logger.info('Performing quality control')
        utils.quality(list(fastqDict.keys()), options)
    logger.info('Done')
    if not options.no_assembly:
        logger.info('Running assembly using SPAdes')
        utils.run_spades(options)
        logger.info('Done')
        logger.info('Running retrieval of assembled genes.')
        retrievedContigs,hits = utils.retrieve_assembled_genes(options)
        if path.isfile(retrievedContigs):
            logger.info('Predicting ORFS.')
            elongatedFasta ='%s/%s-gene-elongated.fasta' %(path.abspath(options.tmp_dir), path.basename(retrievedContigs).rpartition('.')[0])
            orfFile = '%s/%s-long-orfs.fasta' %(options.tmp_dir, path.basename(retrievedContigs).rpartition('.')[0])
            utils.retrieve_surroundings(hits, retrievedContigs, elongatedFasta)
            predict_orfs_orfFinder(elongatedFasta, options.tmp_dir, orfFile, options.min_orf_length) 
            retrievedOrfs = utils.retrieve_predicted_orfs(options, orfFile)
            Results.predictedOrfs = Results.count_contigs(retrievedOrfs)
        Results.retrievedContigs = Results.count_contigs(retrievedContigs)

def pooled_processing_fastq(fastqfile_options):
    # Cannot send logger object to functions run in a multiprocessing Pool.
    logger = logging.getLogger(__name__ + '.pooled_processing_fastq') 
    print((logger.handlers))
    try:
        fastqfile, options = fastqfile_options[0], fastqfile_options[1]
        modelName = path.splitext(path.basename(options.hmm_model))[0]
        fastqBaseName = path.splitext(path.basename(fastqfile))[0]
        fastqFilesBaseName = path.basename(fastqfile)
        fastafile = '%s/%s.fasta' %(path.abspath(options.tmp_dir), fastqBaseName)
        hmmOut = '%s/%s-%s-hmmsearched.out' %(path.abspath(options.hmm_out_dir), fastqBaseName, modelName)
        hitFile = '%s/%s-positives.out' %(path.abspath(options.tmp_dir), fastqBaseName)
        logger.info('Converting fastq to fasta')
        if options.rerun:
            peptideFile ='%s/%s-amino.fasta' %(options.amino_dir, fastqBaseName)
            if path.isfile(peptideFile):
                logger.info('Performing hmmsearch')
                utils.perform_hmmsearch(peptideFile, options.hmm_model, hmmOut, options)
            else:
                fastafile = '%s/%s.fasta' %(options.fasta_dir, fastqBaseName)
                logger.info('Translating and searching')
                utils.translate_and_search(fastafile, options.hmm_model, hmmOut, options)

        elif options.store_peptides:
            logger.info('Translating')
            peptideFile ='%s/%s-amino.fasta' %(path.abspath(options.tmp_dir), fastqBaseName)
            frame = '6'
            if not options.rerun:
                utils.translate_sequence(fastafile, peptideFile, options, frame)
            logger.info('Performing hmmsearch')
            utils.perform_hmmsearch(peptideFile, options.hmm_model, hmmOut, options)
        else:
            logger.info('Translating and searching')
            utils.translate_and_search(fastafile, options.hmm_model, hmmOut, options)

        logger.info('Start to classify')
        utils.classifier(hmmOut, hitFile, options)
        logger.info('Translating, searching, and classification done')
        
        fastqPath = path.dirname(path.abspath(fastqfile)) # Assuming the path is the same to every input fastqfile
        return fastqFilesBaseName, hitFile
    except KeyboardInterrupt:
        raise KeyboardInterruptError()


class KeyboardInterruptError(Exception): pass


if __name__ == '__main__':
    main()
