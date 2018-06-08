#!/usr/bin/env python2.7
import argparse
from sys import argv
from os import path, makedirs
from estimate_sensitivity import estimate_sensitivity
from estimate_specificity import estimate_specificity
from estimate_specificity import create_model
from MetaData import MetaData
from calculate_performance import calculate_performance, summarize_sens_or_spec
import logging

def parse_args(argv):
    desc = 'A program to create and optimize profile hidden Markov models'
    copyright = 'Copyright (c) Fanny Berglund 2018.'
    parser = argparse.ArgumentParser(description=desc + '. ' + copyright)
    parser.add_argument('--reference-sequences','-rin', dest='reference_sequences',
                        required=True,
                        help = 'The sequences that the model should be built of.')
    parser.add_argument('--negative-sequences','-nin', dest='negative_sequences',
                        required=True,
                        help='The sequences that should be used as the negative dataset.\n'\
                             'Should preferable be similar sequences but without the desired phenotype.')
    parser.add_argument('--output','-o',dest='output_dir',
                        help = 'The directory where the output should be saved.')
    parser.add_argument('--modelname',dest = 'modelname',
                        help = 'The name of the new model' )
    
    parser.add_argument('--fragment-lengths','-l',type=int,dest='fragment_lengths',
                        help = 'The length (aa) of the fragments that should be used to determine the threshold score \n'\
                                'for metagenomic input. (default: 33 AA)')
    parser.add_argument('--num-fragments',type=int,dest='num_fragments',
                        help = 'The number of fragments that should be created from each gene. '\
                                '(default: 10 000)')

    parser.add_argument('--only-sens',dest='specificity',action='store_false',
                        help = 'Should be used if only sensitivity of the model should be estimated.')
    parser.add_argument('--only-spec',dest='sensitivity',action='store_false',
                        help = 'Should be used if only the specificity of the model should be estimated.')
    parser.add_argument('--only-full-length',dest='only_full_length',action='store_true',
                        help = 'Should be used if you only want to optimize the threshold score for \n'\
                                'full length genes.')

    parser.add_argument('--only-fragments',dest='only_fragments',action='store_true',
                        help = 'Should be used if you only want to optimize the threshold score for \n'\
                                'full fragmented genes.')
#    parser.add_argument('--fragment-input',dest='create_fragments',action='store_false',
#                        help = 'If you already have fragments you want to test the specificity on.\n'\
#                                'The location of the file containing the fragments should then be '\
#                                'specified using the --negative-sequences argument')

    parser.set_defaults(modelname = 'new_model',
            output_dir = './model_validation',
            create_fragments = True,
            sensitivity = True,
            specificity = True,
            only_full_length = False,
            only_fragments = False,
            fragment_lengths = [33],
            num_fragments = 10000)
    options = parser.parse_args()
    return options


def main():
    
    options = parse_args(argv)
    if not isinstance(options.fragment_lengths,list):
        options.fragment_lengths = [options.fragment_lengths]
    tmpdir,resultsdir,hmmdir = create_dirs(options)
    logfile = '%s/model_creation_optimization.log' %(path.abspath(options.output_dir))
    logging.basicConfig(filename=logfile,filemode='w',
            format='%(asctime)s %(levelname)s %(message)s', level=logging.DEBUG)

    if options.only_full_length:
        full_lengths = [True]
    elif options.only_fragments:
        full_lengths = [False]
    else:
        full_lengths = [True, False]
    
    for i,full_length in enumerate(full_lengths):
        # TODO handle if the input is already fragments
        est_obj = MetaData(full_length,None,None,None,tmpdir,resultsdir,hmmdir)

        if options.sensitivity:
            estimate_sensitivity(options.reference_sequences,est_obj,options)
        if options.specificity:
            estimate_specificity(options, est_obj)

        if not options.sensitivity:
            summarize_sens_or_spec(est_obj,options,False)
        elif not options.specificity:
            summarize_sens_or_spec(est_obj,options,True)
        else:
            calculate_performance(est_obj,1,options)

    ''' If not specificity, create a HMM-model '''
    if not options.specificity:
        create_hmm(options,est_obj)


def create_dirs(options):
    if not path.isdir(path.abspath(options.output_dir)):
        makedirs(options.output_dir)
    tmpdir = '%s/tmpdir/' %(path.abspath(options.output_dir))
    if not path.isdir(path.abspath(tmpdir)):
        makedirs(tmpdir)
    resultsdir = '%s/results/' %(path.abspath(options.output_dir))
    if not path.isdir(path.abspath(resultsdir)):
        makedirs(resultsdir)
    hmmdir = '%s/hmmerModel/' %(path.abspath(options.output_dir))
    if not path.isdir(path.abspath(hmmdir)):
        makedirs(hmmdir)
    return tmpdir,resultsdir,hmmdir    


def create_hmm(args,est_obj):
    modeldir = est_obj.hmmdir
    tmpdir = est_obj.tmpdir
    if not path.isdir(path.abspath(modeldir)):
        makedirs(modeldir)
    fastabasename = path.splitext(path.basename(args.reference_sequences))[0]
    alignfile = '%s/%s-aligned.fasta' %(path.abspath(tmpdir),fastabasename)
    hmmmodel = '%s/%s.hmm' %(path.abspath(modeldir),fastabasename)
    create_model(args.reference_sequences, alignfile, hmmmodel)


if __name__ == '__main__':
#    options = parse_args()
    main()
