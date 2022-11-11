import numpy as np
from plot_cross_validation import plot_cross_validation
#from suggest_cutoffs import corresponding_scores
import time
import argparse
import shlex, subprocess
#from Estimator import Estimator
#import pandas as pd
from os import path


def summarize_sens_or_spec(est_obj,options,only_sensitivity):
    """ Read files """
    if only_sensitivity:
        ref_scores = np.loadtxt(est_obj.sens_score_file)
        num_ref_sequences = int((get_read_information(options.reference_sequences))[0].strip())
        name1 = 'sensitivity'
        name2 = name1
    else:
        ref_scores = np.loadtxt(est_obj.spec_score_file)
        num_ref_sequences = int((get_read_information(options.negative_sequences))[0].strip())
        name1 = 'specificity'
        name2 = '1-specificity'

    full_length = est_obj.full_length
    resultsdir = est_obj.resultsdir
    max_fpr = est_obj.get_max_fpr()
    if full_length:
        results_file = '%s/resulting_%s_%s_full_length.txt' \
                %(path.abspath(resultsdir),name1,options.modelname)
    else:
        results_file = '%s/resulting_%s_%s_%s.txt' \
                %(path.abspath(resultsdir),name1,options.modelname,str(options.fragment_lengths[0]))

    if not full_length:
        num_ref_seq_runs = num_ref_sequences * options.num_fragments
    else:
        num_ref_seq_runs = num_ref_sequences

    ''' Add missed (low scoring) elements to score arrays '''
    min_score = ref_scores.min()
    max_score = ref_scores.max()

    missed_ref = min_score*np.ones(num_ref_seq_runs - len(ref_scores),dtype=float)
    ref_scores = np.append(missed_ref,ref_scores)

    """ Create score array and initialize sensitivity and specificity arrays """
    dt = 0.1
    scores = np.arange(min_score-dt,max_score+dt,dt)
    scores = np.arange(0,max_score+dt,dt)
    size = len(scores)
    ref_fraction = np.empty(size,float)
    score_per_aa = np.empty(size,float)

    ref_max= float(ref_scores.max())
    len_ref_scores = float(len(ref_scores))

    """ Calculate sensitivity and 1-specificity as a function of score """
    for i,score in enumerate(scores):
        ref_fraction[i] = np.divide(np.sum(np.greater(ref_scores,score)),len_ref_scores) 
        score_per_aa[i] = np.divide(score,float(options.fragment_lengths[0]))

    max_ref = np.amax(ref_fraction)

    """ Calculate corresponding sensitivity and specificity for a given score """
    potential_cutoffs = np.arange(0,max_score,1)
    corresponding_sensitivity(potential_cutoffs,scores,ref_fraction,name2,results_file)
    

def calculate_performance(est_obj,preferred_sensitivity,options):
    """ Read files """
    ref_scores = np.loadtxt(est_obj.sens_score_file)
    neg_scores = np.loadtxt(est_obj.spec_score_file)
    full_length = est_obj.full_length
    resultsdir = est_obj.resultsdir
    max_fpr = est_obj.get_max_fpr()
    if full_length:
        results_file = '%s/resulting_sensitivity_specificity_%s_full_length.txt' \
                %(path.abspath(resultsdir),options.modelname)
        figfile = '%s/resulting_sensitivity_specificity_%s_full_length.png' \
                %(path.abspath(resultsdir),options.modelname)
    else:
        results_file = '%s/resulting_sensitivity_specificity_%s_%s.txt' \
                %(path.abspath(resultsdir),options.modelname,str(options.fragment_lengths[0]))
        figfile = '%s/resulting_sensitivity_specificity_%s_%s.png' \
                %(path.abspath(resultsdir),options.modelname,str(options.fragment_lengths[0]))

    num_ref_sequences = int((get_read_information(options.reference_sequences))[0].strip())
    num_neg_sequences = int((get_read_information(options.negative_sequences))[0].strip())
    if not full_length:
        num_ref_seq_runs = num_ref_sequences * options.num_fragments
        num_neg_seq_runs = num_neg_sequences * options.num_fragments
    else:
        num_ref_seq_runs, num_neg_seq_runs = num_ref_sequences, num_neg_sequences

    ''' Add missed (low scoring) elements to score arrays '''
    min_score = min([ref_scores.min(),neg_scores.min()])
    max_score = max([ref_scores.max(),neg_scores.max()])

    missed_neg = min_score*np.ones(num_neg_seq_runs - len(neg_scores),dtype=float)
    neg_scores = np.append(missed_neg,neg_scores)
    missed_ref = min_score*np.ones(num_ref_seq_runs - len(ref_scores),dtype=float)
    ref_scores = np.append(missed_ref,ref_scores)

    """ Create score array and initialize sensitivity and specificity arrays """
    dt = 0.1
    scores = np.arange(min_score-dt,max_score+dt,dt)
    scores = np.arange(0,max_score+dt,dt)
    size = len(scores)
    ref_fraction = np.empty(size,float)
    neg_fraction = np.empty(size,float)
    score_per_aa = np.empty(size,float)

    neg_max = float(neg_scores.max())
    ref_max= float(ref_scores.max())

    len_ref_scores = float(len(ref_scores))
    len_neg_scores = float(len(neg_scores))

    """ Calculate sensitivity and 1-specificity as a function of score """
    for i,score in enumerate(scores):
        ref_fraction[i] = np.divide(np.sum(np.greater(ref_scores,score)),len_ref_scores) 
        neg_fraction[i] = np.divide(np.sum(np.greater(neg_scores,score)),len_neg_scores) 
        score_per_aa[i] = np.divide(score,float(options.fragment_lengths[0]))


    max_ref = np.amax(ref_fraction)
    max_neg = np.amax(neg_fraction)

    """ Calculate corresponding specificity and score for a given sensitivity """
    sensitivity, specificity,suggested_score  = find_score_for_max_fpr(ref_fraction,neg_fraction,scores,max_fpr) 

    """ Calculate corresponding sensitivity and specificity for a given score """
    potential_cutoffs = np.arange(0,max_score,1)
    
    if est_obj.full_length:
        corresponding_sens_spec_whole(potential_cutoffs,scores,ref_fraction,neg_fraction,results_file)
        print('\nFor a preferred sensitivity of %s and max false positive rate of %s,\n'\
        'the suggested minimal score and corresponding sensitivity and specificity are:\n'\
        'score=%s\tsensitivity=%.4f\tspecificity=%.4f'\
        %(preferred_sensitivity,max_fpr,suggested_score,sensitivity,1-specificity))
    else:
        corresponding_sens_spec_frag(potential_cutoffs,scores,ref_fraction,neg_fraction,score_per_aa,results_file)
        print('\nFor a max false positive rate of %s,\n'\
        'the suggested minimal score, score/AA and corresponding sensitivity and specificity are:\n'\
        'score=%.2f\tscore/AA=%.4f\tsensitivity=%.4f\tspecificity=%.4f'\
        %(max_fpr,float(suggested_score),suggested_score/(float(options.fragment_lengths[0])),sensitivity,1-specificity))
    
    plot_cross_validation(scores,ref_fraction,neg_fraction,figfile)

def find_score_for_preferred_sensitivity(ref_fraction,neg_fraction,scores,preferred_sensitivity,max_fpr):
    tpr = float(preferred_sensitivity)
    neg = max_fpr + 1
    dt = 0.1
    score = 0
    while neg >= max_fpr:
        (ref,neg,score) =  corresponding_scores(ref_fraction,neg_fraction,scores,tpr,dt)
        tpr = tpr - 0.01
    return ref,neg,score

def find_score_for_max_fpr(ref_fraction,neg_fraction,scores,max_fpr):
    sens = 0.1
    return fpr_corresponding_scores(ref_fraction, neg_fraction,scores,max_fpr,sens)

def corresponding_sens_spec_frag(potential_cutoffs,scores,ref_fraction,neg_fraction,score_per_aa,results_file):
    with open(results_file,'w') as f:
        f.write("Score\tScore/AA\tSensitivity\tSpecificity\n")
        for i in potential_cutoffs:
            index = np.where(scores==i)[0][0]
            if i < 10:
                f.write("%0.2f\t%.4f\t%.4f\t%.4f\n" %(i,score_per_aa[index],ref_fraction[index],1-neg_fraction[index]))
            else:
                f.write("%0.1f\t%.4f\t%.4f\t%.4f\n" %(i,score_per_aa[index],ref_fraction[index],1-neg_fraction[index]))

def corresponding_sens_spec_whole(potential_cutoffs,scores,ref_fraction,neg_fraction,results_file):
    with open(results_file,'w') as f:
        f.write("Score\tSensitivity\tSpecificity\n")
        for i in potential_cutoffs:
            index = np.where(scores==i)[0][0]
            if i < 10:
                f.write("%0.2f\t%.4f\t%.4f\n" %(i,ref_fraction[index],1-neg_fraction[index]))
            else:
                f.write("%0.1f\t%.4f\t%.4f\n" %(i,ref_fraction[index],1-neg_fraction[index]))

def corresponding_sensitivity(potential_cutoffs,scores,ref_fraction,name,results_file):
    with open(results_file,'w') as f:
        f.write("Score\t%s\t\n" %name)
        for i in potential_cutoffs:
            index = np.where(scores==i)[0][0]
            if i < 10:
                f.write("%0.2f\t%.4f\n" %(i,ref_fraction[index]))
            else:
                f.write("%0.1f\t%.4f\n" %(i,ref_fraction[index]))




def fpr_corresponding_scores(ref_fraction, neg_fraction,scores,fpr,sens):
    i =  np.where(np.logical_and(neg_fraction > fpr-sens, neg_fraction < fpr+sens))[0][:]
    closest_i = (np.absolute(fpr - neg_fraction[i])).argmin()
    final_i = i[closest_i]
    return ref_fraction[final_i], neg_fraction[final_i],scores[final_i]

def corresponding_scores(ref_fraction, neg_fraction,scores,tpr,sens):
    i =  np.where(np.logical_and(ref_fraction > tpr-sens, ref_fraction < tpr+sens))[0][:]
    closest_i = (np.absolute(tpr - ref_fraction[i])).argmin()
    final_i = i[closest_i]
    return ref_fraction[final_i], neg_fraction[final_i],scores[final_i]

def get_read_information(fastafile):
    call_list = ''.join(['grep -c "^>" ',fastafile])
    commands = shlex.split(call_list)
    return subprocess.Popen(commands, stdin=subprocess.PIPE,
            stderr=subprocess.PIPE,stdout=subprocess.PIPE).communicate()    

if __name__ == '__main__':
    parser = argparse.ArgumentParser()                                            
    parser.add_argument('--reference-sequences','-rin',dest='reference_sequences')
    parser.add_argument('--negative-sequences','-nin',dest='negative_sequences')
    parser.add_argument('--modelname',dest = 'modelname')                         
    parser.add_argument('--lengths',dest='fragment_lengths')                      
    parser.add_argument('--num-fragments',dest='num_fragments')                   
                                                                              
#NUM_FRAGMENTS = 10000                                                        
#fragment_length_list = range(33,43,10)                                       
#FRAGMENT_LENGTH = 33                                                         
    full_seq = False                                                              
                                                                              
    parser.set_defaults(fragment_lengths = [33],                                  
        num_fragments = 10000)                                                
    args = parser.parse_args()                                                    

    ref_seq_score_file = '/storage/fannyb/tmp-files/test/test_33_sensitivity_scores.txt'
    neg_seq_score_file = '/storage/fannyb/tmp-files/test/test_33_specificity_scores.txt'
    est = Estimator('name',True,ref_seq_score_file,neg_seq_score_file,None)

    num_ref_seq_runs = 20*10000;
    num_neg_seq_runs = 38*20000;

    preferred_sensitivity = 1.0
    max_fpr = 0.1
    results_file = '/storage/fannyb/tmp-files/test/resulting_sensitivity_and_specificity.txt'
    figfile = '/storage/fannyb/tmp-files/testfig.png'
    calculate_performance(est, preferred_sensitivity,max_fpr,results_file,figfile,args)
