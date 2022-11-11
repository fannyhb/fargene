import shlex, subprocess
from os import path, makedirs
import os
import argparse
#from read_fasta import read_fasta
from random import randint
import logging

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--reference-sequences','-rin',dest='reference_sequences')
    parser.add_argument('--negative-sequences','-nin',dest='negative_sequences')
    parser.add_argument('--modelname',dest = 'modelname')                        
    parser.add_argument('--lengths',dest='fragment_lengths')                     
    parser.add_argument('--num-fragments',dest='num_fragments')                  
    parser.add_argument('--create-fragments',dest = 'create_fragments')

    parser.set_defaults(fragment_lengths = [33],
            num_fragments = 10000,
            create_fragments = True)
    args = parser.parse_args()

    full_seq = False
    estimate_specificity(args,full_seq)


def estimate_specificity(args,est_obj):
    
    full_seq = est_obj.full_length
    #modelpath = './' + args.modelname + '/'    
    tmpdir = est_obj.tmpdir
    resultsdir = path.abspath(est_obj.resultsdir)
    modeldir = est_obj.hmmdir
    #tmpdir = modelpath + 'tmpdir/'             
    #modeldir = modelpath + 'hmmerModel/'
    #if not path.isdir(path.abspath(modelpath)):
    #    makedirs(modelpath)                    
    #if not path.isdir(path.abspath(tmpdir)):   
    #    makedirs(tmpdir)                        
    #if not path.isdir(path.abspath(modeldir)):   
    #    makedirs(modeldir)                        
    
    fastabasename = path.splitext(path.basename(args.reference_sequences))[0]
    alignfile = '%s/%s-aligned.fasta' %(path.abspath(tmpdir),fastabasename)
    hmmmodel = '%s/%s.hmm' %(path.abspath(modeldir),fastabasename)
    summarized_hmmsearchfile = '%s/%s-hmmsearch-negative-sequences-full-length.txt' %(path.abspath(est_obj.resultsdir),args.modelname)
    create_model(args.reference_sequences, alignfile, hmmmodel)

    print('Estimating specificity...')
    
    #if full_seq:
    #    args.create_fragments = False

    
    for fragment_length in args.fragment_lengths:
        #TODO handle full length sequences
        if not full_seq and args.create_fragments:
            fastabasename = path.splitext(path.basename(args.negative_sequences))[0]
            targetfile = '%s/%s-fragmented-%s.fasta' %(path.abspath(tmpdir),fastabasename,str(fragment_length))
            create_subsets(args.negative_sequences,targetfile,args.num_fragments,int(fragment_length))
        else:
            targetfile = args.negative_sequences

        hmmsearchfile = '%s/specificity_%s-hmmsearched.out' %(path.abspath(tmpdir),str(fragment_length))
        if full_seq:                                                                      
            resultsfile = '%s/%s_full_length_specificity_scores.txt' %(tmpdir,args.modelname)
        else:                                                                             
            resultsfile = '%s/%s_%s_specificity_scores.txt' %(tmpdir,args.modelname,str(fragment_length))
        
        run_hmmsearch(hmmmodel, targetfile, hmmsearchfile)
        
        extract_scores(hmmsearchfile,resultsfile)
        sort_scores(resultsfile)
        est_obj.spec_score_file = resultsfile
        if full_seq:
            move_tmp_file(hmmsearchfile,summarized_hmmsearchfile)
            #remove_tmp_files([hmmsearchfile])
        else:
            remove_tmp_files([targetfile,hmmsearchfile])

def extract_scores(hmmsearchfile,resultsfile):
    commands = "grep -v '^#' " + hmmsearchfile + " | awk '{print $1,$14}' > " + resultsfile
    msg = "Running command: " + commands
    logging.info(msg)
    subprocess.call(commands, shell=True)
    logging.info("Done")

def sort_scores(resultsfile):
    name_dic = {}
    with open(resultsfile,'r') as f:
        for line in f:
            line = line.split()
            name,score = line[0],float(line[1])
            if not name in name_dic:
                name_dic[name] = score
            elif name_dic[name] < score:
                name_dic[name] = score
    with open(resultsfile,'w') as f:
        for score in list(name_dic.values()):
            f.write('%f\n' %(score))

def create_fragments(sequence, num_fragments, fragment_length):
    fragments = []                                             
    sequence_length = len(sequence)                            
    max_start = sequence_length - fragment_length              
    for i in range(0,num_fragments):                           
        start = randint(0,max_start)                           
        fragments.append(sequence[start:start+fragment_length])
    return fragments                                           

def create_subsets(fastafile,fragmentfile,num_fragments,fragment_length):
    f = open(fragmentfile,'a')
    for header, seq in read_fasta(fastafile):
        fragments = create_fragments(seq,num_fragments,fragment_length)
        for i,fragment in enumerate(fragments):
            f.write('>%s_%s\n%s\n' %(header.split()[0],str(i),fragment))


def create_model(fastafile, alignfile, hmmfile):                
    devnull = open(os.devnull, 'w')
    call_list = ''.join(['clustalo -quiet -infile=',fastafile, \
            ' -align -outfile=',alignfile, ' -output=fasta'])   
    logging.info('Running command:\n%s' %(call_list))
    commands = shlex.split(call_list)                                                                        
    subprocess.Popen(commands, stdin=subprocess.PIPE,           
            stderr=subprocess.PIPE,stdout=devnull).communicate()               
    logging.info('Done')                                                            
    call_list = ''.join(['hmmbuild ',hmmfile,' ', alignfile])   
    commands = shlex.split(call_list)                           
    logging.info('Running command:\n%s' %(call_list))
    subprocess.Popen(commands, stdin=subprocess.PIPE,           
            stderr=subprocess.PIPE,stdout=devnull).communicate()               
    logging.info('Done')                                                            
                                                                
    call_list = ''.join(['hmmpress ', hmmfile])                 
    commands = shlex.split(call_list)                           
    logging.info('Running command:\n%s' %(call_list))
    subprocess.Popen(commands, stdin=subprocess.PIPE,           
            stderr=subprocess.PIPE,stdout=devnull).communicate()
    logging.info('Done')                                                            
    devnull.close()
    remove_tmp_files([alignfile])

def run_hmmsearch(hmmfile, fragmentfile, outputfile):                          
    call_list = ''.join(['hmmsearch --max -E 1000 --domE 1000 --domtblout ',outputfile,
            ' ', hmmfile, ' ', fragmentfile])                                          
    commands = shlex.split(call_list)                                                  
    with open(os.devnull,'w') as devnull:                                              
        subprocess.Popen(commands, stdin=subprocess.PIPE,                              
            stderr=subprocess.PIPE, stdout=devnull).communicate()

def remove_tmp_files(files_to_be_removed):               
    for tmp_file in files_to_be_removed:                 
        call_list = ''.join(['rm ', tmp_file])           
        commands = shlex.split(call_list)                
        subprocess.Popen(commands, stdin=subprocess.PIPE,
                stderr=subprocess.PIPE).communicate()    

def move_tmp_file(file_to_be_moved,new_name):               
        call_list = ''.join(['mv ', file_to_be_moved,' ', new_name])           
        commands = shlex.split(call_list)                
        subprocess.Popen(commands, stdin=subprocess.PIPE,
                stderr=subprocess.PIPE).communicate()    

def read_fasta(filename, keep_formatting=True):
    """Read sequence entries from FASTA file
    NOTE: This is a generator, it yields after each completed sequence.
    Usage example:
    for header, seq in read_fasta(filename):
        print ">"+header
        print seq
    """

    with open(filename) as fasta:
        line = fasta.readline().rstrip()
        if not line.startswith(">"):
            raise IOError("Not FASTA format? First line didn't start with '>'")
        if keep_formatting:
            sep = "\n"
        else:
            sep = ""
        first = True
        seq = []
        header = ""
        while fasta:
            if line == "": #EOF
                yield header, sep.join(seq)
                break
            elif line.startswith(">") and not first:
                yield header, sep.join(seq)
                header = line.rstrip()[1:]
                seq = []
            elif line.startswith(">") and first:
                header = line.rstrip()[1:]
                first = False
            else:
                seq.append(line.rstrip())
            line = fasta.readline()                                                         
    
if __name__=='__main__':
    main()
