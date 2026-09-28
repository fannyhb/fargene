from os import path
from random import randint
from multiprocessing import Pool, cpu_count
from shutil import which
import glob
import shlex, subprocess
import argparse
import os
import logging


def estimate_sensitivity(reference_sequences, est_obj,args):
    full_seq = est_obj.full_length
    modelpath = './' + args.modelname + '/'
    tmpdir = est_obj.tmpdir
    resultsdir = est_obj.resultsdir

    for FRAGMENT_LENGTH in args.fragment_lengths:
        if full_seq:
            target = 'for full length genes...'
        else:
            target = 'for fragments of lengths %s AA...' %(str(FRAGMENT_LENGTH))
        print('\nEstimating sensitivity %s' %target)

        create_subsets(reference_sequences,tmpdir, args.num_fragments, int(FRAGMENT_LENGTH), full_seq)
        modelfiles = glob.glob(tmpdir + "without*")
        fragmentfiles = glob.glob(tmpdir + "fragments-*")
        jobs = []
        for modelfile in modelfiles:
            modelbasename = path.basename(modelfile).split('-')[1]
            fragment = None
            for fragmentfile in fragmentfiles:
                if path.basename(fragmentfile).split('-')[1] == modelbasename:
                    fragment = fragmentfile
            jobs.append((modelfile, fragment, tmpdir))
        # Each leave-one-out model (align + hmmbuild + hmmpress + hmmsearch)
        # is independent of the others, so build/search them in parallel
        # instead of one at a time.
        with Pool(cpu_count()) as p:
            p.starmap(build_and_search_one_model, jobs)
        truehmm = glob.glob(tmpdir + "*-hmmsearched.out")
        if full_seq:
            resultsfile = '%s%s_full_length_sensitivity_scores.txt' %(tmpdir,args.modelname)
            sum_hmmsearch_file = '%s/%s-hmmsearch-reference-sequences-full-length.txt' %(path.abspath(resultsdir),args.modelname)
            extract_full_seq_hmm_info(truehmm,sum_hmmsearch_file)
        else:
            resultsfile = '%s%s_%s_sensitivity_scores.txt' %(tmpdir,args.modelname,str(FRAGMENT_LENGTH))
        pooled_sort_hmmerfiles(truehmm,resultsfile)
        est_obj.sens_score_file = resultsfile
        remove_tmp_files(truehmm)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--reference-sequences','-rin',dest='reference_sequences')
    parser.add_argument('--modelname',dest = 'modelname')
    parser.add_argument('--lengths',dest='fragment_lengths')
    parser.add_argument('--num-fragments',dest='num_fragments')
    full_seq = False
    parser.set_defaults(fragment_lengths = [33],
            num_fragments = 10000)
    args = parser.parse_args()
    estimate_sensitivity(args.reference_sequences,full_seq,args)

def create_fragments(sequence, num_fragments, fragment_length):
    fragments = []
    sequence_length = len(sequence)
    max_start = sequence_length - fragment_length
    for i in range(0,int(num_fragments)):
        start = randint(0,max_start)
        fragments.append(sequence[start:start+fragment_length])
    return fragments

def create_subsets(fastafile,subsetpath,num_fragments, fragment_length,full_seq):
    header_list = []
    unique_count = 1
    # Read the reference FASTA once and reuse it for every leave-one-out
    # subset, instead of re-reading the whole file from disk for each
    # sequence (was O(n) file reads for n sequences, now O(1)).
    wrapped_records = list(read_fasta(fastafile, True))
    for header, wrapped_seq in wrapped_records:
        seq = wrapped_seq.replace('\n','')
        headername = header.split()[0]
        headername = '_'.join(headername.split('|'))
        headername = ((headername.replace(',','_')).replace('[','_')).replace(']','_')
        headername = (headername.replace('=','_')).replace('"','_')
        if headername in header_list:
            headername = headername + '_' + str(unique_count)
            unique_count = unique_count + 1
        header_list.append(headername)
        with open(subsetpath + "fragments-" + headername + ".fasta",'w') as fragmentfile:
            if not full_seq:
                fragments = create_fragments(seq, num_fragments, fragment_length)
                for i,fragment in enumerate(fragments):
                    fragmentfile.write('>%s_%s\n%s\n' %(headername,str(i),fragment))
            else:
                fragmentfile.write('>%s\n%s\n' %(header,seq))

        with open(subsetpath + "without-" + headername + ".fasta","w") as subsetfile:
            for header2, seq2 in wrapped_records:
                if not header == header2:
                    subsetfile.write(">%s\n%s\n" %(header2,seq2))

def build_and_search_one_model(modelfile, fragment, tmpdir):
    alignfile = tmpdir + modelfile.split('/')[-1] + ".aligned"
    hmmfile = alignfile + ".hmm"
    outputfile = "%s%s-hmmsearched.out" %(tmpdir,path.basename(modelfile))
    create_model(modelfile, alignfile, hmmfile)
    run_hmmsearch(hmmfile, fragment, outputfile, tmpdir)
    remove_tmp_files([alignfile,modelfile,fragment])
    remove_tmp_files(glob.glob(hmmfile + "*"))

def create_model(fastafile, alignfile, hmmfile):
    clustalo_path = which('clustalo')
    if clustalo_path:
        logging.debug('Found clustalo at: %s', clustalo_path)
    else:
        logging.critical('Cannot find Clustal Omega (clustalo)!')
        exit()
    call_list = ''.join(['clustalo -quiet -infile=',fastafile, \
            ' -align -outfile=',alignfile, ' -output=fasta'])
    commands = shlex.split(call_list)
    with open(os.devnull,'w') as devnull:
        subprocess.run(commands, stderr=subprocess.PIPE, stdout=devnull)

        call_list = ''.join(['hmmbuild ',hmmfile,' ', alignfile])
        commands = shlex.split(call_list)
        subprocess.run(commands, stderr=subprocess.PIPE, stdout=devnull)

        call_list = ''.join(['hmmpress ', hmmfile])
        commands = shlex.split(call_list)
        subprocess.run(commands, stderr=subprocess.PIPE, stdout=devnull)


def run_hmmsearch(hmmfile, fragmentfile, outputfile, tmpdir):
    call_list = ''.join(['hmmsearch --max -E 1000 --domE 1000 --domtblout ',outputfile, \
            ' ', hmmfile, ' ', fragmentfile])
    commands = shlex.split(call_list)
    msg =  'Running command:\n%s' %(call_list)
    logging.info(msg)
    with open(os.devnull,'w') as devnull:
        subprocess.run(commands, stderr=subprocess.PIPE, stdout=devnull)
    logging.info('Done')

def sort_hmmerfiles(hmmerfileslist,outfile,full_seq):
    name_dic = {}
    for hmm in hmmerfileslist:
        with open(hmm,'r') as f:
            for line in f:
                if not line.startswith('#'):
                    line = line.split()
                    name = line[0]
                    score = float(line[13])
                    if not name in name_dic:
                        name_dic[name] = score
                    elif name_dic[name] < score:
                        name_dic[name] = score

    with open(outfile,'w') as out:
        for item in name_dic.values():
            out.write('%f\n' %(item))

def pooled_sort_hmmerfiles(hmmerfileslist,outfile):
    with Pool(cpu_count()) as p:
        scores = p.map(sort_one_hmmerfile,hmmerfileslist)
    with open(outfile,'w') as out:
        for score_list in scores:
            for item in score_list:
                out.write('%f\n' %(item))


def sort_one_hmmerfile(hmmfile):
    name_dic = {}
    with open(hmmfile,'r') as hmm:
        for line in hmm:
            if not line.startswith('#'):
                line = line.split()
                name = line[0]
                score = float(line[13])
                if not name in name_dic:
                    name_dic[name] = score
                elif name_dic[name] < score:
                    name_dic[name] = score
    return list(name_dic.values())

def extract_full_seq_hmm_info(hmmerfileslist,outfile):
    with open(outfile,'w') as out:
        with open(hmmerfileslist[0],'r') as f:
            for i,line in enumerate(f):
                if i < 3:
                    out.write(line)
        for hmmerfile in hmmerfileslist:
            with open(hmmerfile,'r') as f:
                for line in f:
                    if not line.startswith('#'):
                        out.write(line)


def sort_hmmerfiles_deprecated(hmmerfileslist,outfile,full_seq):
    score = []
    for hmm in hmmerfileslist:
        with open(hmm,'r') as f:
            for line in f:
                if not line.startswith('#'):
                    score.append(float(line.split()[13]))
                    if full_seq:
                        break
    score.sort()
    with open(outfile,'w') as out:
        for item in score:
            out.write('%f\n' %(item))

def remove_tmp_files(files_to_be_removed):
    for tmp_file in files_to_be_removed:
        call_list = ''.join(['rm ', tmp_file])
        commands = shlex.split(call_list)
        subprocess.run(commands, stderr=subprocess.PIPE)

def read_fasta(filename, keep_formatting=True):
    """Read sequence entries from FASTA file
    NOTE: This is a generator, it yields after each completed sequence.
    Usage example:
    for header, seq in read_fasta(filename):
        print(">"+header)
        print(seq)
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
