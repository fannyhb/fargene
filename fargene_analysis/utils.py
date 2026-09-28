#!/usr/bin/env python3

import shlex
import subprocess as sp
from collections import defaultdict
from os.path import basename, splitext, abspath, isfile, isdir, getsize
from os import makedirs
from multiprocessing import Pool, cpu_count
import itertools
import glob
import logging

def convert_fastq_to_fasta(fastqInfile,fastaOutfile):
    msg = 'seqtk seq -a %s' %(fastqInfile)
    commands = shlex.split(msg)
    with open(fastaOutfile,'w') as fastaOut:
        sp.run(commands, stdout=fastaOut, stderr=sp.PIPE)

def translate_sequence(infile,aminofile,options,frame):
    msg = 'transeq %s %s -frame=%s -table=11 sformat=%s'\
            %(infile,aminofile,frame,options.trans_format)
    commands = shlex.split(msg)
    sp.run(commands, stderr=sp.PIPE)

def perform_hmmsearch(aminofile,hmmModel,hmmOutfile,options):
    if options.sensitive:
        flag = '--max'
    else:
        flag = ''
    msg = ' hmmsearch --domtblout %s -E 1000 --domE 1000 %s %s %s' \
            % (hmmOutfile, flag,hmmModel,aminofile)
    tmpfile = '%s/hmm_tmp.out' %(abspath(options.tmp_dir))
    commands = shlex.split(msg)
    with open(tmpfile,'w') as tmp:
        sp.run(commands, stdout=tmp, stderr=sp.PIPE)
    logging.info('Running command: %s' %(msg))

def translate_and_search(infile,hmmModel,hmmOutfile,options):
    # Equivalent to `cat infile | transeq ... | hmmsearch ... - > tmpout`,
    # but without spawning a shell and a `cat` process just to feed infile
    # into transeq's stdin.
    tmpout = '%s/tmp.out' %abspath(options.tmp_dir)
    if options.sensitive:
        flag = '--max'
    else:
        flag = ''
    transeq_cmd = shlex.split('transeq -filter -frame=6 -table=11 sformat=pearson')
    hmmsearch_cmd = shlex.split('hmmsearch %s -E 1000 --domE 1000 --domtblout %s %s -'
            % (flag, hmmOutfile, hmmModel))
    msg = 'cat %s | %s | %s > %s' % (infile, ' '.join(transeq_cmd), ' '.join(hmmsearch_cmd), tmpout)
    logging.info('Running command: %s' %(msg))
    with open(infile) as inp, open(tmpout,'w') as out:
        transeq_proc = sp.Popen(transeq_cmd, stdin=inp, stdout=sp.PIPE)
        sp.run(hmmsearch_cmd, stdin=transeq_proc.stdout, stdout=out)
        transeq_proc.stdout.close()
        transeq_proc.wait()


def classifier(hmmOutfile,hitFile,options):
    '''
    Returns the sequence id,sequence length (in peptides),
    env_start and env_end of the sequences classified as positives
    in a file
    '''
    # Equivalent to piping hmmOutfile through
    # `grep -v '^#' | awk '<threshold>' | awk '{print $1,$3,$20,$21}'`,
    # done in-process instead of spawning grep/awk/awk for every input file.
    logging.info('Classifying hits from %s into %s' %(hmmOutfile,hitFile))
    with open(hmmOutfile) as infile, open(hitFile,'w') as outfile:
        for line in infile:
            if line.startswith('#'):
                continue
            fields = line.split()
            score = float(fields[13])
            if options.meta:
                passed = score/(float(fields[20])-float(fields[19])) > float(options.meta_score)
            else:
                passed = score > float(options.long_score)
            if passed:
                outfile.write('%s %s %s %s\n' %(fields[0],fields[2],fields[19],fields[20]))


def add_hits_to_fastq_dictionary(hitFile,fastqDict,fastqInfile,options,transformer):
    if not transformer:
        fastqBaseName,sep,readNr = fastqInfile.rpartition('_')
    else:
        fastqBaseName = transformer.get_fastq_basename(fastqInfile)
        logging.debug('\nthe fastqinfile = ' + fastqInfile)
        logging.debug('\nthe basename = ' + fastqBaseName)
    with open(hitFile,'r') as f:
        for line in f:
            readID_tmp, length, start, end = line.split()
            readID,_,indexed_read = readID_tmp.rpartition('/')
            if len(readID) == 0: #For reads that not ends with /1_3
                readID = readID_tmp
                if not options.protein:
                    readID,sep,frame = readID.rpartition('_')
            if not fastqBaseName in fastqDict or not readID in fastqDict[fastqBaseName]:
                fastqDict[fastqBaseName].append(readID)
    return fastqDict

def create_file_with_ids(nameOfIdFile,listOfIds,transformer):
    if not transformer:
        nameOfIdFile = nameOfIdFile + '.txt'
        with open(nameOfIdFile,'w') as f:
            for hit in listOfIds:
                f.write('%s\n' %(hit))
    else:
        nameOfIdFile = ['%s_%s.txt' %(nameOfIdFile,str(i)) for i in range(1,3)]
        with open(nameOfIdFile[0],'w') as f1, open(nameOfIdFile[1],'w') as f2:
            for hit in listOfIds:
                headers = transformer.get_full_read_header(hit)
                f1.write('%s\n' %(headers[0]))
                f2.write('%s\n' %(headers[1]))
    return nameOfIdFile


def extract_fastq_from_file(nameOfIdFile,fastqInfile,fastqOutfile):
    call_list = ''.join(['seqtk subseq ',fastqInfile,' ',nameOfIdFile])
    logging.info('Running command: ' + call_list)
    commands = shlex.split(call_list)
    with open(fastqOutfile,'w') as fastqOut:
        sp.run(commands, stdout=fastqOut, stderr=sp.PIPE)

def retrieve_paired_end_fastq(fastqDict,fastqPath,options,transformer):
    tmpfile = abspath(options.tmp_dir) + '/listOfIds'
    name,_,endsuffix = basename(options.infiles[0]).rpartition('.')
    endsuffix = '.%s' %(endsuffix)
    for key, item in fastqDict.items():
        nameOfIdFile = create_file_with_ids(tmpfile, item,transformer)
        if not transformer:
            fastqBase = abspath(fastqPath) + '/' + key
            fastqInfiles = ['%s_%s%s' %(fastqBase,str(i),endsuffix) for i in range(1,3)]
        else:
            fastqnames = transformer.get_full_fastq_filename(key,endsuffix)
            fastqInfiles = ['%s/%s' %(abspath(fastqPath),fastqnames[i]) for i in range(0,2)]
        fastqOutfiles = ['%s/%s_%s_retrieved.fastq' %(abspath(options.res_dir),key,str(i)) for i in range(1,3)]
        i = 0
        for fastqInfile,fastqOutfile in zip(fastqInfiles,fastqOutfiles):
            if isfile(fastqInfile):
                if not transformer:
                    extract_fastq_from_file(nameOfIdFile,fastqInfile,fastqOutfile)
                else:
                    extract_fastq_from_file(nameOfIdFile[i],fastqInfile,fastqOutfile)
                    i = i+1

def quality(fastqBases,options):
    if options.processes > cpu_count():
        options.processes = cpu_count()
    with Pool(options.processes) as p:
        p.starmap(quality_control_and_adapter_removal, zip(fastqBases,itertools.repeat(options)))

def quality_control_and_adapter_removal(fastqBase,options):
    fastqOutfiles = ['%s/%s_%s_retrieved.fastq' %(abspath(options.res_dir),fastqBase,str(i)) for i in range(1,3)]
    msg = 'trim_galore --paired %s %s -q 30 --output_dir %s' \
            %(fastqOutfiles[0],fastqOutfiles[1],options.trimmed_dir)
    sp.run(msg, shell=True)

def run_spades(options):
    retrievedFastqGrouped = ['%s/all_retrieved_%s.fastq' %(abspath(options.res_dir),str(i))
            for i in range(1,3)]
    if not options.no_quality_filtering:
        sourcepatterns = ['%s/*val_%s.fq*' %(abspath(options.trimmed_dir),str(i)) for i in range(1,3)]
    elif not glob.glob('%s/*_1_retrieved.fq*' %(options.res_dir)):
        sourcepatterns = ['%s/*_%s_retrieved.fastq*' %(abspath(options.res_dir),str(i)) for i in range(1,3)]
    else:
        sourcepatterns = ['%s/*_%s.fq*' %(abspath(options.res_dir),str(i)) for i in range(1,3)]
    msg = ['for f in %s; do case "$f" in *.gz) gunzip -c "$f";; *) cat "$f";; esac; done > %s'
            %(pattern,outfile) for pattern,outfile in zip(sourcepatterns,retrievedFastqGrouped)]
    for command in msg:
        sp.run(command, shell=True)

    tmp_spades_out = '%s/spades_out.txt' %(abspath(options.tmp_dir))
    spades_msg = 'spades.py --meta -1 %s -2 %s -o %s > %s'\
            %(retrievedFastqGrouped[0],retrievedFastqGrouped[1],options.assembly_dir,tmp_spades_out)
    if isfile(retrievedFastqGrouped[0]) and getsize(retrievedFastqGrouped[0]) > 0:
        sp.run(spades_msg,shell=True)
    else:
        msg = 'No retrieved data to assemble'
        print('\n%s\n' %msg)
        logging.error(msg)

def create_dictionary(hitFile,options):
    hitDict = defaultdict(list)
    with open(hitFile,'r') as f:
        for line in f:
            name, length, start, end = line.split()
            if not options.protein:
                name,sep,frame = name.rpartition('_')
            else:
                frame = '-'
            if name in hitDict:
                hitDict[name].append((length,start,end,frame))
            else:
                hitDict[name] = [(length,start,end,frame)]
    return hitDict

def retrieve_fasta(hitDict,fastaInfile,fastaOutfile,options):
    fastaBaseName = splitext(basename(fastaInfile))[0]
    if not hitDict:
        msg = 'No hits in file %s' %(abspath(fastaInfile))
        print('\n%s\n' %msg)
        logging.error(msg)
        return
    with open(fastaOutfile,'a') as outfile:
        for header, seq in read_fasta(fastaInfile,False):
            written = False
            s_id = header.split()[0]
            if s_id in hitDict:
                header = '>' + fastaBaseName + '_' + header
                for i in range(0,len(hitDict[s_id])):
                    info = hitDict[s_id][i]
                    if options.retrieve_whole and not written:
                        outfile.write('%s\n%s\n' %(header,seq))
                        written = True
                    elif not options.retrieve_whole:
                        ali_start, ali_end = int(info[1]),int(info[2])
                        if not options.protein:
                            ali_start, ali_end = translate_position(info[1],info[2],info[3],info[0],len(seq))
                        outfile.write('%s\n%s\n' %(header,seq[ali_start:ali_end]))

def retrieve_peptides(hitDict,aminoInFile,aminoOut,options):
    fastaBaseName = splitext(basename(aminoInFile))[0]
    if not hitDict:
        return
    with open(aminoOut,'a') as outfile:
        for header, seq in read_fasta(aminoInFile,False):
            written = False
            s_id,sep,frame = (header.split()[0]).rpartition('_')
            if s_id in hitDict:
                header = '>' + fastaBaseName + '_' + header
                for i in range(0,len(hitDict[s_id])):
                    info = hitDict[s_id][i]
                    if (not options.protein and info[3]==frame) or options.protein:
                        if options.retrieve_whole and not written:
                            outfile.write('%s\n%s\n' %(header,seq))
                            written = True
                        elif not options.retrieve_whole:
                            ali_start, ali_end = int(info[1]),int(info[2])
                            outfile.write('%s\n%s\n' %(header,seq[ali_start:ali_end]))

def make_fasta_unique(fastaout,options):
    tmp_fastaout = '%s/fastaout_tmp.fasta' %(abspath(options.tmp_dir))
    header_dictionary = {}
    with open(tmp_fastaout,'w') as f:
        for header, seq in read_fasta(fastaout,False):
            header = header.split()[0]
            if header in header_dictionary:
                header_dictionary[header] = header_dictionary[header] + 1
            else:
                header_dictionary[header] = 1
            header = '>%s_seq%s' %(header,str(header_dictionary[header]))
            f.write('%s\n%s\n' %(header,seq))
    return tmp_fastaout

def retrieve_assembled_genes(options):
    options.meta = False
    frame = '6'
    modelName = splitext(basename(options.hmm_model))[0]
    contigFile = '%s/contigs.fasta' %(abspath(options.assembly_dir))
    aminoFile = '%s/contigs-translated.fasta' %(abspath(options.tmp_dir))
    fastaOut = '%s/retrieved-contigs.fasta' %(abspath(options.final_gene_dir))
    aminoOut = '%s/retrieved-contigs-peptides.fasta' %(abspath(options.final_gene_dir))
    hmmOut = '%s/contigs-%s-hmmsearched.out' %(abspath(options.hmm_out_dir),modelName)
    hitFile = '%s/contigs-positives.out' %(abspath(options.tmp_dir))
    if isfile(contigFile):
        translate_sequence(contigFile,aminoFile,options,frame)
        perform_hmmsearch(aminoFile,options.hmm_model,hmmOut,options)
        classifier(hmmOut,hitFile,options)
        hitDict = create_dictionary(hitFile,options)

        retrieve_peptides(hitDict,aminoFile,aminoOut,options)
        options.retrieve_whole = True
        retrieve_fasta(hitDict,contigFile,fastaOut,options)
        return (fastaOut,hitDict)
    else:
        logging.error('The file %s does not exist' %(contigFile))
        return ('', None)

def retrieve_predicted_orfs(options,orfFile):
    options.meta = False
    options.retrieve_whole = True
    frame ='1'
    modelName = splitext(basename(options.hmm_model))[0]
    aminoFile = '%s/orfs-translated.fasta' %(abspath(options.tmp_dir))
    fastaOut = '%s/predicted-orfs.fasta' %(abspath(options.final_gene_dir))
    aminoOut = '%s/predicted-orfs-amino.fasta' %(abspath(options.final_gene_dir))
    hmmOut = '%s/orfs-%s-hmmsearched.out' %(abspath(options.hmm_out_dir),modelName)
    hitFile = '%s/orfs-positives.out' %(abspath(options.tmp_dir))
    if isfile(orfFile) and getsize(orfFile) > 0:
        translate_sequence(orfFile,aminoFile,options,frame)
        perform_hmmsearch(aminoFile,options.hmm_model,hmmOut,options)
        hitDict = orf_classifier(hmmOut,hitFile,options)

        retrieve_peptides(hitDict,aminoFile,aminoOut,options)
        retrieve_fasta(hitDict,orfFile,fastaOut,options)
        return fastaOut
    else:
        logging.error('The file %s does not exist' %(orfFile))
        return fastaOut

def orf_classifier(hmmOut,hitFile,options):
    # Equivalent to piping hmmOut through
    # `grep -v '^#' | awk '$14>threshold' | awk '{print $1,$3,$20,$21,$14}'`,
    # done in-process instead of spawning grep/awk/awk.
    logging.info('Classifying ORF hits from %s into %s' %(hmmOut,hitFile))
    with open(hmmOut) as infile, open(hitFile,'w') as outfile:
        for line in infile:
            if line.startswith('#'):
                continue
            fields = line.split()
            if float(fields[13]) > float(options.long_score):
                outfile.write('%s %s %s %s %s\n' %(fields[0],fields[2],fields[19],fields[20],fields[13]))

    hitDict = defaultdict(list)
    # Tracks, per shortName, the most recently inserted hitDict key with that
    # shortName - equivalent to what a fresh linear rescan of hitDict.keys()
    # would find each time (nothing is ever removed from hitDict), but O(1)
    # instead of O(n) per line.
    last_identifier_by_shortname = {}
    with open(hitFile,'r') as f:
        for line in f:
            name, length, start, end, score = line.split()
            score = float(score)
            shortName = name.split(':')[0]
            if not options.protein:
                name,sep,frame = name.rpartition('_')
            else:
                frame = '-'
            previousOrf = last_identifier_by_shortname.get(shortName)
            if previousOrf is not None:
                if hitDict[previousOrf][0][4] < score:
                    hitDict[name] = [(length,start,end,frame,score)]
                    last_identifier_by_shortname[shortName] = name
            else:
                hitDict[name] = [(length,start,end,frame,score)]
                last_identifier_by_shortname[shortName] = name
    return hitDict



def retrieve_predicted_genes_as_amino(options,retrievedNucFile,aminoOut,frame):
    options.meta = False
    options.retrieve_whole = True
    modelName = splitext(basename(options.hmm_model))[0]
    aminoTmpFile = '%s/retrieved-translated.fasta' %(abspath(options.tmp_dir))
    hmmOut = '%s/retrieved-genes-%s-hmmsearched.out' %(abspath(options.hmm_out_dir),modelName)
    hitFile = '%s/retrieved-genes-positives.out' %(abspath(options.tmp_dir))
    if isfile(retrievedNucFile):
        translate_sequence(retrievedNucFile,aminoTmpFile,options,frame)
        perform_hmmsearch(aminoTmpFile,options.hmm_model,hmmOut,options)
        classifier(hmmOut,hitFile,options)
        hitDict = create_dictionary(hitFile,options)
        retrieve_peptides(hitDict,aminoTmpFile,aminoOut,options)
    else:
        logging.error('The file %s does not exist' %(retrievedNucFile))

def retrieve_surroundings(hitDict,fastaInfile,elongatedFastaOutfile):
    extension = 200
    fastaBaseName = splitext(basename(fastaInfile))[0]
    if not hitDict:
        msg = 'No hits in file %s' %(abspath(fastaInfile))
        print('\n%s\n' %msg)
        logging.error(msg)
        return
    if fastaBaseName == 'retrieved-contigs':
        addition = 'contigs_'
    else:
        addition = ''
    with open(elongatedFastaOutfile,'w') as outfile:
        for header, seq in read_fasta(fastaInfile,False):
            nlen = len(seq)
            s_id = header.split()[0]
            s_id = s_id.lstrip(addition)
            if s_id in hitDict:
                for i in range(0,len(hitDict[s_id])):
                    header = '>' + s_id + '_seq' + str(i+1)
                    info = hitDict[s_id][i]
                    ali_start, ali_end = int(info[1]),int(info[2])
                    ali_start, ali_end = translate_position(info[1],info[2],info[3],info[0],nlen)
                    ali_start, ali_end = include_surroundings(ali_start,ali_end,nlen,extension)
                    outfile.write('%s\n%s\n' %(header,seq[ali_start:ali_end+1]))

def is_fasta(infile):
    with open(infile,'r') as f:
        first = f.readline()
    if first.startswith('>'):
        return True
    return False

def is_fastq(infile):
    with open(infile,'r') as f:
        first = f.readline()
    if first.startswith('@'):
        return True
    return False

def remove_tmp_file(fileToRemove):
    msg = "rm " + fileToRemove
    logging.info("Removing file " + fileToRemove)
    sp.run(msg, shell=True)

def create_dir(directory):
    if not isdir(abspath(directory)):
        try:
            makedirs(directory)
        except OSError as e:
            logging.critical('OS error ({0}): {1}\nCould not create directory {2}.\
                    \nExiting pipeline'.format(e.errno, e.strerror,directory))
            exit()

def decide_min_ORF_length(hmmModel):
    with open(hmmModel,'r') as f:
        for line in f:
            line = line.split()
            if line[0].startswith('LENG'):
                return round((0.9*float(line[1])*3))

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

def translate_position(a_start, a_end, frame, alen, nlen):
    frame = int(frame)
    a_start = int(a_start)
    a_end = int(a_end)
    alen = int(alen)
    nlen = int(nlen)
    if frame < 4:
        start = 3*(a_start-1) + frame
        end = 3*(a_end -1) + frame + 2
        if end > nlen:
            end = nlen
    else:
        frame = frame - 3
        if a_end == alen:
            start = 1
            a_end = alen - a_start + 1
            end = 3*(a_end -1) + frame + 2
        else:
            tmp = alen - a_end + 1
            a_end = alen - a_start + 1
            a_start = tmp
            start = 3*(a_start -1) + frame
            end = 3*(a_end -1) + frame + 2
            if not a_start == 1:
                start = start -3
    if end > nlen:
        end = nlen
    if start < 1 or start == 2:
        start = 1
    return (start-1,end-1)

def include_surroundings(ali_start, ali_end, nlen, extension):
    if ali_start < ali_end:
        start = max(ali_start - extension,1)
        end = min(ali_end + extension,nlen)
    else:
        end = min(ali_start + extension,nlen)
        start = max(ali_end-extension,1)
    return (start-1,end-1)

def remove_files(targetDir, resDir):
    files = glob.glob(targetDir + '/*') + \
    glob.glob(resDir + '/*.fastq') + \
    glob.glob(resDir + '/trimmedReads/*.fq')
    if len(files) == 0:
        return True
    print("\nThe following files will be DELETED!!\n\n {}\n".format(
            '\n'.join(files)))
    ans = str(input("type [y]es to continue or [n]o to abort:"))
    if ans.lower() == 'y' or ans.lower() == 'yes':
        for fileToRemove in files:
            remove_tmp_file(fileToRemove)
        return True
    else:
        return False
