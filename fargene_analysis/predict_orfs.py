from os import path
import shlex
import subprocess as sp
import logging

from .utils import read_fasta

def run_prodigal(infile,outfile):
    'prodigal -i infile -f gff -o genes.gff'
    call_list = ''.join(['prodigal -i ',infile,' -p meta -f gff -o ',outfile])
    commands = shlex.split(call_list)
    try:
        sp.run(commands, stderr=sp.PIPE)
    except OSError as e:
        logging.error("OS error ({0}) : {1}\nCan't find prodigal in path".format(e.errno,e.strerror))
        print("Can't find prodigal in path")

# Minimum prodigal confidence (its own 0-100 "conf=" estimate) required to
# accept a gene call that's truncated at one end of the input sequence, when
# allow_partial=True. This isn't the final filter on the candidate - it just
# needs to be good enough to bother translating and hmmsearching downstream,
# where the real resistance-gene HMM score threshold does the actual
# classification. 50 (better than a coin flip, by prodigal's own estimate)
# favors recall here, consistent with the rest of the pipeline's E-value=1000
# ORF-search thresholds, which are deliberately loose ahead of that final
# HMM classification step.
PARTIAL_ORF_MIN_CONFIDENCE = 50.0

def parse_prodigal(prodigal_file,min_orf_length,allow_partial=False,min_confidence=PARTIAL_ORF_MIN_CONFIDENCE):
    orfs = {}
    if not path.isfile(prodigal_file):
        return
    with open(prodigal_file,'r') as f:
        for line in f:
            if not line.startswith('#'):
                line = line.split()
                attrs = dict(field.split('=',1) for field in line[8].split(';') if '=' in field)
                partial = attrs.get('partial','00')
                if partial != '00':
                    if not allow_partial:
                        continue
                    # Truncated at both ends means neither a start nor a stop
                    # codon was actually observed - too speculative to keep.
                    if partial == '11':
                        continue
                    if float(attrs.get('conf',0)) < min_confidence:
                        continue
                seq_id,start,end,strand = line[0],int(line[3]),int(line[4]),line[6]
                length = int(end) - int(start)
                if length >= min_orf_length:
                    if not seq_id in orfs:
                        orfs[seq_id] = [start,end,strand,length]
                    else:
                        if length > (orfs[seq_id])[3]:
                            orfs[seq_id] = [start,end,strand,length]
    return orfs

def retrieve_orfs(orfs,fastaFile,orfFile):
    with open(orfFile,'w') as orfOut:
        for header,seq in read_fasta(fastaFile):
            header = header.split()[0]
            if header in orfs:
                if orfs[header][2] == '+':
                    orfOut.write('>%s\n%s\n' %(header,seq[orfs[header][0]-1:orfs[header][1]]))
                else:
                    rev_seq = reverse_complement(seq[orfs[header][0]-1:orfs[header][1]])
                    orfOut.write('>%s\n%s\n' %(header,rev_seq))


def reverse_complement(sequence):
    comp = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G',\
            'a': 't', 't': 'a', 'g': 'c', 'c': 'g'}
    return ''.join(reversed([comp.get(base,base) for base in sequence]))

def predict_orfs_prodigal(infile,outdir,orfFile,minLength,allow_partial=False):
    basename = path.basename(infile).rpartition('.')[0]
    outdir = path.abspath(outdir)
    prodigalOut = '%s/%s-predicted-orfs.gff' %(outdir,basename)
    run_prodigal(infile,prodigalOut)
    orfs = parse_prodigal(prodigalOut,minLength,allow_partial)
    retrieve_orfs(orfs,infile,orfFile)

def predict_orfs_orfFinder(infile,tmpdir,orfFile,minLength):
    basename = path.basename(infile).rpartition('.')[0]
    tmpdir = path.abspath(tmpdir)
    orfFinderOut = '%s/%s-predicted-orfs.fasta' %(tmpdir,basename)
    run_ORFFinder(infile,orfFinderOut)
    if path.isfile(orfFinderOut) and path.getsize(orfFinderOut) > 0:
        parse_orfs(orfFinderOut,orfFile,minLength)


def parse_orfs(orfFinderOut,orfFile,minLength):
    with open(orfFile,'w') as orfOut:
        for header, seq in read_fasta(orfFinderOut,False):
            if len(seq) > minLength:
                seq = find_common_startcodon(seq,minLength)
                orfOut.write('>%s\n%s\n' %(header,seq))

def find_common_startcodon(seq,minLength):
    startCodons = ['ATG','GTG','TTG']
    start = 0
    codon_stop = 3
    gotAstart = False
    while not gotAstart:
        if len(seq[start:]) < minLength:
            gotAstart = True
            start = 0
        elif seq[start:codon_stop] in startCodons:
            gotAstart = True
        else:
            start = start + 3
            codon_stop = codon_stop + 3
    return seq[start:]

def run_ORFFinder(infile,orfFile):
    call_list = ''.join(['ORFfinder -in ',infile,
                        ' -outfmt 1 -ml 200 -s 1 -g 11 -out ',orfFile])
    commands = shlex.split(call_list)
    try:
        sp.run(commands, stderr=sp.PIPE)
    except OSError as e:
        logging.error("OS error ({0}) : {1}\nCan't find ORFfinder in path".format(e.errno,e.strerror))
        print("Can't find ORFfinder in path")
