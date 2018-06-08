from utils import read_fasta
#from Bio.Seq import Seq
import shlex
import subprocess as sp
from os import path
import logging

def run_prodigal(infile,outfile):
    'prodigal -i infile -f gff -o genes.gff'
    call_list = ''.join(['prodigal -i ',infile,' -p meta -f gff -o ',outfile])
    commands = shlex.split(call_list)
    try:
        sp.Popen(commands, stdin=sp.PIPE,stderr=sp.PIPE).communicate()
    except OSError as e:
        logging.error("OS error ({0}) : {1}\nCan't find prodigal in path".format(e.errno,e.strerror))
        print "Can't find prodigal in path"

def parse_prodigal(prodigal_file,min_orf_length):
    orfs = {}
    if not path.isfile(prodigal_file):
        return
    with open(prodigal_file,'r') as f:
        for line in f:
            if not line.startswith('#'):
                line = line.split()
                info = line[8].split(';')[1]
                if info[-2:]=='00':
                    seq_id,start,end,strand = line[0],int(line[3]),int(line[4]),line[6]
                    length = int(end) - int(start)
                    if length >= min_orf_length:
                        if not seq_id in orfs:
                            orfs[seq_id] = [start,end,strand,length]
                        else:
                            if length > (orfs[seq_id])[3]:
                                orfs[seq_id] = [start,end,strand,length]
    return orfs;

def retrieve_orfs(orfs,fastaFile,orfFile):
    orfFile = open(orfFile,'w')
    for header,seq in read_fasta(fastaFile):
        header = header.split()[0]
        if header in orfs:
            if orfs[header][2] == '+':
                orfFile.write('>%s\n%s\n' %(header,seq[orfs[header][0]-1:orfs[header][1]]))
            else:
                rev_seq = reverse_complement(seq[orfs[header][0]-1:orfs[header][1]])
                orfFile.write('>%s\n%s\n' %(header,rev_seq))


def reverse_complement(sequence):
    comp = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G',\
            'a': 't', 't': 'a', 'g': 'c', 'c': 'g'}
    return ''.join(reversed([comp.get(base,base) for base in sequence]))
    
def predict_orfs_prodigal(infile,outdir,orfFile,minLength):
    basename = path.basename(infile).rpartition('.')[0]
    outdir = path.abspath(outdir)
    prodigalOut = '%s/%s-predicted-orfs.gff' %(outdir,basename)
    #orfFile = '%s/%s-predicted-orfs.fasta' %(outdir,basename)
    run_prodigal(infile,prodigalOut)
    orfs = parse_prodigal(prodigalOut,minLength)
    retrieve_orfs(orfs,infile,orfFile)

def predict_orfs_orfFinder(infile,tmpdir,orfFile,minLength):
    basename = path.basename(infile).rpartition('.')[0]
    tmpdir = path.abspath(tmpdir)
    orfFinderOut = '%s/%s-predicted-orfs.fasta' %(tmpdir,basename)
    #orfFile = '%s/%s-predicted-orfs.fasta' %(outdir,basename)
    run_ORFFinder(infile,orfFinderOut)
#    orfs = parse_prodigal(prodigalOut,minLength)
    if path.isfile(orfFinderOut) and path.getsize(orfFinderOut) > 0:
        parse_orfs(orfFinderOut,orfFile,minLength)


def parse_orfs(orfFinderOut,orfFile,minLength):
    orfOut = open(orfFile,'w')
    for header, seq in read_fasta(orfFinderOut,False):
        if len(seq) > minLength:
            orfOut.write('>%s\n%s\n' %(header,seq))

def run_ORFFinder(infile,orfFile):
    call_list = ''.join(['ORFfinder -in ',infile,
                        ' -outfmt 1 -ml 200 -s 1 -g 11 -out ',orfFile])
    commands = shlex.split(call_list)
    try:
        sp.Popen(commands, stdin=sp.PIPE,stderr=sp.PIPE).communicate()
    except OSError as e:
        logging.error("OS error ({0}) : {1}\nCan't find ORFfinder in path".format(e.errno,e.strerror))
        print "Can't find ORFfinder in path"


if __name__=='__main__':
    paths = '/storage/fannyb/runs_for_other/nachiket_pathogen_pool/b3/'
    prod_file = paths + 'contigs_prodigal.gff'
    #prod_file = '/storage/fannyb/runs_for_other/nachiket_pathogen_pool/b1_b2/tmpdir/contigs_prodigal.gff'
    min_orf = 200
    fastaFile = '/storage/fannyb/runs_for_other/nachiket_pathogen_pool/b3/retrievedContigs/retrieved-contigs.fasta'
    #fastaFile = '/storage/fannyb/runs_for_other/nachiket_pathogen_pool/b1_b2/retrievedContigs/retrieved-contigs.fasta'
    #
    #fastaFile = '/storage/fannyb/tmp-files/test_prodigal/test2/alistipes.fasta'
    #fastaFile = '/storage/fannyb/tmp-files/test/trylongreads/teredinibacter_extended.fasta'
    orfFile = '/storage/fannyb/runs_for_other/nachiket_pathogen_pool/b3/retrievedContigs/orf_prodigal.fasta'
    outdir = '/storage/fannyb/runs_for_other/nachiket_pathogen_pool/b3/tmpdir/'

    predict_orfs_prodigal(fastaFile,outdir,orfFile,min_orf)
#    run_prodigal(fastaFile,prod_file)
#    orfs = parse_prodigal(prod_file,min_orf)
#    retrieve_orfs(orfs,fastaFile,orfFile)
