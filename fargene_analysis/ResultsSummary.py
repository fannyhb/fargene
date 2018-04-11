import shlex
import subprocess as sp
from os import path


class ResultsSummary(object):

    def __init__(self, summaryFile, numInputFiles,hmmModel):
        self.outfile = summaryFile
        self.hmmerModel = path.basename(hmmModel)
        self.inputFiles = numInputFiles
        self.retrievedSequences = 0
        self.contigs = None
        self.retrievedContigs = None
        self.predictedOrfs = 0


    def count_hits_deprecated(self,infiles,hitDir):
        count = 0
        for fastafile in infiles:
            fastaBaseName = path.basename(fastafile).rpartition('.')[0] 
            hitFile = '%s/%s-positives.out' %(path.abspath(hitDir),fastaBaseName)
            msg = 'wc -l %s' %hitFile
            commands = shlex.split(msg)
            c,_ = sp.Popen(commands, stdin=sp.PIPE,stdout = sp.PIPE).communicate()
            count = count + int(c.split()[0])
        self.retrievedSequences = count

    def count_hits(self,hitFile):
        msg = 'wc -l %s' %hitFile
        commands = shlex.split(msg)
        c,_ = sp.Popen(commands, stdin=sp.PIPE,stdout = sp.PIPE).communicate()
        self.retrievedSequences = self.retrievedSequences + int(c.split()[0])


    def count_contigs(self,contigFile):
        count = 0
        if contigFile==None or not path.isfile(contigFile):
            return count
        with open(contigFile) as f:
            for line in f:
                if line.startswith('>'):
                    count = count + 1
        return count


    def count_orfs_genomes(self,orfFile):
        if not path.isfile(orfFile):
            return
        with open(orfFile) as f:
            for line in f:
                if line.startswith('>'):
                    self.predictedOrfs = self.predictedOrfs + 1

    def write_summary(self,isMeta):
        #outfile = '%s/ResultsSummary.txt' %(path.abspath(outdir))
        outfile = self.outfile
        with open(outfile,'w') as f:
            f.write('Summary Results Novel Gene Finder\n'\
                    'The used HMM-model was:\t%s\n'\
                    'Number of input files:\t%s\n' %(self.hmmerModel,str(self.inputFiles)))
            if isMeta:
                f.write('\nMetagenomic Input\n------------------------\n'\
                        'Number of retrieved contigs:\t%s\n'\
                        'Number of predicted ORFs:\t%s'\
                        %(str(self.retrievedContigs),str(self.predictedOrfs)))
            else:
                f.write('\nGenomic Input\n------------------------\n'\
                        'Number of predicted genes:\t%s\n'\
                        'Number of predicted ORFs:\t%s'\
                        %(str(self.retrievedSequences),str(self.predictedOrfs)))

if __name__ == '__main__':
    infiles = ['/storage/fannyb/tmp-files/test/trylongreads/zobellia_wholegenome.fasta']
    hitDir = '/storage/fannyb/tmp-files/test/trylongreads/test/tmpdir'
    contigFile = '/storage/fannyb/tmp-files/test/oilspill/retrievedContigs/retrieved-contigs.fasta'
    orfFile = '/storage/fannyb/tmp-files/test/oilspill/retrievedContigs/retrieved-contigs-predicted-orfs.fasta'
    orgContigFile = '/storage/fannyb/tmp-files/test/oilspill/spades_assembly/contigs.fasta'
    outdir = '/storage/fannyb/tmp-files/test/oilspill'
    Results = ResultsSummary('txt',20,'B1.hmm')
    Results.count_hits(infiles,hitDir)
    Results.contigs = Results.count_contigs(orgContigFile)
    Results.retrievedContigs = Results.count_contigs(contigFile)
    Results.predictedOrfs = Results.count_contigs(orfFile)
    print Results.contigs
    print Results.retrievedContigs
    print Results.predictedOrfs
    Results.write_summary(outdir,False)
