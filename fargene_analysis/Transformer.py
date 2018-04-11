import difflib
from os import path

class Transformer(object):

    def __init__(self):
        self.fastqFileDiff1 = None
        self.fastqFileDiff2 = None
        self.headerEnd1 = ''
        self.headerEnd2 = ''
        #self.endsuffix = endsuffix
        self.fastqFileDiffPos = None

    def find_file_difference(self,left_reads,right_reads):
        left_reads = path.basename(left_reads)
        right_reads = path.basename(right_reads)
        differ = difflib.Differ()
        res = list(differ.compare(left_reads,right_reads))
        pos = []
        vals = []
        for i,val in enumerate(res):
            if val.startswith('+'):
                pos.append(i-1)
                vals.append(val[2])
        self.fastqFileDiff1 = (left_reads)[pos[0]:]  
        self.fastqFileDiff2 = (right_reads)[pos[0]:]  
        self.fastqFileDiffPos = len(right_reads) - pos[0]

    def find_header_endings(self,left_reads,right_reads):
        with open(left_reads,'r') as f:
            header1 = f.readline().split()[0].strip()
            f.close()
        with open(right_reads,'r') as f:
            header2 = f.readline().split()[0].strip()
            f.close()
        if header1.endswith('/1') and header2.endswith('/2'):
            self.headerEnd1 = '/1'
            self.headerEnd2 = '/2'
        else:
            self.headerEnd1 = ''
            self.headerEnd2 = ''

    def get_fastq_basename(self,fastqname):
        return fastqname[0:-self.fastqFileDiffPos]
        
    def get_full_read_header(self,header):
        return [header + self.headerEnd1, header + self.headerEnd2]

    def get_full_fastq_filename(self,fastqBaseName,endsuffix):
        return [fastqBaseName + self.fastqFileDiff1, fastqBaseName + self.fastqFileDiff2]

    def verify_transform_is_working(self,left_read,right_read):
        fastqfile1 = path.basename(left_read)
        fastqfile2 = path.basename(right_read)
        fastqPath = path.dirname(path.abspath(left_read))

        base = self.get_fastq_basename(fastqfile1)
        transformed_fastqfiles = self.get_full_fastq_filename(base,'fq')
        fullFastqFile = fastqPath + '/' + transformed_fastqfiles[0]

        try:
            with open(fullFastqFile,'r') as f:
                header = f.readline().strip()
        except IOError as e:
            errorMsg = 'Transformer.py: Transformation of fastq file names is not working\n'\
                    'Infile: {0}\nTransformed infile: {1}\n'\
                    'I/O error({2}): {3}'.format(left_read,fullFastqFile,e.errno,e.strerror)
            print errorMsg
            exit()
        trans_header = header.split()[0]
        readID,_,indexedRead = trans_header.rpartition('/')
        if len(readID) == 0:
            readID = trans_header
        transformedHeader = self.get_full_read_header(readID)
        transformedHeader[0] = transformedHeader[0]
        if not transformedHeader[0] == header.split()[0]:
            errorMsg = 'Transformer.py: Tranformation of fastq header is not working\n'\
                    'Original fastq header: {0}\n'\
                    'Transformed fastq header: {1}'.format(header.split()[0],transformedHeader[0])
            print errorMsg
            exit()

