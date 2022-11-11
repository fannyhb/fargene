import argparse
import glob

from utils import read_fasta

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i','--input_files',required=True)
    parser.add_argument('--length',type=int,default=600)
    parser.add_argument('-o','--output',default='longsequences.fasta')
    parser.add_argument('--cut_stars',default=False,action='store_true')

    args = parser.parse_args()
    files = glob.glob(args.input_files)
    outputfile = open(args.output,'w')
    pick_long_reads(files,args.length,outputfile,args)

def pick_long_reads(files,length,outputfile,args):
    totalcount = 0
    passedcount = 0

    for fastafile in files:
        for header, seq in read_fasta(fastafile, False):
            totalcount = totalcount + 1
            if len(seq) >= length:
                longread = True
                if args.cut_stars:
                    seq = seq.split('*')
                    tmp_length = 0
                    for i in range(0,len(seq)):
                        if len(seq[i]) > tmp_length:
                            tmp_length = len(seq[i])
                            j = i
                    if tmp_length < length:
                        print(header, tmp_length)
                        longread = False
                    seq = seq[j]
                if longread:
                    passedcount = passedcount + 1
                    if header.startswith("Contig"):
                        filename = (fastafile.split("/")[-1]).split(".")[0]
                        header = filename + "_" + header
                    outputfile.write('>%s\n%s\n' %(header,seq))

    outputfile.close()
    print('searched %s sequences.\n \
            %s of them were longer than %s' \
            %(totalcount,passedcount,length))

if __name__=='__main__':
    main()
