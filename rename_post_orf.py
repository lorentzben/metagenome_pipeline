#!/usr/bin/env python3 
import argparse
from Bio import SeqIO

def rename(infile, outfile):
   
    with open (infile) as original, open(outfile, 'w') as renamed:
        records = SeqIO.parse(infile, 'fasta')
        for record in records:
            seq_id = record.description.split("#")[4].split(';')[0].split('=')[1].strip() 
            record.id = seq_id +"_"+record.id
            record.description = seq_id
            SeqIO.write(record, renamed, 'fasta')

def main(args):
    if args.outseq:
        outname = args.outseq
    else:
        outname = args.inseq[:-4]+"_renamed.fna"
    if args.inseq: 
        rename(args.inseq, outname)
    
    

if __name__ == "__main__":
    # Build Argument Parser in order to facilitate ease of use for user
    parser = argparse.ArgumentParser(
        description="Checks dependencies are installed, and validates a manifest file for qiime 2")
    parser.add_argument('-i', '--in', action='store', required=True,
                        help="name of the sequence file to renamep", dest='inseq')
    parser.add_argument('-o', "--out", action='store', required=False,
                        help="name of the fasta outfile once renamed", dest="outseq")
    parser.add_argument('-v', '--version', action='version',
                        version='%(prog)s 1.0')

    args = parser.parse_args()
    # print(args)
    main(args)
