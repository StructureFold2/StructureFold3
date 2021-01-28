#!/usr/bin/env python3

#Imports
import argparse
import sf3libs.sf3io as sfio

#Functions
def fasta_stepwise(sequences,size,step):
    '''Breaks sequences into steps'''
    bits = {}
    for name, seq in sequences.items():
        mini = {'_'.join([name,str(i+1),str(i+size)]):seq[i:i+size] for i in range(0,len(seq)-(size-1),step)}
        bits.update(mini)
    return bits

#Workflow
def main():
    parser = argparse.ArgumentParser(description='Splits sequences into sliding window segments')
    in_files = parser.add_argument_group('Input')
    in_files.add_argument('fasta',type=str,help='Input Fasta File')
    settings = parser.add_argument_group('Settings')
    settings.add_argument('-size',default=10,type=int,help='[default = 10] Size of the Windows')
    settings.add_argument('-step',default=5,type=int,help='[default = 5] Step of the Windows')
    out_files = parser.add_argument_group('Output')
    out_files.add_argument('-name',default=None,help='Specify output file name')
    out_files.add_argument('-sort_file',action='store_true',help='Sort the output by name')
    args = parser.parse_args()

    #Get Sequences
    seqs = sfio.read_fasta(args.fasta)

    #Create Segments
    segements = fasta_stepwise(seqs,args.size,args.step)

    #Output
    base = sfio.rm_ext(args.fasta,'.fasta','.fas','.fa')
    default_name = '_'.join([base,str(args.size)+'window',str(args.step)+'step'])+'.fasta'
    out_name = sfio.check_extension(args.name,'.fasta') if args.name else default_name
    sfio.write_fasta(segements,out_name,args.sort_file)

if __name__ == '__main__': 
    main()