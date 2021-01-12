#!/usr/bin/env python3

'''
Reformats <.react> files into <.csv> files where that may be useful.
Requires the base <.fasta> the <.react> files were generated with.
'''

#Imports
import os
import argparse
import sf3libs.sf3io as sfio

#Functions
def write_out_csv(sequence,values,outfile='stats.csv'):
    '''Writes out the data'''
    with open(outfile,'w') as g:
        g.write(','.join(['Position','Nucleotide','Reactivity'])+'\n')
        lines = zip(range(1,len(sequence)+1),sequence,values)
        for line in lines:
            g.write(','.join([str(x) for x in line])+'\n')

def batch_write_out_csv(sequence_dict,values_dict,name_suffix,outdyr):
    '''batch writes a lot of csvs'''
    q_keys = set(sequence_dict.keys()).intersection(set(values_dict.keys()))
    for key in q_keys:
        new_fyle = os.path.join(outdyr,'_'.join([key,name_suffix])+'.csv')
        sequence, values = sequence_dict[key],values_dict[key]
        if len(sequence) == len(values):
            with open(new_fyle,'w') as g:
                g.write(','.join(['Position','Nucleotide','Reactivity'])+'\n')
                lines = zip(range(1,len(sequence)+1),sequence,values)
                for line in lines:
                    g.write(','.join([str(x) for x in line])+'\n')

def main():
    parser = argparse.ArgumentParser(description='Reformats reactivity values to <.csv> format')
    in_files = parser.add_argument_group('Input')
    in_files.add_argument('react',type=str,help = 'File to pull values from')
    in_files.add_argument('fasta',type=str,help = 'File used to generate the <.react>')
    settings = parser.add_argument_group('Settings')
    settings.add_argument('-mode',type=str.upper,default=None,choices=['S','M'],help='Single or Multi-Transcript mode')
    settings.add_argument('-transcript',default = None,type=str,help='[S] Specific Transcript to reformat',dest='name')
    settings.add_argument('-restrict',default=None,metavar='<.txt>',help='[M] Limit analysis to these specific transcripts')
    settings.add_argument('-outdir',type=str,default = None, help='[M] Name of the out directory, overrides default')
    args = parser.parse_args()

    #Read in files
    sequences,reactivities = sfio.read_fasta(args.fasta),sfio.read_react(args.react)

    #Single Transcript Mode
    if args.mode == 'S':
        #Get data
        outname= '_'.join([args.name,args.react.replace('.react','')])+'.csv'
        seq = sequences[args.name] if args.name in sequences else None
        reacts = reactivities[args.name] if args.name in reactivities else None
        #Write
        if seq and reacts:
            if len(seq) == len(reacts):
                write_out_csv(seq,reacts,outname)
            else:
                print('Sequence length does not match number of reactivities for {}'.format(args.transcript))
        else:
            print('One or more of the files did not contain entry {}'.format(args.transcript))

    #Multi Mode
    elif args.mode == 'M':
        #Directory
        new_dir = args.outdir if args.outdir else args.react.replace('.react','')+'_'+'all_csvs'
        if args.restrict:
            covered = sfio.read_restrict(args.restrict)
            sequences = {n:s for n,s in sequences.items() if n in covered}
            reactivities = {n:s for n,s in reactivities.items() if n in covered}
        if new_dir not in os.listdir('.'):
            os.mkdir(new_dir)
            batch_write_out_csv(sequences,reactivities,args.react.replace('.react',''),new_dir)
        else:
            print('{} folder already exists. Remove or rename and try again.'.format(new_dir))

    else:
        print('Check options! Must have a valid mode!')

if __name__ == '__main__':
    main()
