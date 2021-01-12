#!/usr/bin/env python3

'''
This script will generate a convenient <.csv> detailing the nucleotide or amino-acid
composition of sequences in fasta format. This <.csv> may be imported and merged with
other <.csv>s to build a comprehensive data sheet.
'''

#Imports
import argparse
import collections
import sf3libs.sf3io as sfio

#Functions
def generate_composition(fasta):
    '''Creates a frequency and length dictionary from a fasta'''
    return {k:collections.Counter(v) for k,v in sfio.read_fasta(fasta).items()}

def write_comp_csv(data,outfyle,thres,colname,rnumb=7):
    '''Writes out a composition report'''
    overall = sum(data.values(),collections.Counter())
    total,thres2 = sum(overall.values()),thres/100
    fracts = {k:float(v)/total for k, v in overall.items()}
    passing = sorted(dict(filter(lambda x: x[1]>=thres2, fracts.items())).keys())
    header = ','.join([colname]+[key+'_percent' for key in passing]+['length'])
    with open(outfyle,'w') as g:
        g.write(header+'\n')
        for seq_name,seq_data in data.items():
            sub_total = sum(seq_data.values())
            fracts = [round(seq_data[r]/sub_total,rnumb) for r in passing]
            outline = [str(x) for x in [seq_name]+fracts+[sub_total]]
            outline = ','.join(outline)
            g.write(outline+'\n')

#Workflow
def main():
    parser = argparse.ArgumentParser(description='Creates a composition report for a  <.fasta> file')
    in_files = parser.add_argument_group('Input')
    in_files.add_argument('fasta',default = None, help = 'Input Fasta File')
    settings = parser.add_argument_group('Settings')
    settings.add_argument('-threshold',type=float,default=1.0,metavar='<number>',help='[default = 1.0] Minimum percent to report')
    settings.add_argument('-description',type=str,default='transcript',metavar='<name>',help='[default = transcript] Sequence column heading')
    out_files = parser.add_argument_group('Output')
    out_files.add_argument('-affix',type=str,default='composition', help='[default = composition] <.csv> file affix')
    args = parser.parse_args()

    #Generate Frequencies
    seq_data = generate_composition(args.fasta)

    #Generate Outname
    base_name = '.'.join(args.fasta.split('.')[:-1])
    affixes = [args.affix,str(args.threshold)]
    outname = '_'.join([base_name]+affixes)+'.csv'

    #Write Out
    write_comp_csv(seq_data,outname,args.threshold,args.description)


if __name__ == '__main__':
    main()
