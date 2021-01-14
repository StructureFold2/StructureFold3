#!/usr/bin/env python3

'''
Calculates reverse transcriptase stop specificity on one or more <.rtsc> files
'''

#Imports
import argparse
import collections
import sf3libs.sf3io as sfio

#Functions
def collect_specificity(fyle_lyst,fasta):
    '''Generates a nested dictionary of specificities'''
    ref = sfio.read_fasta(fasta)
    all_data = {f.replace('.rtsc',''):rtsc_specificity(f,ref) for f in fyle_lyst}
    return all_data

def rtsc_specificity(fyle,fasa_dict):
    '''Generates a file specific specificity dictionary'''
    data,counts = sfio.read_rtsc(fyle),collections.Counter()
    for transcript, stops in data.items():
        seq = fasa_dict[transcript].upper()[:-1]
        matched_values = zip(list(seq),stops[1:])
        for pair in matched_values:
            counts[pair[0]]+=int(pair[1])
    return counts

def write_specificity(adict,outfyle,allowed,rnumb=5):
    '''Writes all the specificity data to a <.csv>'''
    cols = ['file','base','count','specificity']
    with open(outfyle,'w') as g:
        g.write(','.join(cols)+'\n')
        for fyle,sub in sorted(adict.items()):
            for letter,value in sorted(sub.items()):
                if letter in allowed:
                    line = [fyle,letter,value,round(value/sum(sub.values()),rnumb)]
                    g.write(','.join([str(x) for x in line])+'\n')

#Workflow
def main():
    parser = argparse.ArgumentParser(description='Analyzes native/reagent nucleotide RT stop specificity')
    in_files = parser.add_argument_group('Input')
    in_files.add_argument('index',type=str,help="Fasta used to generate the <.rtsc>")
    in_files.add_argument('rtsc',default = None, help='Input <.rtsc>',nargs='+')
    settings = parser.add_argument_group('Settings')
    settings.add_argument('-report',default='ACGT',help='Include these nucelotides in report')
    settings.add_argument('-round',type=int,default=5,help='[default = 5] Decimal places to report',dest='digits')
    out_files = parser.add_argument_group('Output')
    out_files.add_argument('-name',default = None, help='Specify output file name')
    args = parser.parse_args()

    #Outfile Nomenclature
    base_name = sorted([f.replace('.rtsc','') for f in args.rtsc])
    suffixes = [args.report,'specificity']
    default_name = '_'.join(base_name+suffixes)+'.csv'
    out_fyle = sfio.check_extension(args.name,'.csv') if args.name else default_name

    #Calculate Specificity
    specificity_data = collect_specificity(args.rtsc,args.index)

    #Writeout
    write_specificity(specificity_data,out_fyle,args.report,args.digits)

if __name__ == '__main__': 
    main()
