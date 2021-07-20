#!/usr/bin/env python3

#Imports
import glob
import argparse
import collections
import sf3libs.sf3io as sfio

#Functions
def rtsc_coverage(rtsc_file,fasta_index,specificity='AC'):
    '''Generates a coverage dictionary for an <.rtsc> file with a given specificity'''
    rtsc_data,coverage = sfio.read_rtsc(rtsc_file), {}
    for transcript,stops in rtsc_data.items():
        effective_sequence = fasta_index[transcript].upper()[:-1]
        effective_stops = stops[1:]
        matched = zip(effective_sequence,effective_stops)
        specific_bases_seq = sum([v for k, v in collections.Counter(effective_sequence).items() if k in specificity])
        specific_bases_stops = sum([item[1] for item in filter(lambda x: x[0] in specificity, matched)])
        try:
            coverage[transcript] = float(specific_bases_stops)/specific_bases_seq
        except ZeroDivisionError:
            coverage[transcript] = 0
    return coverage

def collect_coverages(fyle_list,fasta_fyle,specificity='AC'):
    '''Applies rtsc_coverage to files'''
    fasta_index,coverage_data = sfio.read_fasta(fasta_fyle),{}
    for fyle in fyle_list:
        coverage_data[fyle.replace('.rtsc','')] = rtsc_coverage(fyle,fasta_index,specificity)
    return coverage_data

def write_coverage(data,out_fyle='derp.csv'):
    '''Writes out coverages'''
    h_keys = sorted(data.keys())
    v_keys = sorted(list(set.union(*map(set, data.values()))))
    header = ','.join(['transcript']+[derp+'_coverage' for derp in h_keys])
    with open(out_fyle,'w') as g:
        g.write(header+'\n')
        for transcript in v_keys:
            entry = [str(data[h_key].get(transcript,'NA')) for h_key in h_keys]
            out_line = ','.join([transcript]+entry)
            g.write(out_line+'\n')

def write_ol(data,out_fyle='derp.txt',threshold=1.0):
    '''Writes out an overlap file'''
    h_keys = sorted(data.keys())
    shared_transcripts = sorted(list(set.intersection(*map(set, data.values()))))
    with open(out_fyle,'w') as g:
        for transcript in shared_transcripts:
            passing = all([data[h_key][transcript] >= threshold for h_key in h_keys])
            if passing:
                g.write(transcript+'\n')

def main():
    parser = argparse.ArgumentParser(description='Calculates RT stop coverage from <.rtsc> file(s)')
    in_files = parser.add_argument_group('Input')
    in_files.add_argument('fasta',type=str,metavar='fasta',help='Reference Fasta')
    in_files.add_argument('f',type=str,help='Input <.rtsc> files', nargs='+')
    settings = parser.add_argument_group('Settings')
    settings.add_argument('-bases',type=str,default='AC',metavar='ACGT',help='[default = AC] Coverage Specificity')
    settings.add_argument('-ot',type=float,default=1.0,help='[default = 1.0] Overlap file threshold')
    out_files = parser.add_argument_group('Output')
    out_files.add_argument('-name',type=str,default = None, help='Output file name')
    out_files.add_argument('-ol',action='store_true',help='Create an overlap file')
    out_files.add_argument('-on',type=str,metavar='',default=None,help='Overlap file name')
    args = parser.parse_args()
    
    #Generate or assign name
    default_name = '_'.join(sorted([fyle.replace('.rtsc','') for fyle in args.f])+['coverage'])+'.csv'
    out_name = sfio.check_extension(args.name,'.csv') if args.name else default_name
    
    #Collect Data
    coverage_data = collect_coverages(args.f,args.fasta,args.bases)
    
    #Write Data
    write_coverage(coverage_data,out_name)
    
    #Create overlap file
    if args.ol:
        default_ol = '_'.join(sorted([fyle.replace('.rtsc','') for fyle in args.f])+['overlap',str(args.ot)])+'.txt'
        out_ol = default_ol if args.on == None else check_extension(args.on,'.txt')
        write_ol(coverage_data,out_ol,args.ot)

if __name__ == '__main__':
    main()
