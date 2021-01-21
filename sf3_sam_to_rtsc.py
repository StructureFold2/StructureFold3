#!/usr/bin/env python3

#Imports
import argparse
import sf3libs.sf3io as sfio
import sf3libs.sf3sam as sf3sam
from collections import Counter

#Functions
def sam_to_rtsc(in_sam,seq_lims,mm,fst,keyflag):
    '''Takes a SAM file and writes it as an rtsc file, returns log'''
    #Read in SAM, gather passing stops and metrics
    report,stops,rtsc = Counter(),{},{}
    with open(in_sam,'r') as f:
        for line in f:
            if line.startswith('@'):
                continue
            else:
                item = sf3sam.MappedRead(line)
                if item.pass_bitflag(keyflag):
                    tests = {'mismatches':item.pass_mismatches(mm),
                             'firstbase':item.pass_first(fst)}
                    if all(tests.values()):
                        report['passing']+=1
                        stops.setdefault(item.rname,Counter())[int(item.position)-1]+=1
                    else:
                        failed_tests = sorted([k for k,v in tests.items() if not v])
                        failed_key = '_and_'.join(failed_tests)
                        report[failed_key]+=1
                else:
                    report['bitflag']+=1
    #Reorganize the stops
    for seq,limit in seq_lims.items():
        if seq in stops:
            rtsc[seq] = [str(stops[seq][index]) for index in range(0,limit)]
        else:
            rtsc[seq] = ['0']*limit
    #Write the stops
    out_name = in_sam.replace('.sam','.rtsc')
    sfio.write_rtsc(rtsc,out_name)
    #Return file specific metrics
    return report

def write_sam_filter_report(data,out_name):
    '''Writes a report on all the filtering metrics'''
    keyring = sorted(set.union(*map(set,data.values())))
    header = ','.join(['sam_file'] + keyring)
    with open(out_name,'w') as g:
        g.write(header+'\n')
        for fyle, sub in sorted(data.items()):
            line = ','.join([fyle]+[str(sub.get(key,0)) for key in keyring])
            g.write(line+'\n')

#Workflow
def main():
    parser = argparse.ArgumentParser(description='Converts <.sam> into reverse transcriptase stop files <.rtsc>')
    in_files = parser.add_argument_group('Input')
    in_files.add_argument('fasta',default=None,help='Index Fasta File')
    in_files.add_argument('sam',default=None,help='Input SAM file(s)',nargs='+')
    settings = parser.add_argument_group('Settings')
    settings.add_argument('-mismatches',type=int,default=3,metavar='<number>',help='[default = 3] Maximum allowed mismatches/indels')
    settings.add_argument('-firstmm',action='store_true',default=False,help='Accept alignments with first base mismatches')
    settings.add_argument('-reverse',action='store_true',default=False,help='Accept alignments to the reverse strand')
    settings.add_argument('-rm_secondary',action='store_false',default=True,help='Remove secondary alignments',dest='secondary')
    out_files = parser.add_argument_group('Output')
    out_files.add_argument('-logname',type=str,default='filter_log.csv',help='[default = filter_log.csv] Log file name')
    parser.set_defaults(r1_unmap=False)
    args = parser.parse_args()

    #Generate Seq Limits
    reference = sfio.read_fasta(args.fasta)
    limits = {name:len(seq) for name, seq in reference.items()}

    #Generate SAM filterflag
    keys = {'r1_unmap':args.r1_unmap,'r1_reverse':args.reverse,'secondary':args.secondary}
    passing_keys = [k for k,v in keys.items() if not v]
    keyflag = sum([sf3sam.flag_values[R] for R in passing_keys])

    #Iterate through files
    log_data = {}
    for fyle in args.sam:
        entry = sam_to_rtsc(fyle,limits,args.mismatches,args.firstmm,keyflag)
        log_data[fyle] = entry

    #Write out log
    write_sam_filter_report(log_data,sfio.check_extension(args.logname,'.csv'))

if __name__ == '__main__':
    main()
