#!/usr/bin/env python3

'''
Converts <.rtsc> to <.csv> of relative transcript abundances. 
Traditional RNA-seq may be more precise at this task than Structure-Seq libraries.
'''

#Imports
import argparse
import sf3libs.sf3io as sfio

#Functions
def read_in_total_stops(afile):
    '''Reads in <.rtsc>, returns the total stops and length for each transcript'''
    base_data = sfio.read_rtsc(afile)
    new = {k:(sum(v),len(v)) for k,v in base_data.items()}
    return new

def values_to_TPM(data):
    '''Calculates TPM or Transcripts Per Kilobase Million reads'''
    reads_per_kb = {k:v[0]/(float(v[1])/1000) for k, v in data.items()}
    normalize = sum(reads_per_kb.values())/1000000
    TPM_values = {k:v/normalize for k,v in reads_per_kb.items()}
    return TPM_values

def values_to_RPKM(data):
    '''Calculates RPKM or Reads Per Kilobase per Million reads'''
    norm = sum(v[0] for v in data.values())
    RPKM_values = {k:(float(v[0])*1000*1000000)/(v[1]*norm) for k, v in data.items()}
    return RPKM_values

def populate_dictionary(fylelyst,funct):
    '''Generates a nested dictionary of abundances'''
    data,new = {fyle.replace('.rtsc',''):funct(fyle) for fyle in fylelyst},{}
    for f_name, sub_dict in rx_data.items():
        for transcript, value in sub_dict.items():
            new.setdefault(transcript,{})[f_name] = value
    return new

def write_data(adict,outfyle,data_unit):
    '''Writes the data to a <.csv>'''
    header = ','.join(['transcript',outfyle.replace('.csv','')])
    with open(outfyle,'w') as g:
        g.write(header+'\n')
        for transcript,abundance_stat in adict.items():
            g.write(','.join([transcript,str(abundance_stat)])+'\n')

def write_data_batch(data,outname,suffix,missing):
    '''writes the data to a <.csv>'''
    f_keys = sorted(list(set.union(*map(set, data.values()))))
    transcript_keys = sorted(data.keys())
    header = ['transcript']+['_'.join([fyle,suffix]) for fyle in f_keys]
    with open(outname,'w') as g:
        g.write(','.join(header)+'\n')
        for transcript in transcript_keys:
            data = [data[transcript].get(fyle,missing) for fyle in f_keys]
            g.write(','.join([transcript]+data)+'\n')

def main():
    parser = argparse.ArgumentParser(description='Determines approximate transcript abundance using <.rtsc> files')
    in_files = parser.add_argument_group('Input')
    in_files.add_argument('rtsc',default=None,nargs='+',help='Input <.rtsc> file(s)')
    settings = parser.add_argument_group('Settings')
    settings.add_argument('-mode',type=str.upper,choices=['RPKM','TPM'],help='Abundance Metric to Calculate')
    settings.add_argument('-zero',action='store_true',help='Set missing values to zero')
    out_files = parser.add_argument_group('Output')
    out_files.add_argument('-name',default=None,help='Specify output file name')
    args = parser.parse_args()

    #Files to operate on, dictionary of functions
    abundance_methods = {'RPKM':values_to_RPKM,'TPM':values_to_TPM}

    #Nomenclature
    name = sorted([x.replace('.rtsc','') for x in args.rtsc])+[args.mode]
    default_name = '_'.join(name_1)+'.csv'
    out_name = default_name if not args.name else sfio.check_extension(args.name,'.csv')

    #Build data set
    abundances = populate_dictionary(args.rtsc,abundance_methods[args.mode])

    #Write out
    blank = 0.0 if zero else 'NA'
    write_data_batch(abundances,out_name,args.mode,blank)

if __name__ == '__main__':
    main()
