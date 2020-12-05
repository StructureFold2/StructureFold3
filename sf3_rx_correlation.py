#!/usr/bin/env python3 

'''
Reformats either <.rtsc> or <.react> files for correlation analysis.
Using a coverage overlap file or any other list of specific transcripts
will truncate the output only to those transcripts. By default, 
this reports all bases for both file types, requiring an acompanying 
<.fasta> file. Note that the stops in <.rtsc> files have not been offset 
to reference the actual stop one base immediately five prime of the last 
sequenced base of the reads; when the specificity option is enabled, it 
is this base that is referenced, wheras <.react> files are precisely 1:1 
and require no offset.

The specificity option references the exact letters input by the user,
so if the <.fasta> has mixed or lower case letters, both need to be entered;
if the intention was to probe A and C, and there were lower case variants
in the fasta, the input would be ACac.
'''

#Imports
import argparse
import sf3libs.sf3io as sfio

#Functions
def write_react_repeatability(react_data,sequences,out_fyle,specificity='ACGT'):
   '''Writes out react correlation data.'''
   all_keys = sorted(set.union(*map(set,react_data.values())))
   header=','.join(['transcript','position','base']+all_keys)
   with open(out_fyle,'w') as g:
       g.write(header+'\n')
       for transcript, data in react_data.items():
           temp_data = [data.get(key,'NA'*len(sequences[transcript])) for key in all_keys]
           for pos, base in enumerate(sequences[transcript],1):
               if base in specificity:
                   info_cols = [transcript,str(pos),base]
                   data_cols = [str(values[pos-1]) for values in temp_data]
                   g.write(','.join(info_cols+data_cols)+'\n')

def write_rtsc_repeatability(rtsc_data,sequences,out_fyle,specificity='ACGT'):
   '''Writes out the rtsc correlation data. The offset (1) is built in'''
   all_keys = sorted(set.union(*map(set,react_data.values())))
   header=','.join(['transcript','position','base']+all_keys)
   with open(out_fyle,'w') as g:
       g.write(header+'\n')
       for transcript, data in react_data.items():
           temp_data = [data.get(key,'NA'*len(sequences[transcript]))[1:] for key in all_keys]
           for pos, base in enumerate(sequences[transcript][:-1],1):
               if base in specificity:
                   info_cols = [transcript,str(pos),base]
                   data_cols = [str(values[pos-1]) for values in temp_data]
                   g.write(','.join(info_cols+data_cols)+'\n')

#Workflow
def main():
    parser = argparse.ArgumentParser(description='Reformats <.rtsc>/<.react> for easy correlation analysis')
    parser.add_argument('fasta',help='Reference Fasta')
    parser.add_argument('rx',help='Input <.rx> files', nargs='+')
    parser.add_argument('-restrict',default=None, help='Filter to these transcripts via coverage file')
    parser.add_argument('-name',default=None, help='Specify output file name')
    parser.add_argument('-bases',default='AGCT', help='[ACGT] Nucleotide specifictiy')
    parser.add_argument('-verbose',action='store_true', help='Display metrics')
    args = parser.parse_args()

    #Set up
    covered = sfio.read_restrict(args.restrict) if args.restrict else None
    fasta_dict = sfio.read_fasta(args.fasta)
    #desc = [x.replace('.rtsc','').replace('.react','') for x in args.rx]
    desc = [sfio.rm_ext(x,'.rtsc','.react') for x in args.rx)
    
    
    #All files are <.rtsc>
    if all(fyle.split('.')[-1] == 'rtsc' for fyle in args.rx):
        rx_data = sfio.read_rx_files(args.rx,mode='rtsc',verbose=args.verbose)
        default_name = '_'.join(desc+['rtsc']+['correlation.csv'])
        out_name = sfio.check_extension(args.name,'.csv') if args.name else default_name
        if args.restrict:
            rx_data = {n:s for n,s in rx_data.items() if n in covered}
        if args.verbose:
            print('Remaining after filtering',len(rx_data),sep=',')
            print('Writing File:',out_name,sep=',')
        write_rtsc_repeatability(rx_data,fasta_dict,out_name,args.bases)

    #All files are <.react>
    elif all(fyle.split('.')[-1] == 'react' for fyle in args.rx):
        rx_data = sfio.read_rx_files(args.rx,mode='react',verbose=args.verbose)
        default_name = '_'.join(desc+['react']+['correlation.csv'])
        out_name = sfio.check_extension(args.name,'.csv') if args.name else default_name
        if args.restrict:
            rx_data = {n:s for n,s in rx_data.items() if n in covered}
        if args.verbose:
            print('Remaining after filtering',len(rx_data),sep=',')
            print('Writing File:',out_name,sep=',')
        write_react_repeatability(rx_data,fasta_dict,out_name,args.bases)

    #Files do not make sense
    else:
        print('All Files must be the same type, and either .rtsc or .react!')

if __name__ == '__main__':
    main()
