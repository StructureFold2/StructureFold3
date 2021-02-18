#!/usr/bin/python3

#Imports
import argparse
import sf3libs.sf3io as sfio
from sf3libs.sf3reactutils import ReactStats

#Functions
def write_out_stats(data,outfyle,trim,minlen):
    '''Write out the data in <.csv> format'''
    f_keys = sorted(list(set.union(*map(set, data.values()))))
    transcript_keys = sorted(data.keys())
    types = ['_max','_average','_std','_gini']
    header = ['transcript']+[fyle+mod for fyle in f_keys for mod in types]
    with open(outfyle,'w') as g:
        g.write(','.join(header)+'\n')
        for transcript in transcript_keys:
            sub,new = data[transcript],[]
            for fyle in f_keys:
                if fyle in sub:
                    new.extend(ReactStats(sub[fyle],trim,minlen).as_list())
                else:
                    new.extend(ReactStats().as_list())
            g.write(','.join(new)+'\n')

def main():
    parser = argparse.ArgumentParser(description='Generates a simple statistical summary for <.react> files.')
    in_files = parser.add_argument_group('Input')
    in_files.add_argument('react',default=None,help='Input <.react> files',nargs='+')
    settings = parser.add_argument_group('Settings')
    settings.add_argument('-restrict',default=None,metavar='<.txt>',help='Limit analysis to these specific transcripts')
    settings.add_argument('-trim',type=int,default=20,metavar='<number>',help='[default = 20] ignore n last bp of reactivity')
    settings.add_argument('-minlen',type=int,default=10,metavar='<number>',help='[default = 10] minimum effective length of transcript')
    out_files = parser.add_argument_group('Output')
    out_files.add_argument('-name',default=None,help='Specify output file name')
    args = parser.parse_args()

    #File Input
    rx_data = sfio.read_rx_files(sorted(args.react),'.react',verbose=False)
    
    #Filter by coverage
    restrict = sfio.read_restrict(args.restrict) if args.restrict else None
    if restrict:
            rx_data = {n:s for n,s in rx_data.items() if n in restrict}

    #Nomenclature
    name_1 = sorted([x.replace('.react','') for x in in_files])
    name_2 = [str(qq)+q for qq,q in zip([args.trim,args.minlen],['trim','minlen'])]
    default_name = '_'.join(name_1+name_2+['statistics'])+'.csv'
    out_name = sfio.check_extension(args.name,'.csv') if args.name else default_name

    #Write Out File
    write_out_stats(data,out_name,args.trim,args.minlen)

if __name__ == '__main__':
    main()
