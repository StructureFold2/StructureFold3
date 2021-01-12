#!/usr/bin/env python3

#Imports
import argparse
import sf3libs.sf3io as sfio

#Functions
def merge_rtsc(rtsc_lyst):
    '''Takes a list of <.rtsc> files, returns a dictionary that is sum of the RT stops'''
    all_stops= {}
    for rtsc in sorted(rtsc_lyst):
        data = sfio.read_in_rtsc(rtsc)
        for transcript, values in data.items():
            if transcript in all_stops:
                current_values = all_stops[transcript]
                all_stops[transcript] = [sum(spot) for spot in zip(current_values,values)]
            else:
                all_stops[transcript] = values
    return all_stops

def main():
    parser = argparse.ArgumentParser(description='Combines <.rtsc> files, typically replicates of the same sample')
    in_files = parser.add_argument_group('Input')
    in_files.add_argument('rtsc',help='Input <.rtsc> files',nargs='+')
    out_files = parser.add_argument_group('Output')
    out_files.add_argument('-sort',action='store_true',default=False,help='Sort output by transcript name')
    out_files.add_argument('-name',default=None,help='Specify output file name')
    args = parser.parse_args()

    #Generate name or assign the user provided name
    default_name = '_'.join(sorted([x.replace('.rtsc','') for x in args.rtsc]))+'.rtsc'
    out_name = default_name if args.name == None else sfio.check_extension(args.name,'.rtsc')

    #Pool all <.rtsc> into a dictionary
    all_stops = merge_rtsc(args.rtsc)

    #Write out the dictionary
    sfio.write_rtsc(all_stops,args.sort,out_name)


if __name__ == '__main__': 
    main()
