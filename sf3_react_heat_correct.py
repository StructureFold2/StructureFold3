#!/usr/bin/env python3

'''
Must be used on <.react> files that have already been filtered to only 
include the transcripts above the accepted coverage threshold, and must be the exact same transcripts!
'''

#Imports
import sys
import argparse
import sf3libs.sf3io as sfio

#Functions
def sum_react(react_dict):
    '''Sum of all the reactivities'''
    return sum([sum(filter(lambda x: isinstance(x, float), v)) for v in react_dict.values()])

def apply_correction(react_dict,correction):
    '''Applies a correction to all values'''
    return {k:[x*correction if x!= 'NA' else 'NA' for x in v] for k, v in react_dict.items()}

#Main Function
def main():
    parser = argparse.ArgumentParser(description='Corrects two <.react>s for differential temperature')
    in_files = parser.add_argument_group('Input')
    in_files.add_argument('lower',type=str,help='lower temp <.react> file')
    in_files.add_argument('higher',type=str,help='higher temp <.react> file')
    out_files = parser.add_argument_group('Output')
    out_files.add_argument('-suffix',type=str,default='corrected',help='[default = corrected] Suffix for out files')
    args = parser.parse_args()

    #Sum all reactivities
    cold_react,hot_react = map(sfio.read_react,[args.lower,args.higher])
    cold_sum,hot_sum = map(sum_react,[cold_react,hot_react])

    #Check files
    if not cold_react.keys() == hot_react.keys():
        print('Warning! Non-parallel transcript sets between reacts! Quitting...')
        sys.exit() 

    #Calculate Corrections
    heat_correction = (hot_sum+cold_sum)/(2*hot_sum)
    cold_correction = (hot_sum+cold_sum)/(2*cold_sum)
    
    #Show Corrections
    print('Higher Temp values to be downscaled by factor: {}'.format(heat_correction))
    print('Lower Temp values to be upscaled by factor: {}'.format(cold_correction))
  
    #Apply corrections
    new_hot = apply_correction(hot_react,heat_correction)
    new_cold = apply_correction(cold_react,cold_correction)
  
    #Write Out
    cold_name = args.lower.replace('.react','_'+args.suffix+'.react')
    hot_name = args.higher.replace('.react','_'+args.suffix+'.react')
    sfio.write_react(new_cold,l_name)
    sfio.write_react(new_hot,h_name)

if __name__ == '__main__': 
    main()
