#!/usr/bin/env python3

'''
This script takes two <.rtsc> along with the base <.fasta> to calculate the reactivity
scores of each nucelotide of each transcript. If used without a <.scale> file, the script will generate 
a <.scale> to be used with this script when calculating any additinonal reactvities from other samples that
are to be compared with the original.
Transcripts which produce a 0 average on the normalization scale will not have their final reactivity calculated, 
and these transcripts may be logged to a file (optional). If given a <.scale> file, it will apply it as the normalization scale. 
Transcripts not included in the <.scale> will not be calculated and can be logged to a file (optional).
Normalization (2-8%) may be turned off for the reactivity calculation, as can taking the natural log of the reactivity values, 
via options.
'''

#Imports
import math
import argparse
import os
import sf3libs.sf3io as sfio

#Functions
def calculate_raw_reactivity(reagent_minus,reagent_plus,nlog_off=False):
    '''Calculates raw reactivity, with or without the natural log'''
    key_set,data_out = set(reagent_plus.keys()).intersection(set(reagent_minus.keys())),{}
    for key in key_set:
        plus_vector = [math.log(value+1,math.e) for value in reagent_plus[key]] if nlog_off == False else reagent_plus[key]
        minus_vector = [math.log(value+1,math.e) for value in reagent_minus[key]] if nlog_off == False else reagent_minus[key]
        sum_plus,sum_minus,length = sum(plus_vector),sum(minus_vector),len(plus_vector)
        if sum_plus != 0 and sum_minus != 0:
            nrm_plus_vector = [float(y)/float(sum_plus)*length for y in plus_vector]
            nrm_minus_vector = [float(z)/float(sum_minus)*length for z in minus_vector]
            plus_minus = [max(0,a-b) for a, b in zip(nrm_plus_vector,nrm_minus_vector)]
            data_out[key] = plus_minus
    return data_out

def generate_normalization_scale(derived_reactivities,transcript_seqs,specificity):
    '''Generates the 2-8% scale to normalize against'''
    data = {}
    for transcript, reactivities in derived_reactivities.items():
        sequence = transcript_seqs[transcript]
        accepted = sorted([reactivities[k] for k in range(1,len(reactivities)) if sequence[k-1] in specificity],reverse=True)
        top = accepted[int(len(accepted)*0.02):int(len(accepted)*0.1)]
        top_average = sum(top)/len(top) if len(top) > 0 else 0
        if top_average > 0:
            data[transcript] = top_average
    return data

def write_norm_scale(scale_dictionary,outfile):
    '''Writes out a normalization scale file'''
    with open(outfile,'w') as g:
        g.write(','.join(['transcript','value'])+'\n')
        for transcript, value in sorted(scale_dictionary.items()):
            g.write(','.join([transcript,str(value)])+'\n')

def read_norm_scale(normalization_file):
    '''Reads in a normalization scale file'''
    info = {}
    with open(normalization_file, 'r') as f:
        for line in f:
            if line.startswith('transcript'):
                continue
            else:
                transcript,value = line.strip().split(',')
                info[transcript] = float(value)
    return info

def calculate_final_reactivity(derived_reactivities,sequences,specificity,threshold,nrm_scale,norm_off=False):
    '''Calculates the final reactivity'''
    data_out,missing_transcripts = {},{}
    for transcript, reactivities in derived_reactivities.items():
        if transcript in nrm_scale:
            normalizer = nrm_scale[transcript] if norm_off == False else 1
            seq = sequences[transcript]
            normalized_values =[str(float('%.3f'%min((reactivities[x]/normalizer), threshold))) if seq[x-1] in specificity else 'NA' for x in range(1,len(reactivities))]+['NA']
            data_out[transcript] = normalized_values
        else:
            missing_transcripts[transcript] = None
    return data_out,missing_transcripts

def main():
    parser = argparse.ArgumentParser(description='Generates a <.react> file from two <.rtsc> files')
    in_files = parser.add_argument_group('Input')
    in_files.add_argument('control',type=str,help='Control <.rtsc> file')
    in_files.add_argument('treatment',type=str,help='Reagent <.rtsc> file')
    in_files.add_argument('fasta',type=str,help='Transcript <.fasta> file')
    settings = parser.add_argument_group('Settings')
    settings.add_argument('-restrict',default = None,metavar='<.txt>',help='Limit analysis to these specific transcripts')
    settings.add_argument('-scale',type=str,default= None,metavar='<.scale>',help='Provide a normalizaiton file for calculation')
    settings.add_argument('-threshold',type=float,default=7.0,help='[default = 7.0] Maximum Reactivity Cap')
    settings.add_argument('-ln_off',action='store_true', help='Do not take the natural log of the stop counts')
    settings.add_argument('-nrm_off',action='store_true',help='Do not perform final 2-8'+u"\uFF05"+' reactivity normalization')
    settings.add_argument('-bases',type=str,default = 'AC',metavar='ACGT',help='[default = AC] Reaction Specificity')
    out_files = parser.add_argument_group('Output')
    out_files.add_argument('-save_fails',action='store_true',help='Log transcripts with zero or missing scales')
    out_files.add_argument('-name',type=str,default = None, help='Specify output file name')
    args = parser.parse_args()
    
    #Create output name
    base_tag = [x.split(os.sep)[-1].replace('.rtsc','') for x in [args.control,args.treatment]]
    log_tag = ['ln'] if args.ln_off == False else []
    nrm_tag = ['nrm'] if args.ln_off == False else []
    base_name = '_'.join(base_tag+log_tag+nrm_tag)
    out_name = base_name+'.react' if args.name == None else sfio.check_extension(args.name,'.react')

    #Read in data
    control_data,treatment_data = map(sfio.read_rtsc,[args.control,args.treatment])

    #Apply filter if enabled
    if args.restrict:
        covered = sfio.read_restrict(args.restrict)
        control_data = {n:s for n,s in control_data.items() if n in covered}
        treatment_data = {n:s for n,s in treatment_data.items() if n in covered}

    #Calculate Derived Reactivity
    data = calculate_raw_reactivity(control_data,treatment_data,args.ln_off)
    
    #Read in transcript sequences
    seqs = sfio.read_fasta(args.fasta)
    
    #Generate and write scale, or read a <.scale> file in
    normalizaiton_scale = generate_normalization_scale(data,seqs,args.bases) if args.scale == None else read_norm_scale(args.scale)
    if args.scale == None:
        write_norm_scale(normalizaiton_scale,out_name.replace('.react','.scale'))
    
    #Calculate Final Reactivity
    out_reactivity,out_missing = calculate_final_reactivity(data,seqs,args.bases,args.threshold,normalizaiton_scale,args.nrm_off)
    
    #Write Out
    sfio.write_react(out_reactivity,out_name,sort_flag=True)

    #Write Out Fails
    if args.save_fails:
        sfio.write_keys(out_missing,out_name.replace('.react','_unresolvable_transcripts.txt'))

if __name__ == '__main__':
    main()
