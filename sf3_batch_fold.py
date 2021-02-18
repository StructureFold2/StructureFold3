#!/usr/bin/env python3

#Imports
import os
import argparse
import subprocess
import sf3libs.sf3io as sfio

#Functions
def gen_length_errors(error_dict,seqs,s=10,l=5000):
    '''Populates Length Errors: too short or long'''
    length_errors = {k:len(v) for k,v in seqs.items() if len(v)<s or len(v)>l}
    for k, v in length_errors.items():
        error_dict.setdefault(k,{})['length_error'] = v

def gen_xstraint_errors(error_dict,seqs,cons):
    '''Populates xstraint Errors: missing sequence or length desynch'''
    if cons:
        missing = {k:'missing' for k in cons.keys() if k not in seqs}
        s_keys = seqs.keys() & cons.keys()
        desynch = {k:'desynch' for k in s_keys if len(cons[k]) != len(seqs[k])}
        merged = {**missing,**desynch}
        for k, v in merged.items():
            error_dict.setdefault(k,{})['xstraint_error'] = v
    else:
        pass

def write_error_log(error_dict,out_name):
    '''Writes the Error Log'''
    try:
        fields = sorted(list(set.union(*map(set, error_dict.values()))))
        header = ['transcript']+fields
        with open(out_name,'w') as g:
            g.write(','.join(header)+'\n')
            for name, sub in error_dict.items():
                values = [name]+[str(sub.get(f,'NA')) for f in fields]
                g.write(','.join(values)+'\n')
    except TypeError:
        pass

def write_params_log(params,out_name):
    '''Writes parameter log'''
    with open(out_name,'w') as g:
        g.write(','.join(['Parameter','Value'])+'\n')
        for p,v in sorted(vars(params).items()):
            g.write(','.join([p,str(v)])+'\n')

def gen_temp_sequence(name,seq,path):
    '''Generates a temporary fasta file'''
    out = os.path.join(path,'_'.join([name,'temp','sequence.fa']))
    with open(out,'w') as g:
        g.write('>'+name+'\n'+seq+'\n')
    return out

def gen_rest(name,react_entry,path,trim=0):
    '''Generates a temporary restraint file'''
    entry = react_entry[:-trim] if trim else react_entry
    out = os.path.join(path,'_'.join([name,'temp','restraint.txt']))
    with open(os.path.join(out),'w') as g:
        for pos,value in enumerate(entry,1):
            if value != 'NA':
                g.write('\t'.join([str(pos),str(value)])+'\n')
    return out

def gen_cons(name,react_entry,path,trim=0,threshold=7):#WARNING Function not yet done
    '''Generates a temporary constraint file'''
    entry = react_entry[:-trim] if trim else react_entry
    out = os.path.join(path,'_'.join([name,'temp','constraint.txt']))
    cons = [(pos, value) for pos, value in enumerate(entry,1)]
    constraints = filter(lambda x: x[1] >= threshold,cons)
    with open(out,'w') as g:
        g.write('DS:\n-1\nSS:\n')
        for constraint in constraints:
            g.write(str(contstraint[0])+'\n')
        g.write('-1\nMod:\n-1\nPairs:-1 -1\nFMN:\n-1\nForbids:\n-1 -1\n')
    return out

def gen_paths(base,*sub_paths):
    '''Generate paths'''
    empty = []
    for sub in sub_paths:
        sub_dir = os.path.join(base,sub)
        os.mkdir(sub_dir)
        empty.append(sub_dir)
    return empty

def gen_file_paths(file_name,file_paths,file_suffixes):
    '''Generates file paths'''
    empty = []
    for f_path,f_sfx in zip(file_paths,file_suffixes):
        out = os.path.join(f_path,'_'.join([file_name]+f_sfx))
        empty.append(out)
    return empty

def remove_temp(*fyles):
    '''Remove files'''
    for fyle in fyles:
        subprocess.run(['rm',fyle]) 

def file_len(fname):
    '''Gives file length'''
    with open(fname,'r') as f:
        for i, l in enumerate(f):
            pass
    return i + 1

#Workflow
def main():
    parser = argparse.ArgumentParser(description='Batch runs folding programs')
    in_files = parser.add_argument_group('Input')
    in_files.add_argument('fasta',default=None,help='Reference Fasta File')
    settings = parser.add_argument_group('Settings')
    settings.add_argument('-mode',type=str.upper,default='R',choices=['R','V'],help='RNAStructure or Vienna')
    settings.add_argument('-react',default=None,metavar='<.react>',help='React file to use as restraints/constraints')
    settings.add_argument('-restrict',default=None,metavar='<.txt>',help='Limit folding to these specific transcripts')
    settings.add_argument('-temp',type=float,default=310.15,help='[default = 310.15] Temperature to use for folding')
    settings.add_argument('-cores',default='4',type=str,help="[default = 4] Number of cores to use")
    settings.add_argument('-distance',default=99999,type=int,help='[default = 99999] Maximum pairing distance')
    settings.add_argument('-minlen',default=10,type=int,help='[default = 10] Minimum length to fold')
    settings.add_argument('-maxlen',default=5000,type=int,help='[default = 5000] Maximum length to fold')
    settings.add_argument('-truncate',default=0,type=int,help='[default = 0] Ignore <.react> for last N nucleotides')
    settings.add_argument('-threshold',default=None,type=float,help='Apply hard constraints using this threshold',dest='th')
    rna = parser.add_argument_group('RNAStructure Settings')
    rna.add_argument('-slope',default=1.8,type=float,help='[default = 1.8] Parameter for RNAstructure')
    rna.add_argument('-intercept',default=-0.6,type=float,help='[default = -0.6] Parameter for RNAstructure')
    rna.add_argument('-partition',action='store_true',help='Use the partition function')
    rna.add_argument('-multiple',action='store_true',help='Output all structures rather than just MFE')
    out_files = parser.add_argument_group('Output')
    out_files.add_argument('-errorname',default='errors.csv',help='Name the error log')
    out_files.add_argument('-paraname',default='parameters.csv',help='Name the parameter log')
    out_files.add_argument('-outdir',default='folded',help='Name the out directory')
    out_files.add_argument('-ct',default='CT',help='Name the ct folder')
    out_files.add_argument('-ps',default='PS',help='Name the ps folder')
    out_files.add_argument('-bp',default='PS',help='Name the bp folder (Partition Only)')
    out_files.add_argument('-pfs',default='PFS',help='Name the pfs folder (Partition Only)')
    args = parser.parse_args()

    #Prepare out directory
    if os.path.isdir(args.outdir):
        print('Out directory already exists, quitting...')
        quit()
    else:
        os.mkdir(args.outdir)

    #Input
    sequences = sfio.read_fasta(args.fasta)
    xstraints = sfio.read_react(args.react) if args.react else None
    restrict = sfio.read_restrict(args.restrict) if args.restrict else None

    #Populate and filter errors
    errors = {}
    gen_length_errors(errors,sequences,args.minlen,args.maxlen)
    gen_xstraint_errors(errors,sequences,xstraints)
    if errors:
        sl = len(sequences)
        sequences = {z:y for z,y in sequences.items() if z not in errors}
        el = len(sequences)
        print('After screening errors, {} of {} folds remain'.format(el,sl))
    if restrict:
        sl = len(sequences)
        sequences = {n:s for n,s in sequences.items() if n in restrict}
        el = len(sequences)
        print('After restricting analysis, {} of {} folds remain'.format(el,sl))

    #Housekeeping in base directory
    write_params_log(args,os.path.join(args.outdir,sfio.check_extension(args.paraname,'.csv')))
    write_error_log(errors,os.path.join(args.outdir,sfio.check_extension(args.errorname,'.csv')))

    #RNAStructure Suite
    if args.mode == 'R':
        rna_env = {**os.environ,'OMP_NUM_THREADS':args.cores}
        options = {'-t':args.temp,'-md':args.distance}

        #partition-smp Suite
        if args.partition:
            paths = gen_paths(args.outdir,args.pfs,args.bp,args.ps,args.ct)
            extensions = ['.pfs','.bp','.ps','.ct']

            #Restrained/Constrained
            if xstraints:
                opts = {} if args.th else {'-si':args.intercept,'-sm':args.slope}
                options.update(opts)
                out_type = 'constraint' if args.th else 'restraint'
                suffixes = [[str(args.temp),out_type+X] for X in extensions]
                params = sfio.flatten_list([[k,str(v)] for k, v in options.items()])
                #Iterate through sequences
                for name,seqeunce in sequences.items():
                    temp_seq = gen_temp_sequence(name,seq,args.outdir)
                    gargs = [name,xstraints[name],args.outdir,args.truncate]
                    guide = gen_cons(*gargs,args.th) if args.th else gen_rest(*gargs)
                    #Out Files
                    pf_out,bp_out,ps_out,ct_out = gen_file_paths(name,paths,suffixes)
                    #Run Commands
                    gflag = '-c' if args.th else '-sh'
                    pf_command = ['partition-smp',temp_seq,pf_out,gflag,guide]+params
                    bp_command = ['ProbabilityPlot',pf_out,bp_out,'-t']
                    ps_command = ['ProbabilityPlot',pf_out,ps_out]
                    subprocess.run(pf_command,stdout=subprocess.DEVNULL,env=rna_env)
                    subprocess.run(bp_command,env=rna_env,stdout=subprocess.DEVNULL)
                    subprocess.run(ps_command,env=rna_env,stdout=subprocess.DEVNULL)
                    #Generate CT file
                    if file_len(bp_out) > 2:
                        ct_command = ['MaxExpect',pf_out,ct_out]
                        subprocess.run(ct_command,stdout=subprocess.DEVNULL)
                    remove_temp(temp_seq,guide)

            #in-silico
            else:
                suffixes = [[str(args.temp),'silico'+X] for X in extensions]
                params = sfio.flatten_list([[k,str(v)] for k, v in options.items()])
                #Iterate through sequences
                for name,seq in sequences.items():
                    temp_seq = gen_temp_sequence(name,seq,args.outdir)
                    #Out Files
                    pf_out,bp_out,ps_out,ct_out = gen_file_paths(name,paths,suffixes)
                    #Run Commands
                    pf_command = ['partition-smp',temp_seq,out_pfs]+params
                    bp_command = ['ProbabilityPlot',pf_out,bp_out,'-t']
                    ps_command = ['ProbabilityPlot',pf_out,ps_out]
                    subprocess.run(pf_command,stdout=subprocess.DEVNULL,env=rna_env)
                    subprocess.run(bp_command,stdout=subprocess.DEVNULL)
                    subprocess.run(ps_command,stdout=subprocess.DEVNULL)
                    #Generate CT file
                    if file_len(out_bp) > 2:
                        ct_command = ['MaxExpect',pf_out,out_ct]
                        subprocess.run(ct_command,stdout=subprocess.DEVNULL)
                    remove_temp(temp_seq)

        #Fold-smp Suite
        else:
            paths = gen_paths(args.outdir,args.ct,args.ps)
            extensions = ['.ct','.ps']
            flags = {'-mfe': not args.multiple}
            applied_flags = [k for k,v in flags.items() if v ]

            #Restrained/Constrained
            if xstraints:
                out_type = 'constraint' if args.th else 'restraint'
                suffixes= [[str(args.temp),out_type+X] for X in extensions]
                opts = {} if args.th else {'-si':args.intercept,'-sm':args.slope}
                options.update(opts)
                params = sfio.flatten_list([[k,str(v)] for k, v in options.items()])
                params.extend(applied_flags)
                #Iterate through sequences
                for name,seq in sequences.items():
                    temp_seq = gen_temp_sequence(name,seq,args.outdir)
                    gargs = [name,xstraints[name],args.outdir,args.truncate]
                    guide = gen_cons(*gargs,args.th) if args.th else gen_rest(*gargs)
                    gflag = '-c' if args.th else '-sh'
                    #Out Files
                    ct_out,ps_out = gen_file_paths(name,paths,suffixes)
                    #Run Commands
                    fold_command = ['Fold-smp',temp_seq,ct_out,gflag,guide]+params
                    draw_command = ['draw',ct_out,ps_out]
                    subprocess.run(fold_command,env=rna_env,stdout=subprocess.DEVNULL)
                    subprocess.run(draw_command,stdout=subprocess.DEVNULL)
                    #Remove Temporary Files
                    remove_temp(temp_seq,guide)

            #in-silico
            else:
                out_type = 'silico'
                suffixes = [[str(args.temp),out_type+X] for X in extensions]
                params = sfio.flatten_list([[k,str(v)] for k, v in options.items()])
                params.extend(applied_flags)
                #Iterate through sequences
                for name,seq in sequences.items():
                    #Generate Temporary Files
                    temp_seq = gen_temp_sequence(name,seq,args.outdir)
                    #Out Files
                    ct_out,ps_out = gen_file_paths(name,paths,suffixes)
                    #Run Commands
                    fold_command = ['Fold-smp',temp_seq,ct_out]+params
                    draw_command = ['draw',ct_out,ps_out]
                    subprocess.run(fold_command,env=rna_env,stdout=subprocess.DEVNULL)
                    subprocess.run(draw_command,stdout=subprocess.DEVNULL)
                    #Remove Temporary Files
                    remove_temp(temp_seq)

    #Vienna Package Suite
    if args.mode == 'V':
        print('Vienna Package not yet supported')

if __name__ == '__main__':
    main()
