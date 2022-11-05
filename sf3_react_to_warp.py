#!/usr/bin/python3

#Imports
import argparse
import sf3libs.sf3io as sfio

#Functions
def write_warp(seqz,rctz,out='file.warp'):
    '''Writes out a warp file'''
    key_set = seqz.keys() & rctz.keys()
    subs = {'0.0':'NaN','Na':'NaN'}
    with open(out,'w') as g:
        for key in key_set:
            new = [str(x) for x in rctz[key]]
            CLEAN = ','.join([subs[q] if q in subs else q for q in new])
            g.write(key+'\n')
            g.write(seqz[key]+'\n')
            g.write(CLEAN+'\n')

def main():
    parser = argparse.ArgumentParser(description='Convert .react to ShapeWarp')
    in_files = parser.add_argument_group('Input')
    in_files.add_argument('fasta',metavar='<.fasta>',help='Input <.fasta> file')
    in_files.add_argument('react',metavar='<.react>',help='Input <.react> file')
    settings = parser.add_argument_group('Settings')
    settings.add_argument('-restrict',default=None,metavar='<.txt>'
                          ,help='Limit analysis to these specific transcripts')
    out_files = parser.add_argument_group('Output')
    out_files.add_argument('-name',type=str,default=None,help='Output ShapeWarp file')
    args = parser.parse_args()

    #File Input
    seqs,rct = sfio.read_fasta(args.fasta),sfio.read_react(args.react)

    #Filter by coverage
    restrict = sfio.read_restrict(args.restrict) if args.restrict else None
    if restrict:
            seqs = {n:s for n,s in seqs.items() if n in restrict}
            rcts = {a:b for a,b in rct.items() if a in restrict}

    #Nomenclature
    default_name = args.react.replace('.react','.warp')
    out_name = sfio.check_extension(args.name,'.warp') if args.name else default_name

    #Write Out
    write_warp(seqs,rct,out=out_name)

if __name__ == '__main__':
    main()
