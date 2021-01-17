#Imports
import re

#Classes
class MappedRead(object):
    '''Mapped Read in a SAM file'''
    def __init__(self,sam_line):
        #Get Basic Attributes
        attributes = ('qname','bitflag','rname','position','mapq',
            'cigar','rnext','pnext','tlen','sequence','quality')
        fields = sam_line.strip().split()
        for attribute,field in zip(attributes,fields):
            setattr(self,attribute,field)
        optionals = {o.split(':')[0]:o.split(':')[2] for o in fields[11:]}
        #Get Derived Attributes
        self.mismatches = int(optionals.get('NM',0))
        #Attributes based on MD
        MD = optionals.get('MD','Null')
        self.first_match = not MD.startswith('0')
        self.snps = re.findall("\d+[A-Z]|\d+\^[A-Z]+|\d+",MD)
        #Attributes based on bitflag
        for attr, val in bitflag_to_attributes(self.bitflag).items():
            setattr(self,attr,val)

#Functions
def flag_to_bin(flag):
    '''Converts a SAM bitflag into binary'''
    return bin(int(flag)).replace('0b','').zfill(12)

def bitflag_to_attributes(flag):
    '''Creates a dictionary of attributes from a SAM bitflag'''
    bits = [int(x) for x in list(flag_to_bin(flag))]
    bits.reverse()
    types = ['read_paired','proper_pair','r1_unmap','r2_unmap',
             'r1_reverse','r2_reverse','1_in_pair','2_in_pair',
             'secondary_alignment','quality_fail','PCR_duplicate',
             'supplementary']
    return {k: True if v else False for k, v in zip(types,bits)}



