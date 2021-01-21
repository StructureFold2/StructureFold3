#Imports
import re

#Classes
class MappedRead(object):
    '''Mapped Read in a SAM file'''
    def __init__(self,sam_line):
        fields = sam_line.strip().split()
        #Assign Base Attributes
        attributes = ('qname','bitflag','rname','position','mapq',
            'cigar','rnext','pnext','tlen','sequence','quality')
        for attribute,field in zip(attributes,fields):
            setattr(self,attribute,field)
        #Assign derived and optional attributes
        optionals = {o.split(':')[0]:o.split(':')[2] for o in fields[11:]}
        #Get Derived Attributes
        self.mismatches = int(optionals.get('NM',0))
        MD = optionals.get('MD','Null')
        self.first = MD.startswith('0')
        self.snps = re.findall("\d+[A-Z]|\d+\^[A-Z]+|\d+",MD)

    def pass_bitflag(self,bitflag):
        '''Checks for unaccepted bits in a bitflag'''
        A,B = map(flag_to_bin,[self.bitflag,bitflag])
        bad_flag = any([a+b=='11' for a,b in zip(A,B)])
        result = False if bad_flag else True
        return result

    def pass_mismatches(self,number):
        '''Checks if this read exceeds number of mismatches'''
        return number >= self.mismatches
    
    def pass_first(self,accept):
        '''Given False/True, checks status of first base'''
        if not accept and not self.first:
            return True
        elif not accept and self.first:
            return False
        else:
            return True

#Functions
def flag_to_bin(flag):
    '''Converts a SAM bitflag into binary'''
    return bin(int(flag)).replace('0b','').zfill(12)

def cigar_to_tuples(cigarstring):
    '''Returns a tuple version of a CIGAR'''
    return re.findall("(\d+)([A-Z])", cigarstring)

def bitflag_to_attributes(flag):
    '''Creates a dictionary of attributes from a SAM bitflag'''
    bits = [int(x) for x in list(flag_to_bin(flag))]
    bits.reverse()
    types = ['read_paired','proper_pair','r1_unmap','r2_unmap',
             'r1_reverse','r2_reverse','1_in_pair','2_in_pair',
             'secondary','quality_fail','PCR_duplicate','supplementary']
    return {k: True if v else False for k, v in zip(types,bits)}

#Variables
flag_values={'read_paired':1,'proper_pair':2,'r1_unmap':4,'r2_unmap':8,
             'r1_reverse':16,'r2_reverse':32,'1_in_pair':64,'2_in_pair':128,
             'secondary':256,'quality_fail':512,'PCR_duplicate':1024,
             'supplementary':2048}
