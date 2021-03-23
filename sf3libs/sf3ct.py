#Imports
import re
import os
import glob

class CT_Report(object):
    '''MFE <.ct> file'''
    def __init__(self,ct_file,DeltaG='NA'):
        names = ['index','bases','left','right','paired']
        with open(ct_file,'r') as f:
            #Header
            header = next(f)
            header_bits = header.strip().split()
            self.length = int(header_bits[0])
            self.name = header_bits[-1]
            energy_info = re.findall(r"ENERGY = [-+]?\d*\.\d+",header)
            self.energy = float(energy_info[0].split()[2]) if energy_info else DeltaG
            #Body
            lines = [line.strip().split() for line in f]
            info = [[q[i] for q in lines] for i in range(0,5)]
            for name,data in zip(names,info):
                setattr(self,name,data)
            #Derived Attributes
            self.single_st = self.paired.count('0')/float(self.length)
            self.double_st = 1 - self.single_st

    def stats(self,digits=5):
        '''Returns a list of basic attritubes to report'''
        strands = [round(s,digits) for s in [self.double_st,self.single_st]]
        basic = [self.energy]+strands
        return [str(x) for x in basic]

def collect_cts(directory,dg_val='NA'):
    '''Snags a directory worth of CT files'''
    data = {}
    for fyle in glob.glob(os.path.join(directory,'*.ct')):
        data[os.path.split(fyle)[-1]] = CT_Report(fyle,dg_val)
    return data
