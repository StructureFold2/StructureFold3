'''
Gini function modified from 
https://planspacedotorg.wordpress.com/2013/06/21/how-to-calculate-gini-coefficient-from-raw-data-in-python/
'''

#Imports
import numpy

#Classes
class ReactStats(object):
    '''Just holds all the components of an entry in a clean way'''
    def __init__(self,react_vector=[],trim=0,minlen=0):
        if react_vector:
            vektor = react_vector[:-trim] if trim else react_vector
            vektor = filter(lambda x: isinstance(x, float),vektor)
            if len(vektor) >= minlen:
                try:
                    self.mahx = max(vektor)
                    self.average = numpy.average(vektor)
                    self.std = numpy.std(vektor)
                    self.gini = gini(vektor)
                except ValueError:
                    self.mahx,self.average,self.std,self.gini = ['NA']*4
            else:
                self.mahx,self.average,self.std,self.gini = ['NA']*4
        else:
            self.mahx,self.average,self.std,self.gini = ['NA']*4
            
    def as_list(self):
        return [self.mahx,self.average,self.std,self.gini]

#Functions
def gini(list_of_values):
    '''Returns the Gini value of a list of values.'''
    try:
        sorted_list = sorted(list_of_values)
        height, area = 0, 0
        for value in sorted_list:
            height += value
            area += height - value / 2.
        fair_area = height * len(list_of_values) / 2.
        return (fair_area - area) / fair_area
    except ZeroDivisionError:
        return 'NA'