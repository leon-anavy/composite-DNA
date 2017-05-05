# -*- coding: utf-8 -*-
"""
Created on Wed Mar 01 15:57:56 2017

@author: t-leonan
"""
import composite_alphabet as cAB
import itertools, pickle

def _generate_ratios(basic_alphabet,resolution):
    alphabet_ratios = []
    candidates = range(resolution+1)
    for grp in itertools.combinations_with_replacement(candidates,len(basic_alphabet)):
        if sum(grp) == resolution:
            for g in set(prmt for prmt in itertools.permutations(grp)):
                alphabet_ratios.append(g)
    return alphabet_ratios

class composite_DNA_alphabet(cAB.composite_alphabet):
    def __repr__(self):
        return repr(self.letters)
    
    def __str__(self):
        return '[{}]'.format(','.join(str(l) for l in self.letters))
        
    def __len__(self):
        return len(self.letters)

    def __init__(self,resolution = None, ratios = None, dist = None):
        basic_alphabet = ('A','C','G','T')
        if ratios is None:
            ratios = _generate_ratios(basic_alphabet,resolution)
        ratio_dicts = [{l:v for l,v in zip(basic_alphabet,ratio)} for ratio in ratios ]
        super(composite_DNA_alphabet,self).__init__(basic_alphabet,ratio_dicts,dist)
    
    def ratios_from_payloads(self, payloads, oligo_len):
        counts = [{l:0 for l in self.basic_alphabet} for _ in range(oligo_len)]
        total = 0
        for idx,payload in enumerate(payloads):
            payload = payload.rstrip()
            if 'N' not in payload:
                if len(payload) == oligo_len:
                    total += 1
                    for i in range(oligo_len):
                        counts[i][payload[i]] += 1
        if total == 0:
            return total, None
        ratios = [{l:1.0*counts[i][l]/total for l in self.basic_alphabet} for i in range(oligo_len)]
        return total,ratios
        
def from_file(filename):
    return pickle.load(filename)



