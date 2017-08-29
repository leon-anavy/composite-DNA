# -*- coding: utf-8 -*-
"""
Created on Wed Jul 06 14:51:30 2016

@author: t-leonan
"""
import pickle, operator, re
import numpy as np

def L1_dist(ratio1,ratio2):
    return sum(abs(r1-r2) for (r1,r2) in zip(ratio1,ratio2))

def KL_dist(ratio1,ratio2):
    eps = 0.001
    tweaked_ratio2 = [r + eps for r in ratio2]
    tweaked_ratio2 = [r/sum(tweaked_ratio2) for r in tweaked_ratio2]
    return sum(r1*np.log(r1/r2) if r1 > 0.0 else 0 for (r1,r2) in zip(ratio1,tweaked_ratio2))

class composite_letter(object):
    def __repr__(self):
        return repr(self.ratio_dict)
    
    def __str__(self):
        return '{}'.format(''.join((str(self.ratio_dict[l])+l for l in self.basic_alphabet)))

    def stacked_str(self):
        return '{}'.format('\n'.join((str(self.ratio_dict[l]) + l for l in self.basic_alphabet)))
        
    def __init__(self, basic_alphabet, ratio_dict):
        self.basic_alphabet = basic_alphabet
        self.ratio_dict = ratio_dict

    basic_alphabet = property(operator.attrgetter('_basic_alphabet'))

    @basic_alphabet.setter
    def basic_alphabet(self,a):
        if len(a) == 0:
            raise Exception('basic alphabet must be non-empty')
        self._basic_alphabet = a

    ratio_dict = property(operator.attrgetter('_ratio_dict'))

    @ratio_dict.setter
    def ratio_dict(self,r):
        if not all(l in r.keys() for l in self.basic_alphabet):
            raise Exception('ratio must match basic alphabet')
        self._ratio_dict = r

    def __add__(self,other):
        return composite_word([self, other])

    def __hash__(self):
        return hash((self.basic_alphabet,tuple(self.ratio_dict[l] for l in self.basic_alphabet)))

    def __eq__(self, other):
        return type(self) == type(other) and (self.basic_alphabet,self.ratio_dict) == (other.basic_alphabet,other.ratio_dict)

    def __ne__(self, other):
        # Not strictly necessary, but to avoid having both x==y and x!=y
        # True at the same time
        return not(self == other)

    def freqs(self):
        ratio = self.ratio_dict
        alp = self.basic_alphabet
        return tuple(float(ratio[l])/sum(ratio.values()) for l in alp)

    # Expects a string in the format #A#C#G#T where # stands for the number of each letter in the ratio (e.g 1A2C0G0T)
    @classmethod
    def parse(cls,letter_as_string):
        mm = re.finditer(pattern='([0-9]*)([A-Z])',string=letter_as_string)
        basic_alphabet = []
        ratio = {}
        for m in mm:
            r,l = m.groups()
            basic_alphabet.append(l)
            ratio[l] = int(r)
        return cls(basic_alphabet,ratio)
    
class composite_word(object):
    def __repr__(self):
        return repr(self.letters)
    
    def __str__(self):
        return '|'.join(str(l) for l in self.letters)

    def __len__(self):
        return len(self.letters)
        
    def __init__(self, letters):
        self.letters = letters
        
    letters = property(operator.attrgetter('_letters'))
    
    @letters.setter
    def letters(self,l):
        l1 = l[0]
        if not all(ll.basic_alphabet == l1.basic_alphabet for ll in l):
            raise Exception('All letters must have the same basic alphabet')
        # disabled for IUPAC
#        if not all(sum(ll.ratio_dict.values()) == sum(l1.ratio_dict.values()) for ll in l):
#            raise Exception('All letters must have the same resolution')   
        self._letters = l
            
    def __add__(self,other):
        if other is '':
            return self
        if type(other) is composite_word:
            return composite_word(self.letters + other.letters)
        if type(other) is composite_letter:
            return composite_word(self.letters + [other])

    def __hash__(self):
        return hash(tuple(self.letters))

    def __eq__(self, other):
        return tuple(self.letters) == tuple(other.letters)

    def __ne__(self, other):
        # Not strictly necessary, but to avoid having both x==y and x!=y
        # True at the same time
        return not(self == other)

    def __getitem__(self, item):
        return self.letters.__getitem__(item)
    
    @classmethod
    def parse(cls,word_as_string):
        letters_as_strings = word_as_string.split('|')
        return cls([composite_letter.parse(l) for l in letters_as_strings])

        
class composite_alphabet(object):
    def __repr__(self):
        return repr(self.letters)
    
    def __str__(self):
        return '[{}]'.format(','.join(str(l) for l in self.letters))
        
    def __len__(self):
        return len(self.letters)
        
    def __init__(self, basic_alphabet, ratios, dist = None):
        self.basic_alphabet = basic_alphabet
        self.letters = [composite_letter(basic_alphabet,r) for r in ratios]
        self.dist = dist

    dist = property(operator.attrgetter('_dist'))

    @dist.setter
    def dist(self,dist_str):
        if dist_str is None or dist_str == 'L1':
            self._dist = L1_dist
        elif dist_str == 'KL':
            self._dist = KL_dist
        else:
            raise Exception('invalid distance measure')

    def get_dist(self):
        if self.dist == L1_dist:
            return 'L1'
        elif self.dist == KL_dist:
            return 'KL'
        else:
            raise Exception('something went terribly wrong')

    @classmethod
    def from_letters(cls, letters):
        basic_alphabet = letters[0].basic_alphabet
        ratios = [l.ratio_dict for l in letters]
        return cls(basic_alphabet, ratios)

    # expects observed_ratio as a dictionary of "letter:count" items
    def dist_from_letter(self,observed_ratio,letter):
        observed_ratio_p = [float(observed_ratio[l]) / sum(observed_ratio.values()) for l in self.basic_alphabet]
        letter_ratio_p = [float(letter.ratio_dict[l]) / sum(letter.ratio_dict.values()) for l in self.basic_alphabet]
        try:
            return self.dist(observed_ratio_p,letter_ratio_p)
        except:
            return L1_dist(observed_ratio_p,letter_ratio_p)
    
    def identify_letter(self,observed_ratio,letters = None):
        if letters is None:
            letters = self.letters
        # change observed to sum 1
        # set alphabet resolution and possible values
        #   for value in observed:
        #       if close enough to possible value:
        #           set as constant
        # for letter in letters:
        #   change letter to sum 1
        #   if constants not equal to letter values:
        #       continue
        #   check distance and keep if minimal so far
        return(min(self.letters,key = lambda l: self.dist_from_letter(observed_ratio,l)))
    
    def get_two_closest_letters(self, observed_ratio, letters = None):
        if letters is None:
            letters = self.letters
        letters_and_dists = [(l,self.dist_from_letter(observed_ratio,l)) for l in letters]
        return(sorted(letters_and_dists, key = lambda x:x[1])[:2])
    
    def as_words(self):
        return [composite_word([l]) for l in self.letters]

    def to_file(self,filename):
        pickle.dump(self,filename)
        
def from_file(filename):
    return pickle.load(filename)