# -*- coding: utf-8 -*-
"""
Created on Thu Jul 07 11:14:59 2016
Encodes a bitstream to given alphabet

@author: t-leonan
"""
import math

class encoder:
    def __repr__(self):
        return repr(self.encoding_dict)
    
    def __str__(self):
        return str(self.encoding_dict)
        
    def __init__(self, output_alphabet):
        self.output_alphabet = output_alphabet
        get_bin = lambda x, n: format(x, 'b').zfill(n)
        self.block_size = int(math.floor(math.log(len(output_alphabet),2)))
        self.encoding_dict = {get_bin(v,self.block_size):output_alphabet.letters[v] for v in range(2**self.block_size)}
        self.decoding_dict = {v: k for k, v in self.encoding_dict.iteritems()}

    def encode(self,bit_stream):
        # pad with zeros
        extra = len(bit_stream) % self.block_size
        if extra:
            pad = self.block_size - extra
            padded_stream = bit_stream + '0'*pad
        else:
            padded_stream = bit_stream
        # break to blocks
        blocks = [padded_stream[i:i+self.block_size] for i in range(0,len(padded_stream),self.block_size)]
        return (self.encoding_dict[b] for b in blocks)

    def decode(self,message):
        out_message = ''
        for idx,l in enumerate(message):
            try:
                out_l = self.decoding_dict[l]
            except:
                # just picj a random letter and hope for the best
                out_l = self.decoding_dict.values()[1]
            out_message += out_l
        return out_message
        