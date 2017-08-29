import re
from Bio import SeqIO
import numpy as np
import composite_alphabet as cAB
from collections import defaultdict

IUPAC_dict = {tuple((1, 0, 0, 0)): 'A', tuple((0, 1, 0, 0)): 'C', tuple((0, 0, 1, 0)): 'G', tuple((0, 0, 0, 1)): 'T', \
              tuple((1, 0, 1, 0)): 'R', tuple((0, 1, 0, 1)): 'Y', tuple((0, 1, 1, 0)): 'S', tuple((1, 0, 0, 1)): 'W', \
              tuple((0, 0, 1, 1)): 'K', tuple((1, 1, 0, 0)): 'M', tuple((0, 1, 1, 1)): 'B', tuple((1, 0, 1, 1)): 'D', \
              tuple((1, 1, 0, 1)): 'H', tuple((1, 1, 1, 0)): 'V', tuple((1, 1, 1, 1)): 'N'}
IUPAC_dict_rev = {v: k for k, v in IUPAC_dict.items()}
eps = 1E-3
np.random.seed(1)

def to_IUPAC(letter):
    ratio = [letter.ratio_dict[l] for l in letter.basic_alphabet]
    k = float(sum(ratio))
    for i in range(1, 21):
        if (k / i).is_integer():
            ratio_i = tuple([(float(r) / i) for r in ratio])
            if all([r.is_integer() for r in ratio_i]):
                ratio_i = tuple([int(r) for r in ratio_i])
                try:
                    return IUPAC_dict[ratio_i]
                except:
                    pass
    ratio_100 = [int(r * 100.0 / sum(ratio)) for r in ratio]
    ratio_100[ratio_100.index(max(ratio_100))] += 100 - sum(ratio_100)
    return '(N:{:02}{:02}{:02}{:02})'.format(ratio_100[0], ratio_100[1], ratio_100[2], ratio_100[3])


def change_res(ratio, target_res):
    orig_res = sum(ratio.values())
    if all((target_res * float(v) / orig_res).is_integer() for v in ratio.values()):
        return {l: int(target_res * float(v) / orig_res) for l, v in ratio.items()}


# Expects a string which is one of the IUPAC letters or an (N:########) string where #### sum up to 100 and represent A,C,G,T
def parse_IUPAC(letter_as_IUPAC_string, resolution):
    try:
        ratio = {l: v for l, v in zip('ACGT', IUPAC_dict_rev[letter_as_IUPAC_string])}
        return cAB.composite_letter(('A', 'C', 'G', 'T'), change_res(ratio, resolution))
    except Exception as e:
        a, c, g, t = map(int, re.match(pattern='\(N:([0-9]{2})([0-9]{2})([0-9]{2})([0-9]{2})\)',
                                       string=letter_as_IUPAC_string).groups())
        ratio = {'A': a, 'C': c, 'G': g, 'T': t}
        return cAB.composite_letter(('A', 'C', 'G', 'T'), change_res(ratio, resolution))


def parse_IUPAC_sequence(seq, resolution):
    letter_matches = []
    # IUPAC letters
    pat = '([NHARTCSKDVMBWGY])'
    res = re.finditer(pattern=pat, string=seq)
    for r in res:
        grp1 = r.group()
        pos = r.start()
        if pos == len(seq)-1 or seq[pos+1] != ':':
            letter_matches.append((grp1, pos))
    # "real" composite
    pat = '\(N:[0-9]{2}[0-9]{2}[0-9]{2}[0-9]{2}\)'
    res = re.finditer(pattern=pat, string=seq)
    for r in res:
        letter_matches.append((r.group(), r.start()))
    letter_matches.sort(key=lambda x: x[1])
    return [parse_IUPAC(l, resolution) for l, _ in letter_matches]


def read_composite_fasta(filename, resolution):
    with open(filename, 'r') as f:
        for title, seq in SeqIO.FastaIO.SimpleFastaParser(f):
            yield (title, parse_IUPAC_sequence(seq, resolution))


def simulate_composite_letter(letter, depth, eps=eps):
    freqs = letter.freqs()
    freqs = [(f + eps) / (1 + 4 * eps) for f in freqs]
    return np.random.multinomial(depth, freqs)


def simulate_composite_reads(seq, depth, seq_name):
    letters = ('A', 'C', 'G', 'T')
    positions = []
    freqs = [simulate_composite_letter(l, depth) for l in seq]
    for pos in range(len(seq)):
        pos_letters = [ll for sublist in [l * f for l, f in zip(letters, freqs[pos])] for ll in sublist]
        positions.append(np.random.permutation(pos_letters))
    read_seqs = zip(*positions)
    temp = 'sequence:{};read:{}'
    reads = []
    for idx, r in enumerate(read_seqs):
        reads.append(SeqIO.SeqRecord(r, id=temp.format(seq_name, idx), name=temp.format(seq_name, idx),
                                     description=temp.format(seq_name, idx)
                                     , letter_annotations={'phred_quality': [40, ] * len(r)}))
    return reads

def composite_reads_to_ratios((bc, payloads), alphabet, expected_len):
    ratios = alphabet.ratios_from_payloads(payloads, expected_len)[1]
    return ratios

def ratios_to_oligo((bc, ratios), alphabet, expected_len):
    if ratios is not None:
        oligo = cAB.composite_word([alphabet.identify_letter(ratios[i]) for i in range(expected_len)])
        return bc,oligo

# reads a fastq file and split the reads according to the barcode located at the first comp['BC_bytes'] bytes
def demultiplex_reads(reads, comp):
    BC_pos = comp['BC_bytes']*4
    sample_reads = defaultdict(list)
    for idx, r in enumerate(reads):
        bc = str(r[0:BC_pos])
        sample_reads[bc].append(r)
    return sample_reads

def simulate_fasta_reads(fasta_file, depth, out_file,k):
    seqs = read_composite_fasta(fasta_file,k)
    with open(out_file, 'w') as outf:
        outf.write('')
    for idx,seq in enumerate(seqs):
        seq_reads = simulate_composite_reads(seq[1],depth,idx)
        with open(out_file,'a') as outf:
            outf.writelines([''.join(r.seq)+'\n' for r in seq_reads])