# -*- coding: utf-8 -*-

def read_seq_from_file(filename):
    """ Reads a sequence from a multi-line text file. """
    fh = open(filename, "r")
    lines = fh.readlines()
    seq = ""
    for l in lines:
        seq += l.replace("\n","")
    fh.close()
    return seq

def write_seq_to_file(seq, filename):
    """ Writes a sequence to file. """
    with open(filename) as f:
        f.write(seq)


def read_genetic_code_from_file(filename):
    """ Reads the genetic code to a dictionary from a multi-line text file. """
    with open(filename) as f:
        lines = f.readlines()
    gc = {}
    for l in lines:
        gc[l[1:4]] = l[7]

def validate_dna (dna_seq):
    """ Checks if DNA sequence is valid. Returns True is sequence is valid, or False otherwise. """
    seqm = dna_seq.upper()
    valid = seqm.count("A") + seqm.count("C") + seqm.count("G") + seqm.count("T")
    return valid == len(seqm)

def frequency (seq):
    """ Calculates the frequency of each symbol in the sequence. Returns a dictionary. """
    dic = {}
    for s in seq.upper():
        if s in dic: dic[s] += 1
        else: dic[s] = 1
    return dic

def gc_content (dna_seq):
    """ Returns the percentage of G and C nucleotides in a DNA sequence. """
    gc_count = 0
    for s in dna_seq:
        if s in "GCgc": gc_count += 1
    return gc_count / len(dna_seq)

def gc_content_subseq (dna_seq, k=100):
    """ Returns GC content of non-overlapping sub-sequences of size k. """
    res = []
    for i in range(0, len(dna_seq)-k+1, k):
        subseq = dna_seq[i:i+k]
        res.append(gc_content(subseq))
    return res


def transcription (dna_seq):
    """ Function that computes the RNA corresponding to the transcription of the DNA sequence provided. """
    assert validate_dna(dna_seq), "Invalid DNA sequence"
    return dna_seq.upper().replace("T","U")


def reverse_complement (dna_seq):
    """ Computes the reverse complement of the inputted DNA sequence. """
    assert validate_dna(dna_seq), "Invalid DNA sequence"
    comp_rev = ""
    for s in dna_seq.upper()[::-1]:
        if s == "A": comp_rev += "T"
        elif s == "T": comp_rev += "A"
        elif s == "C": comp_rev += "G"
        elif s == "G": comp_rev += "C"
    return comp_rev


def translate_codon (cod):
    """Translates a codon into an aminoacid using an internal dictionary with the standard genetic code."""
    tc = {"GCT":"A", "GCC":"A", "GCA":"A", "GCG":"A",
      "TGT":"C", "TGC":"C",
      "GAT":"D", "GAC":"D",
      "GAA":"E", "GAG":"E",
      "TTT":"F", "TTC":"F",
      "GGT":"G", "GGC":"G", "GGA":"G", "GGG":"G",
      "CAT":"H", "CAC":"H",
      "ATA":"I", "ATT":"I", "ATC":"I",
      "AAA":"K", "AAG":"K",
      "TTA":"L", "TTG":"L", "CTT":"L", "CTC":"L", "CTA":"L", "CTG":"L",
      "ATG":"M", "AAT":"N", "AAC":"N",
      "CCT":"P", "CCC":"P", "CCA":"P", "CCG":"P",
      "CAA":"Q", "CAG":"Q",
      "CGT":"R", "CGC":"R", "CGA":"R", "CGG":"R", "AGA":"R", "AGG":"R",
      "TCT":"S", "TCC":"S", "TCA":"S", "TCG":"S", "AGT":"S", "AGC":"S",
      "ACT":"T", "ACC":"T", "ACA":"T", "ACG":"T",
      "GTT":"V", "GTC":"V", "GTA":"V", "GTG":"V",
      "TGG":"W",
      "TAT":"Y", "TAC":"Y",
      "TAA":"_", "TAG":"_", "TGA":"_"}
    if cod in tc: return tc[cod]
    else: return None


def translate_seq (dna_seq, ini_pos = 0):
    """ Translates a DNA sequence into an aminoacid sequence. """
    assert validate_dna(dna_seq), "Invalid DNA sequence"
    seqm = dna_seq.upper()
    seq_aa = ""
    for i in range(ini_pos, len(seqm)-2, 3):
        cod = seqm[i:i+3]
        aa = translate_codon(cod)
        if aa: seq_aa += aa
        else: seq_aa += "X"
    return seq_aa


def codon_usage(dna_seq, aa):
    """Provides the frequency of each codon encoding a given aminoacid, in a DNA sequence ."""
    assert validate_dna(dna_seq), "Invalid DNA sequence"
    seqm = dna_seq.upper()
    dic = {}
    total = 0
    for i in range(0, len(seqm)-2, 3):
        cod = seqm[i:i+3]
        if translate_codon(cod) == aa:
            if cod in dic:
                dic[cod] += 1
            else: dic[cod] = 1
            total += 1
    if total >0:
        for k in dic:
            dic[k] /= total
    return dic


def reading_frames (dna_seq):
    """Computes the six reading frames of a DNA sequence (including the reverse complement."""
    assert validate_dna(dna_seq), "Invalid DNA sequence"
    res = []
    res.append(translate_seq(dna_seq,0))
    res.append(translate_seq(dna_seq,1))
    res.append(translate_seq(dna_seq,2))
    rc = reverse_complement(dna_seq)
    res.append(translate_seq(rc,0))
    res.append(translate_seq(rc,1))
    res.append(translate_seq(rc,2))
    return res


def all_proteins_rf (aa_seq):
    """Computes all posible proteins in an aminoacid sequence."""
    aa_seq = aa_seq.upper()
    current_prot = []
    proteins = []
    for aa in aa_seq:
        if aa == "_":
            if current_prot:
                for p in current_prot:
                    proteins.append(p)
                current_prot = []
        else:
            if aa == "M":
                current_prot.append("")
            for i in range(len(current_prot)):
                current_prot[i] += aa
    return proteins

def all_orfs (dna_seq):
    """Computes all possible proteins for all open reading frames."""
    assert validate_dna(dna_seq), "Invalid DNA sequence"
    res = []
    rfs = reading_frames(dna_seq)
    for rf in rfs:
        prots = all_proteins_rf(rf)
        for p in prots:
            res.append(p)
    return res

def all_orfs_ord (dna_seq, minsize = 0):
    """Computes all possible proteins for all open reading frames. Returns ordered list of proteins with minimum size."""
    assert validate_dna(dna_seq), "Invalid DNA sequence"
    rfs = reading_frames (dna_seq)
    res = []
    for rf in rfs:
        prots = all_proteins_rf(rf)
        for p in prots:
            if len(p) >= minsize:
                insert_prot_ord(p, res)
    return res

def insert_prot_ord (prot, list_prots):
    ''' inserts prot in list_prots in a sorted way '''
    i = 0
    while i < len(list_prots) and len(prot) < len(list_prots[i]):
        i += 1
    list_prots.insert(i, prot)
    

def longest_protein_sequence(seq):
    return all_orfs_ord(seq)[0]


def read_fasta_2dictionary(filename):
    with open(filename) as f:
        data = f.read()
    data = data.split('>')[1:]
    res = {}
    for x in data:
        key, value = x.split('\n', 1)
        key = key.split(' ', 1)[0]
        res[key] = value.replace('\n', '')
    return res


def test_frequency():
    seq_aa = input("Protein sequence:")
    freq_aa = frequency(seq_aa)
    list_f = sorted(freq_aa.items(), key=lambda x: x[1], reverse = True)
    for (k,v) in list_f:
        print("Aminoacid:", k, ":", v)


def test_sequence_from_file():
    filename = input("File name:")
    seq = read_seq_from_file(filename)
    transcribed_seq = transcription(seq)
    reverse_compl = reverse_complement(seq)
    gc_cont = gc_content(seq)
    translated_seq = translate_seq(seq)
    for x in [seq, transcribed_seq, reverse_compl, gc_cont, translated_seq]:
        print(x)


def test_longest_prot_seq():
    filename = input("File name:")
    seq = read_seq_from_file(filename)
    print(longest_protein_sequence(seq))   


if __name__ == "__main__":
    print(read_fasta_2dictionary('code_and_data\PS00727.fasta'))
