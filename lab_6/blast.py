# -*- coding: utf-8 -*-
from collections import namedtuple
from sequence_alignments import create_submat, get_global_alignment, get_local_alignment

alignment = namedtuple('alignment', ['start', 'start_seq', 'size', 'score', 'seq_index'])
hit_extended = namedtuple('hit_extended', ['start', 'start_seq', 'size', 'score'])
hit_base = namedtuple('hit_base', ['start', 'start_seq'])

class MyBlast:
    '''
    Class for the point matrices.
    '''

    def __init__(self, filename: str = None, w: int = 3, only_seq: bool = True):
        '''
        Constructor
        '''
        if filename is not None:
            self.read_db(filename, only_seq)
        else:
            self.db = []
        self.w = w
        self.map = None

    def read_db(self, filename: str, only_seq: bool):
        """From file with sequences line by line read the sequences to a list"""
        if only_seq:
            with open(filename, 'r') as f:
                lines = f.readlines()
            self.db = [line.strip() for line in lines]
        else:
            self.db = read_fas_file(filename, only_seq = True)
        
    def add_sequence_to_db(self, seq: str):
        """Add an extra sequence to DB"""
        self.db.append(seq)

    def build_map(self, query: str):
        self.map = {}
        for i in range(0, len(query) - self.w):
            subseq = query[i: i + self.w]
            if subseq in self.map:
                self.map[subseq].append(i)
            else:
                self.map[subseq] = [i]

    def get_hits(self, seq: str, mismatch: int = 0) -> list[hit_base]:
        hits = []
        if mismatch > self.w:
            raise ValueError("Mismatch must be less than window size")
        for i in range(0, len(seq) - self.w):
            subseq = seq[i: i + self.w]
            if mismatch > 0:
                for k in self.map:
                    if hamming_distance(k, subseq) <= mismatch:
                        for j in self.map[k]:
                            hits.append(hit_base(j, i))
            else:
                if subseq in self.map:
                    for j in self.map[subseq]:
                        hits.append(hit_base(j, i))
        return hits
    
    def extends_hit(self, seq: str, hit: hit_base, query: str) -> hit_extended:
        # move forward
        matfw = 0
        k = 0
        bestk = 0
        while 2 * matfw >= k and hit.start + self.w + k < len(query) and hit.start_seq + self.w + k < len(seq):
            if query[hit.start + self.w + k] == seq[hit.start_seq + self.w + k]:
                matfw += 1
                bestk = k + 1
            k += 1
        size = self.w + bestk

        # move backward
        k = 0
        matbw = 0
        bestk = 0
        while 2 * matbw >= k and hit.start > k and hit.start_seq > k:
            if query[hit.start - k - 1] == seq[hit.start_seq - k - 1]:
                matbw += 1
                bestk = k + 1
            k += 1
        size += bestk

        return hit_extended(hit.start - bestk, hit.start_seq - bestk, size, self.w + matfw + matbw)

    def hit_best_score(self, seq: str, query: str, mismatch: int = 0) -> hit_extended:
        hits = self.get_hits(seq, mismatch)
        best_score = -1.0
        best = None
        for h in hits:
            ext = self.extends_hit(seq, h, query)
            score = ext.size
            if score > best_score or (score == best_score and ext.score < best.score):
                best_score = score
                best = ext
        return best


    def best_alignment(self, query: str, n: int = 1, mismatch: int = 0) -> list[alignment]:
        self.build_map(query)
        res = []
        for i, item in enumerate(self.db):
            if (hit := self.hit_best_score(item, query, mismatch)):
                res.append(alignment(hit.start, hit.start_seq, hit.size, hit.score, i))
        return sorted(res, key = lambda x: x.score, reverse = True)[:min(n, len(res))]

def read_fas_file(filename: str, only_seq: bool = False) -> dict[str, str] | list[str]:
    with open(filename, 'r') as f:
        content = f.readlines()
    content_dict = {}
    for line in content:
        if line[0] == '>':
            key = line[1:]
            content_dict[key] = ''
        else:
            content_dict[key] += line
        
    return content_dict if not only_seq else list(content_dict.values())

def hamming_distance(seq1: str, seq2: str) -> int:
    return sum([1 for i in range(len(seq1)) if seq1[i] != seq2[i]])

# Tests
def test_1():
    mb = MyBlast("data\\seqBlast.txt", 11)
    query = "gacgcctcgcgctcgcgcgctgaggcaaaaaaaaaaaaaaaaaaaatcggatagctagctgagcgctcgatagcgcgttcgctgcatcgcgtatagcgctgaagctcccggcgactgtctgtaaatcggatctcatctcgctctatcct"
    r = mb.best_alignment(query)
    
    print('test_1')
    print(r)
    print()

def test_2():
    mb = MyBlast("data\\seqBlast.txt", 11)
    query2 = "cgacgacgacgacgaatgatg"
    r = mb.best_alignment(query2)
    
    print('test_2')
    print(r)
    print()

def test_3():
    mb = MyBlast("data\\seqBlast.txt", 11)
    with open("data\\query1.fasta", 'r') as f:
        query1 = f.read()
    with open("data\\query2.fasta", 'r') as f:
        query2 = f.read()
    r1 = mb.best_alignment(query1)
    r2 = mb.best_alignment(query2)
    
    print('test_3')
    print(r1)
    print(r2)
    print()

def test_4():
    mb = MyBlast("data\\db.fas", 11, True)
    query = read_fas_file("data\\query.fas", only_seq = True)[0]
    r = mb.best_alignment(query)
    seq = mb.db[r[0].seq_index]
    # todo: change the implementation of the functions implementeing alignment recovery because they exceed the maximum recursion depth
    # sm = create_submat(1, -1, set(query + seq))
    # global_alignment = get_global_alignment(query, seq, sm, -8)
    # local_alignment = get_local_alignment(query, seq, sm, -8)
    
    print('test_4')
    print(r)
    # print(global_alignment)
    # print(local_alignment)
    print()

def test_5():
    query = "ATTGTTG"
    target = "ATGAAATGTGCCAGTCGATC"
    mb = MyBlast()
    mb.add_sequence_to_db(target)
    r = mb.best_alignment(query, mismatch = 1)

    print('test_5')
    print(r)
    print()

def test_6():
    mb = MyBlast("data\\db.fas", 5, True)
    query = read_fas_file("data\\query.fas", only_seq = True)[0]
    r = mb.best_alignment(query, n = 10)

    print('test_6')
    for i in r:
        print(i)

if __name__ == "__main__":
    test_1()
    test_2()
    test_3()
    test_4()
    test_5()
    test_6()