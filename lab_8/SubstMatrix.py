class SubstMatrix:
    def __init__(self):
        self.alphabet = ""
        self.sm = {}

    def __getitem__(self, ij):
        i, j = ij
        return self.score_pair(i, j)

    def score_pair(self, c1, c2):
        if c1 not in self.alphabet or c2 not in self.alphabet:
            return None
        return self.sm[c1 + c2]

    def read_submat_file(self, filename, sep):
        with open(filename, "r", encoding="utf-8") as f:
            data = f.read()
        header, *body = data.splitlines()
        self.alphabet = "".join([x[0] for x in header.split(sep)])
        for i, line in body:
            for j, token in enumerate(line.split(sep)):
                k = self.alphabet[i]+self.alphabet[j]
                self.sm[k] = int(token)

    def create_submat(self, match, mismatch, alphabet):
        self.alphabet = alphabet
        for c1 in alphabet:
            for c2 in alphabet:
                if c1 == c2:
                    self.sm[c1 + c2] = match
                else:
                    self.sm[c1 + c2] = mismatch


def test1():
    sm = SubstMatrix()
    sm.read_submat_file("blosum62.mat", "\t")
    print(sm.alphabet)
    print(sm.score_pair("G", "M"))
    print(sm.score_pair("W", "W"))
    print(sm.score_pair("A", "S"))
    print(sm.score_pair("X", "X"))
    print(sm["G","K"])
    print(sm["T","T"])

def test2():
    sm = SubstMatrix()
    sm.create_submat(3, -1, "ACGT")
    print(sm.alphabet)
    print(sm.score_pair("A", "A"))
    print(sm.score_pair("A", "T"))
    print(sm.score_pair("T", "T"))
    print(sm["G","G"])


if __name__ == "__main__":
    test1()
    print()
    test2()
