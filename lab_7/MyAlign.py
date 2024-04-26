class MyAlign:
    def __init__(self, lseqs, al_type = "protein"):
        self.listseqs = lseqs
        self.al_type = al_type

    def __len__(self): # number of columns
        return len(self.listseqs[0])

    def __getitem__(self, n):
        if isinstance(n, tuple) and len(n) == 2:
            i, j = n
            return self.listseqs[i][j]
        if isinstance(n, int):
            return self.listseqs[n]
        return None

    def __str__(self):
        return "\n".join(self.listseqs)

    def num_seqs(self):
        return len(self.listseqs)

    def column(self, indice):
        return [seq[indice] for seq in self.listseqs]

    def consensus(self):
        res = ""
        for i in range(len(self)):
            col = self.column(i)
            counts = {}
            for c in col:
                counts[c] = counts.get(c, 0) + 1
            most_common = [k for k, v in counts.items() if v == max(counts.values())]
            most_common.sort()
            if len(most_common) == 1:
                res += most_common[0]
            else:
                res += most_common[0] if most_common[0] != "-" else most_common[1]
        return res

if __name__ == "__main__":
    alig = MyAlign(["ATGA-A","AA-AT-"], "dna")
    print(alig)
    print(len(alig))
    print(alig.column(2))
    print(alig[1,1])
    print(alig[0])
    print(alig.consensus())
