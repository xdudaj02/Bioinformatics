class Gene:
    name = ""
    chromosome = ""
    start = 0
    end = 0
    sequence = ""
    numberExons = 0

    def __init__(self, name, chromosome, start, end, sequence, numberExons):
        self.name = name
        self.chromosome = chromosome
        self.start = start
        self.end = end
        self.sequence = sequence
        self.numberExons = numberExons

    def set_name(self, name):
        self.name = name

    def set_chromosome(self, chromosome):
        self.chromosome = chromosome

    def set_start(self, start):
        self.start = start

    def set_end(self, end):
        self.end = end

    def set_sequence(self, sequence):
        self.sequence = sequence

    def set_numberExons(self, numberExons):
        self.numberExons = numberExons

    def get_name(self):
        return self.name
    
    def get_chromosome(self):
        return self.chromosome
    
    def get_start(self):
        return self.start
    
    def get_end(self):
        return self.end
    
    def get_sequence(self):
        return self.sequence
    
    def get_numberExons(self):
        return self.numberExons
    
    def get_gene_length(self):
        return self.end - self.start + 1
    
    def get_gc_content(self):
        return (self.sequence.count("G") + self.sequence.count("C")) / len(self.sequence)
    
    def to_upper(self):
        self.sequence = self.sequence.upper()