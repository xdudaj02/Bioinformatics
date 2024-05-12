"""
Filename: pwm.py

This file implements the Position Weight Matrix (PWM) representation of a set of motifs.
If run as a script, an illustrative test is performed.

Author: Jakub Duda
Email: up202311235@edu.fc.up.pt

Usage:
    python pwm.py
"""

ALPHABET_DNA = 'ACGT'
ALPHABET_RNA = 'ACGU'
ALPHABET_PROT = 'ACDEFGHIKLMNPQRSTVWY'

class NumMatrix:
    '''Class for a matrix of numbers.'''
    # pylint: disable=too-few-public-methods
    def __init__(self, rows, cols):
        self.mat = []
        for i in range(rows):
            self.mat.append([])
            for _ in range(cols):
                self.mat[i].append(0)

    def __getitem__(self, n):
        return self.mat[n]


class PWM:
    '''Class for Position Weight Matrix (PWM) representation of a set of motifs.'''
    decimals = 2

    def __init__(self, lst_of_words: list[str] = None):
        # Note: There is a discrepancy between Task 1 and Task 3 where the 1st task assumes that
        #   a list of sequences is always supplied to the PWM class constructor and it may be
        #   empty, while the 3rd task assumes that the list of sequences may be omitted completely.
        #   Both approaches are supported.
        self.lst_of_words = lst_of_words if lst_of_words is not None else []
        if not lst_of_words:
            self.bio_type = 'DNA'
            self.alphabet = 'ACGT'
            self.word_length = 0
            self._initialize_empty_pwm()
        else:
            self.word_length = len(lst_of_words[0])
            if not all(len(word) == self.word_length for word in lst_of_words):
                raise ValueError('All words must be of the same length')
            self._infer_bio_type()
            self._infer_alphabet()
            self._calculate_pwm()

    def __str__(self):
        header = 'Base\t'
        for i in range(self.word_length):
            header += f'Pos {i + 1}\t'
        body = ''
        for i, base in enumerate(self.alphabet):
            body += f'\n{base}\t'
            for j in range(self.word_length):
                body += f'{self.pwm[i][j]:.{self.decimals}f}\t'
        return header + body

    def _infer_bio_type(self):
        # Infer the sequence type from the first word in lst_of_words
        first_word = self.lst_of_words[0]
        if 'U' in first_word:
            self.bio_type = 'RNA'
        elif set(first_word).issubset(ALPHABET_DNA):
            self.bio_type = 'DNA'
        elif set(first_word).issubset(ALPHABET_PROT):
            self.bio_type = 'Protein'
        else:
            raise ValueError('Unable to infer sequence type from motif words')

    def _infer_alphabet(self):
        # Infer the alphabet from the bio_type
        match self.bio_type:
            case 'DNA':
                self.alphabet = ALPHABET_DNA
            case 'RNA':
                self.alphabet = ALPHABET_RNA
            case 'Protein':
                self.alphabet = ALPHABET_PROT
            case _:
                raise ValueError('Unknown bio_type')

    def _initialize_empty_pwm(self):
        # Initialize the frequency matrix with zeros
        self.pwm = NumMatrix(len(self.alphabet), self.word_length)

    def _calculate_pwm(self):
        # Note: From the task description, it is not clear whether the 'pwm' matrix should hold the
        #   frequency (absolute frequency/count) or the probability (relative frequency) values.
        #   This implementation calculates the relative frequencies.
        # Calculate PWM matrix
        self._initialize_empty_pwm()
        # Calculate the frequency matrix
        for word in self.lst_of_words:
            for j, base in enumerate(word):
                i = self.alphabet.index(base)
                self.pwm[i][j] += 1
        # Calculate the probability matrix from the frequency matrix
        for i in range(len(self.alphabet)):
            for j in range(self.word_length):
                self.pwm[i][j] /= len(self.lst_of_words)

    def updateSequenceInPWM(self, seq: str):
        '''Update the PWM matrix with a new sequence'''
        # pylint: disable=invalid-name
        if len(seq) != self.word_length:
            raise ValueError('Trying to add a sequence of different length')

        # Add the new sequence to the list of sequences
        self.lst_of_words.append(seq)

        # Recalculate the frequence matrix
        self._calculate_pwm()

        return self.pwm


def test_pwm():
    '''Test PWM class'''
    # Create an empty PWM instance
    pwm = PWM()
    print(pwm)
    print()

    # Create a PWM instance with a list of sequences
    lst_of_words = ['CCTGCCY', 'CCTTCCY', 'AGCGACA', 'AGAGTCA', 'AGCGACC']
    print(f'Motif words: {lst_of_words}')
    pwm = PWM(lst_of_words)
    print(pwm)
    print()

    # Update the PWM matrix with a new sequence
    new_seq = 'AGCGACA'
    print(f'New sequence: {new_seq}')
    pwm.updateSequenceInPWM(new_seq)
    print(pwm)


if __name__ == '__main__':
    test_pwm()
