'''PWM'''
import math
from NumMatrix import NumMatrix

class PWM:
    '''Class for Position Weight Matrix (PWM) representation of a set of words.'''
    def __init__(self, lst_of_words: list[str]):
        if not lst_of_words:
            self.bio_type = 'DNA'
            self.alphabet = 'ACGT'
            self.pwm = None
            self.lst_of_words = None
            self.word_length = None
        else:
            self.word_length = len(lst_of_words[0])
            if not all(len(word) == self.word_length for word in lst_of_words):
                raise ValueError('All words must be of the same length')
            self.lst_of_words = lst_of_words
            self._infer_bio_type()
            self._infer_alphabet()
            self._calculate_pwm()

    def _infer_bio_type(self):
        # Infer the sequence type from the first word in lst_of_words
        first_word = self.lst_of_words[0]
        if 'U' in first_word:
            self.bio_type = 'RNA'
        elif set(first_word).issubset('ACGTU'):
            self.bio_type = 'DNA'
        elif set(first_word).issubset('ABCDEFGHIKLMNPQRSTVWXYZ'):
            self.bio_type = 'Protein'
        else:
            raise ValueError('Unable to infer sequence type from motif words')

    def _infer_alphabet(self):
        # Infer the alphabet from the bio_type
        if self.bio_type == 'DNA':
            self.alphabet = 'ACGT'
        elif self.bio_type == 'RNA':
            self.alphabet = 'ACGU'
        elif self.bio_type == 'Protein':
            self.alphabet = 'ACDEFGHIKLMNPQRSTVWY'
        else:
            raise ValueError('Unknown bio_type')

    def _initialize_empty_pwm(self):
        # Initialize the frequency matrix with zeros
        self.pwm = NumMatrix(len(self.alphabet), self.word_length)

    def _calculate_pwm(self):
        # Calculate PWM matrix
        # todo: should this be pfm, ppm or pwm?
        self._initialize_empty_pwm()
        for word in self.lst_of_words:
            for i, base in enumerate(word):
                # self.pwm[base][i] += 1
                letter_index = self.alphabet.index(base)
                self.pwm[letter_index][i] += 1
        self.add_pseudocounts()
        for i, base in enumerate(self.alphabet):
            for j in range(self.word_length):
                # self.pwm[base][j] /= len(self.lst_of_words)
                self.pwm[i][j] /= len(self.lst_of_words)

    def add_pseudocounts(self, pseudo: float = 0.01):
        '''Add pseudocounts to the PWM matrix.'''
        # for base in self.alphabet:
        #     self.pwm[base] = [freq + pseudo for freq in self.pwm[base]]
        for i, _ in enumerate(self.alphabet):
            for j in range(self.word_length):
                # todo: add to all or only to zeros?
                # if not self.pwm[i][j]:
                self.pwm[i][j] += pseudo

    def print_pwm(self):
        '''Print PWM matrix with labels'''
        print('PWM Matrix:')
        print('Base\t', end='')
        for i in range(self.word_length):
            print(f'Pos {i + 1}\t', end='')
        print()
        for i, base in enumerate(self.alphabet):
            print(base, end='\t')
            # for freq in self.pwm[base]:
            #     print(f'{freq:.2f}\t', end='')
            for j in range(self.word_length):
                print(f'{self.pwm[i][j]:.3f}\t', end='')
                # print(f'{self.pwm[i][j]:g}\t', end='')
            print()

    def information_content(self):
        '''Calculate information content for each position in the motif'''
        info_content = []
        for j in range(self.word_length):
            # todo: whats the correct formula and what to calculate from (ppm vs pwm)
            # pos_ic = 0
            pos_ic = 2
            for i, _ in enumerate(self.alphabet):
                # freq = self.pwm[base][j]
                freq = self.pwm[i][j]
                if freq > 0:
                    pos_ic += freq * math.log2(freq)
            info_content.append(pos_ic)
        return info_content

    def consensus(self):
        '''Derive the consensus word from the PWM matrix'''
        consensus_word = ''
        for i in range(self.word_length):
            max_freq = 0
            max_base = ''
            for j, base in enumerate(self.alphabet):
                freq = self.pwm[j][i]
                if freq > max_freq:
                    max_freq = freq
                    max_base = base
            consensus_word += max_base
        return consensus_word

    def seq2freqnorm(self, lst_of_seqs: list[str]):
        '''Calculate PWM matrix from a list of sequences'''
        word_length = len(lst_of_seqs[0])
        pwm = NumMatrix(len(self.alphabet), word_length)
        for word in lst_of_seqs:
            for i, base in enumerate(word):
                letter_index = self.alphabet.index(base)
                pwm[letter_index][i] = pwm[letter_index][i] + 1
        for i, base in enumerate(self.alphabet):
            for j in range(word_length):
                relative_value = pwm[i][j] / len(lst_of_seqs)
                pwm[i][j] = relative_value

    def score_sequence_by_motif(self, target_seq: str):
        '''Score a target sequence using the PWM matrix'''
        seq_length = len(target_seq)

        scores = []
        max_score = -1
        max_score_pos = None

        # Slide the PWM along the target sequence and calculate scores
        for i in range(seq_length - self.word_length + 1):
            score = 1
            for j in range(self.word_length):
                base = target_seq[i + j]
                if base in self.alphabet:
                    # todo: when calculating from ppm, multiplication should be used
                    score *= self.pwm[self.alphabet.index(base)][j]
            scores.append(score)

            # Update max_score and max_score_pos if needed
            if score > max_score:
                max_score = score
                max_score_pos = i

        return (scores, max_score_pos)

    def update_sequence_in_pwm(self, seq: str):
        '''Update the PWM matrix with a new sequence'''
        # todo: pseudocounts ?
        if len(seq) != self.word_length:
            raise ValueError("Input sequence length doesn't match the PWM length")

        # Add the new sequence to the list of sequences
        self.lst_of_words.append(seq)

        # Recalculate the frequence matrix
        self._calculate_pwm()

        return self.pwm

    def log_odds(self, bckgrd_freq: dict[str, float] = None):
        '''Calculate log odds scores for the PWM matrix'''
        if bckgrd_freq and not set(self.alphabet).issubset(set(bckgrd_freq)):
            raise ValueError('Background frequencies not provided for all symbols in the alphabet')

        # If background frequencies are not provided, assume equal frequencies for all symbols
        if bckgrd_freq is None:
            # todo: this or always 0,25 ?
            freq = 1 / len(self.alphabet)
            bckgrd_freq = {base: freq for base in self.alphabet}

        # Initialize matrix to store log odds scores
        log_odds_matrix = self.pwm.copy()

        # Calculate log odds scores for each symbol in the PWM matrix
        for i, base in enumerate(self.alphabet):
            for j in range(self.word_length):
                log_odds_matrix[i][j] = math.log2(self.pwm[i][j] / bckgrd_freq[base])

        return log_odds_matrix


def test_pwm():
    '''Test PWM class with example motif words'''
    # Example motif words
    lst_of_words = ['CCTGCCY', 'CCTTCCY', 'AGCGACA', 'AGAGTCA', 'AGCGACC']

    # Create PWM instance
    pwm = PWM(lst_of_words)

    # Add pseudocounts
    # pwm.add_pseudocounts()

    # Print PWM matrix
    pwm.print_pwm()

    # Calculate information content
    info_content = pwm.information_content()
    print(f'Information content for each position: {info_content}')

    # Get consensus word
    consensus_word = pwm.consensus()
    print(f'Consensus word: {consensus_word}')

    # Score a target sequence
    target_seq = 'AGCGACAGTCCY'
    print(f'Target sequence: {target_seq}')
    scores, max_score_pos = pwm.score_sequence_by_motif(target_seq)
    print(f'Scores: {scores}')
    print(f'Max score position: {max_score_pos}')

    # Update the PWM matrix with a new sequence
    new_seq = 'AGCGACA'
    print(f'New sequence: {new_seq}')
    pwm.update_sequence_in_pwm(new_seq)
    pwm.print_pwm()

    # Calculate log odds scores
    log_odds_matrix = pwm.log_odds()
    log_odds_matrix.print_mat()


if __name__ == '__main__':
    test_pwm()

# todo: triple check everything is in accordance with the assignment
