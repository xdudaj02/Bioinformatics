'''
Filename: classify_proteins.py

This script implements the classification of protein sequences based on k-mer counts.
It reads two input files in FASTA format, extracts the protein sequences and for each sequence
counts the number of k-mers of a given length k. The script then runs multiple classifiers (SVM,
Random Forest, Naive Bayes) to classify the protein sequences and outputs the classification
results. The classification results include the mean and standard deviation of the accuracy,
precision, recall and F1-score for each classifier.

Authors:
    - Mariana Lourenço <up201906985@edu.fc.up.pt>
    - Jakub Duda <up202311235@edu.fc.up.pt>
    - Soulaimane Salehddine <up202312271@edu.fc.up.pt>

Usage:
    python classify_proteins.py -a <file_a> -b <file_b> -k <k>

Arguments:
    -a <file_a>: The first input file in FASTA format
    -b <file_b>: The second input file in FASTA format
    -k <k>: The length of the k-mer

Output:
    A table with the classification results for the protein sequences in the input files.
'''

import argparse
import itertools
import pandas as pd
from sklearn.model_selection import StratifiedKFold, cross_validate
from sklearn.svm import SVC
from sklearn.ensemble import RandomForestClassifier
from sklearn.naive_bayes import GaussianNB

PROTEIN_ALPHABET = 'ACDEFGHIKLMNPQRSTVWY'


def parse_args() -> argparse.Namespace:
    '''Parse command line arguments

    Returns:
        argparse.Namespace: The parsed arguments
    '''
    def fasta_filename(file: str) -> str:
        '''Check if file is a FASTA file

        Args:
            file (str): The filename

        Returns:
            str: The filename
        '''
        if not file.endswith('.fasta'):
            raise argparse.ArgumentTypeError('Not a FASTA file')
        return file.split('.fasta')[0]
    parser = argparse.ArgumentParser(description='Classify protein sequences')
    parser.add_argument('-a', help='First input file', required=True, type=fasta_filename)
    parser.add_argument('-b', help='Second input file', required=True, type=fasta_filename)
    parser.add_argument('-k', help='Length of the k-mer', required=True, type=int)
    return parser.parse_args()


def fasta_to_dict(filename: str) -> dict[str, str]:
    '''Reads a FASTA file and returns a dictionary

    Args:
        filename (str): Path to the FASTA file

    Returns:
        dict: A dictionary with the id as the key and the sequence as the value
    '''
    try:
        with open(f'{filename}.fasta', 'r', encoding='utf-8') as f:
            lines = f.readlines()
    except FileNotFoundError as e:
        raise FileNotFoundError(f'File not found: {filename}') from e
    fasta_dict = {}
    while lines:
        try:
            header = lines.pop(0).split('|', 2)[1]
        except IndexError as e:
            raise ValueError('Invalid header in FASTA file') from e
        seq = ''
        while lines and not lines[0].startswith('>'):
            seq += lines.pop(0).strip()
        if header in fasta_dict:
            raise ValueError('Duplicate entry in FASTA file')
        if not seq:
            raise ValueError('Empty sequence in FASTA file')
        fasta_dict[header] = seq
    if not fasta_dict:
        raise ValueError('Empty FASTA file')
    return fasta_dict


def generate_kmers(k: int, char_set: str) -> list[str]:
    '''Generate all possible k-mers

    Args:
        k (int): Length of the k-mer
        char_set (str): The character set to use

    Returns:
        list: A list of all possible k-mers
    '''
    if k < 1:
        raise ValueError('Invalid k-mer length')
    if not char_set:
        raise ValueError('Invalid character set')
    return [''.join(p) for p in itertools.product(char_set, repeat=k)]


def kmer_count(seq: str, k: int, kmers: list[str]) -> dict[str, int]:
    '''Count the number of k-mers in a sequence

    Args:
        seq (str): The protein sequence
        k (int): Length of the k-mer
        kmers (list): List of all possible k-mers

    Returns:
        dict: A dictionary with the k-mer as the key and the count as the value
    '''
    if not seq:
        raise ValueError('Empty sequence')
    if k < 1:
        raise ValueError('Invalid k-mer length')

    kmer_dict = {x: 0 for x in kmers} if kmers else {}
    for i in range(len(seq) - k + 1):
        kmer = seq[i:i + k]
        try:
            kmer_dict[kmer] += 1
        except KeyError:
            # ignore invalid k-mers (e.g. containing X, B, Z, J, O or U)
            continue
    return kmer_dict


def create_kmer_df(
    fasta_dict: dict[str, str],
    class_name: str,
    k: int,
    alphabet: str = PROTEIN_ALPHABET
) -> pd.DataFrame:
    '''Populate a DataFrame with k-mer counts

    Args:
        fasta_dict (dict): Dictionary with the protein sequences
        class_name (str): The class name (of the protein family)
        k (int): Length of the k-mer
        alphabet (str): The character set to use

    Returns:
        pd.DataFrame: A DataFrame with the k-mer counts
    '''
    if not fasta_dict:
        raise ValueError('Empty FASTA dictionary')
    if k < 1:
        raise ValueError('Invalid k-mer length')

    kmers = generate_kmers(k, alphabet)

    kmer_cnts = {seq_id: kmer_count(seq, k, kmers) for seq_id, seq in fasta_dict.items()}

    df = pd.DataFrame.from_dict(kmer_cnts, orient='index')
    df.loc[:, 'class'] = pd.Series(class_name, index=df.index, dtype='category')
    return df


def get_classification_results(df: pd.DataFrame) -> pd.DataFrame:
    '''Classify protein sequences and return the results

    Args:
        df (pd.DataFrame): DataFrame with the k-mer counts

    Returns:
        pd.DataFrame: DataFrame with the classification results
    '''
    classifiers = [
        ('SVM', SVC()),
        ('Random Forest', RandomForestClassifier()),
        ('Naive Bayes', GaussianNB())
    ]
    cv = StratifiedKFold(n_splits=10, shuffle=True, random_state=42)
    scoring = ['accuracy', 'precision_weighted', 'recall_weighted', 'f1_weighted']

    x = df.drop('class', axis=1)
    y = df['class']
    results = {}
    for clf_name, clf in classifiers:
        out = cross_validate(clf, x, y, cv=cv, scoring=scoring)
        results[(clf_name, 'mean')] = {
            x.replace('_weighted', ''): out[f'test_{x}'].mean() for x in scoring
        }
        results[(clf_name, 'std')] = {
            x.replace('_weighted', ''): out[f'test_{x}'].std() for x in scoring
        }

    results_df = pd.DataFrame.from_dict(results, orient='index')
    return results_df


def main():
    '''Main function'''
    args = parse_args()

    fasta_dict_a = fasta_to_dict(args.a)
    fasta_dict_b = fasta_to_dict(args.b)

    if not set(fasta_dict_a.keys()).isdisjoint(fasta_dict_b.keys()):
        raise ValueError('Duplicate entries in input files')

    df_a = create_kmer_df(fasta_dict_a, args.a, args.k)
    df_b = create_kmer_df(fasta_dict_b, args.b, args.k)

    df = pd.concat([df_a, df_b])

    results = get_classification_results(df)
    print(results)


if __name__ == '__main__':
    main()
