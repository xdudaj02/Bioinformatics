# Protein sequence classification
This script explores the possibility of classifying protein sequences based solely on their nucleotide n-mer composition. It takes two FASTA files as input, both containing protein sequences of a certain protein family and performs binary classification on them. Multiple classification algorithms are tested and compared using multiple metrics. The results are outputted in the form of a table. The script is written in Python and uses the scikit-learn library for machine learning algorithms.

### Usage
The script can be run using the following command:
```
python protein_classification.py <sequences_1.fasta> <sequences_2.fasta> <k-mer_length>
```
where `<sequences_1.fasta>` and `<sequences_2.fasta>` are the paths to the FASTA files containing the protein sequences of the two classes and `<k-mer_length>` is the length of the nucleotide n-mers to be used for feature extraction.

### Models
The script uses the following classification algorithms:
- Support Vector Machine (SVM)
- Random Forest
- Naive Bayes

### Evaluation
The performance of the models is evaluated using the following metrics:
- Accuracy
- Precision
- Recall
- F1 score

For each metric, the mean and standard deviation are calculated over 10-fold stratified cross-validation.

### Example output
```bash
$ python protein_classification.py data/sequences_1.fasta data/sequences_2.fasta 2
```
```
                    accuracy  precision    recall        f1
SVM           mean  0.993769   0.993933  0.993769  0.993754
              std   0.004845   0.004659  0.004845  0.004848
Random Forest mean  0.995148   0.995312  0.995148  0.995132
              std   0.006236   0.006029  0.006236  0.006265
Naive Bayes   mean  0.961231   0.967829  0.961231  0.962566
              std   0.013602   0.010131  0.013602  0.012873
```
