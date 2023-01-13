# Given a file with multiple FASTA sequences, calculate the log odds ratio for each pair of amno acids

import argparse
import math

parser = argparse.ArgumentParser(description = 'Calculate the log odds ratio for each pair of amino acids')
parser.add_argument('filename', help = 'FASTA file with multiple sequences')
args = parser.parse_args()

# Simple sequence class
class Sequence:
    def __init__(self, header, sequence):
        self.header = header
        self.sequence = sequence

    def getHeader(self):
        return self.header

    def getSequence(self):
        return self.sequence


# FASTA reader that tracks sequence as Sequence objects
def readFASTA(filename):
    # Create a dictionary to store the sequences
    sequences = {}
    with open(filename, 'r') as f:
        lines = f.readlines()
        # Loop through the lines
        for line in lines:
            # If the line is a header, create a new sequence object
            if line[0] == '>':
                header = line[1:].strip()
                sequences[header] = Sequence(header, '')
            # If the line is a sequence, add it to the current sequence object
            else:
                sequences[header].sequence += line.strip()
    print('Read', len(sequences), 'sequences')
    return sequences

# Calculate log odds ratio
def calcLogOdds(sequences):
    # Get counts of each amino acid across all sequences 
    aaCount = {}
    for sequence in sequences.values():
        for aa in sequence.sequence:
            if aa not in aaCount:
                aaCount[aa] = 0
            aaCount[aa] += 1

    # Get the length of the aligned sequences by looking at the first sequence 
    seqLength = len(list(sequences.values())[0].sequence)


    # 
    sequenceList = sequences.values()

    ### Calculating the "joint probability" of each pair of amino acids
    # Loop through each position in the sequence 
    pairCount = {}
    for i in range(seqLength):
        # Between each pair of sequences, count the number of times each pair of amino acids occurs
        for seq1 in sequenceList:
            for seq2 in sequenceList:
                # Make sure we're not comparing the same sequence 
                if seq1 != seq2:
                    # Get the pair of amino acids at this position between the two sequences
                    pair = (seq1.sequence[i], seq2.sequence[i])
                    # If the pair is not in the dictionary, add it
                    if pair not in pairCount:
                        pairCount[pair] = 0
                    # Increment the count for this pair
                    pairCount[pair] += 1

    # Calculate log odds ratio 
    logOddsRatio = {}

    for pair in pairCount:
        # Calculate the log odds ratio = log (P(A,B) / (P(A) * P(B)))
        logOddsRatio[pair] = math.log(pairCount[pair] / (aaCount[pair[0]] * aaCount[pair[1]]))


    # print(logOddsRatio)
    return logOddsRatio


# Print the log odds ratio between every amino acid in a nice matrix
def printLogOddsRatio(logOddsRatio):
    # Get the list of amino acids contained from our fasta files
    aaList = []
    for pair in logOddsRatio:
        if pair[0] not in aaList:
            aaList.append(pair[0])
        if pair[1] not in aaList:
            aaList.append(pair[1])
    # Print the header
    print(' ', end = '')
    for aa in aaList:
        print('\t', aa, end = '')
    print()
    # Print the matrix
    for aa1 in aaList:
        print(aa1, end = '')
        for aa2 in aaList:
            # If the pair is in the dictionary, print the log odds ratio
            if (aa1, aa2) in logOddsRatio:
                print('\t', logOddsRatio[(aa1, aa2)], end = '')
            else:
                print('\t', 0, end = '')
            # If the pair is not in the dictionary, print 0
        print()



# Main function
def main():
    # Read the FASTA file 
    sequences = readFASTA(args.filename)
    # Calculate the log odds ratio for each pair of amino acids
    logOddsRatio = calcLogOdds(sequences)
    # Print log odds ratio between each pair of amino acids as a matrix
    printLogOddsRatio(logOddsRatio)

if __name__ == '__main__':
    main()
