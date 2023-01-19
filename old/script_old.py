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
    # Get the length of the aligned sequences by looking at the first sequence 
    seqLength = len(list(sequences.values())[0].sequence)


    # Sequences as list
    sequenceList = sequences.values()

    ### Calculating the "joint probability" of each pair of amino acids
    # Loop through each position in the sequence 
    pairCount = {}
    for i in range(seqLength):
        # Between each pair of non-same sequences, count the number of times each pair of amino acids occurs
        for seq1 in sequenceList:
            for seq2 in sequenceList:
                if seq1 != seq2:
                    # Keep track of count of the (sorted)pair of amino acids in dictionary, 
                    pair = (seq1.sequence[i], seq2.sequence[i])
                    pair = tuple(sorted(pair))
                    if pair not in pairCount:
                        pairCount[pair] = 0
                    pairCount[pair] += 1


    # Calculate total number of posible pairs, used later
    totalPairs = len(sequenceList) * (len(sequenceList) - 1) / 2

    # Normalise frequency of each pair of amino acids i.e (num. seq choose 2) 
    pairFreq = {}
    for pair in pairCount:
        pairFreq[pair] = pairCount[pair] / totalPairs
    
    print('Pair count:', pairCount)

    logOddsRatios = {}

    # Get probability of occurence of each amino acid in the pairs 
    aaFreq = {}
    for pair in pairFreq:
        for aa in pair:
            if aa not in aaFreq:
                aaFreq[aa] = 0
            aaFreq[aa] += pairFreq[pair]


    # print(logOddsRatio)
    return logOddsRatios


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
