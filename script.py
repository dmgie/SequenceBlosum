# Given a file with multiple FASTA sequences, calculate the log odds ratio for each pair of amno acids

import argparse
import math
from itertools import combinations

parser = argparse.ArgumentParser(description = 'Calculate the log odds ratio for each pair of amino acids')
parser.add_argument('filename', help = 'FASTA file')
parser.add_argument('temp-scale', help = 'Whether to enable the temperature scale', action = 'store_false')
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

    ### Joint probability of amino acids pairs and counts of each amino acids
    aaCount = {}
    pairCount = {}
    for i in range(seqLength):
        for seq in sequenceList:
            aa = seq.sequence[i]
            if aa not in aaCount:
                aaCount[aa] = 0
            aaCount[aa] += 1
        # Go through each pair of sequences and count pairs occurence
        for seq1, seq2 in combinations(sequenceList, 2):
            # Count the pairs
            aa1 = seq1.sequence[i]
            aa2 = seq2.sequence[i]
            pair = (aa1, aa2)
            pair = tuple(sorted(pair))
            if pair not in pairCount:
                pairCount[pair] = 0
            pairCount[pair] += 1

    
    totalAAs = sum(aaCount.values())
    aaFreq = {}
    for aa in aaCount:
        aaFreq[aa] = aaCount[aa] / totalAAs
    
    # Calculate total number of posible pairs => len * n(n-1)/2 (or rather n choose 2)
    # Normalise frequency of each pair of amino acids i.e (num. seq choose 2) = observed frequency
    totalPairs = seqLength * len(sequenceList) * (len(sequenceList) - 1) / 2
    pairFreq = {}
    for pair in pairCount:
        pairFreq[pair] = pairCount[pair] / totalPairs

    # print('AA count:', aaCount)
    # print('AA freq:', aaFreq)
    # print('Pair count:', pairCount)
    # print('Pair freq:', pairFreq)

    logOddsRatios = {}
    # Calculate log odds ratio for each pair of amino acids
    for pair in pairFreq:
        aa1 = pair[0]
        aa2 = pair[1]
        expectedFreq = 0
        if aa1 == aa2:
            expectedFreq = aaFreq[aa1] * aaFreq[aa2]
        else:
            expectedFreq = 2 * aaFreq[aa1] * aaFreq[aa2]
        logOddsRatio = 2 * math.log(pairFreq[pair] / (expectedFreq),2) 
        logOddsRatios[pair] = logOddsRatio

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
    aaList.sort()
    # Print the header
    print(' ', end = '')
    for aa in aaList:
        print('\t', aa, end = '')
    print()
    # Print the matrix
    for aa1 in aaList:
        print(aa1, end = '')
        for aa2 in aaList:
            # If the pair is in the dictionary, print the log odds ratio, rounded to an int
            if (aa1, aa2) in logOddsRatio:
                print('\t', int(round(logOddsRatio[(aa1, aa2)], 0)), end = '')
            else:
                print('\t', 0, end = '')
            # If the pair is not in the dictionary, print 0
        print()



# Main function
def main():
    sequences = readFASTA(args.filename)
    logOddsRatio = calcLogOdds(sequences)
    printLogOddsRatio(logOddsRatio)

if __name__ == '__main__':
    main()
